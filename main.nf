nextflow.enable.dsl = 2

include { PANGROWTH } from './modules/pangrowth'
include { PLOT } from './modules/plotting'
include { PLOT_ALL } from './modules/plotting_all'


def datasetStem(inputPath, inputType) {
    def stem = inputPath.name
    if (inputType != 'directory') {
        stem = stem.replaceFirst(/(?i)\.(tar\.gz|tgz|zip|txt|list)$/, '')
    }
    stem = stem.replaceAll(/[^A-Za-z0-9._-]+/, '_')
        .replaceAll(/^[-_.]+/, '')
        .replaceAll(/[-_.]+$/, '')
    return stem ?: 'pangenome'
}

def isFastaPath(path) {
    return path.name ==~ /(?i).+\.(fa|fasta|fna|ffn)(\.gz)?/
}

def readFastaList(inputList) {
    def inputUri = inputList.toUri()
    def genomes = []
    def seen = new HashSet()
    inputList.readLines().eachWithIndex { rawLine, index ->
        def entry = rawLine.trim()
        if (entry && !entry.startsWith('#')) {
            if (entry.contains('://') && !entry.startsWith('s3://')) {
                error("FASTA list line ${index + 1}: use a local path or an s3:// path: '${entry}'.")
            }
            if (inputUri.scheme == 's3' && entry.startsWith('/')) {
                error("FASTA list line ${index + 1}: an S3 list must reference S3 files, not local paths: '${entry}'.")
            }
            if (!(entry ==~ /(?i).+\.(fa|fasta|fna|ffn)(\.gz)?/)) {
                error("FASTA list line ${index + 1}: unsupported FASTA filename: '${entry}'.")
            }
            def candidates
            if (entry.startsWith('/') || entry.startsWith('s3://')) {
                candidates = [file(entry, glob: false)]
            } else {
                candidates = [inputList.parent.resolve(entry)]
                if (inputUri.scheme == 's3') {
                    // S3Path.root is the bucket, independent of how toUri()
                    // represents the endpoint and bucket (including s3:///...).
                    candidates.add(inputList.root.resolve(entry))
                }
            }
            candidates = candidates.collect { it.normalize() }.unique()
            def matches = candidates.findAll { it.exists() && !it.isDirectory() }
            if (matches.isEmpty()) {
                error("FASTA list line ${index + 1}: file not found or not readable: '${entry}'. Checked: ${candidates.collect { it.toUri() }.join(', ')}")
            }
            if (matches.size() > 1) {
                error("FASTA list line ${index + 1}: ambiguous path '${entry}' exists relative to both the list and bucket root. Use a full s3://BUCKET/key path in the list.")
            }
            def genome = matches[0]
            def uri = genome.toUri()
            if (!seen.add(uri.toString())) {
                error("FASTA list line ${index + 1}: duplicate genome path: '${entry}'.")
            }
            genomes.add(genome)
        }
    }
    if (genomes.size() < 3) {
        error("The FASTA list must contain at least three genomes; found ${genomes.size()}.")
    }
    return genomes
}

def readFastaDirectory(inputDir) {
    // A single '*' deliberately limits directory inputs to immediate children.
    def resolved = file(inputDir.resolve('*').toString())
    def candidates = resolved instanceof Collection ? resolved : [resolved]
    def genomes = candidates
        .findAll { it.exists() && !it.isDirectory() && isFastaPath(it) }
        .unique { it.toUri().toString() }
        .sort { left, right -> left.toUri().toString() <=> right.toUri().toString() }
    if (genomes.size() < 3) {
        error("Input directory '${inputDir}' must contain at least three FASTA files directly inside it; found ${genomes.size()}. Subdirectories are not searched.")
    }
    return genomes
}


workflow {
    if (!params.input) {
        error('Provide --input with a directory, archive, FASTA list, or a pattern matching several such inputs.')
    }

    def rawInputs = params.input instanceof Collection ? params.input : [params.input]
    def inputPatterns = rawInputs.collectMany { rawInput ->
        rawInput.toString().split(/[\s,]+/).findAll { it }
    }
    def selectedInputs = inputPatterns.collectMany { pattern ->
        def resolved = file(pattern.toString())
        def candidates = resolved instanceof Collection ? resolved : [resolved]
        return candidates.findAll { it.exists() }
    }.unique { it.toUri().toString() }
        .sort { left, right -> left.toUri().toString() <=> right.toUri().toString() }
    if (selectedInputs.isEmpty()) {
        error("Input did not match any readable files: '${params.input}'.")
    }

    def usedIds = [:]
    def datasets = selectedInputs.collect { inputPath ->
        def lowerName = inputPath.name.toLowerCase()
        def inputType
        if (inputPath.isDirectory()) {
            inputType = 'directory'
        } else if (lowerName.endsWith('.txt') || lowerName.endsWith('.list')) {
            inputType = 'list'
        } else if (lowerName.endsWith('.zip') || lowerName.endsWith('.tar.gz') || lowerName.endsWith('.tgz')) {
            inputType = 'archive'
        } else {
            error("Unsupported input '${inputPath.name}': use a directory, a .txt or .list FASTA list, or a .zip, .tar.gz, or .tgz archive.")
        }

        def baseId = datasetStem(inputPath, inputType)
        if (baseId == 'all') {
            baseId = 'all_input'
        }
        def occurrence = (usedIds[baseId] ?: 0) + 1
        usedIds[baseId] = occurrence
        def datasetId = occurrence == 1 ? baseId : "${baseId}_${occurrence}"
        def inputFiles
        if (inputType == 'list') {
            inputFiles = readFastaList(inputPath)
        } else if (inputType == 'directory') {
            inputFiles = readFastaDirectory(inputPath)
        } else {
            inputFiles = [inputPath]
        }
        return tuple(datasetId, inputPath.name, inputType, inputFiles)
    }

    PANGROWTH(Channel.fromList(datasets))

    PLOT(PANGROWTH.out.results)

    if (datasets.size() > 1) {
        def combinedPlotInput = PANGROWTH.out.results
            .collect(flat: false)
            .map { resultRows ->
                def rows = resultRows.sort { left, right -> left[0] <=> right[0] }
                tuple(
                    rows.collect { it[0] },
                    rows.collect { it[1] },
                    rows.collect { it[2] },
                    rows.collect { it[3] },
                    rows.collect { it[4] },
                    rows.collect { it[5] }
                )
            }
        PLOT_ALL(combinedPlotInput)
    }
}
