nextflow.enable.dsl = 2

include { PANGROWTH } from './modules/pangrowth'
include { PLOT } from './modules/plotting'


def datasetStem(inputFile) {
    def stem = inputFile.name.replaceFirst(/(?i)\.(tar\.gz|tgz|zip|txt|list)$/, '')
    stem = stem.replaceAll(/[^A-Za-z0-9._-]+/, '_')
        .replaceAll(/^[-_.]+/, '')
        .replaceAll(/[-_.]+$/, '')
    return stem ?: 'pangenome'
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


workflow {
    if (!params.input) {
        error('Provide --input with an archive, FASTA list, or a pattern matching several such files.')
    }

    def rawInputs = params.input instanceof Collection ? params.input : [params.input]
    def inputPatterns = rawInputs.collectMany { rawInput ->
        rawInput.toString().split(/[\s,]+/).findAll { it }
    }
    def selectedInputFiles = inputPatterns.collectMany { pattern ->
        def resolved = file(pattern.toString())
        def candidates = resolved instanceof Collection ? resolved : [resolved]
        return candidates.findAll { it.exists() && !it.isDirectory() }
    }.unique { it.toUri().toString() }
        .sort { left, right -> left.toUri().toString() <=> right.toUri().toString() }
    if (selectedInputFiles.isEmpty()) {
        error("Input did not match any readable files: '${params.input}'.")
    }

    def usedIds = [:]
    def datasets = selectedInputFiles.collect { inputFile ->
        def lowerName = inputFile.name.toLowerCase()
        def inputType
        if (lowerName.endsWith('.txt') || lowerName.endsWith('.list')) {
            inputType = 'list'
        } else if (lowerName.endsWith('.zip') || lowerName.endsWith('.tar.gz') || lowerName.endsWith('.tgz')) {
            inputType = 'archive'
        } else {
            error("Unsupported input '${inputFile.name}': use a .txt or .list FASTA list, or a .zip, .tar.gz, or .tgz archive.")
        }

        def baseId = datasetStem(inputFile)
        if (baseId == 'all') {
            baseId = 'all_input'
        }
        def occurrence = (usedIds[baseId] ?: 0) + 1
        usedIds[baseId] = occurrence
        def datasetId = occurrence == 1 ? baseId : "${baseId}_${occurrence}"
        def inputFiles = inputType == 'list' ? readFastaList(inputFile) : [inputFile]
        return tuple(datasetId, inputFile.name, inputType, inputFiles)
    }

    PANGROWTH(Channel.fromList(datasets))

    // PLOT accepts lists of labels and result files, allowing the same process
    // to create both per-pangenome and combined visualisations.
    def plotInputs = PANGROWTH.out.results.map { datasetId, histFile, growthFile, coreFile, quorumFile, hillFile ->
        tuple(
            datasetId,
            [datasetId],
            [histFile],
            [growthFile],
            [coreFile],
            [quorumFile],
            [hillFile]
        )
    }

    if (datasets.size() > 1) {
        def combinedPlotInput = PANGROWTH.out.results
            .collect(flat: false)
            .map { resultRows ->
                def rows = resultRows.sort { left, right -> left[0] <=> right[0] }
                tuple(
                    'all',
                    rows.collect { it[0] },
                    rows.collect { it[1] },
                    rows.collect { it[2] },
                    rows.collect { it[3] },
                    rows.collect { it[4] },
                    rows.collect { it[5] }
                )
            }
        plotInputs = plotInputs.concat(combinedPlotInput)
    }

    PLOT(plotInputs)
}
