nextflow.enable.dsl = 2

def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

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

process PANGROWTH {
    tag "${datasetId}"
    label 'highmemMedium'
    container 'ghcr.io/gi-bielefeld/pangrowth:clowm'

    publishDir params.outdir, mode: 'copy', \
        saveAs: { filename -> "${datasetId}/${filename}" }

    input:
    tuple val(datasetId), val(inputName), val(inputType), path(inputFiles, stageAs: 'source????/*')

    output:
    tuple val(datasetId), path('pangrowth_hist.txt'), path('pangrowth_growth.txt'), \
        path('pangrowth_core.txt'), path('pangrowth_quorum.txt'), \
        path('pangrowth_hill.tsv'), emit: comparison
    path 'pangrowth_hill.pdf', optional: true
    path 'pangrowth_hist.pdf', optional: true
    path 'pangrowth_hist_percentage.pdf', optional: true
    path 'pangrowth_growth.pdf', optional: true
    path 'pangrowth_core.pdf', optional: true
    path 'pangrowth_quorum.pdf', optional: true
    path 'pangrowth_growth_fit.txt', optional: true
    path 'pangrowth_core_fit.txt', optional: true
    path 'pangrowth_quorum_fit.txt', optional: true
    path 'pangrowth.log'

    script:
    def stagedFiles = inputFiles instanceof List ? inputFiles : [inputFiles]
    def prepareInput
    if (inputType == 'list') {
        // Use Nextflow's staged paths, including the numbered directories that
        // prevent collisions when multiple genomes have the same filename.
        prepareInput = "printf '%s\\n' ${stagedFiles.collect { shellQuote(it) }.join(' ')} > fasta_files.list"
    } else {
        def archive = stagedFiles[0]
        prepareInput = """
        mkdir input_files
        case ${shellQuote(archive.name.toLowerCase())} in
            *.zip)
                python -m zipfile -e ${shellQuote(archive)} input_files
                ;;
            *.tar.gz|*.tgz)
                tar -xzf ${shellQuote(archive)} -C input_files
                ;;
            *)
                echo "ERROR: Archive input must be a .zip, .tar.gz, or .tgz file." >&2
                exit 1
                ;;
        esac
        find input_files -type f \\( -iname '*.fa' -o -iname '*.fa.gz' -o -iname '*.fasta' -o -iname '*.fasta.gz' -o -iname '*.fna' -o -iname '*.fna.gz' -o -iname '*.ffn' -o -iname '*.ffn.gz' \\) -print | LC_ALL=C sort > fasta_files.list
        """
    }
    def canonicalArg = params.canonical ? '' : '-b'
    def plotLabelArg = "--label ${shellQuote(datasetId)}"
    def quorumCalculation = (params.quorum as double) == 1.0d ? """
    echo "Quorum is 1.0; reusing the strict-core curve" >> pangrowth.log
    cp pangrowth_core.txt pangrowth_quorum.txt
    """ : """
    echo "Running quorum-core calculation" >> pangrowth.log
    pangrowth core -q ${params.quorum} -h pangrowth_hist.txt 2>> pangrowth.log \\
        | awk '{ print \$NF }' > pangrowth_quorum.txt
    """
    """
    set -euo pipefail

    ${prepareInput}

    genome_count=\$(wc -l < fasta_files.list)
    if [ "\${genome_count}" -lt 3 ]; then
        echo "ERROR: The input must contain at least three FASTA files; found \${genome_count}." >&2
        exit 1
    fi

    {
        echo "pangrowth CloWM workflow"
        printf 'Input: %s (%s)\\n' ${shellQuote(inputName)} ${shellQuote(inputType)}
        echo "Genome files: \${genome_count}"
        echo "k-mer size: ${params.kmer}"
        echo "Minimum within-genome count: ${params.minimum_count}"
        echo "Canonical k-mers: ${params.canonical}"
        echo "Worker threads: ${task.cpus}"
        echo "Hashtable suffix bits: ${params.suffix_size}"
        echo "Quorum-curve threshold: ${params.quorum}"
        echo
        echo "Input files:"
        sed 's/^/  /' fasta_files.list
        echo
    } > pangrowth.log

    echo "Running histogram calculation" >> pangrowth.log
    pangrowth hist \
        -k ${params.kmer} \
        -t ${task.cpus} \
        -s ${params.suffix_size} \
        -c ${params.minimum_count} \
        ${canonicalArg} \
        -i fasta_files.list > pangrowth_hist.txt 2>> pangrowth.log

    echo "Running pangenome growth calculation" >> pangrowth.log
    pangrowth growth -h pangrowth_hist.txt 2>> pangrowth.log \
        | awk '{ print \$NF }' > pangrowth_growth.txt

    echo "Running strict core calculation" >> pangrowth.log
    pangrowth core -h pangrowth_hist.txt 2>> pangrowth.log \
        | awk '{ print \$NF }' > pangrowth_core.txt

    ${quorumCalculation}

    echo "Running Hill-number calculation" >> pangrowth.log
    pangrowth hill \
        -p ${params.hill_points} \
        pangrowth_hist.txt > pangrowth_hill.tsv 2>> pangrowth.log

    echo "Creating plots" >> pangrowth.log
    if ! plot_hill.py ${plotLabelArg} pangrowth_hill.tsv pangrowth_hill.pdf >> pangrowth.log 2>&1; then
        echo "WARNING: Hill-number plot generation failed." >> pangrowth.log
        rm -f pangrowth_hill.pdf
    fi
    if ! plot_hist.py ${plotLabelArg} \
        pangrowth_hist.txt pangrowth_hist.pdf >> pangrowth.log 2>&1; then
        echo "WARNING: Histogram plot generation failed." >> pangrowth.log
        rm -f pangrowth_hist.pdf
    fi
    if ! plot_hist.py --norm_x ${plotLabelArg} \
        pangrowth_hist.txt pangrowth_hist_percentage.pdf >> pangrowth.log 2>&1; then
        echo "WARNING: Percentage histogram plot generation failed." >> pangrowth.log
        rm -f pangrowth_hist_percentage.pdf
    fi
    if ! plot_growth.py ${plotLabelArg} pangrowth_growth.txt pangrowth_growth.pdf \
        > pangrowth_growth_fit.txt 2>> pangrowth.log; then
        echo "WARNING: Growth plot generation failed." >> pangrowth.log
        rm -f pangrowth_growth.pdf pangrowth_growth_fit.txt
    fi
    if ! plot_core.py ${plotLabelArg} pangrowth_core.txt pangrowth_core.pdf \
        > pangrowth_core_fit.txt 2>> pangrowth.log; then
        echo "WARNING: Core plot generation failed." >> pangrowth.log
        rm -f pangrowth_core.pdf pangrowth_core_fit.txt
    fi
    if ! plot_core.py ${plotLabelArg} pangrowth_quorum.txt pangrowth_quorum.pdf \
        > pangrowth_quorum_fit.txt 2>> pangrowth.log; then
        echo "WARNING: Quorum plot generation failed." >> pangrowth.log
        rm -f pangrowth_quorum.pdf pangrowth_quorum_fit.txt
    fi
    """
}

process COMPARE_PLOTS {
    tag 'all pangenomes'
    container 'ghcr.io/gi-bielefeld/pangrowth:clowm'

    publishDir "${params.outdir}/all", mode: 'copy'

    input:
    tuple val(datasetIds), path(histFiles, stageAs: 'hist????/*'), \
        path(growthFiles, stageAs: 'growth????/*'), \
        path(coreFiles, stageAs: 'core????/*'), \
        path(quorumFiles, stageAs: 'quorum????/*'), \
        path(hillFiles, stageAs: 'hill????/*')

    output:
    path 'pangrowth_hist.pdf', optional: true
    path 'pangrowth_hist_percentage.pdf', optional: true
    path 'pangrowth_growth.pdf', optional: true
    path 'pangrowth_core.pdf', optional: true
    path 'pangrowth_quorum.pdf', optional: true
    path 'pangrowth_hill.pdf', optional: true
    path 'pangrowth_growth_fit.txt', optional: true
    path 'pangrowth_core_fit.txt', optional: true
    path 'pangrowth_quorum_fit.txt', optional: true
    path 'pangrowth.log'

    script:
    def labelArgs = datasetIds.collect { "--label ${shellQuote(it)}" }.join(' ')
    def histArgs = histFiles.collect { shellQuote(it) }.join(' ')
    def growthArgs = growthFiles.collect { shellQuote(it) }.join(' ')
    def coreArgs = coreFiles.collect { shellQuote(it) }.join(' ')
    def quorumArgs = quorumFiles.collect { shellQuote(it) }.join(' ')
    def hillArgs = hillFiles.collect { shellQuote(it) }.join(' ')
    """
    set -euo pipefail

    echo "Combined plots for: ${datasetIds.join(', ')}" > pangrowth.log

    if ! plot_hist.py ${labelArgs} ${histArgs} pangrowth_hist.pdf >> pangrowth.log 2>&1; then
        echo "WARNING: Combined histogram plot generation failed." >> pangrowth.log
        rm -f pangrowth_hist.pdf
    fi
    if ! plot_hist.py --norm_x ${labelArgs} ${histArgs} pangrowth_hist_percentage.pdf >> pangrowth.log 2>&1; then
        echo "WARNING: Combined percentage histogram plot generation failed." >> pangrowth.log
        rm -f pangrowth_hist_percentage.pdf
    fi
    if ! plot_growth.py ${labelArgs} ${growthArgs} pangrowth_growth.pdf \
        > pangrowth_growth_fit.txt 2>> pangrowth.log; then
        echo "WARNING: Combined growth plot generation failed." >> pangrowth.log
        rm -f pangrowth_growth.pdf pangrowth_growth_fit.txt
    fi
    if ! plot_core.py ${labelArgs} ${coreArgs} pangrowth_core.pdf \
        > pangrowth_core_fit.txt 2>> pangrowth.log; then
        echo "WARNING: Combined core plot generation failed." >> pangrowth.log
        rm -f pangrowth_core.pdf pangrowth_core_fit.txt
    fi
    if ! plot_core.py ${labelArgs} ${quorumArgs} pangrowth_quorum.pdf \
        > pangrowth_quorum_fit.txt 2>> pangrowth.log; then
        echo "WARNING: Combined quorum plot generation failed." >> pangrowth.log
        rm -f pangrowth_quorum.pdf pangrowth_quorum_fit.txt
    fi
    if ! plot_hill.py ${labelArgs} ${hillArgs} pangrowth_hill.pdf >> pangrowth.log 2>&1; then
        echo "WARNING: Combined Hill-number plot generation failed." >> pangrowth.log
        rm -f pangrowth_hill.pdf
    fi
    """
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

    if (datasets.size() > 1) {
        def comparisonInputs = PANGROWTH.out.comparison
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
        COMPARE_PLOTS(comparisonInputs)
    }
}
