nextflow.enable.dsl = 2

def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
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
    tag "${inputName}"
    label 'highmemMedium'
    container 'ghcr.io/gi-bielefeld/pangrowth:clowm-v0.1.0'

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(inputName), val(inputType), path(inputFiles, stageAs: 'source????/*')

    output:
    path 'pangrowth_hist.txt'
    path 'pangrowth_hist_infix.txt', optional: true
    path 'pangrowth_growth.txt'
    path 'pangrowth_core.txt'
    path 'pangrowth_hill.tsv'
    path 'pangrowth_hill.pdf', optional: true
    path 'pangrowth_hist.pdf', optional: true
    path 'pangrowth_growth.pdf', optional: true
    path 'pangrowth_core.pdf', optional: true
    path 'pangrowth_growth_fit.txt', optional: true
    path 'pangrowth_core_fit.txt', optional: true
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
    def telomereArg = params.cdbg && params.account_telomeres ? '-T' : ''
    def cdbgArg = params.cdbg ? '--cdbg -o pangrowth' : ''
    def histRedirect = params.cdbg ? '' : '> pangrowth_hist.txt'
    def hillInputs = params.cdbg ? 'pangrowth_hist.txt pangrowth_hist_infix.txt' : 'pangrowth_hist.txt'
    def hillTailArg = params.cdbg && !params.hill_force_exact ? "-B ${params.hill_exact_tail}" : ''
    def hillToleranceArg = params.cdbg && !params.hill_force_exact && params.hill_adaptive_tolerance != null \
        ? "-A ${params.hill_adaptive_tolerance}" : ''
    def hillExactArg = params.cdbg && params.hill_force_exact ? '-E' : ''
    def plotNormXArg = params.histogram_normalize_x ? '--norm_x' : ''
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
        echo "CDBG diversity: ${params.cdbg}"
        echo "Core quorum: ${params.quorum}"
        echo
        echo "Input files:"
        sed 's/^/  /' fasta_files.list
        echo
    } > pangrowth.log

    echo "Running histogram calculation" >> pangrowth.log
    pangrowth hist ${cdbgArg} \
        -k ${params.kmer} \
        -t ${task.cpus} \
        -s ${params.suffix_size} \
        -c ${params.minimum_count} \
        ${canonicalArg} ${telomereArg} \
        -i fasta_files.list ${histRedirect} 2>> pangrowth.log

    echo "Running pangenome growth calculation" >> pangrowth.log
    pangrowth growth -h pangrowth_hist.txt 2>> pangrowth.log \
        | awk '{ print \$NF }' > pangrowth_growth.txt

    echo "Running pangenome core calculation" >> pangrowth.log
    pangrowth core -q ${params.quorum} -h pangrowth_hist.txt 2>> pangrowth.log \
        | awk '{ print \$NF }' > pangrowth_core.txt

    echo "Running Hill-number calculation" >> pangrowth.log
    pangrowth hill \
        -p ${params.hill_points} \
        ${hillTailArg} ${hillToleranceArg} ${hillExactArg} \
        ${hillInputs} > pangrowth_hill.tsv 2>> pangrowth.log

    if ${params.create_plots}; then
        echo "Creating plots" >> pangrowth.log
        if ! plot_hill.py pangrowth_hill.tsv pangrowth_hill.pdf >> pangrowth.log 2>&1; then
            echo "WARNING: Hill-number plot generation failed." >> pangrowth.log
            rm -f pangrowth_hill.pdf
        fi
        if ! plot_hist.py ${plotNormXArg} --norm_y=${params.histogram_normalize_y} \
            pangrowth_hist.txt pangrowth_hist.pdf >> pangrowth.log 2>&1; then
            echo "WARNING: Histogram plot generation failed." >> pangrowth.log
        fi
        if ! plot_growth.py pangrowth_growth.txt pangrowth_growth.pdf \
            > pangrowth_growth_fit.txt 2>> pangrowth.log; then
            echo "WARNING: Growth plot generation failed." >> pangrowth.log
            rm -f pangrowth_growth.pdf pangrowth_growth_fit.txt
        fi
        if ! plot_core.py pangrowth_core.txt pangrowth_core.pdf \
            > pangrowth_core_fit.txt 2>> pangrowth.log; then
            echo "WARNING: Core plot generation failed." >> pangrowth.log
            rm -f pangrowth_core.pdf pangrowth_core_fit.txt
        fi
    fi
    """
}


workflow {
    if (!params.input) {
        error('Provide --input with an archive or FASTA list file.')
    }
    def inputFile = file(params.input, glob: false)
    if (!inputFile.exists() || inputFile.isDirectory()) {
        error("Invalid input: '${params.input}' is not a file.")
    }

    def inputName = inputFile.name.toLowerCase()
    def inputType
    if (inputName.endsWith('.txt') || inputName.endsWith('.list')) {
        inputType = 'list'
    } else if (inputName.endsWith('.zip') || inputName.endsWith('.tar.gz') || inputName.endsWith('.tgz')) {
        inputType = 'archive'
    } else {
        error("Unsupported input '${inputFile.name}': use a .txt or .list FASTA list, or a .zip, .tar.gz, or .tgz archive.")
    }

    def inputFiles = inputType == 'list' \
        ? readFastaList(inputFile) : [inputFile]
    PANGROWTH(Channel.value(tuple(inputFile.name, inputType, inputFiles)))
}
