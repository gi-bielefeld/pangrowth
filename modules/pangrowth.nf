def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
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
        path('pangrowth_hill.tsv'), emit: results
    path 'pangrowth.log'

    script:
    def stagedFiles = inputFiles instanceof List ? inputFiles : [inputFiles]
    def prepareInput
    if (inputType == 'list' || inputType == 'directory') {
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
    """
}
