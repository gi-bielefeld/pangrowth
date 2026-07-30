nextflow.enable.dsl = 2

process PANGROWTH {
    tag "${archive.simpleName}"
    label 'highmemMedium'
    container 'ghcr.io/gi-bielefeld/pangrowth:clowm-v0.1.0'

    publishDir params.outdir, mode: 'copy'

    input:
    path archive

    output:
    path 'pangrowth_hist.txt'
    path 'pangrowth_hist_infix.txt', optional: true
    path 'pangrowth_growth.txt'
    path 'pangrowth_core.txt'
    path 'pangrowth_hill.tsv'
    path 'pangrowth_hist.pdf', optional: true
    path 'pangrowth_growth.pdf', optional: true
    path 'pangrowth_core.pdf', optional: true
    path 'pangrowth_growth_fit.txt', optional: true
    path 'pangrowth_core_fit.txt', optional: true
    path 'pangrowth.log'

    script:
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

    mkdir input_files
    case "${archive.name}" in
        *.zip)
            python -m zipfile -e "${archive}" input_files
            ;;
        *.tar.gz|*.tgz)
            tar -xzf "${archive}" -C input_files
            ;;
        *)
            echo "ERROR: --input must be a .zip, .tar.gz, or .tgz archive." >&2
            exit 1
            ;;
    esac

    find input_files -type f \\( -iname '*.fa' -o -iname '*.fa.gz' -o -iname '*.fasta' -o -iname '*.fasta.gz' -o -iname '*.fna' -o -iname '*.fna.gz' -o -iname '*.ffn' -o -iname '*.ffn.gz' \\) -print | LC_ALL=C sort > fasta_files.list

    genome_count=\$(wc -l < fasta_files.list)
    if [ "\${genome_count}" -lt 3 ]; then
        echo "ERROR: The archive must contain at least three FASTA files; found \${genome_count}." >&2
        exit 1
    fi

    {
        echo "pangrowth CloWM workflow"
        echo "Input archive: ${archive.name}"
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
    // CloWM inputs must be real files rather than S3 directory prefixes.
    // Requiring an archive also lets Nextflow stage the complete genome
    // collection safely.
    def inputArchive = file(params.input)
    if (!inputArchive.exists() || inputArchive.isDirectory()) {
        error("Invalid input: '${params.input}' is not a valid archive file for --input.")
    }

    input_archive_ch = Channel.value(inputArchive)
    PANGROWTH(input_archive_ch)
}
