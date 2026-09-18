def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

process PLOT {
    tag "${outputId == 'all' ? 'all pangenomes' : outputId}"
    container 'ghcr.io/gi-bielefeld/pangrowth:clowm'

    publishDir params.outdir, mode: 'copy', \
        saveAs: { filename -> "${outputId}/${filename}" }

    input:
    tuple val(outputId), val(datasetIds), path(histFiles, stageAs: 'hist????/*'), \
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
    path 'pangrowth_plot.log'

    script:
    def labelArgs = datasetIds.collect { "--label ${shellQuote(it)}" }.join(' ')
    def histArgs = histFiles.collect { shellQuote(it) }.join(' ')
    def growthArgs = growthFiles.collect { shellQuote(it) }.join(' ')
    def coreArgs = coreFiles.collect { shellQuote(it) }.join(' ')
    def quorumArgs = quorumFiles.collect { shellQuote(it) }.join(' ')
    def hillArgs = hillFiles.collect { shellQuote(it) }.join(' ')
    def scope = outputId == 'all' ? 'Combined plots' : 'Plots'
    """
    set -euo pipefail

    echo "${scope} for: ${datasetIds.join(', ')}" > pangrowth_plot.log

    if ! plot_hist.py ${labelArgs} ${histArgs} pangrowth_hist.pdf >> pangrowth_plot.log 2>&1; then
        echo "WARNING: Histogram plot generation failed." >> pangrowth_plot.log
        rm -f pangrowth_hist.pdf
    fi
    if ! plot_hist.py --norm_x ${labelArgs} ${histArgs} pangrowth_hist_percentage.pdf >> pangrowth_plot.log 2>&1; then
        echo "WARNING: Percentage histogram plot generation failed." >> pangrowth_plot.log
        rm -f pangrowth_hist_percentage.pdf
    fi
    if ! plot_growth.py ${labelArgs} ${growthArgs} pangrowth_growth.pdf \
        > pangrowth_growth_fit.txt 2>> pangrowth_plot.log; then
        echo "WARNING: Growth plot generation failed." >> pangrowth_plot.log
        rm -f pangrowth_growth.pdf pangrowth_growth_fit.txt
    fi
    if ! plot_core.py ${labelArgs} ${coreArgs} pangrowth_core.pdf \
        > pangrowth_core_fit.txt 2>> pangrowth_plot.log; then
        echo "WARNING: Core plot generation failed." >> pangrowth_plot.log
        rm -f pangrowth_core.pdf pangrowth_core_fit.txt
    fi
    if ! plot_core.py ${labelArgs} ${quorumArgs} pangrowth_quorum.pdf \
        > pangrowth_quorum_fit.txt 2>> pangrowth_plot.log; then
        echo "WARNING: Quorum plot generation failed." >> pangrowth_plot.log
        rm -f pangrowth_quorum.pdf pangrowth_quorum_fit.txt
    fi
    if ! plot_hill.py ${labelArgs} ${hillArgs} pangrowth_hill.pdf >> pangrowth_plot.log 2>&1; then
        echo "WARNING: Hill-number plot generation failed." >> pangrowth_plot.log
        rm -f pangrowth_hill.pdf
    fi
    """
}
