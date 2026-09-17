# Output files

Each input collection has its own output folder, named after its list or
archive file. If several collections are supplied, `all/` contains combined
versions of the plots and fit summaries below. The `all/` folder does not
duplicate numerical result tables.

## Numerical results

* `pangrowth_hist.txt` — one line per genome frequency. Line *i* is the number
  of distinct *k*-mers present in exactly *i* genomes.
* `pangrowth_growth.txt` — exact expected pangenome size for sample sizes from
  one through the total number of genomes. Values are written one per line.
* `pangrowth_core.txt` — expected strict-core size over the same range of sample
  sizes, written one value per line.
* `pangrowth_quorum.txt` — expected quorum-core size using the configured
  quorum threshold, written one value per line. At a quorum of 1.0, this reuses
  the strict-core result.
* `pangrowth_hill.tsv` — tab-separated Hill-number estimates. Columns are
  `fit`, `m`, `richness`, `exp_entropy`, and `inv_simpson`. `fit` identifies
  interpolation (`int`), the observed sample (`obs`), or extrapolation (`ext`),
  and `m` is the number of genomes.

## Visualisations and fits

These files are present when the automatically run plot generation and the
corresponding fit succeed:

* `pangrowth_hist.pdf` — distribution of *k*-mers by absolute genome
  frequency. Tick labels are thinned automatically while retaining the first
  and final genome counts.
* `pangrowth_hist_percentage.pdf` — the same raw *k*-mer counts grouped into
  5%, 10%, ..., 100% genome-frequency bins.
* `pangrowth_growth.pdf` — expected pangenome growth and new-item curves with
  power-law fits.
* `pangrowth_core.pdf` — expected core-genome curve with an exponential-decay
  fit.
* `pangrowth_quorum.pdf` — expected quorum-core curve with an
  exponential-decay fit.
* `pangrowth_hill.pdf` — Hill-number richness, exponential entropy, and inverse
  Simpson curves. Solid lines show interpolation, points mark the observed
  sample, and dashed lines show extrapolation. All panels use the richness
  range on the y-axis.
* `pangrowth_growth_fit.txt` — fitted growth and average-new-item equations.
* `pangrowth_core_fit.txt` — fitted asymptotic core equation and the predicted
  core fraction of an average genome.
* `pangrowth_quorum_fit.txt` — fitted asymptotic quorum-core equation and the
  predicted quorum-core fraction of an average genome.

## Log

`pangrowth.log` records the effective workflow parameters, discovered FASTA
files, progress messages, diagnostic output from pangrowth, and any plotting
warnings.
