# Output files

## Numerical results

* `pangrowth_hist.txt` — one line per genome frequency. Line *i* is the number
  of distinct *k*-mers present in exactly *i* genomes.
* `pangrowth_hist_infix.txt` — infix-equivalent frequency histogram used for
  compacted de Bruijn graph diversity. It is generated only when CDBG diversity
  is enabled.
* `pangrowth_growth.txt` — exact expected pangenome size for sample sizes from
  one through the total number of genomes. Values are written one per line.
* `pangrowth_core.txt` — expected strict or quorum-based core size over the same
  range of sample sizes, written one value per line.
* `pangrowth_hill.tsv` — tab-separated Hill-number estimates. Columns are
  `fit`, `m`, `richness`, `exp_entropy`, and `inv_gini_simp`. `fit` identifies
  interpolation (`int`), the observed sample (`obs`), or extrapolation (`ext`),
  and `m` is the number of genomes.

## Visualisations and fits

These files are present when PDF plotting is enabled and the corresponding fit
succeeds:

* `pangrowth_hist.pdf` — distribution of *k*-mers by genome frequency.
* `pangrowth_growth.pdf` — expected pangenome growth and new-item curves with
  power-law fits.
* `pangrowth_core.pdf` — expected core-genome curve with an exponential-decay
  fit.
* `pangrowth_growth_fit.txt` — fitted growth and average-new-item equations.
* `pangrowth_core_fit.txt` — fitted asymptotic core equation and the predicted
  core fraction of an average genome.

## Log

`pangrowth.log` records the effective workflow parameters, discovered FASTA
files, progress messages, diagnostic output from pangrowth, and any plotting
warnings.
