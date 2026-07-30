# Usage

## Input / Output

### Input archive

The input must be a `.zip`, `.tar.gz`, or `.tgz` archive containing at least
three genomes. Each genome must be stored in a separate FASTA file. Supported
file-name endings are `.fa`, `.fasta`, `.fna`, and `.ffn`; each may additionally
end in `.gz`.

Files may contain multiple records, such as the contigs of one assembly, and
may be arranged in subdirectories. Non-FASTA files are ignored. A FASTA file is
the sampling unit: do not split one genome over several files, because those
files would be interpreted as separate genomes.

### Output directory

All numerical tables, optional PDF plots, fit summaries, and `pangrowth.log`
are written to the selected output directory.

## Analysis parameters

### *k*-mer length

The default length is 17. Short or heterogeneous sequences may benefit from a
smaller value. Values from 5 through 31 are accepted by this workflow.

### Minimum within-genome count

By default, a *k*-mer occurring at least once in a genome is considered present.
Increasing this threshold can filter low-abundance *k*-mers introduced by
sequencing errors when read data are used.

### Core-genome quorum

The quorum is the minimum fraction of genomes in which an item must occur to be
part of the core. The default `1.0` computes a strict core. For example, `0.9`
computes a soft core whose items occur in at least 90% of the sampled genomes.

### Canonical *k*-mers

Canonical counting identifies each *k*-mer with its reverse complement and is
enabled by default. Disable it only when strand orientation is meaningful for
the analysis.

### Compacted de Bruijn graph diversity

When enabled, pangrowth calculates both the *k*-mer histogram and an
infix-equivalent histogram. It combines them to estimate Hill-number diversity
for the colored compacted de Bruijn graph. This requires an additional input
scan and therefore more runtime.

### PDF plots

Plots are enabled by default. Plotting failures do not discard numerical
results; details are recorded in `pangrowth.log`. The histogram x-axis can be
shown either as an absolute genome count or normalized to the interval `(0,1]`.
The y-axis can be left unchanged, weighted by multiplicity, converted to a
percentage, or both weighted and converted to a percentage.

## Expert parameters

### Hashtable suffix size

The suffix size partitions the counting hashtable. The default is 4. Larger
values create more partitions and can impose substantial memory overhead, so
the workflow limits this setting to 12.

### Telomeres

When compacted de Bruijn graph diversity is enabled, sequence ends can be
treated as telomeres that break unitigs.

### Hill-number sampling

By default, 30 interpolation and extrapolation sample points are emitted. Use 0
to emit every point.

### CDBG interpolation

The default Bernoulli-hybrid interpolation calculates five right-tail bins
exactly. The limit can be changed, and an adaptive relative tolerance can stop
the corrections once three consecutive corrections change richness,
exponential entropy, and inverse Gini–Simpson diversity by no more than that
tolerance. Exact interpolation is available but may be considerably slower.
