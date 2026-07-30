# pangrowth

`pangrowth` estimates how a pangenome changes as genomes are sampled. From an
archive containing one nucleotide FASTA file per genome, this workflow:

* counts the frequency of each *k*-mer across genomes;
* calculates the exact expected pangenome growth curve;
* calculates the expected strict or quorum-based core-genome curve; and
* reports Hill-number diversity (richness, exponential entropy, and inverse
  Gini–Simpson diversity).

It can optionally estimate diversity on the colored compacted de Bruijn graph
by combining *k*-mer and infix-equivalent frequency histograms. PDF plots and
plain-text numerical results are produced for downstream use.

## Input at a glance

Upload a `.zip`, `.tar.gz`, or `.tgz` archive containing at least three FASTA
files. Each FASTA file is treated as one genome. Files may be gzip-compressed
and may be placed in subdirectories inside the archive.

See [Usage](usage.md) for parameter details and [Output](output.md) for a
description of every generated file.

## Local execution

With Nextflow and Docker installed, run the same workflow locally with:

```bash
nextflow run . --input genomes.tar.gz --outdir results
```

## Contact

For questions, feedback, or problems, contact
[pangenomics-service@cebitec.uni-bielefeld.de](mailto:pangenomics-service@cebitec.uni-bielefeld.de)
or open an issue in the
[pangrowth repository](https://github.com/gi-bielefeld/pangrowth/issues).

Pangrowth is provided as a service of the
[German Network for Bioinformatics Infrastructure (de.NBI)](https://www.denbi.de/).
