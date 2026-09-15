# pangrowth

`pangrowth` estimates how a pangenome changes as genomes are sampled. From a
genome archive or a list of nucleotide FASTA paths, this workflow:

* counts the frequency of each *k*-mer across genomes;
* calculates the exact expected pangenome growth curve;
* calculates the expected strict or quorum-based core-genome curve; and
* reports Hill-number diversity (richness, exponential entropy, and inverse
  Simpson diversity).

It can optionally estimate diversity on the colored compacted de Bruijn graph
by combining *k*-mer and infix-equivalent frequency histograms. PDF plots and
plain-text numerical results are produced for downstream use.

## Input at a glance

Upload a `.zip`, `.tar.gz`, or `.tgz` archive containing at least three FASTA
files. Each FASTA file is treated as one genome. Files may be gzip-compressed
and may be placed in subdirectories inside the archive.

Alternatively, select a `.txt` or `.list` file such as `all_list.txt`; the
workflow detects list input automatically from its filename.
Paths in the list are resolved automatically relative to the list's directory
or its S3 bucket root. Full `s3://BUCKET/key` paths are also accepted.
Nextflow stages the listed files automatically using the selected S3 provider
and your bucket permissions.

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
