# pangrowth

`pangrowth` estimates how a pangenome changes as genomes are sampled. From one
or more genome archives or lists of nucleotide FASTA paths, this workflow:

* counts the frequency of each *k*-mer across genomes;
* calculates the exact expected pangenome growth curve;
* calculates separate expected strict-core and quorum-core curves; and
* reports Hill-number diversity (richness, exponential entropy, and inverse
  Simpson diversity).

PDF plots and plain-text numerical results are produced automatically for
downstream use.

The workflow records computation and plotting as separate Nextflow processes:
`PANGROWTH` creates the numerical results, and `PLOT` consumes those results to
create the individual and combined visualisations.

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

To compare several pangenomes, enter comma-separated file paths using the
Input field's Raw mode, or enter a wildcard that matches several list files or
archives. Each collection is analysed separately. Its results are placed in a
folder named after the collection file, and combined plots are placed in
`all/`.

See [Usage](usage.md) for parameter details and [Output](output.md) for a
description of every generated file.

## Local execution

With Nextflow and Docker installed, run the same workflow locally with:

```bash
nextflow run . --input genomes.tar.gz --outdir results
```

Quote wildcard input so that Nextflow, rather than the shell, resolves it:

```bash
nextflow run . --input 'data/pangenomes/*.txt' --outdir results
```

## Contact

For questions, feedback, or problems, contact
[pangenomics-service@cebitec.uni-bielefeld.de](mailto:pangenomics-service@cebitec.uni-bielefeld.de)
or open an issue in the
[pangrowth repository](https://github.com/gi-bielefeld/pangrowth/issues).

Pangrowth is provided as a service of the
[German Network for Bioinformatics Infrastructure (de.NBI)](https://www.denbi.de/).
