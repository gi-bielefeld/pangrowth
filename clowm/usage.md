# Usage

## Input / Output

### Input archive

Select a `.zip`, `.tar.gz`, or `.tgz` archive containing at least
three genomes. Each genome must be stored in a separate FASTA file. Supported
file-name endings are `.fa`, `.fasta`, `.fna`, and `.ffn`; each may additionally
end in `.gz`.

Files may contain multiple records, such as the contigs of one assembly, and
may be arranged in subdirectories. Non-FASTA files are ignored. A FASTA file is
the sampling unit: do not split one genome over several files, because those
files would be interpreted as separate genomes.

### FASTA list input

Select a `.txt` or `.list` file such as `all_list.txt` as **Input** (`input`).
The workflow detects the input type from the filename, ignoring case; no
separate input-type parameter is needed. Unsupported extensions are rejected
before a process is submitted. Each nonempty line contains one FASTA
path, without quotes. Blank lines and lines starting with `#` are ignored.
The same FASTA extensions as archive input are supported, including `.gz`.
At least three genomes are required; duplicate paths and missing files are
reported before analysis starts.

The workflow finds each relative entry beside the list file or, for an S3
list, at the root of the same bucket. No base-folder setting is needed.
Full `s3://BUCKET/key` paths are used unchanged. If a relative entry matches
different files in both locations, the workflow reports the ambiguity; use
a full S3 path for that entry.
For CloWM, all referenced genomes must be readable using your bucket
permissions and the same S3 provider as the list and output directory.
Local absolute paths are supported only for local list files; laptop paths
are not accessible from CloWM.

For the `yeast/cerevisiae-n81/` layout, keep the existing list entries:

```text
yeast/cerevisiae-n81/GCA_039880145.1/GCA_039880145.1_ASM3988014v1_genomic.fna
yeast/cerevisiae-n81/GCA_025727145.1/GCA_025727145.1_ASM2572714v1_genomic.fna
yeast/cerevisiae-n81/GCA_025727285.1/GCA_025727285.1_ASM2572728v1_genomic.fna
```

For your bucket, select:

| Parameter | Value |
| --- | --- |
| `input` | `s3://initial-bucket-4853e34c/yeast/cerevisiae-n81/all_list.txt` |
| `outdir` | `s3://initial-bucket-4853e34c/results/cerevisiae-n81/` |

The list above works unchanged: its entries are found relative to the bucket
root. For a portable list that also works locally, you can instead write
paths relative to `all_list.txt`, such as:

```text
GCA_039880145.1/GCA_039880145.1_ASM3988014v1_genomic.fna
GCA_025727145.1/GCA_025727145.1_ASM2572714v1_genomic.fna
GCA_025727285.1/GCA_025727285.1_ASM2572728v1_genomic.fna
```

Use just a filename only when the FASTA is in the same directory as the list;
genomes in subdirectories need the subdirectory in their relative path.

Use the actual bucket name shown in the S3 path copied from CloWM's file
browser, which may differ from its display name. Do not use the HTTPS storage
endpoint as the input path. CloWM supplies the endpoint and credentials to
Nextflow; the list and genomes must already be uploaded to the bucket.
See [Nextflow S3 paths](https://docs.seqera.io/nextflow/amazons3) and
[CloWM's S3 provider selection](https://wiki.clowm.de/04_wf_parameter_form/).

Nextflow reads the list and stages each referenced FASTA as a process input.
Inside the task, paths look like `source0001/assembly.fna`. The workflow writes
a new `fasta_files.list` using these staged paths for `pangrowth hist -i`.
Numbered directories keep identical filenames from different source folders
separate. Genome data do not need to be downloaded to your laptop for a CloWM
run.

For a local run from the repository root, use the portable list with paths
relative to `all_list.txt`:

```bash
nextflow run . --input data/yeast/cerevisiae-n81/all_list.txt \
    --outdir results/cerevisiae-n81
```

This command runs the analysis and requires the configured compute resources;
it is not a lightweight input check.

### Comparing pangenomes

One selected list or archive defines one pangenome. To analyse and compare
specific collections in CloWM, switch **Input** to **Raw** and separate their
paths with commas (recommended) or whitespace:

```text
s3://initial-bucket-4853e34c/yeast/cerevisiae.txt,s3://initial-bucket-4853e34c/yeast/paradoxus.txt
```

Input paths must not themselves contain commas or whitespace. Alternatively,
use a wildcard that matches the collection files:

```text
s3://initial-bucket-4853e34c/yeast/*/all_list.txt
```

CloWM's current parameter form has no multi-file picker for a single workflow
parameter, so multiple explicit paths use Raw mode. You can still select one
exact file normally. A wildcard must match only supported collection files;
each matched file is validated independently. Direct Nextflow use additionally
accepts an array of paths if parameters are supplied through a JSON or YAML
file.

The collection filename, without `.txt`, `.list`, `.zip`, `.tar.gz`, or
`.tgz`, becomes its output folder and plot label. Thus `cerevisiae.txt` becomes
`cerevisiae/`. Unsafe characters are replaced with underscores, duplicate
names receive `_2`, `_3`, and so on, and the filename `all` becomes
`all_input`. With at least two collections, the workflow also creates `all/`
containing plots with every pangenome.

### Output directory

Numerical tables, PDF plots, fit summaries, and `pangrowth.log` are
written below one named folder per pangenome. Combined plots are written below
`all/` when at least two collections are supplied.

## Analysis parameters

### *k*-mer length

The default length is 31. Short or heterogeneous sequences may benefit from a
smaller value. Values from 5 through 31 are accepted by this workflow.

### Minimum within-genome count

By default, a *k*-mer occurring at least once in a genome is considered present.
Increasing this threshold can filter low-abundance *k*-mers introduced by
sequencing errors when read data are used.

### Core-genome quorum

The workflow always computes the strict-core curve. It also computes a separate
quorum-core curve for items occurring in a minimum fraction of genomes. The
default quorum is `0.9`, meaning at least 90% of the sampled genomes.
If the quorum is set to `1.0`, the workflow reuses the strict-core result
instead of running the same calculation twice.

### Canonical *k*-mers

Canonical counting identifies each *k*-mer with its reverse complement and is
enabled by default. Disable it only when strand orientation is meaningful for
the analysis.

### PDF plots

Plots are generated automatically. Plotting failures do not discard numerical
results; details are recorded in `pangrowth.log`. Two raw-count histogram plots
are attempted: one uses absolute genome frequencies, and the other groups
frequencies into five-percentage-point bins.

## Expert parameters

### Hashtable suffix size

The suffix size is measured in bits, not nucleotide characters. It partitions
the counting hashtable into `2^s` parts. The default is 10, creating 1,024
partitions. Larger values create more partitions and can impose substantial
memory overhead, so the workflow limits this setting to 12.

### Worker threads

Pangrowth receives the full CPU allocation of its Nextflow task through
`-t ${task.cpus}`. The bundled `highmemMedium` profile requests eight CPUs;
CloWM or another server configuration can override that allocation, and
Pangrowth will automatically use the resulting number.

### Hill-number sampling

By default, 30 interpolation and extrapolation sample points are emitted. Use 0
to emit every point.
