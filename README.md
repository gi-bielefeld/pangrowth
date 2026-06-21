# Pangrowth

<p align="center">
    <img src="png/logo.png" width="512" alt="logo_pangrowth">
</p>

`pangrowth` is an efficient tool designed for genomic researchers to predict the [openness of a pangenome](https://doi.org/10.1101/2022.11.15.516472 ), 
estimate the **core** genome size and the pangenome diversiy using [Hill numbers](https://en.wikipedia.org/wiki/Diversity_index ). 

This tool is capable of analyzing fasta sequences using **_k_-mers**, as well
as any other genomic elements such as genes, CDS, ORFs, as long as it is provided as either a **frequency histogram** or a **pan-matrix** (with columns representing genomes and rows representing items; see `panmatrix_ecoli_n50.txt` for an example).

### Key features

- **_k_-mer counting**: utilizes a modified version of [yak](https://github.com/lh3/yak) to count _k_-mers
- **growth/core calculation**: computes the _exact_ expected genomic growth/core size quadratically in the number of genomes
- **hill numbers**: compute Hill numbers from frequency list
- **colored compacted de Bruijn graph (cdbg)**: estimates pangenome diversity of the ccdbg using Hill numbers, by combining _k_-mer and infix equivalents histograms

## Table of Contents
<!--ts-->
   * [Pangrowth](#pangrowth)
      * [Table of Contents](#table-of-contents)
      * [Install](#install)
      * [Usage](#usage)
         * [Histogram from fasta files](#histogram-from-fasta-files)
         * [Pangenome growth (from histogram or pan-matrix)](#pangenome-growth-from-histogram-or-pan-matrix)
         * [Pangenome core (from histogram or pan-matrix)](#pangenome-core-from-histogram-or-pan-matrix)
         * [Hill numbers from k-mer histogram](#hill-numbers-from-k-mer-histogram)
         * [Colored compacted de Bruijn graph](#colored-compacted-de-bruijn-graph)
      * [Contact](#contact)
<!--te-->

## Install

Default build:

```bash
git clone https://github.com/gi-bielefeld/pangrowth
cd pangrowth
cmake -B build
cmake --build build -j 8
```

This builds the main `pangrowth` executable. The default build does not link
the ggcat C++ API. The executable is `build/pangrowth`; either call it by that
path, add it to `PATH`, or run examples from inside the build directory.

Optional ggcat-enabled build:

```bash
./scripts/build_ggcat_cpp_api.sh

cmake -B build-ggcat -DPANGROWTH_WITH_GGCAT=ON
cmake --build build-ggcat -j 8
```

The helper script downloads ggcat and builds the static C++ API libraries with
`cargo`. The ggcat-enabled executable is `build-ggcat/pangrowth`.

To plot the results we need the following python libraries: numpy, pandas, matplotlib, scipy and searbon. You can install them with:
```
pip install -r requirements.txt
```

## Usage

### Histogram from fasta files

```bash
./pangrowth hist -k 17 -t 12 data/fa/*.fna.gz > hist.txt
```
- `pangrowth` also accepts a file containing a list of fasta files (each one on a single line) 
    passed with the paremeter `-i fasta_list.txt`

To visualize the histogram:
```bash
python scripts/plot_hist.py hist.txt hist.pdf
```

<p align="center">
    <img src="png/hist.png" width="750" alt="k-mer frequency histogram of 12 ecoli">
</p>

If you have multiple histograms that you want to compare with different number
of genomes you can use:

```bash
python scripts/plot_hist.py --norm_x --norm_y=both hist.txt data/hist_ecoli_n50.txt data/hist_ecoli_n200.txt hist_multiple.pdf
```

- The flag `--norm_x` normalize the x-axis to be between (0,1].
- The flag `--norm_y` allows two types of normalization: 
    - `multiplicity` which adjusts each histogram value h[i] multiplying it by its
      index i (i.e., h[i] * i, this means that values appearing once remain the same,
      values appearing twice are doubled, and so on)
    - `percentage` which divides the values of h[i] by the total sum of h (its total
      sum equals 1) 
    The `--norm_y=both` applies both in series.

<p align="center">
    <img src="png/hist_multiple.png" width="800" alt="k-mer frequency histogram of multiple ecoli">
</p>

### Pangenome growth from histogram (or pan-matrix)

```bash
./pangrowth growth -h data/hist_ecoli_n50.txt > growth.txt
#./pangrowth growth -p data/panmatrix_ecoli_n50.txt > growth.txt
```

To fit the openness and visualize the growth:

```bash
python scripts/plot_growth.py growth.txt growth.pdf
```

<p align="center">
    <img src="png/growth.png" width="850" alt="k-mer growth of ecoli">
</p>

We can again pass multiple growth files to `scripts/plot_growth.py` to
compare with other species.

```bash
python scripts/plot_growth.py growth.txt data/growth_ecoli_n200.txt growth_multiple.pdf
```

<p align="center">
    <img src="png/growth_multiple.png" width="850" alt="k-mer growth of multiple ecoli">
</p>

### Pangenome core from histogram (or pan-matrix)

```bash
./pangrowth core -h data/hist_ecoli_n50.txt > core.txt
#./pangrowth core -p data/panmatrix_ecoli_n50.txt > core.txt
./pangrowth core -h data/hist_ecoli_n50.txt -q 0.9 > core_q90.txt
```

- The `-q` takes a quorum to considered the item in the core (default 1.0).

To fit the core genome and report the percentage of core item over the expected
genome size:


```bash
python scripts/plot_core.py core_q90.txt data/core_q90_ecoli_n200.txt core.pdf
```  

The expected genome size is calculated as the total sum of the histogram divided by
the number of genomes.

<p align="center">
    <img src="png/core.png" width="700" alt="k-mer core size of multiple ecoli">
</p>

### Hill numbers from k-mer histogram

Hill numbers measure pangenome diversity (species richness, exponential entropy, inverse Simpson index) from a _k_-mer frequency histogram:

```bash
./pangrowth hill -p 30 data/hist_ecoli_n50.txt
```

- `-p INT` sets the number of sample points (default: 30); use `-p 0` to output all points
- `-f FILE` reads sample points from a file (one integer per line), overriding `-p`

The output is a tab-separated table with columns: `fit`, `m`, `richness`, `exp_entropy`, `inv_gini_simp`, where `fit` is `int` (interpolation), `obs` (observed), or `ext` (extrapolation) and `m` is the number of genomes.

### Colored compacted de Bruijn graph

The colored compacted de Bruijn graph (cdbg) compacts non-branching path of _k_-mers into **unitigs**. Its diversity can be estimated by combining a _k_-mer histogram with an **infix equivalents histogram**.

**Step 1**: generate the _k_-mer and infix-equivalent histograms.

From fasta files:

```bash
./pangrowth hist --cdbg -k 17 -t 12 -T -o ecoli data/fa/*.fna.gz
```

This writes `ecoli_hist.txt` and `ecoli_hist_infix.txt`. Without `-o`, the
default names are `out_hist.txt` and `out_hist_infix.txt`.

Options:
- `-k INT` _k_-mer size (default: 17)
- `-t INT` number of worker threads (default: 4)
- `-i PATH` file containing a list of fasta files (one per line)
- `-b` turn off canonical _k_-mer transformation
- `-T` account for telomeres breaking unitigs
- `-c INT` minimum _k_-mer count to consider a (_k_+1)-mer (default: 1)

From a colored ggcat _k_-mer graph:

```bash
ggcat build -c -k 17 -d inputs.tsv -o graph_k17.fa -t ggcat_tmp -j 8 --min-multiplicity 1
./pangrowth hist --cdbg --ggcat -k 17 -t 8 -o ecoli graph_k17.fa
```

The `inputs.tsv` file must contain one color name and one fasta path per line,
separated by a tab. The ggcat build used here currently requires `-j` to be a
power of two (eg, `1`, `2`, `4`, `8`...).

The ggcat mode computes both histograms in one graph scan. For an edge candidate
$axb$, its color set is $C(ax) \cap C(xb)$ from the adjacent _k_-mer color
sets.

To compute only the _k_-mer histogram, retain stdout:

```bash
./pangrowth hist --ggcat -k 17 -t 8 graph_k17.fa > hist.txt
```

**Step 2**: compute Hill numbers for the cdbg using both histograms:

```bash
./pangrowth hill ecoli_hist.txt ecoli_hist_infix.txt
```

By default, interpolation uses the Bernoulli-hybrid approximation with exact
values for the last 5 bins (`-B 5`). Use `-B INT` to change the exact right
tail, or `-E` to force the original exact interpolation.

Adaptive interpolation stops before the `-B` limit after three consecutive
exact tail corrections each change richness, exponential entropy, and inverse
Gini--Simpson by at most the `-A` relative tolerance. For example:

```bash
./pangrowth hill -B 20 -A 0.001 ecoli_hist.txt ecoli_hist_infix.txt
```

If `-B` is omitted with `-A`, its limit is 20. The adaptive criterion controls
the effect of each additional exact bin; it is not a bound on the final error
relative to fully exact interpolation.

The output format is identical to `hill`: a tab-separated table with columns `fit`, `m`, `richness`, `exp_entropy`, `inv_gini_simp`.

The previous `hist_infix`, `hist_infix_ggcat`, and `hill_cdbg` commands remain
available as compatibility aliases.

## Publication

Parmigiani, L., Wittler, R., Stoye, J.,: [Revisiting pangenome openness with _k_-mers](https://peercommunityjournal.org/articles/10.24072/pcjournal.415/  ). PCI Comp & Biol.  (2024).

## Contact

For any question, feedback or problem, please feel free to file an [issue](https://github.com/gi-bielefeld/pangrowth/issues/new) on Github or contact me [here](mailto:pangenomics-service@cebitec.uni-bielefeld.de) and I will get back to you as soon as possible.

Pangrowth is provided as a service of the [German Network for Bioinformatics Infrastructure (de.NBI)](https://www.denbi.de/). We would appriciate if you would participate in the evaluation of Pangrwoth by completing this [very short survey](https://www.surveymonkey.de/r/denbi-service?sc=bigi&tool=pangrowth).
