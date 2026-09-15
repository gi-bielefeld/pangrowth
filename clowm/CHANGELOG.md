# Changelog

All notable changes to the pangrowth CloWM workflow are documented in this
file. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and workflow versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Fixed

* Build the development container from the actual `clowm` branch and use its
  `clowm` image tag, keeping workflow code and bundled plotting scripts in
  sync.
* Defer per-pangenome output-folder naming until the process input name is
  available, avoiding an early `No such variable: datasetId` error.

### Changed

* Always create absolute and five-percentage-point histogram plots with raw
  *k*-mer counts. Histogram normalization controls are no longer exposed in
  CloWM, and absolute-axis tick labels are thinned automatically.
* Detect list (`.txt`, `.list`) and archive (`.zip`, `.tar.gz`, `.tgz`) input
  automatically from the filename, ignoring case.
* Resolve list entries automatically relative to the list or its S3 bucket
  root, with an error for ambiguous matches. Removed the `input_type` and
  `list_base` parameters; users only select their input file.

### Added

* Multi-pangenome execution from comma- or whitespace-separated input paths,
  an input wildcard, or parameter arrays for direct Nextflow use, with one
  result folder per collection and combined comparison plots in `all/`.
* Hill-number PDF with interpolation, observed, and extrapolation styling.
* FASTA list input with automatic resolution of relative entries, including
  bucket-root paths in CloWM. Referenced genomes are staged by Nextflow into
  separate directories to avoid filename collisions.
* Validation of list entries for missing files, duplicate paths, FASTA
  extensions, minimum genome count, and accidental local paths in S3 lists.

## [0.1.0] - 2026-07-30

### Added

* Initial CloWM Nextflow workflow for archived nucleotide FASTA collections.
* Frequency histogram, pangenome growth, core-genome, and Hill-number outputs.
* Optional colored compacted de Bruijn graph diversity calculation.
* Optional PDF plots and fit summaries.
* CloWM parameter schema, metadata, documentation, and container publishing
  workflow.

### Fixed

* Include `ps` in the runtime image so Nextflow can collect task metrics on
  CloWM.
