# Changelog

All notable changes to the pangrowth CloWM workflow are documented in this
file. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and workflow versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

* FASTA list input with an optional base path for relative entries, including
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
