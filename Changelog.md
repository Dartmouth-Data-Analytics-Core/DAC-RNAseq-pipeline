## v1.3

- PCA is calculated using a dynamic number of HVGs (specified by asymptote detection) as well as using ALL genes
- Switched plotting over to plotnine (updated pcaplot_env.yaml)
- Added fast clustering (update pcaplot_env.yaml)
- Added log file output that informs how many HVGs were used for PCA calculation
- Updates to multiqc_env.yaml to help github tests pass

## v2.0

Modular refactor of the GDSC Bulk RNA-Seq Pipeline: alignment, optional-feature,
and reference-building rules now live in `additional_rules/` and are loaded via
conditional `include:` based on config. Adds Singularity support and config validation.

### Added
- **Singularity support** — containers hosted on ghcr.io, driven by `job.script.sh`.
  `job.script.conda.sh` retained for Conda runs.
- **Config validation** — `schemas/config.schema.yaml` (JSON Schema draft-07) covering
  all parameters with type and enum constraints, plus a conditional requiring
  `read_length` when `remove_rRNA: true`. Validated at startup via `snakemake.utils.validate()`.

### Changed
- **Modular rules** — optional and aligner-specific rules extracted from the Snakefile:

  | File | Loaded when |
  |------|-------------|
  | `additional_rules/alignment/star.smk` | `aligner_name: "star"` |
  | `additional_rules/alignment/hisat.smk` | `aligner_name: "hisat2"` |
  | `additional_rules/ribodetector/ribodetector.smk` | `remove_rRNA: true` |
  | `additional_rules/RustQC/rustqc.smk` | `run_rustqc: true` |
  | `additional_rules/rsem/rsem.smk` | `run_rsem: true` |
  | `additional_rules/check_refs/check_refs.smk` | always (target-only) |
  | `additional_rules/build_refs/build_refs.smk` | always (target-only) |

  `check_refs` and `build_refs` declare no `output:`, so they never enter the default
  DAG and run only when targeted by name (e.g. `snakemake build_refs`) — CI/CD behavior unchanged.
- **`rule all`** — rebuilt as an `all_inputs = []` list with labeled sections (trimming,
  ribodetector, alignment, deduplication, QC, RSEM, featureCounts, PCA); optional blocks
  appended by config flag.
- **`run_rsem`** — migrated from `"yes"`/`"no"` strings to booleans across all 18 config
  files (`config.yaml`, `prebuilt_configs/*.yaml`, `tests/test_config_*.yaml`).
- **`aligner_name`** — `"hisat"` renamed to `"hisat2"`; index paths standardized on
  `hisat2_index/`.
- **README** — rewritten to mirror the GDSC miRNA-Seq pipeline: Optional Features
  (RiboDetector, RustQC, RSEM), Aligner Selection, Example Data (with a warning that
  the sample data is a chr5/6/7 subset), Job Submission (Singularity vs Conda),
  Prebuilt Configs, Utilities, and a schema reference in Configuration.

### Fixed
- Wrapped `sns.clustermap` in `pca_plotting.py` in `try/except RecursionError` — scipy's
  recursive dendrogram algorithm blows the default recursion limit past ~10k genes.
  PCA plots are written before the clustermap step, so declared outputs now survive the failure.

## v2.0.1

### Fixed

- Added singularity support for conda job script since RustQC requires contanerization

- Fixed `hisat` path and aligner mismatch, use `hisat2` instead. 

- Added `tmpdir` directive to `ribodetector.smk` to avoid internal clean-up error