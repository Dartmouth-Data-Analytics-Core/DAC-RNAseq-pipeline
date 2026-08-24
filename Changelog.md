## v2.1 (unreleased)

### Added
- **Interactive HTML dashboard** — every run now emits `RNAseq_Dashboard.html`, a single
  self-contained page (no CDN, no network, no plugins) collating the run's QC and
  quantification results: headline QC tiles, per-sample metric bar panels, Picard
  gene-body coverage, a sortable/filterable metrics table, the top 100 genes by mean TPM
  with mitochondrial genes flagged, and interactive PCA, scree, sample-distance and
  gene-variance panels. Hovering a sample in any panel highlights it in all the others.
  Light/dark aware, with CSV export from both tables, and the CQB mark in the top-right
  corner (a 4.6 KB pre-scaled copy of `img/cqb_logo.jpg`, since it is embedded in every
  report; `--logo` overrides it and `--logo ""` omits it).

  | File | Purpose |
  |------|---------|
  | `scripts/build_dashboard.py` | Builder — parses the pipeline outputs, computes PCA/correlation, writes the HTML |
  | `scripts/dashboard_template.html` | Page shell and stylesheet |
  | `scripts/dashboard.js` | Chart and table runtime (inlined at build time) |
  | `img/cqb_logo_small.png` | CQB mark for the header, inlined as a data URI |

  Metrics are parsed from the upstream tool outputs (STAR/HISAT2 logs, Picard
  `CollectRnaSeqMetrics` and `MarkDuplicates`, the featureCounts `.summary`, Cutadapt
  reports, RiboDetector) rather than from `multiqc_report.html`, so gene-body coverage
  and strand correctness come from the authoritative files instead of MultiQC's
  compressed plotdata. The script therefore also runs standalone on any finished run
  directory — see the README's *Interactive Dashboard* section.

- **`rule dashboard`** — added to the Snakefile and to `rule all`; runs in the
  `rna_pcaplot` environment/container after featureCounts and Picard.
- **`dashboard_script` config key** — optional (not in the schema's `required` list), so
  configs written before v2.1 remain valid; the Snakefile falls back to
  `scripts/build_dashboard.py` when it is absent.

- **Tabular PCA outputs** — `scripts/pca_plotting.py` now writes the numbers behind its
  figures, not just the PNGs:

  | File in `plots/` | Contents |
  |------------------|----------|
  | `PCA_top_PCs.csv` | Sample x PC scores, high-variance gene set |
  | `PCA_top_variance_explained.csv` | Percent variance per PC, high-variance gene set |
  | `PCA_all_PCs.csv`, `PCA_all_variance_explained.csv` | The same for the all-gene PCA |
  | `PCA_top_normalized_counts.tsv` | The log2 median-of-ratios matrix fed to the HVG PCA |
  | `PCA_all_normalized_counts.tsv` | The same, before HVG selection |
  | `PCA_gene_variance.csv` | Per-gene variance and rank, i.e. the HVG cut |

  All are declared as outputs of `rule pca_plots` and as inputs to `rule dashboard`.

- **Optional `group` column in the sample sheet** — when present, `pca_plotting.py`
  colours every PCA figure by it (with a `Group` legend) and writes the grouping into
  `PCA_*_PCs.csv`; the dashboard reads it back and colours its interactive PCA the same
  way, with legend hover highlighting a whole group across every panel. Colours come
  from the previously unused `sample_ref/sample_colors_hex.tsv`, matched by name when
  the group names are the ones it lists and otherwise ordered for visual separation
  (its 3rd and 4th entries are near-identical purples, so file order would draw two
  groups the same). No new dependencies: plotnine and pandas were already in
  `rna_pcaplot`. Sheets without a `group` column produce byte-identical output.
  `sample_fastq_list_paired.csv` now carries the column so CI exercises both paths.

### Changed
- **The dashboard reads the pipeline's PCA instead of recomputing it** — it plots
  `PCA_top_PCs.csv` (the high-variance feature set), so the interactive PCA and the
  delivered `PCA_top_*.png` cannot disagree. `--feature-set all` switches to the
  all-gene PCA; `--recompute-pca` forces the in-process path. That fallback still exists
  for runs predating these outputs, and now mirrors `pca_plotting.py` step for step
  (all-zero genes are no longer pre-filtered) so it selects the same gene set.
- **Sample-sample correlation replaced by a distance matrix** — computed as Euclidean
  distance on `PCA_top_normalized_counts.tsv`, the same matrix the PCA used, so the two
  panels are directly comparable. `--dist-method` also accepts `spearman` / `pearson`,
  applied as `1 - r`. Replaces the previous `--corr-method` flag.
- **`pca_plotting.py` plotting imports are now guarded** — plotnine/seaborn/matplotlib
  are imported behind a `try`, so a missing plotting backend no longer blocks the
  tabular outputs. `rule pca_plots` still declares the PNGs as outputs, so a genuinely
  broken plotting environment fails the rule as before.

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

- **Optional `group` column in the sample sheet** — when present, `pca_plotting.py`
  colours every PCA figure by it (with a `Group` legend) and writes the grouping into
  `PCA_*_PCs.csv`; the dashboard reads it back and colours its interactive PCA the same
  way, with legend hover highlighting a whole group across every panel. Colours come
  from the previously unused `sample_ref/sample_colors_hex.tsv`, matched by name when
  the group names are the ones it lists and otherwise ordered for visual separation
  (its 3rd and 4th entries are near-identical purples, so file order would draw two
  groups the same). No new dependencies: plotnine and pandas were already in
  `rna_pcaplot`. Sheets without a `group` column produce byte-identical output.
  `sample_fastq_list_paired.csv` now carries the column so CI exercises both paths.

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