## v2.1 (unreleased)

### Added
- **Interactive HTML dashboard** — every run emits `RNAseq_Dashboard.html`, a single
  self-contained page (no CDN, no network, no plugins) collating the run's QC and
  quantification results: headline tiles, per-sample metric bar panels, Picard gene-body
  coverage, a sortable/filterable metrics table, the top 100 genes by mean TPM with
  mitochondrial genes flagged, and interactive PCA, scree, sample-distance and
  gene-variance panels. Hovering a sample in any panel highlights it in all the others.
  Light/dark aware, CSV export from the sample-metrics table, CQB mark top-right.

  | File | Purpose |
  |------|---------|
  | `scripts/build_dashboard.py` | Builder — parses pipeline outputs, computes what is missing, writes the HTML |
  | `scripts/dashboard_template.html` | Page shell and stylesheet |
  | `scripts/dashboard.js` | Chart and table runtime (inlined at build time) |
  | `scripts/dashboard_report.py` | Verdict-section renderer, standard library only |
  | `scripts/add_qc_report.py` | Adds the verdict table to an already-built dashboard |
  | `img/cqb_logo_small.png` | Header logo, inlined as a data URI |

  Metrics are parsed from the upstream tool outputs (STAR/HISAT2 logs, Picard
  `CollectRnaSeqMetrics` and `MarkDuplicates`, the featureCounts `.summary`, Cutadapt
  reports, RiboDetector) rather than from `multiqc_report.html`, so gene-body coverage
  and strand correctness come from the authoritative files. The script also runs
  standalone on any finished run directory.

- **`rule dashboard`** — in the Snakefile and in `rule all`; runs in the `rna_pcaplot`
  environment/container after featureCounts and Picard.

- **QC verdict section** — a Markdown `qc_assessment.md` (kit context, per-sample
  PASS / PASS WITH CAVEATS / FAIL table, cohort paragraph) renders as a section at the
  bottom of the dashboard, with verdicts as coloured badges. Picked up automatically on
  rebuild, or inserted into an existing dashboard with `scripts/add_qc_report.py`
  (idempotent, marker-bounded, `--remove` to strip). The Markdown subset escapes all
  authored text, so a report can never break the page. The client results email stays a
  separate deliverable and is not rendered into the dashboard.

- **Tabular PCA outputs** — `scripts/pca_plotting.py` writes the numbers behind its
  figures: `PCA_top_PCs.csv`, `PCA_top_variance_explained.csv`, the `PCA_all_*`
  equivalents, `PCA_top_normalized_counts.tsv`, `PCA_all_normalized_counts.tsv`, and
  `PCA_gene_variance.csv`. Declared as outputs of `rule pca_plots` and inputs to
  `rule dashboard`.

- **Optional `group` column in the sample sheet** — when present, `pca_plotting.py`
  colours every PCA figure by it (with a `Group` legend) and writes the grouping into
  `PCA_*_PCs.csv`; the dashboard reads it back and colours its interactive PCA the same
  way, with legend hover highlighting a whole group across every panel. Colours come from
  the previously unused `sample_ref/sample_colors_hex.tsv`, matched by name when the group
  names are the ones it lists and otherwise ordered for visual separation. Optional for
  both layouts, matched by name rather than position; blank entries render as *ungrouped*,
  and sheets without the column produce byte-identical output.

- **`dashboard_script` config key** — optional (not in the schema's `required` list), so
  configs written before v2.1 remain valid; the Snakefile falls back to
  `scripts/build_dashboard.py` when it is absent. Added to the shipping configs and to
  `tests/*.yaml`, which previously inherited script paths from the repo's `config.yaml`
  rather than stating them.

### Changed
- **The dashboard reads the pipeline's PCA instead of recomputing it** — it plots
  `PCA_top_PCs.csv`, so the interactive PCA and the delivered PNGs cannot disagree.
  `--feature-set all` switches to the all-gene PCA; `--recompute-pca` forces the
  in-process path, which now mirrors `pca_plotting.py` step for step.
- **Sample-sample correlation replaced by a distance matrix** — Euclidean distance on
  `PCA_top_normalized_counts.tsv`, the same matrix the PCA used. `--dist-method` also
  accepts `spearman` / `pearson`, applied as `1 - r`.
- **`pca_plotting.py` plotting imports are guarded** — a missing plotting backend no
  longer blocks the tabular outputs. `rule pca_plots` still declares the PNGs as outputs.
- **Trimming split into a count and a rate** — `Reads post-trim` (absolute, in M) plus
  `% Reads retained`, replacing the single `% reads post-trim` column.
- **No CSV export on the top-genes table** — it is a top-100 subset, and a download
  button invites treating it as the gene list. A note names the full matrices shipped
  with the results. The sample-metrics table keeps its export, since that one is complete.
- **The top-genes table no longer sorts by gene name** — the table's premise is rank by
  mean TPM. Chromosome, mean, and per-sample columns still sort.
- **CI runs on `dev-dashboard`** — added to the workflow's push and pull_request branches.

### Fixed
- **Top-genes header covered the first row** — the sticky column header carried a 31px
  offset meant for the sample-metrics table's group-header row. Now scoped to a header
  row that actually follows a group row.

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