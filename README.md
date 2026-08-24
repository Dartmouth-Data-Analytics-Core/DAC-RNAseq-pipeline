# Dartmouth GDSC Bulk RNA-Seq Pipeline
<img src="img/cqb_logo.jpg" alt="CQB Logo" width="200" align="right"/>

![Version](https://img.shields.io/badge/version-2.0.1-blue)

The GDSC Bulk RNA-Seq pipeline provides preprocessing, alignment, and quantification of bulk RNA sequencing data with robust quality control and data visualization, implemented through [Snakemake](https://snakemake.readthedocs.io/en/stable/) for use on the [Dartmouth Discovery HPC](https://rc.dartmouth.edu/discoveryhpc/). The pipeline supports both single- and paired-end libraries, a choice of [HISAT2](https://daehwankimlab.github.io/hisat2/) or [STAR](https://github.com/alexdobin/STAR) for alignment, and is compatible with human (hg38) and mouse (GRCm39). Software dependencies are managed via Singularity containers hosted on GitHub Container Registry (GHCR) or Conda environment yaml files.

**What we return:**

An example pre-processing report can be found [here](https://github.com/Dartmouth-Data-Analytics-Core/Example_Preprocessing_Reports/blob/main/Bulk-RNA-Seq/Example_Bulk_RNA.md). In addition to a quality assessment, deliverables include a standard set of MultiQC reports (raw and alignment-level), counts (raw and normalized), PCA figures, and an interactive HTML dashboard (`RNAseq_Dashboard.html`) that collates all of the above into one page, viewable via a Dartmouth WebShare link.

## Documentation

- [Summary](#summary)
- [Installation](#installation)
- [Example Data](#example-data)
- [Configuration](#configuration)
- [Aligner Selection](#aligner-selection)
- [Optional Features](#optional-features)
- [Job Submission](#job-submission)
- [Running the Pipeline](#running-the-pipeline)
- [Prebuilt Configs](#prebuilt-configs)
- [Utilities](#utilities)
- [Contact](#contact)

## Summary

The pipeline supports Singularity containers for all software dependencies (`job.script.sh`).

To run this pipeline:
1. Populate [`sample_fastq_list_paired.csv`](sample_fastq_list_paired.csv) or [`sample_fastq_list_single.csv`](sample_fastq_list_single.csv) with your sample information
2. Select a config from [`prebuilt_configs/`](prebuilt_configs/) or edit [`config.yaml`](config.yaml) directly
3. Build reference files with `snakemake build_refs` and verify them with `snakemake check_refs`
4. Submit [`job.script.sh`](job.script.sh) to the SLURM scheduler

Currently the pipeline performs the following steps:

- Adapter and quality trimming with [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- Alignment to the reference genome using [HISAT2](https://daehwankimlab.github.io/hisat2/) or [STAR](https://github.com/alexdobin/STAR)
- Duplicate marking with [Picard MarkDuplicates](https://broadinstitute.github.io/picard/)
- RNA-seq alignment metrics with [Picard CollectRnaSeqMetrics](https://broadinstitute.github.io/picard/)
- BAM flagstat and idxstats with [Samtools](http://www.htslib.org/)
- Read quantification with [featureCounts](http://subread.sourceforge.net/), normalized to TPM, RPKM, and FPKM
- PCA, variance plots, and the underlying PC/normalised-count tables from read count matrices using a custom Python script
- Aggregated QC reporting with [MultiQC](https://multiqc.info/)
- An interactive, self-contained HTML QC dashboard (`RNAseq_Dashboard.html`)
- **(optional)** Ribosomal RNA filtering with [RiboDetector](https://github.com/hzi-bifo/RiboDetector)
- **(optional)** Comprehensive QC with [RustQC](https://github.com/seqeralabs/rustqc)
- **(optional)** Isoform quantification with [RSEM](https://deweylab.github.io/RSEM/)

## Installation

Clone the repository:

```shell
git clone https://github.com/Dartmouth-Data-Analytics-Core/DAC-RNAseq-pipeline
cd DAC-RNAseq-pipeline
```

Activate an environment containing Snakemake:

```shell
conda activate /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/snakemake
```

## Example Data

The repository includes small example datasets in [`sample_data/`](sample_data/) and [`sample_ref/`](sample_ref/) for testing the pipeline without access to full reference files or sequencing data.

>[!WARNING]
> The reference genome and annotation in `sample_ref/` are a **heavily subsetted subset of hg38 covering only chromosomes 5, 6, and 7 (first 100,000 bp each)**. They are intended solely for pipeline testing and CI/CD validation. Do not use these files for real analyses — results will be incomplete and biologically meaningless.

## Configuration

**1. Sample sheet**

Populate a sample CSV with your sample information. This is a comma-separated file with the following columns:

| Column | Description |
|--------|-------------|
| `sample_id` | Short sample identifier used to name all output files |
| `fastq_1` | Path to the R1 FASTQ file |
| `fastq_2` | Path to the R2 FASTQ file (paired-end only) |
| `group` | *(optional)* Experimental group. When present, the PCA figures and the dashboard's interactive PCA are coloured by it. |

Adding `group` changes nothing else about the run — it is read only for plotting, and omitting it leaves every output exactly as before:

```
sample_id,fastq_1,fastq_2,group
sample1,/path/s1_R1.fq.gz,/path/s1_R2.fq.gz,control
sample2,/path/s2_R1.fq.gz,/path/s2_R2.fq.gz,control
sample3,/path/s3_R1.fq.gz,/path/s3_R2.fq.gz,treated
```

Group colours come from [`sample_ref/sample_colors_hex.tsv`](sample_ref/sample_colors_hex.tsv). If your group names match that file's `group_name` values they are matched by name; otherwise its colours are used in an order chosen for visual separation, so a four-group run does not draw two near-identical purples. Samples with a blank `group` are drawn as *ungrouped*.

Then set `sample_csv` in your config to point to this file.

**2. Core settings**

| Parameter | Description | Values |
|-----------|-------------|--------|
| `layout` | Library layout | `"single"` or `"paired"` |
| `aligner_name` | Aligner to use | `"hisat"` or `"star"` |
| `featurecounts_strand` | Strandedness for featureCounts | `"0"` (unstranded), `"1"` (stranded), `"2"` (reverse) |
| `picard_strand` | Strandedness for Picard metrics | `"FIRST_READ_TRANSCRIPTION_STRAND"`, `"SECOND_READ_TRANSCRIPTION_STRAND"`, `"NONE"` |
| `rsem_strandedness` | Strandedness for RSEM | `"forward"`, `"reverse"`, `"none"` |

For a full description of every parameter and its accepted values, see [`schemas/config.schema.yaml`](schemas/config.schema.yaml).

**3. Reference files**

Reference files can be built automatically using `snakemake build_refs` (see [Running the Pipeline](#running-the-pipeline)), or provided directly:

| Parameter | Description |
|-----------|-------------|
| `reference_fa` | Path to reference genome FASTA |
| `annotation_gtf` | Path to gene annotation GTF |
| `aligner_index` | Path to HISAT2 index prefix or STAR index directory |
| `picard_refflat` | Path to Picard RefFlat annotation file |
| `picard_rrna_list` | Path to Picard rRNA interval list |

Prebuilt references for human (hg38) and mouse (GRCm39) are available to the Dartmouth community on Discovery/DartFS. See [DAC Genome References](https://github.com/Dartmouth-Data-Analytics-Core/DAC-Genome-References) for details.

**4. Pipeline parameters**

Each organism and configuration has a prebuilt config in [`prebuilt_configs/`](prebuilt_configs/). These files contain all tunable settings. For custom runs, edit [`config.yaml`](config.yaml) directly. All config fields are validated against [`schemas/config.schema.yaml`](schemas/config.schema.yaml) at startup.

## Aligner Selection

The pipeline supports two aligners, selected via `aligner_name` in the config. Only the chosen aligner's rules are loaded at runtime.

### HISAT2

```yaml
aligner_name: "hisat"
aligner_path: "hisat2"
aligner_index: "path/to/hisat2_index/genome"
```

HISAT2 is recommended for most standard RNA-seq experiments. It is fast, memory-efficient, and does not require a pre-alignment index-building step within the pipeline.

### STAR

```yaml
aligner_name: "star"
aligner_path: "STAR"
aligner_index: "path/to/star_index/"
```

STAR supports `--quantMode TranscriptomeSAM`, which is required when running RSEM for isoform quantification. If `aligner_index` does not exist, the pipeline will generate the index automatically before alignment.

## Optional Features

Optional rules are conditionally loaded at runtime based on config flags. Setting a flag to `false` (the default) means the associated rule is never registered and adds no overhead.

### rRNA Filtering — RiboDetector

>[!IMPORTANT]
> `read_length` is required when `remove_rRNA` is enabled. Set it to the read length of your library in base pairs.

```yaml
remove_rRNA: true
read_length: 150
```

When enabled, [RiboDetector](https://github.com/hzi-bifo/RiboDetector) filters ribosomal RNA reads from trimmed FASTQs before alignment. A per-sample rRNA percentage summary is generated and included in the MultiQC report.

### Comprehensive QC — RustQC

```yaml
run_rustqc: true
```

When enabled, [RustQC](https://github.com/seqeralabs/rustqc) runs a comprehensive 14-tool RNA-seq QC analysis on deduplicated BAMs in a single pass, replacing the default Samtools flagstat/idxstats step. Results are included in the MultiQC report.

### Isoform Quantification — RSEM

>[!IMPORTANT]
> RSEM requires STAR as the aligner (`aligner_name: "star"`), as it depends on transcriptome-aligned BAMs produced by STAR's `--quantMode TranscriptomeSAM` flag.

```yaml
run_rsem: true
rsem_strandedness: "reverse"
rsem_ref_path: "path/to/rsem_ref/genome"
```

When enabled, [RSEM](https://deweylab.github.io/RSEM/) quantifies transcript isoform expression in addition to gene-level featureCounts. An RSEM reference can be built automatically using `snakemake build_refs`.

## Job Submission

The pipeline can be run using either Singularity containers or Conda environments.

| Script | Method | Notes |
|--------|--------|-------|
| [`job.script.sh`](job.script.sh) | Singularity | Recommended. Pulls pre-built containers from GHCR — no environment setup required. |
| [`job.script.conda.sh`](job.script.conda.sh) | Conda | Builds environments from [`env_config/`](env_config/) YAML files on first run. Slower to start and still requires Singularity/Apptainer if running with `RustQC` enabled. |

Both scripts submit to SLURM via `sbatch`. Open the relevant script and confirm the `--configfile` path and any cluster resource settings before submitting.

## Running the Pipeline

**Build and verify reference files:**

```shell
# Build aligner index, Picard flat reference, and rRNA interval list
snakemake -s Snakefile build_refs --cores 4 --use-singularity --configfile prebuilt_configs/human_config_paired_hisat.yaml

# Append the generated reference paths to your config
cat ref/pipeline_refs/hg38.entries.yaml >> prebuilt_configs/human_config_paired_hisat.yaml

# Verify all reference paths exist and are correctly formatted
snakemake -s Snakefile check_refs --cores 4 --use-singularity --configfile prebuilt_configs/human_config_paired_hisat.yaml
```

**Submit to the SLURM scheduler:**

```shell
sbatch job.script.sh
```

**Run on a single machine:**

```shell
snakemake -s Snakefile --use-singularity --cores 40 --configfile prebuilt_configs/human_config_paired_hisat.yaml
```

**Run with a cluster profile:**

```shell
snakemake -s Snakefile --use-singularity --profile cluster_profile --configfile prebuilt_configs/human_config_paired_hisat.yaml
```

## Prebuilt Configs

Prebuilt configs are available in [`prebuilt_configs/`](prebuilt_configs/) for common combinations of organism, library layout, and aligner:

| Config | Organism | Layout | Aligner | RSEM |
|--------|----------|--------|---------|------|
| `human_config_paired_hisat.yaml` | hg38 | paired | HISAT2 | no |
| `human_config_single_hisat.yaml` | hg38 | single | HISAT2 | no |
| `human_config_paired_star.yaml` | hg38 | paired | STAR | no |
| `human_config_single_star.yaml` | hg38 | single | STAR | no |
| `human_config_paired_star_rsem.yaml` | hg38 | paired | STAR | yes |
| `mouse_config_paired_hisat.yaml` | GRCm39 | paired | HISAT2 | no |
| `mouse_config_single_hisat.yaml` | GRCm39 | single | HISAT2 | no |
| `mouse_config_paired_star.yaml` | GRCm39 | paired | STAR | no |
| `mouse_config_single_star.yaml` | GRCm39 | single | STAR | no |
| `mouse_config_paired_star_rsem.yaml` | GRCm39 | paired | STAR | yes |

When using a prebuilt config, you still need to create a sample CSV for your specific run and set the `sample_csv` field accordingly.

## Interactive Dashboard

Every run produces `RNAseq_Dashboard.html`, a single self-contained page that collates the run's QC and quantification results. It has no external dependencies — no CDN, no network, no plugins — so it can be dropped straight into a WebShare link and opened offline.

The dashboard contains:

| Section | Contents |
|---------|----------|
| Run summary | Headline tiles — mean reads obtained / uniquely aligned / assigned, then median alignment, mRNA %, strand, rRNA, duplication and PC1 variance — plus the run configuration. Read depths fall back to the featureCounts summary when the aligner logs are not in the folder, and any tile whose metric is missing is dropped rather than shown blank. |
| QC metrics | Per-sample bar panels for the key Cutadapt / aligner / Picard / featureCounts metrics, each with the cohort median marked |
| Gene-body coverage | Picard normalised 5'→3' coverage, one line per sample |
| Sample metrics table | Every parsed metric, sortable and filterable, with CSV export |
| Top genes | The 100 genes with the highest mean TPM, mitochondrial genes flagged, switchable between TPM and raw counts, with CSV export |
| PCA and distances | Interactive PCA with a selectable PC pair, scree plot, sample-to-sample distance heatmap, and the gene-variance curve |

Hovering a sample anywhere on the page highlights it in every other panel; clicking pins the highlight (`Esc` clears it). The report follows the reader's light/dark preference and remembers the toggle.

### Colouring the PCA by group

If the sample sheet has a `group` column, `pca_plotting.py` colours every PCA figure by it, adds a `Group` legend, and writes the grouping into `PCA_*_PCs.csv`. The dashboard reads that column back, so its interactive PCA uses the same grouping — hovering a group in the legend highlights all of its samples across every panel on the page. With no `group` column the plots are byte-for-byte what they were before.

Running the script by hand, point `-m` at the sample sheet:

```bash
python scripts/pca_plotting.py featurecounts/featurecounts.readcounts.tsv plots -m sample_fastq_list_paired.csv
```

### Where the PCA comes from

The dashboard does **not** recompute the PCA. [`scripts/pca_plotting.py`](scripts/pca_plotting.py) writes the numbers behind its figures, and the dashboard reads those back, so the interactive PCA and the delivered `PCA_top_*.png` are guaranteed to show the same components:

| File in `plots/` | Contents |
|------------------|----------|
| `PCA_top_PCs.csv` | Sample × PC scores, high-variance (HVG) gene set — what the dashboard plots |
| `PCA_top_variance_explained.csv` | Percent variance per PC, HVG set |
| `PCA_all_PCs.csv`, `PCA_all_variance_explained.csv` | The same for the all-gene PCA |
| `PCA_top_normalized_counts.tsv` | The log2 median-of-ratios matrix fed to the HVG PCA |
| `PCA_all_normalized_counts.tsv` | The same, before HVG selection |
| `PCA_gene_variance.csv` | Per-gene variance and rank, i.e. the HVG cut |

The **sample-to-sample distance** heatmap is Euclidean distance computed on that same normalised matrix, so samples that sit together in the PCA are the ones that read as close in the heatmap. Pass `--feature-set all` to show the all-gene PCA and its matrix instead.

If those files are missing — an older run, or a project processed before this existed — the dashboard falls back to computing the PCA in-process. The fallback mirrors `pca_plotting.py` step for step (median-of-ratios → log2 → variance-plateau HVG cut → per-gene standardisation), so it lands on the same gene set and the same components.

These outputs are also useful on their own: `PCA_top_normalized_counts.tsv` is a ready-to-use normalised matrix, and `PCA_top_PCs.csv` can be joined to sample metadata to colour the PCA by experimental group.

### What the script needs

Three files, which must stay in the same directory:

```
scripts/build_dashboard.py        # the builder
scripts/dashboard_template.html   # page shell + stylesheet
scripts/dashboard.js              # chart and table runtime
img/cqb_logo_small.png            # header logo, inlined as a data URI
```

The logo is a pre-scaled, palette-reduced copy of `img/cqb_logo.jpg` (4.6 KB rather than 227 KB) because it is embedded in every report; the full-size original would add ~296 KB to each one. Override it with `--logo <path>`, or pass `--logo ""` to omit it.

Python requirements are `numpy` and `pandas`. `scipy` is used when present (PCA variance-plateau fit, distance-heatmap ordering) and `PyYAML` when present (config parsing); both have fallbacks, so neither is required.

Nothing else is needed — no MultiQC, no R, no Quarto, no network. Point it at a directory and it uses whatever it finds:

| Present in the run directory | What you get |
|------------------------------|--------------|
| `featurecounts/*_tpm.ann.tsv` | Top-genes table, mitochondrial flags |
| `featurecounts/*.ann.tsv` | PCA + sample distances (when `plots/` is absent), counts toggle |
| `featurecounts/*.summary` | Assigned-read metrics; read-depth tiles when the aligner logs are absent |
| `plots/PCA_*` | PCA and distances taken from the pipeline's own numbers |
| `alignment/`, `markdup/`, `metrics/picard/`, `trimming/` | The QC metric panels, gene-body coverage, and the full sample table |
| `config.yaml` | Organism, aligner, layout and strand in the run-configuration card |

Sample IDs are taken from the sample sheet if one is given, then from the BAMs, the Picard metrics, the featureCounts summary, the counts header, and finally `PCA_top_PCs.csv` — so a trimmed-down delivered folder holding only `featurecounts/` and `plots/` still works without a sample sheet.

The directory's **name** is irrelevant, only the subdirectory layout matters (`featurecounts/`, `plots/`, `alignment/`, …). If you point at the folder one level above the run and exactly one child looks like a run directory, it descends into it and says so.

### Which config it reports

The run-configuration card is only as trustworthy as the config it reads, and the default is `<run-dir>/config.yaml`. That matters when a run was driven by a prebuilt config: the repo's untouched `config.yaml` is often still sitting in the directory, and reading it would describe the CI test reference as though it were the run's genome.

- **Under Snakemake** this is already handled. `rule dashboard` passes `config.yaml` *and* any `--configfile` override, merged in Snakemake's own order (later wins), so a prebuilt config is reported correctly.
- **Standalone**, pass the config the run actually used: `-c prebuilt_configs/mouse_config_single_star.yaml`. The flag is repeatable and later files win, mirroring the same merge.

As a backstop, organism is taken from the **gene ID prefixes in the count matrix** (`ENSG` → human, `ENSMUSG` → mouse, and so on) rather than from the config, since that comes from the data itself. If the config disagrees, or points at `sample_ref/`, the run-configuration card gains a `⚠ Check config` row naming the discrepancy and the script warns on stderr.

### Rebuilding by hand

The dashboard is built by [`scripts/build_dashboard.py`](scripts/build_dashboard.py), which reads the pipeline's own outputs rather than `multiqc_report.html`. That means it can be re-run at any time on a finished run directory — useful after re-sequencing a sample, or for a project processed before the dashboard existed:

```shell
# from inside a finished pipeline run directory
python scripts/build_dashboard.py

# or point it anywhere
python scripts/build_dashboard.py \
    --run-dir /path/to/run \
    --output QC_Dashboard.html
```

Useful options:

| Option | Description |
|--------|-------------|
| `-d`, `--run-dir` | Pipeline run directory (default: current directory) |
| `-o`, `--output` | Output HTML path (default: `RNAseq_Dashboard.html`) |
| `-c`, `--config` | Pipeline config to read run context from (default: `<run-dir>/config.yaml`) |
| `--fastqc-dir` | Extra directory to search for a raw-read MultiQC, e.g. `fastqc_results/` (repeatable) |
| `--top-n` | Number of genes in the expression table (default: 100) |
| `--pcs` | Principal components to compute (default: 10) |
| `--hvg` | Fixed number of variable genes for PCA (default: variance-plateau auto-detect, matching `pca_plotting.py`) |
| `--plots-dir` | Where `pca_plotting.py` wrote its outputs (default: `<run-dir>/plots`) |
| `--feature-set` | `top` (default, the high-variance gene set) or `all` |
| `--recompute-pca` | Ignore the precomputed PCA and recompute from the counts matrix |
| `--dist-method` | `euclidean` (default), or `spearman` / `pearson` as `1 - r` |
| `--title`, `--subtitle` | Header text |
| `--logo` | Logo inlined top-right (default `img/cqb_logo_small.png`; `""` for none) |

Every panel degrades on its own: if a run has no Picard metrics, no TPM matrix, or only one sample, the affected panels say so and the rest of the report still builds. Raw-read FastQC metrics (% GC, per-read duplication) are picked up automatically when a `fastqc_results/` MultiQC from [`Utilities/run_fastqc.sh`](Utilities/run_fastqc.sh) sits alongside the run.

The script needs only `numpy` and `pandas`; `scipy` is used when present for the PCA variance-plateau fit and for ordering the distance heatmap, and both fall back cleanly when it is not. It runs in the `rna_pcaplot` environment/container.

## Utilities

The [`Utilities/`](Utilities/) folder contains helper scripts for common pre-pipeline tasks. See [`Utilities/README.md`](Utilities/README.md) for full usage instructions. Scripts include:

- **Sample sheet generation** — automatically links GSR sequencing metadata to external sample IDs and generates a pipeline-ready CSV (`make_sample_sheet.sh` + `linkMeta.R`)
- **Raw read QC** — runs FastQC on raw FASTQ files and aggregates results into a MultiQC report with sample-name remapping (`run_fastqc.sh`)
- **Contamination screening** — screens raw reads against common contaminants using FastQ Screen (`run_fastq_screen.sh`)
- **TPM annotation** — generates a color-annotated Excel spreadsheet of TPM expression values for post-pipeline QC (`Annotate_TPMs.R`)

## Contact

**Contact and questions:** Please address questions to *DataAnalyticsCore@groups.dartmouth.edu* or submit an issue in the GitHub repository.

**This pipeline was created with funds from the COBRE grant 1P20GM130454. If you use the pipeline in your own work, please acknowledge the pipeline by citing the grant number in your manuscript.**
