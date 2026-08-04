# Dartmouth GDSC Bulk RNA-Seq Pipeline
<img src="img/cqb_logo.jpg" alt="CQB Logo" width="200" align="right"/>

![Version](https://img.shields.io/badge/version-2.0.1-blue)

The GDSC Bulk RNA-Seq pipeline provides preprocessing, alignment, and quantification of bulk RNA sequencing data with robust quality control and data visualization, implemented through [Snakemake](https://snakemake.readthedocs.io/en/stable/) for use on the [Dartmouth Discovery HPC](https://rc.dartmouth.edu/discoveryhpc/). The pipeline supports both single- and paired-end libraries, a choice of [HISAT2](https://daehwankimlab.github.io/hisat2/) or [STAR](https://github.com/alexdobin/STAR) for alignment, and is compatible with human (hg38) and mouse (GRCm39). Software dependencies are managed via Singularity containers hosted on GitHub Container Registry (GHCR) or Conda environment yaml files.

**What we return:**

An example pre-processing report can be found [here](https://github.com/Dartmouth-Data-Analytics-Core/Example_Preprocessing_Reports/blob/main/Bulk-RNA-Seq/Example_Bulk_RNA.md). In addition to a quality assessment, deliverables include a standard set of MultiQC reports (raw and alignment-level), counts (raw and normalized) and PCA figures, viewable via a Dartmouth WebShare link.

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
- PCA and variance plots from read count matrices using a custom Python script
- Aggregated QC reporting with [MultiQC](https://multiqc.info/)
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

## Utilities

The [`Utilities/`](Utilities/) folder contains helper scripts for common pre-pipeline tasks. See [`Utilities/README.md`](Utilities/README.md) for full usage instructions. Scripts include:

- **Sample sheet generation** — automatically links GSR sequencing metadata to external sample IDs and generates a pipeline-ready CSV (`make_sample_sheet.sh` + `linkMeta.R`)
- **Raw read QC** — runs FastQC on raw FASTQ files and aggregates results into a MultiQC report with sample-name remapping (`run_fastqc.sh`)
- **Contamination screening** — screens raw reads against common contaminants using FastQ Screen (`run_fastq_screen.sh`)
- **TPM annotation** — generates a color-annotated Excel spreadsheet of TPM expression values for post-pipeline QC (`Annotate_TPMs.R`)

## Contact

**Contact and questions:** Please address questions to *DataAnalyticsCore@groups.dartmouth.edu* or submit an issue in the GitHub repository.

**This pipeline was created with funds from the COBRE grant 1P20GM130454. If you use the pipeline in your own work, please acknowledge the pipeline by citing the grant number in your manuscript.**
