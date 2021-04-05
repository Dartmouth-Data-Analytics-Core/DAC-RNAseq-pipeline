# Dartmouth CQB RNA-seq analysis pipeline

## Introduction 
The pipeline is designed to provide efficient pre-processing and quality control of bulk RNA-sequencing (RNA-seq) data on high performance computing clusters (HPCs) using the Torque/PBS scheduler or using a single high CPU, high RAM machine, and has been made available by the *Data Analytics Core (DAC)* of the *Center for Quantitative Biology (CQB)*, located at Dartmouth College. Both single- and paired-end datasets are supported, in addition to both library preparation methods for full-length or 3'-only analysis. The pipeline has been built and tested using human and mouse data sets. Required software can be installed using Conda with the enrionment file (environment.yml), or specified as paths in the config file.

<img src="logo.jpg" width="250" height="140" >

## Pipeline summary:
The major steps implmented in the pipeline include: 

- FASTQ quality control assesment using [*FASTQC*](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- Read trimming for Poly-A tails, specified adapters, and read quality using [*cutadapt*](https://cutadapt.readthedocs.io/en/stable/)
- Alignment using [*Hisat*] or [*STAR*](https://github.com/alexdobin/STAR)
- Quantification with [*Featurecouts*](http://subread.sourceforge.net/) and [*RSEM*](https://deweylab.github.io/RSEM/)

All of these tools can be installed in the [conda environment](https://docs.conda.io/en/latest/). As input, the pipeline takes raw data in FASTQ format, and produces quantified read counts (using *HTSeq-Count* or *RSEM*) as well as a detailed quality control report (including pre- and post-alignment QC metrics) for all processed samples. Quality control reports are aggregated into HTML files using *MultiQC*. 


## Implementation
The pipeline uses Snakemake to submit jobs to the scheduler, or spawn processes on a single machine, and requires several variables to be configured by the user when running the pipeline: 

* **sample_tsv** - A TSV file containing sample names and paths to fastq paths.  See example in this repository for formatting.
* **layout** - Either "single" or "paired" library construction.

* **aligner_name** - 'hisat' or 'star'
* **aligner_index** - Path to Hisat or STAR genome reference index

* **picard_rrna_list** - Absolute path to coordinates of ribosomal RNA sequences in reference genome, in [interval-list format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists)
* **picard_refflat** - Absolute path to genome annotation in [RefFlat format](https://gatk.broadinstitute.org/hc/en-us/articles/360040509431-CollectRnaSeqMetrics-Picard-)

* **annotation_gtf** - Absolute path to genome annotation file (.gtf) of [*Featurecouts*](http://subread.sourceforge.net/) or [*RSEM*](https://deweylab.github.io/RSEM/)
* **picard_strand** - "FIRST_READ_TRANSCRIPTION_STRAND" "SECOND_READ_TRANSCRIPTION_STRAND"
* **featurecounts_strand** - "1" or "2" #1 for first read transcription strand, 2 for second.*

### Command Line Examples



> **Contact & questions:** 
> Please address questions to *DataAnalyticsCore@groups.dartmouth.edu* or generate a issue in the GitHub repository. 

> **This pipeline was created with funds from the COBRE grant **1P20GM130454**. 
> If you use the pipeline in your own work, please acknowledge the pipeline by citing the grant number in your manuscript.**

