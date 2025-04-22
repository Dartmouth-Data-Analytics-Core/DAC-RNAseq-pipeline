# Dartmouth CQB RNA-seq analysis pipeline

## Introduction 
The pipeline is designed to provide efficient pre-processing and quality control of bulk RNA-sequencing (RNA-seq) data on high performance computing clusters (HPCs) using the Torque/PBS scheduler or using a single high CPU, high RAM machine, and has been made available by the *Data Analytics Core (DAC)* of the *Center for Quantitative Biology (CQB)*, located at Dartmouth College. Both single- and paired-end datasets are supported, in addition to both library preparation methods for full-length or 3'-only analysis. The pipeline has been built and tested using human and mouse data sets. Required software can be installed using Conda with the enrionment file (environment.yml), or specified as paths in the config.yaml file.

<img src="img/cqb_logo.jpg" width="250" height="140" >

## Pipeline summary:
The major steps implmented in the pipeline include: 

- FASTQ quality control assesment using [*FASTQC*](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- Read trimming for Poly-A tails, specified adapters, and read quality using [*cutadapt*](https://cutadapt.readthedocs.io/en/stable/)
- Alignment using [*Hisat*](https://daehwankimlab.github.io/hisat2/) or [*STAR*](https://github.com/alexdobin/STAR)
- Quantification with [*Featurecouts*](http://subread.sourceforge.net/) and [*RSEM*](https://deweylab.github.io/RSEM/)

All of these tools can be installed in a [conda environment](https://docs.conda.io/en/latest/) or on paths available to a computing server. As input, the pipeline takes raw data in FASTQ format, and produces quantified read counts (using *HTSeq-Count* or *RSEM*) as well as a detailed quality control report (including pre- and post-alignment QC metrics) for all processed samples. Quality control reports are aggregated into HTML files using *MultiQC*. 

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

  
## Running tests using pre-built environments on Discovery
Clone this repository:
```shell
git clone -b dev-Mike https://github.com/Dartmouth-Data-Analytics-Core/DAC-RNAseq-pipeline.git
cd DAC-RNAseq-pipeline
```
Activate an environment containing Snakemake:
```shell
conda activate /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/snakemake
```

Build, configure, and check reference files:
```shell

# INDEX SEQUENCE/ANNOTATION FILES FOR MAPPING 
snakemake -s Snakefile  --use-conda -j 6 --conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/DAC-RNAseq-pipeline build_refs

# ADD REFERENCE DETAILS TO CONFIG FILE
cat ref/pipeline_refs/hg38_chr567_100k.entries.yaml >> config.yaml

# CHECK REFERENCE IS CORRECTLY FORMATED
snakemake -s Snakefile  --use-conda -j 6 --conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/DAC-RNAseq-pipeline check_refs
```
Run the pipeline:
```shell
snakemake -s Snakefile  --use-conda -j 6 --conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/DAC-RNAseq-pipeline
```
  
## Running the pipeline using pre-built references and config files on Discovery
The DAC has made public references and their corresponding aligner index and annotation files available to the Dartmouth community on Discovery/DartFS.  Additional documentation on the public references can be found [in thieir repository](https://github.com/Dartmouth-Data-Analytics-Core/DAC-Genome-References).  Pre-built config.yaml files for this RNA-Seq pipeline have also been added to the prebuilt_configs directory of this repository.  As of 4/29/24, there are configs for any combination of human/mouse, single/paired reads, and Hisat2/STAR/RSEM.  An example of using a prebuilt config for human, Hisat2, paired-end reads is as follows: 
```shell
snakemake -s Snakefile --configfile prebuilt_configs/human_config_paired_hisat.yaml  --use-conda -j 6 --conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/DAC-RNAseq-pipeline
```
When using a pre-built config, one will still have to create a sample_fastq_list.txt for each specific run, and ensure this file is specified correctly at the top of the config file.


  
## More Command Line Examples
Submit the pipeline to a single machine, allowing usage of 40 cores:
```shell
snakemake --use-conda -s Snakefile -j 40
```

Submit the pipeline to a computing cluster using the profile defined in cluster_profile/config.yaml, and allow jobs to be re-run twice in case of failure:
```shell
snakemake --use-conda -s Snakefile --profile cluster_profile -T 2
```

### Snakemake job graph example for three samples:
<img src="img/dag.svg" width="1024" height="300" >

**Contact & questions:** 
Please address questions to *DataAnalyticsCore@groups.dartmouth.edu* or submit an issue in the GitHub repository. 

**This pipeline was created with funds from the COBRE grant **1P20GM130454**. 
If you use the pipeline in your own work, please acknowledge the pipeline by citing the grant number in your manuscript.**

