# Utilitiy_Scripts
🛠️ Helpful day to day utilities that tie into GDSC-Pipelines and make life easier!


## Sample Sheet Generation

Utility scripts to automatically generate sample sheets for the pipeline that assign external sample names provided to the Genomics Shared Resources. There are wo scripts that work in conjunction with one another:

1. `make_sample_sheet.sh`
2. `linkMeta.R`

`make_sample_sheet.sh` is a driver script which you will call to generate a temporary sample sheet. `linkMeta.R` is called internally to open the slims metadata sheet from the sequencing facility and link file names to external IDs, preventing human error in making the sample sheet. By default, this pipeline outputs comma-separated files (.csv) to seamlessly tie into the pipeline.

### Implementation

To view usage menu for the code, run the following:

```shell
bash make_sample_sheet.sh
```
The script takes two main arguments, the path to the raw data files on `GSR_Active` and a library layout (one of single or paired).

**By deafult, in the GSR folder, there should be an xlsx file named `metadata.xlsx`. Ensure this file is present before running the script!**

1. Run from **within** the `Utilities` folder by activating `sampleSheets` conda environment
```shell

# Mv into utility fodler
cd Utilities

# Activate conda environment
conda activate /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/sampleSheets

# Run the code to make the sample sheet
bash make_sample_sheet.sh /dartfs-hpc/rc/lab/G/GSR_Active/Labs/YourLab/YourProject paired

# Move newly generated sample sheet out of utility folder
mv sample_fastq_list_<layout>.csv ../

```

This will generate either `sample_fastq_list_single.csv` or `sample_fastq_list_paired.csv` **in the utilities folder** depending on your layout within your cloned DAC-RNASeq-Pipeline folder.
You can now specify this file name in `config.yaml` or in any of the configs in `prebuilt_configs`.

**Be sure to manually check the sample sheet before running!**

## Fastqc runs
Contains two scripts that work in conjunction:

1. `run_fastqc.sh`

2. `make_fastqc_config.sh`

`run_fastqc.sh` is an `SBATCH` driver script that you will interact with that internally calls `make_fastqc_config.sh`. By default, after running `multiqc` on a folder of fastqc reports, the sample names reflect the fastq file names. `make_fastqc_config.sh` maps sample names from your single or paired end sample file (generated in the sample sheet generation code above) and automatically generated a `fastqc_multiqc_config.yaml`.
`Fastqc` is then run on every raw file and then `multiqc` is ran on the resulting log files using the generated config in verbose mode. 

The script will automatically generate a folder called `fastqc_results` where all the outputs will live.

**Do not change the name of the sample sheet generated in the sample sheet generation script, as the fastqc script relies on those names**

To run this code, you will need to edit the following variables

1. `DATA_DIR`: path to the raw data folder

2. `LAYOUT`: one of single or paired

### Implementation

To run this code:

1. Run from the utilies folder after pointing to the proper paths in `run_fastqc.sh`

```shell
cd Utilities
sbatch run_fastqc.sh
```

Similarly, to check for contamination, you can run FastqScreen using the `run_fastq_screen` script. Modify the following variables: `DATA_DIR` and `OUTPUT_DIR`.

## TPM Annotation

The `Annotate_TPMs.R` script generates a visually annotated Excel spreadsheet of transcript-per-million (TPM) expression values for genes. This script is designed to help quickly assess gene expression patterns as a post-pipeline QC check. The key features of this script are:

- **Average TPM calculation:** Computes the mean TPM for each gene across all samples

- **Sorting:** Orders genes by decreasing average TPM to show highly expressed genes

- **Heatmap-style annotation:** Applies a gradient color scale to TPM values across all samples to enable the visualization of highly expressed genes. This also enables the assessment of whether top genes are mitochondrial in nature which could indicate sample stress.

This script requires a TPM tsv file as output by the pipeline (`featurecounts.readcounts_tpm.tsv` or `featurecounts.readcounts_tpm.ann.tsv`) and an output path/filename that **must** end in .xlsx. 


Example usage:

```
Rscript Annotate_TPMs.R <input_tpm_file> <output_excel_path>
```
