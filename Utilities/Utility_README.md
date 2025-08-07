# Utilitiy_Scripts
🛠️ Helpful day to day utilities that tie into GDSC-Pipelines and make life easier!


## Sample Sheet Generation
There are wo scripts that work in conjunction with one another:

1. `make_sample_sheet.sh`
2. `linkMeta.R`

The first script is a driver which generates a temporary sample sheet. The R-script makes use of `openxlsx` and `stringr` to open the slims metadata sheet from the sequencing facility and link file names to external IDs, preventing human error in making the sample sheet. By default, this pipeline outputs comma-separated files (.csv) but can be easily modified to output tab-separated files if desired in the `linkMeta.R` file. 

### Implementation
To view usage menu for the code, run the following:
```shell
bash make_sample_sheet.sh
```
The script takes two main arguments, the path to the raw data files on `GSR_Active` and a library layout (one of single or paired).
**By deafult, in the GSR folder, there should be an xlsx file named `metadata.xlsx`. Ensure this file is present before running the script!**

1. Run from within the `Utilities` folder by activating `sampleSheets` conda environment
```shell

cd Utilities

conda activate /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/sampleSheets

bash make_sample_sheet.sh /dartfs-hpc/rc/lab/G/GSR_Active/Labs/YourLab/YourProject paired

```

This will generate either `sample_fastq_list_single.csv` or `sample_fastq_list_paired.csv` depending on your layout within your cloned DAC-RNASeq-Pipeline folder.
You can now specify this file name in `config.yaml` or in any of the configs in `prebuilt_configs`.

## Fastqc runs
Contains two scripts that work in conjunction:

1. run_fastqc.sh

2. make_fastqc_config.sh

The first script is an `SBATCH` driver that internally calls `make_fastqc_config.sh`. By default, after running `multiqc` on a folder of fastqc reports, the sample names are just the regular file names. `make_fastqc_config.sh` maps sample names from your single or paired end sample file (generated in the sample sheet generation code above) and automatically generated a `fastqc_multiqc_config.yaml`.
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

