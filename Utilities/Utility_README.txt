# Utilitiy_Scripts
🛠️ Helpful day to day utilities that tie into GDSC-Pipelines and make life easier!

## Table of Contents
- [Sample Sheet Generation](#sample-sheet-generation)
- [Fastqc Runs](#fastqc-runs)


## Sample Sheet Generation
The [sample sheet generation folder](https://github.com/mikemartinez99/Utilitiy_Scripts/tree/main/sample_sheet_generation) contains two scripts that work in conjunction with one another:

1. `make_sample_sheet.sh`
2. `linkMeta.R`

The first script is a driver which generates a temporary sample sheet. The R-script makes use of `openxlsx` and `stringr` to open the slims metadata sheet from the sequencing facility and link file names to external IDs, preventing human error in making the sample sheet. By default, this pipeline outputs comma-separated files (.csv) but can be easily modified to output tab-separated files if desired in the `linkMeta.R` file. 

### Implementation
To view usage menu for the code, run the following:
```shell
bash make_sample_sheet.sh
```

1. Copy the scripts from the utilities folder to your working directory
```shell
cd <workingDir>
cp /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/utilities/general_pipeline_utilities/sample_sheet_generation/* .
```

2. Ensure that you have the slims metadata file in your working directory. It should be named `metadata.xlsx` by default. If it is not, change the file name to represent this
```shell
mv <some_file.xlsx` metadata.xlsx
```

**If on discovery, polaris, or andes, follow steps 3A. If local, skip to step 4**

3A. If on discovery, activate the conda environment with all R dependencies
```shell
conda activate sampleSheets
```

4. Run the script directly on data path or on a symlink of data
```shell
# If data is stored in paired end
bash make_sample_sheet.sh ./data paired

# If data is single end
bash make_sample_sheet.sh ./data single
```

**NOTE** Depending on the pipeline you are running, you may need to modify the sample sheet header to include extra fields. This can be easily achieved in `bbedit`, `sublime`, `nano` or your favorite text editor. 

## Fastqc runs
The [fastqc runs folder](https://github.com/mikemartinez99/Utilitiy_Scripts/tree/main/fastqc_runs) contains two scripts that work in conjunction:

1. run_fastqc.sh

2. make_fastqc_config.sh

The first script is an `SBATCH` driver that internally calls `make_fastqc_config.sh`. By default, after running `multiqc` on a folder of fastqc reports, the sample names are just the regular file names. `make_fastqc_config.sh` maps sample names from your single or paired end sample file (generated in the sample sheet generation code above) and automatically generated a `fastqc_multiqc_config.yaml` file (*which lives in your working directory, NOT the `fastqc_results` directory*)
`Fastqc` is then run on every raw file and then `multiqc` is ran on the resulting log files using the generated config in verbose mode. 

This code assumes that your fastqc results will always be stored in a folder called `fastqc_results` and that your config will always be called `fastqc_multiqc_config.yaml`.

To run this code, you will need to edit the following variables

1. `DATA_DIR`: path to the raw data folder

2. `OUTPUT_DIR`: path to `fastqc_results` folder

3. `SAMPLESHEET`: Path to sample sheet (either single or paired)

4. `LAYOUT`: one of single or paired

### Implementation

To run this code:

1. copy these two files into your working directory

```shell
cd <workingDir>
cp /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/utilities/general_pipeline_utilities/fastqc_runs/* .
```

2. Modify the variables in the driver script `run_fastq.sh`

3. Submit the script

```shell
sbatch run_fastqc.sh
```




