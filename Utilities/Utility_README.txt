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

1. Run from within the `Utilities` folder by activating `sampleSheets` conda environment
```shell

cd Utilities

conda activate sampleSheets

bash make_sample_sheet.sh /full/path/to/your/data paired

```

## Fastqc runs
Contains two scripts that work in conjunction:

1. run_fastqc.sh

2. make_fastqc_config.sh

The first script is an `SBATCH` driver that internally calls `make_fastqc_config.sh`. By default, after running `multiqc` on a folder of fastqc reports, the sample names are just the regular file names. `make_fastqc_config.sh` maps sample names from your single or paired end sample file (generated in the sample sheet generation code above) and automatically generated a `fastqc_multiqc_config.yaml` file (*which lives in your working directory, NOT the `fastqc_results` directory*)
`Fastqc` is then run on every raw file and then `multiqc` is ran on the resulting log files using the generated config in verbose mode. 

This code assumes that your fastqc results will always be stored in a folder called `fastqc_results` and that your config will always be called `fastqc_multiqc_config.yaml`. The script will automatically generate this output folder for you. 

To run this code, you will need to edit the following variables

1. `DATA_DIR`: path to the raw data folder

2. `SAMPLESHEET`: Path to sample sheet (either single or paired)

3. `LAYOUT`: one of single or paired

### Implementation

To run this code:

1. Run from the utilies folder after pointing to the proper paths in `run_fastqc.sh`

```shell
cd Utilities
sbatch run_fastqc.sh
```

