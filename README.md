# Dartmouth Analytics core RNA-seq pipeline
**This pipeline was created with funds from the COBRE grant number 1P20GM130454. 
If you use this code and find it helpful please cite us!**

- Both full-length transcript & 3'-end counting assays are supported
- Handles both single- and paired-end data sets 



This pipeline runs through a quality assessment of your raw data (fastQC), trimming of low quality bases (cutadapt), read alignment (STAR), and quantifying reads mapped (picard and RSEM).



This R code is written using R functions. The functions require some variables to be defined by the user when calling the script. The following variables must be defined by the user. 

**Lab** - The name of the lab, this is used as a prefix for some of the files generated,

**FastqRaw** - The location (absolute path) of the raw fastq files.

**SamNames** - The sample names used on the raw fastq files - everything that comes before "_R1_001.fastq.gz or _R2_001.fastq.gz" should be in your SamNames file this variable is generally going to be a list and should be defined as a list SamName<- c("SamName1", "SamName2", "SamName3", etc.).

**SeqMethod** - This is either "fullLength" or "3Prime" depending on the type of data you are analyzing.

**AlignInd** - This is the reference index that you would like to use during the alignment step please give an absolute path(*_index).

**AlignRef** - This is the reference that you would like to use during the alignment step, please give an absolute path (*.gtf).

**PicardInt** - This is the interval list that picard will use to generate raw count data, please give an absolute path(*.rRNA.interval_list).

**PicardRef** - This is the reference file that picard will use to generate raw count data, please give an absolute path (*.refFlat.txt).

**QuantRef** - This is the file used by RSEM to generate raw counts, please give an absolute path(RSEMref).

**CondaEnv**- This is the environment that includes all of the dependencies needed to run this pipeline, the yml file to create this environment is included in this directory (environment.yml).

**OutputFolder** - This is the folder where you would like the output files to be directed to, you should create the following directories in this folder before running this script:
tmp/
fastqc/
trim/
alignment/
rawcounts/

**This pipeline was created with funds from the COBRE grant number 1P20GM130454. 
If you use this code and find it helpful please cite us!**
