#!/bin/bash

# Name of the job
#SBATCH --job-name=RNAseq.preprocess
#SBATCH --nodes=1
#SBATCH --partition=standard
#SBATCH --time=60:00:00
#SBATCH --mail-user=XXXXXXXXX@dartmouth.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=RNAseq.preprocess_%j.out

source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/snakemake
snakemake -s Snakefile  \
	--conda-frontend conda \
	--use-conda \
	--conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/DAC-RNAseq-pipeline \
	--profile cluster_profile \
	--rerun-incomplete \
	--keep-going
