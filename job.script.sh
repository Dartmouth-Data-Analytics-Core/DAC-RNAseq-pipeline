#!/bin/bash

# Name of the job
#SBATCH --job-name=RNAseq_UMI

# Number of compute nodes
#SBATCH --nodes=1

#SBATCH --partition=standard

#SBATCH --account=dac

# Walltime (job duration)
#SBATCH --time=60:00:00

# Email address
#SBATCH --mail-user=fwk@dartmouth.edu

# Email notifications (comma-separated options: BEGIN,END,FAIL)
#SBATCH --mail-type=FAIL

source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq 
snakemake -T 10 --profile /dartfs-hpc/rc/lab/G/GSR_Active/Labs/Jakubzick/RNAseq_4-18-22/combined/cluster_profile/
