#!/bin/bash

# Name of the job
#SBATCH --job-name=human_snakemake

# Number of compute nodes
#SBATCH --nodes=1

# specify que to submit to (preempt1 only contains q10 DAC node)
#SBATCH --partition=preempt1

# specify account you are submitting from
#SBATCH --account=dac

# Walltime (job duration)
#SBATCH --time=60:00:00

# Email address
#SBATCH --mail-user=omw@dartmouth.edu

# Email notifications (comma-separated options: BEGIN,END,FAIL)
#SBATCH --mail-type=FAIL

source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate snakemake-tutorial
snakemake -T 10 --profile cluster_profile --rerun-incomplete
