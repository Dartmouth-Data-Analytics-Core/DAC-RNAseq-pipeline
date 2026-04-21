#!/bin/bash

#SBATCH --job-name=Hixon_RNA
#SBATCH --nodes=1
#SBATCH --partition=preempt1
#SBATCH --account=dac
#SBATCH --time=60:00:00
#SBATCH --mail-user=f007qps@dartmouth.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=_Hixon_RNA_%j.log

#----- Source conda environment
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/snakemake

#----- Make slurm logs folder
mkdir -p slurm_logs

#! If using 'RUN_RUSTQC: True', include the following line in the snakemake call
# --use-singularity \
# --singularity-args "--bind /dartfs,/optnfs" \

#----- Call Snakemake
snakemake -s Snakefile \
	--use-conda \
	--conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/DAC-RNAseq-pipeline \
	--profile cluster_profile \
	--rerun-incomplete \
	--keep-going
