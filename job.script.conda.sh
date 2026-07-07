#!/bin/bash

#SBATCH --job-name=rnaseq                          
#SBATCH --nodes=1
#SBATCH --partition=standard
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16  
#SBATCH --time=60:00:00
#SBATCH --mail-user=f007qps@dartmouth.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=%x_%j.log
#========================================================#

#----- Environment information
CONDA_BASE="/optnfs/common/miniconda3"
SNAKEMAKE_ENV="/dartfs/rc/nosnapshots/G/GMBSR_refs/envs/snakemake"
CONDA_PREFIX_PATH="/dartfs/rc/nosnapshots/G/GMBSR_refs/envs/DAC-RNAseq-pipeline"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${SNAKEMAKE_ENV}"

#----- LOGGER
cat <<EOF
#───────────────────────── Initialization ──────────────────────────#
Running GDSC-bulkRNA-Seq Pipeline v2.0 with Snakemake $(snakemake --version)

Job:        $SLURM_JOB_NAME
Job ID:     $SLURM_JOB_ID
Node:       $(hostname)
Start time: $(date)
Work dir:   $(pwd)
Conda base: $CONDA_BASE
Snakemake:  $SNAKEMAKE_ENV
Binary:     $(which snakemake)
#───────────────────────────────────────────────────────────────────#

SNAKEMAKE LOG:
EOF

#----- Make slurm logs
mkdir -p slurm_logs

#----- Invoke Snakemake
snakemake -s \
    Snakefile \
    --profile cluster_profile \
    -T 2 \
    --use-conda \
    --conda-frontend conda \
    --conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/DAC-RNAseq-pipeline

#----- Capture exit status
SNAKEMAKE_EXIT=$?

#----- Final status
echo ""
echo "#────────────────────────── Job Complete ──────────────────────────#"
echo "End time: $(date)"
if [ $SNAKEMAKE_EXIT -eq 0 ]; then
    echo "Status: SUCCESS"
else
    echo "Status: FAILED (exit code: $SNAKEMAKE_EXIT)"
fi

exit $SNAKEMAKE_EXIT
echo "End time: $(date)"
