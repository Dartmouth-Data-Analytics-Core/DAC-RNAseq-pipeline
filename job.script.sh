#!/bin/bash

#SBATCH --job-name=rnaseq                          
#SBATCH --nodes=1
#SBATCH --partition=standard
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16  
#SBATCH --time=60:00:00
#SBATCH --mail-type=FAIL
#SBATCH --output=%x_%j.log
#========================================================#

#----- Environment information
CONDA_BASE="/optnfs/common/miniconda3"
SNAKEMAKE_ENV="/dartfs/rc/nosnapshots/G/GMBSR_refs/envs/snakemake"
CONDA_PREFIX_PATH="/dartfs/rc/nosnapshots/G/GMBSR_refs/envs/DAC-RNAseq-pipeline"
CONTAINERS_PATH="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/shared-software/singularity-containers/bulk-RNASeq/"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${SNAKEMAKE_ENV}"

#----- LOGGER
cat <<EOF
#───────────────────────── Initialization ──────────────────────────#
Running GDSC-bulkRNA-Seq Pipeline v2.0.2 with Snakemake $(snakemake --version)

Job:        $SLURM_JOB_NAME
Job ID:     $SLURM_JOB_ID
Node:       $(hostname)
Start time: $(date)
Work dir:   $(pwd)
Conda base: $CONDA_BASE
Snakemake:  $SNAKEMAKE_ENV
Binary:     $(which snakemake)
SIFs:       $CONTAINERS_PATH
#───────────────────────────────────────────────────────────────────#

SNAKEMAKE LOG:
EOF

#----- Make slurm logs
mkdir -p slurm_logs

#----- Singularity arguments
SINGARGS="--bind /dartfs-hpc --bind /dartfs --bind /scratch"

#----- Temp dir for converting images (only matters when an image must be downloaded)
# Images are converted once and cached in CONTAINERS_PATH, shared by everyone. The first
# run after a new image/tag downloads and converts it; on /dartfs scratch that takes ~17 min,
# on the node's local /tmp ~30 s. Uncomment if this run has to download an image.
# export APPTAINER_TMPDIR=/tmp/$USER/apptainer_tmp SINGULARITY_TMPDIR=/tmp/$USER/apptainer_tmp
# mkdir -p "$APPTAINER_TMPDIR"

#----- Invoke Snakemake
snakemake -s \
    Snakefile \
    --profile cluster_profile \
    -T 2 \
    --use-singularity \
    --singularity-args "${SINGARGS}" \
    --singularity-prefix "${CONTAINERS_PATH}"

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
