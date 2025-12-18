#!/bin/bash

#SBATCH --job-name=fastQ_Screen
#SBATCH --nodes=1
#SBATCH --partition=standard
#SBATCH --time=60:00:00
#SBATCH --mail-user=XXXX@dartmouth.edu
#SBATCH --mail-type=FAIL

# Define the data directory
DATA_DIR="/dartfs-hpc/rc/lab/G/GSR_Active/Labs/"

# Output directory for FastQC results
OUTPUT_DIR=""

# Path to FastQC (if not in PATH, specify full path here)
FASTQ_CMD="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/sullivan/tools/fastq-screen/FastQ-Screen-0.14.1/fastq_screen"
CONF_FILE="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/fastq_screen/fastq_screen.conf"

# Check if the data directory exists
if [[ ! -d "$DATA_DIR" ]]; then
  echo "Error: Data directory $DATA_DIR not found."
  exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Navigate to the data directory
cd "$DATA_DIR" || exit

# Run FastQC for every FASTQ file in the directory
for fastq_file in *.fastq.gz; do
  # Check if there are any FASTQ files in the directory
  if [[ ! -e "$fastq_file" ]]; then
    echo "No FASTQ files found in $DATA_DIR."
    exit 0
  fi

  # Run FastQC
  echo "Running FastQC on $fastq_file..."
  "$FASTQ_CMD" --conf "$CONF_FILE" --outdir "$OUTPUT_DIR" --threads 4 "$fastq_file"
done

echo "FastQ Screen processing completed. Results are in $OUTPUT_DIR."

# Run Multiqc
multiqc "$OUTPUT_DIR"
