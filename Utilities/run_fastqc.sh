#!/bin/bash

#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#SBATCH --partition=preempt1
#SBATCH --account=dac
#SBATCH --time=60:00:00
#SBATCH --mail-user=f007qps@dartmouth.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=RNAseq_raw_fastqc%j.out

# Define the symlinked data directory
DATA_DIR="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/Labs/lacefield/250619_Zach_hyperTribe/GDSC-HyperTRIBE-pipeline/data/run2"

# Output directory for FastQC results
OUTPUT_DIR="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/Labs/lacefield/250619_Zach_hyperTribe/GDSC-HyperTRIBE-pipeline/fastqc_results/"

# Path to FastQC (if not in PATH, specify full path here)
FASTQC_CMD="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/sullivan/tools/fastqc/FastQC/fastqc"

# DEFINE THE PATH TO THE SAMPLE SHEET
SAMPLESHEET="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/Labs/lacefield/250619_Zach_hyperTribe/GDSC-HyperTRIBE-pipeline/sample_fastq_list_paired.csv"

# DEFINE THE LAYOUT
LAYOUT="paired"
#-----------------------------------------#
#-----------------------------------------#

# Check if the data directory exists
if [[ ! -d "$DATA_DIR" ]]; then
  echo "Error: Data directory $DATA_DIR not found."
  exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

#----- Create yaml file for sample renaming
bash make_fastqc_config.sh "$SAMPLESHEET" "$LAYOUT"

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
  "$FASTQC_CMD" -o "$OUTPUT_DIR" "$fastq_file"
done

#echo "FastQC processing completed. Results are in $OUTPUT_DIR."
cd "$OUTPUT_DIR"
multiqc . -c ../fastqc_multiqc_config.yaml --verbose
