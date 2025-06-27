#!/bin/bash

#----- Quick helper script to make sample or paired-end sample sheets linked with 
# SLIMS metadata tables generated from GSR-Active. 

#----- Usage menu
usage() {
    echo "Usage: $0 <fastq_dir> <layout>"
    echo ""
    echo "Arguments:"
    echo "  <fastq_dir>    Path to directory containing FASTQ files"
    echo "  <layout>       Either 'paired' or 'single'"
    echo ""
    echo "Example:"
    echo "  $0 ./fastqs paired"
    echo "  $0 ./fastqs single"
    exit 1
}

#----- Check arguments
if [[ $# -ne 2 ]]; then
    echo "Error: Invalid number of arguments." >&2
    usage
fi

#----- Source conda and activate
#source /optnfs/common/miniconda3/etc/profile.d/conda.sh
#conda activate /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/sampleSheets

#----- Set directory containing fastq.gz files
FASTQ_DIR=$1  
OUTPUT_CSV="samples.csv"
LAYOUT=$2

if [[ "$LAYOUT" == "paired" ]]; then
	#----- Generate header
	echo "fastq_1,fastq_2" > "$OUTPUT_CSV"
	#----- Loop over all R1 files and find matching R2
	for r1 in "$FASTQ_DIR"/*_R1*.fastq.gz; do
	echo "$r1"	
#----- Skip if not match 
		[ -e "$r1" ] || continue
		
		#----- Derive R2 filename with sed
		r2="${r1/_R1/_R2}"	
		
		#----- Check that the R2 file exists
		if [[ -f "$r2" ]]; then
			echo "$r1,$r2" >> "$OUTPUT_CSV"
		else
			echo "Warning, No R2 file found for $r1" >&2
		fi
	done
fi
	
if [[ "$LAYOUT" == "single" ]]; then
	#----- Generate header
	echo "fastq_1" > "$OUTPUT_CSV"
	
	for f in "$FASTQ_DIR"/*.fastq.gz; do
  		echo "$f" >> "$OUTPUT_CSV"
	done
fi
	

#----- Run RScript to generate sample sheet
Rscript linkMeta.R "$PWD"  metadata.xlsx samples.csv





