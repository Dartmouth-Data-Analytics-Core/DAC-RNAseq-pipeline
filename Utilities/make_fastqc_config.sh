#!/bin/bash

# Usage check
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <sample_sheet.csv> <layout: paired|single>"
  exit 1
fi

INPUT="$1"
LAYOUT="$2"
OUTPUT="fastqc_results/fastqc_multiqc_config.yaml"

# Validate layout
if [[ "$LAYOUT" != "paired" && "$LAYOUT" != "single" ]]; then
  echo "Error: layout must be either 'paired' or 'single'"
  exit 1
fi

# Start YAML file
echo "sample_names_replace:" > "$OUTPUT"

# Skip header and process lines
tail -n +2 "$INPUT" | while IFS=',' read -r sample fastq1 fastq2; do
  base1=$(basename "$fastq1" .fastq.gz)

  if [[ "$LAYOUT" == "paired" ]]; then
    base2=$(basename "$fastq2" .fastq.gz)
    echo "  \"$base1\": \"${sample}_Forward\"" >> "$OUTPUT"
    echo "  \"$base2\": \"${sample}_Reverse\"" >> "$OUTPUT"
  elif [[ "$LAYOUT" == "single" ]]; then
    echo "  \"$base1\": \"${sample}\"" >> "$OUTPUT"
  fi
done

echo "✅ multiqc config written to $OUTPUT"
