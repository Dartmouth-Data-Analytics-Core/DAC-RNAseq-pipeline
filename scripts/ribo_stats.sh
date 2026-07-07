#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <ribodetector_root_dir>" >&2
  exit 1
fi

RESULTS_DIR="$1"

if [[ ! -d "$RESULTS_DIR" ]]; then
  echo "Error: '$RESULTS_DIR' is not a directory" >&2
  exit 1
fi

out="ribodetector/rrna_norrna_pct_mqc.tsv"
mkdir -p "$(dirname "$out")"

# write MultiQC custom-content header
cat > "$out" <<'EOF'
# id: 'ribodetector'
# section_name: 'Ribodetector rRNA Detection'
# description: 'Percentage of reads classified as rRNA by ribodetector. High values indicate poor rRNA depletion during library prep.'
# format: 'tsv'
# plot_type: 'table'
# pconfig:
#     id: 'ribodetector_table'
#     title: 'Ribodetector: rRNA contamination'
# headers:
#     total_reads:
#         title: 'Total reads'
#         description: 'Total input reads processed by ribodetector'
#         format: '{:,.0f}'
#         scale: 'Blues'
#     rrna_reads:
#         title: 'rRNA reads'
#         description: 'Reads classified as rRNA (filtered out)'
#         format: '{:,.0f}'
#         scale: 'Reds'
#     pct_rrna:
#         title: '% rRNA'
#         description: 'Percentage of input reads classified as rRNA'
#         format: '{:,.2f}'
#         suffix: '%'
#         min: 0
#         max: 100
#         scale: 'RdYlGn-rev'
EOF

# append column header and data
printf "sample\ttotal_reads\trrna_reads\tpct_rrna\n" >> "$out"

for d in "$RESULTS_DIR"/*/; do
  sample=$(basename "$d")
  log="$d/stderr.log"

  [[ -f "$log" ]] || { echo "skip: no stderr.log in $d" >&2; continue; }

  total=$(sed 's/\x1B\[[0-9;]*[A-Za-z]//g' "$log" \
    | grep -oP '[0-9]+(?=\s+sequences loaded)' | head -1)

  rrna=$(sed 's/\x1B\[[0-9;]*[A-Za-z]//g' "$log" \
    | grep -oP '(?<=Detected\s)[0-9]+(?=\s+rRNA sequences)' | head -1)

  [[ -n "${total:-}" && -n "${rrna:-}" ]] || { echo "skip: missing counts in $log" >&2; continue; }
  [[ "$total" =~ ^[0-9]+$ && "$total" -gt 0 ]] || { echo "skip: bad total in $log" >&2; continue; }
  [[ "$rrna"  =~ ^[0-9]+$ ]] || { echo "skip: bad rrna in $log" >&2; continue; }

  pct=$(awk -v r="$rrna" -v t="$total" 'BEGIN{printf "%.4f", 100*r/t}')

  printf "%s\t%s\t%s\t%s\n" "$sample" "$total" "$rrna" "$pct" >> "$out"
done

echo "Saved: $out"