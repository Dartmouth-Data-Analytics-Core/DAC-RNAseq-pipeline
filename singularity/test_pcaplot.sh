#!/usr/bin/env bash
# Build the slim rna_pcaplot image and check it runs rule pca_plots end to end.
#
# Run from the repo root on the cluster (an interactive node, not the login node):
#   bash singularity/test_pcaplot.sh                      # build from the .def, synthetic counts
#   COUNTS=/path/featurecounts.readcounts.tsv bash singularity/test_pcaplot.sh   # real counts
#   IMG=/path/rna_pcaplot.sif bash singularity/test_pcaplot.sh                   # skip the build
#   IMG=docker://ghcr.io/dartmouth-data-analytics-core/rna_pcaplot:2.2 bash ...  # test the pushed image
set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
OUT=${OUT:-$REPO/pcaplot_test}
SING=$(command -v apptainer || command -v singularity) || { echo "apptainer/singularity not found"; exit 1; }
# build/extract on node-local disk: on /dartfs the extraction alone takes minutes
export APPTAINER_TMPDIR=${APPTAINER_TMPDIR:-/tmp/$USER/apptainer_tmp}
mkdir -p "$APPTAINER_TMPDIR" "$OUT"
BINDS="--bind /dartfs-hpc --bind /dartfs --bind /scratch --bind $REPO --bind $OUT"
FAIL=0
ok()  { echo "  OK   $*"; }
bad() { echo "  FAIL $*"; FAIL=1; }

# ---- 1. build
if [[ -z "${IMG:-}" ]]; then
    IMG=$OUT/rna_pcaplot.sif
    echo "== Building $IMG (APPTAINER_TMPDIR=$APPTAINER_TMPDIR)"
    start=$SECONDS
    # --ignore-fakeroot-command: without /etc/subuid, Apptainer wraps %post in the `fakeroot`
    # LD_PRELOAD shim, which crashes micromamba (std::bad_alloc); the user namespace is enough here
    (cd "$REPO/singularity/defs" && "$SING" build --fakeroot --ignore-fakeroot-command --force "$IMG" rna_pcaplot.def) \
        || { echo "Build failed. If --fakeroot is not allowed here, build with Docker and pass IMG=docker://..."; exit 1; }
    echo "  built in $((SECONDS - start)) s"
fi
[[ -f "$IMG" ]] && echo "  image size: $(du -h "$IMG" | cut -f1)"
run() { "$SING" exec $BINDS "$IMG" "$@"; }

# ---- 2. what is inside
echo "== Contents"
run python -c "import pandas, numpy, scipy, sklearn, plotnine, seaborn, matplotlib, fastcluster; print('python imports ok')" \
    && ok "python imports" || bad "python imports"
run Rscript -e 'suppressPackageStartupMessages(library(DESeq2)); cat("DESeq2", as.character(packageVersion("DESeq2")), "\n")' \
    && ok "R + DESeq2" || bad "R + DESeq2"
for tool in gcc g++ gfortran x86_64-conda-linux-gnu-gcc; do
    run sh -c "command -v $tool" >/dev/null && bad "$tool still in the image" || ok "$tool removed"
done
run sh -c 'ls /opt/env/lib | grep -qi "libQt6\|libclang"' && bad "Qt/LLVM GUI stack still in the image" || ok "no Qt/LLVM"
echo "  env size: $(run du -sh /opt/env | cut -f1)"

# ---- 3. input
if [[ -z "${COUNTS:-}" ]]; then
    COUNTS=$OUT/synthetic.readcounts.tsv
    echo "== Writing synthetic counts: $COUNTS"
    run python - "$COUNTS" "$OUT/samples.csv" <<'EOF'
import sys, numpy as np, pandas as pd
rng = np.random.default_rng(1)
genes, samples = 3000, [f"S{i}" for i in range(1, 7)]
mu = rng.lognormal(4, 1.5, genes)
fold = np.ones((genes, 6)); fold[:300, 3:] = rng.lognormal(1, 0.5, (300, 1))  # group B shift
counts = rng.negative_binomial(5, 5 / (5 + mu[:, None] * fold))
df = pd.DataFrame(counts, columns=[f"bams/{s}.srt.bam" for s in samples],
                  index=pd.Index([f"G{i}" for i in range(genes)], name="Geneid"))
for c in ["Chr", "Start", "End", "Strand"]: df.insert(0, c, "1")
df.insert(4, "Length", 1000)
with open(sys.argv[1], "w") as fh:
    fh.write("# featureCounts-like synthetic table\n")
    df.to_csv(fh, sep="\t")
pd.DataFrame({"sample_id": samples, "group": ["A"] * 3 + ["B"] * 3}).to_csv(sys.argv[2], index=False)
EOF
    META=$OUT/samples.csv
else
    META=${META:-$REPO/sample_fastq_list_single.csv}
fi

# ---- 4. run the rule's command exactly as the Snakefile does
echo "== Running pca_plotting.py"
rm -rf "$OUT/plots"
start=$SECONDS
run python "$REPO/scripts/pca_plotting.py" "$COUNTS" "$OUT/plots" -m "$META" > "$OUT/pca.log" 2>&1 \
    && ok "script exited 0 ($((SECONDS - start)) s)" || { bad "script failed, see $OUT/pca.log"; tail -20 "$OUT/pca.log"; }

# ---- 5. outputs declared by rule pca_plots
echo "== Outputs"
for f in PCA_top_PC1_vs_PC2.png PCA_top_PCA_variance_bar.png PCA_top_PCs.csv \
         PCA_top_variance_explained.csv PCA_all_PCs.csv PCA_all_variance_explained.csv \
         PCA_top_normalized_counts.tsv PCA_all_normalized_counts.tsv PCA_gene_variance.csv; do
    [[ -s "$OUT/plots/$f" ]] && ok "$f" || bad "$f missing or empty"
done
grep -q "^method rlog" "$OUT/plots/deseq2_version.txt" 2>/dev/null \
    && ok "DESeq2 rlog ran ($(head -1 "$OUT/plots/deseq2_version.txt"))" \
    || bad "rlog did not run: $(cat "$OUT/plots/deseq2_version.txt" 2>/dev/null)"
if grep -q "Plotting libraries unavailable" "$OUT/pca.log"; then bad "plotting libs failed to import"
elif [[ -s "$OUT/plots/PCA_top_PC1_vs_PC2.png" ]]; then ok "plotting libs loaded"
else bad "no plots written"; fi
grep -iE "error|traceback" "$OUT/pca.log" && bad "errors in $OUT/pca.log" || ok "no errors in log"

echo
[[ $FAIL == 0 ]] && echo "ALL CHECKS PASSED  (outputs in $OUT/plots)" || echo "SOME CHECKS FAILED"
exit $FAIL
