#!/usr/bin/env python3
"""
GDSC Bulk RNA-Seq — interactive HTML dashboard builder.

Collates the QC and quantification outputs of the DAC/GDSC-RNAseq pipeline into a
single self-contained HTML file: no CDN, no network, no runtime dependencies for
the reader. Everything (CSS, JS, data) is inlined, so the report can be dropped on
a WebShare link and opened offline.

Inputs are read straight from the pipeline's own outputs rather than from the
rendered MultiQC HTML, so metrics that MultiQC only stores in its compressed
plotdata (gene-body coverage, strand correctness) come from the authoritative
Picard files instead of a decode step.

Usage
-----
    python scripts/build_dashboard.py                     # run in the pipeline dir
    python scripts/build_dashboard.py -d /path/to/run -o QC_Dashboard.html

See --help for the full option list.
"""

import argparse
import base64
import csv
import glob
import html
import json
import os
import re
import sys
from datetime import datetime

import numpy as np
import pandas as pd

VERSION = "1.0.0"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SMALL HELPERS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def log(msg):
    print(f"[dashboard] {msg}", file=sys.stderr)


def warn(msg):
    print(f"[dashboard] WARNING: {msg}", file=sys.stderr)


def to_float(x):
    """Best-effort numeric coercion; returns None for blanks, '?', 'NA', etc."""
    if x is None:
        return None
    if isinstance(x, (int, float, np.integer, np.floating)):
        return None if (isinstance(x, float) and np.isnan(x)) else float(x)
    s = str(x).strip().replace(",", "").replace("%", "")
    if s in ("", "NA", "N/A", "?", "nan", "NaN", "-", "."):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def first_existing(*paths):
    for p in paths:
        if p and os.path.exists(p):
            return p
    return None


def clean_sample_name(name):
    """Strip the path/extension decoration featureCounts and friends leave behind."""
    n = os.path.basename(str(name).strip())
    for ext in (
        ".srt.bam", ".mkdup.bam", ".bam", ".Log.final.out",
        ".hisat.summary.txt", ".picard.rna.metrics.txt",
        ".mkdup.log.txt", ".cutadapt.report",
    ):
        if n.endswith(ext):
            n = n[: -len(ext)]
    return n


def safe_div(num, den):
    if num is None or den in (None, 0):
        return None
    return num / den


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PARSERS - one per upstream tool
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def parse_star_log(path):
    """STAR Log.final.out -> dict of the fields we surface."""
    out = {}
    keymap = {
        "Number of input reads": "total_reads",
        "Uniquely mapped reads number": "uniquely_mapped",
        "Uniquely mapped reads %": "pct_unique",
        "% of reads mapped to multiple loci": "pct_multi",
        "% of reads mapped to too many loci": "pct_multi_many",
        "% of reads unmapped: too short": "pct_too_short",
        "% of reads unmapped: too many mismatches": "pct_unmapped_mm",
        "% of reads unmapped: other": "pct_unmapped_other",
        "Average input read length": "avg_read_length",
        "Average mapped length": "avg_mapped_length",
        "Number of splices: Total": "n_splices",
        "Mismatch rate per base, %": "mismatch_rate",
    }
    with open(path, errors="ignore") as fh:
        for line in fh:
            if "|" not in line:
                continue
            k, v = line.split("|", 1)
            k, v = k.strip(), v.strip()
            if k in keymap:
                out[keymap[k]] = to_float(v)

    # STAR reports "too many loci" separately; the QC read is the combined figure.
    multi = (out.get("pct_multi") or 0) + (out.get("pct_multi_many") or 0)
    out["pct_multimapped"] = multi if ("pct_multi" in out or "pct_multi_many" in out) else None
    unmapped = sum(
        out.get(k) or 0 for k in ("pct_too_short", "pct_unmapped_mm", "pct_unmapped_other")
    )
    out["pct_unmapped"] = unmapped if out.get("pct_too_short") is not None else None
    out["aligner"] = "STAR"
    return out


def parse_hisat_summary(path):
    """HISAT2 --summary-file -> the same normalized field names as STAR."""
    out = {"aligner": "HISAT2"}
    txt = open(path, errors="ignore").read()

    m = re.search(r"^\s*(\d+)\s+reads; of these:", txt, re.M)
    if m:
        out["total_reads"] = float(m.group(1))

    # "aligned exactly 1 time" appears for both SE and PE blocks; sum them.
    exactly1 = sum(int(x) for x in re.findall(r"^\s*(\d+)\s+\([\d.]+%\)\s+aligned (?:concordantly )?exactly 1 time", txt, re.M))
    gt1 = sum(int(x) for x in re.findall(r"^\s*(\d+)\s+\([\d.]+%\)\s+aligned (?:concordantly )?>1 times", txt, re.M))
    if exactly1:
        out["uniquely_mapped"] = float(exactly1)

    total = out.get("total_reads")
    if total:
        if exactly1:
            out["pct_unique"] = 100.0 * exactly1 / total
        if gt1:
            out["pct_multimapped"] = 100.0 * gt1 / total

    m = re.search(r"([\d.]+)%\s+overall alignment rate", txt)
    if m:
        rate = float(m.group(1))
        out["overall_alignment_rate"] = rate
        out["pct_unmapped"] = 100.0 - rate
    return out


def _parse_picard_metrics_block(path):
    """Return (metrics dict, histogram DataFrame or None) from a Picard metrics file."""
    metrics, hist_rows, hist_header = {}, [], None
    state = None
    with open(path, errors="ignore") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("## METRICS CLASS"):
                state = "metrics_header"
                continue
            if line.startswith("## HISTOGRAM"):
                state = "hist_header"
                continue
            if not line.strip():
                if state in ("metrics_values", "hist"):
                    state = None if state == "metrics_values" else state
                continue
            if line.startswith("#"):
                continue

            fields = line.split("\t")
            if state == "metrics_header":
                hdr = fields
                state = ("metrics_values", hdr)
                continue
            if isinstance(state, tuple) and state[0] == "metrics_values":
                hdr = state[1]
                for k, v in zip(hdr, fields):
                    metrics[k] = v
                state = None
                continue
            if state == "hist_header":
                hist_header = fields
                state = "hist"
                continue
            if state == "hist":
                hist_rows.append(fields)

    hist = None
    if hist_header and hist_rows:
        width = len(hist_header)
        hist = pd.DataFrame(
            [r[:width] + [None] * (width - len(r)) for r in hist_rows],
            columns=hist_header,
        )
    return metrics, hist


def parse_picard_rnaseq(path):
    """Picard CollectRnaSeqMetrics -> percentages plus the gene-body coverage curve."""
    raw, hist = _parse_picard_metrics_block(path)
    g = lambda k: to_float(raw.get(k))
    out = {
        "pf_bases": g("PF_BASES"),
        "pf_aligned_bases": g("PF_ALIGNED_BASES"),
        "pct_coding": (g("PCT_CODING_BASES") or 0) * 100 if g("PCT_CODING_BASES") is not None else None,
        "pct_utr": (g("PCT_UTR_BASES") or 0) * 100 if g("PCT_UTR_BASES") is not None else None,
        "pct_intronic": (g("PCT_INTRONIC_BASES") or 0) * 100 if g("PCT_INTRONIC_BASES") is not None else None,
        "pct_intergenic": (g("PCT_INTERGENIC_BASES") or 0) * 100 if g("PCT_INTERGENIC_BASES") is not None else None,
        "pct_ribosomal": (g("PCT_RIBOSOMAL_BASES") or 0) * 100 if g("PCT_RIBOSOMAL_BASES") is not None else None,
        "pct_mrna": (g("PCT_MRNA_BASES") or 0) * 100 if g("PCT_MRNA_BASES") is not None else None,
        "pct_usable": (g("PCT_USABLE_BASES") or 0) * 100 if g("PCT_USABLE_BASES") is not None else None,
        "pct_correct_strand": (g("PCT_CORRECT_STRAND_READS") or 0) * 100
        if g("PCT_CORRECT_STRAND_READS") is not None else None,
        "median_cv_coverage": g("MEDIAN_CV_COVERAGE"),
        "median_5prime_bias": g("MEDIAN_5PRIME_BIAS"),
        "median_3prime_bias": g("MEDIAN_3PRIME_BIAS"),
        "median_5prime_to_3prime_bias": g("MEDIAN_5PRIME_TO_3PRIME_BIAS"),
    }

    coverage = None
    if hist is not None:
        pos_col = "normalized_position" if "normalized_position" in hist.columns else hist.columns[0]
        val_col = next(
            (c for c in hist.columns if "normalized_coverage" in c), None
        )
        if val_col:
            pos = pd.to_numeric(hist[pos_col], errors="coerce")
            val = pd.to_numeric(hist[val_col], errors="coerce")
            ok = pos.notna() & val.notna()
            if ok.any():
                coverage = [[float(p), round(float(v), 4)] for p, v in zip(pos[ok], val[ok])]
    return out, coverage


def parse_picard_markdup(path):
    raw, _ = _parse_picard_metrics_block(path)
    dup = to_float(raw.get("PERCENT_DUPLICATION"))
    return {
        "pct_duplication": dup * 100 if dup is not None else None,
        "estimated_library_size": to_float(raw.get("ESTIMATED_LIBRARY_SIZE")),
        "read_pairs_examined": to_float(raw.get("READ_PAIRS_EXAMINED")),
        "unpaired_reads_examined": to_float(raw.get("UNPAIRED_READS_EXAMINED")),
    }


def parse_featurecounts_summary(path):
    """featureCounts .summary -> {sample: {assigned, pct_assigned, ...}}."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.columns = [clean_sample_name(c) for c in df.columns]
    per_sample = {}
    totals = df.sum(axis=0)
    for s in df.columns:
        assigned = to_float(df.loc["Assigned", s]) if "Assigned" in df.index else None
        total = to_float(totals[s])
        rec = {
            "assigned_reads": assigned,
            "pct_assigned": safe_div(assigned, total) * 100 if safe_div(assigned, total) is not None else None,
            "fc_total": total,
        }
        # Delivered result folders often keep only featurecounts/ and plots/, with the
        # STAR logs left behind on the cluster. The summary still totals every record in
        # the BAM, so read depth and unique alignments can be recovered from it. Kept in
        # separate keys from the STAR-derived ones so the two are never conflated.
        unmapped = to_float(df.loc["Unassigned_Unmapped", s]) if "Unassigned_Unmapped" in df.index else None
        multi = to_float(df.loc["Unassigned_MultiMapping", s]) if "Unassigned_MultiMapping" in df.index else None
        if total is not None and unmapped is not None and multi is not None:
            rec["fc_unique_aligned"] = max(0.0, total - unmapped - multi)
        for status, key in (
            ("Unassigned_NoFeatures", "unassigned_nofeatures"),
            ("Unassigned_Ambiguity", "unassigned_ambiguity"),
            ("Unassigned_MultiMapping", "unassigned_multimapping"),
        ):
            if status in df.index:
                v = to_float(df.loc[status, s])
                rec[key] = v
                rec["pct_" + key] = safe_div(v, total) * 100 if safe_div(v, total) is not None else None
        per_sample[s] = rec
    return per_sample


def parse_cutadapt_report(path):
    """Cutadapt stdout -> input/output read counts and the fraction trimmed away."""
    txt = open(path, errors="ignore").read()
    out = {}

    def grab(pattern):
        m = re.search(pattern, txt)
        return to_float(m.group(1)) if m else None

    out["raw_reads"] = grab(r"Total reads processed:\s+([\d,]+)") or grab(
        r"Total read pairs processed:\s+([\d,]+)"
    )
    out["trimmed_reads"] = grab(r"Reads written \(passing filters\):\s+([\d,]+)") or grab(
        r"Pairs written \(passing filters\):\s+([\d,]+)"
    )
    out["reads_with_adapter"] = grab(r"Reads with adapters:\s+([\d,]+)")
    out["bp_in"] = grab(r"Total basepairs processed:\s+([\d,]+)")
    out["bp_out"] = grab(r"Total written \(filtered\):\s+([\d,]+)")

    if out.get("raw_reads") and out.get("trimmed_reads"):
        out["pct_reads_retained"] = 100.0 * out["trimmed_reads"] / out["raw_reads"]
    if out.get("bp_in") and out.get("bp_out"):
        out["pct_bp_surviving"] = 100.0 * out["bp_out"] / out["bp_in"]
    return out


def parse_ribodetector_tsv(path):
    """Ribodetector custom-content TSV (comment-prefixed MultiQC header)."""
    per_sample = {}
    with open(path, errors="ignore") as fh:
        rows = [l for l in fh if not l.startswith("#")]
    if not rows:
        return per_sample
    reader = csv.DictReader(rows, delimiter="\t")
    for r in reader:
        s = clean_sample_name(r.get("sample", ""))
        if not s:
            continue
        per_sample[s] = {
            "ribo_total_reads": to_float(r.get("total_reads")),
            "ribo_rrna_reads": to_float(r.get("rrna_reads")),
            "pct_rrna_ribodetector": to_float(r.get("pct_rrna")),
        }
    return per_sample


def parse_flagstat(path):
    txt = open(path, errors="ignore").read()
    out = {}
    m = re.search(r"^(\d+)\s+\+\s+\d+\s+in total", txt, re.M)
    if m:
        out["flagstat_total"] = float(m.group(1))
    m = re.search(r"^(\d+)\s+\+\s+\d+\s+mapped\s+\(", txt, re.M)
    if m:
        out["flagstat_mapped"] = float(m.group(1))
    if out.get("flagstat_total"):
        out["flagstat_pct_mapped"] = 100.0 * out.get("flagstat_mapped", 0) / out["flagstat_total"]
    return out


def parse_multiqc_fastqc(run_dir, extra_dirs=()):
    """
    FastQC metrics, when a raw-read MultiQC lives alongside the run.

    FastQC is not part of rule all - it is run by Utilities/run_fastqc.sh - so this
    is opportunistic: we look for any multiqc_fastqc.txt under the usual spots and
    collapse R1/R2 rows down to one row per sample.
    """
    candidates = []
    for d in list(extra_dirs) + [run_dir, os.path.join(run_dir, ".."), os.path.dirname(run_dir.rstrip("/"))]:
        if not d:
            continue
        candidates += glob.glob(os.path.join(d, "**", "multiqc_fastqc.txt"), recursive=True)
    candidates = [c for c in dict.fromkeys(candidates) if os.path.exists(c)]
    if not candidates:
        return {}, None

    path = candidates[0]
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as e:  # malformed or truncated file - not worth failing the report
        warn(f"could not read {path}: {e}")
        return {}, None

    if "Sample" not in df.columns:
        return {}, None

    df["_sample"] = [
        re.sub(r"[._-](R?[12])(_001)?$", "", clean_sample_name(s)) for s in df["Sample"]
    ]
    per_sample = {}
    for s, grp in df.groupby("_sample"):
        rec = {}
        for col, key in (
            ("%GC", "pct_gc"),
            ("total_deduplicated_percentage", "fastqc_pct_unique"),
            ("Total Sequences", "fastqc_total_sequences"),
            ("avg_sequence_length", "fastqc_avg_length"),
            ("percent_duplicates", "fastqc_pct_duplicates"),
        ):
            if col in grp.columns:
                vals = [to_float(v) for v in grp[col] if to_float(v) is not None]
                if vals:
                    rec[key] = float(np.mean(vals))
        if rec:
            per_sample[s] = rec
    return per_sample, path


def parse_multiqc_general_stats(run_dir):
    """
    MultiQC's general stats table, used only to fill gaps.

    Everything here is also derivable from the raw tool outputs, so this is a
    convenience for runs where an intermediate file was cleaned up.
    """
    hits = glob.glob(os.path.join(run_dir, "**", "multiqc_general_stats.txt"), recursive=True)
    if not hits:
        return {}, None
    try:
        df = pd.read_csv(hits[0], sep="\t")
    except Exception:
        return {}, None
    if "Sample" not in df.columns:
        return {}, None
    out = {}
    for _, row in df.iterrows():
        s = clean_sample_name(row["Sample"])
        rec = {k: to_float(v) for k, v in row.items() if k != "Sample"}
        out.setdefault(s, {}).update({k: v for k, v in rec.items() if v is not None})
    return out, hits[0]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COLLECTION - walk the run directory and build one record per sample
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def discover_samples(run_dir, sample_csv=None):
    """
    Sample IDs, from whichever source the run directory still has.

    Ordered most to least authoritative. The later entries matter for delivered result
    folders, which are often trimmed down to featurecounts/ and plots/ - the sample
    names are in those headers, so needing a sample sheet to read them would be a
    pointless hurdle.
    """
    if sample_csv and os.path.exists(sample_csv):
        try:
            df = pd.read_csv(sample_csv)
            if "sample_id" in df.columns:
                return [str(s) for s in df["sample_id"]]
        except Exception as e:
            warn(f"could not read sample sheet {sample_csv}: {e}")

    bams = sorted(glob.glob(os.path.join(run_dir, "alignment", "*.srt.bam")))
    if bams:
        return [clean_sample_name(b) for b in bams]

    picard = sorted(glob.glob(os.path.join(run_dir, "metrics", "picard", "*.picard.rna.metrics.txt")))
    if picard:
        return [clean_sample_name(p) for p in picard]

    fc_summary = first_existing(
        os.path.join(run_dir, "featurecounts", "featurecounts.readcounts.raw.tsv.summary"),
        *glob.glob(os.path.join(run_dir, "featurecounts", "*.summary")),
    )
    if fc_summary:
        try:
            with open(fc_summary, errors="ignore") as fh:
                cols = fh.readline().rstrip("\n").split("\t")[1:]
            if cols:
                return [clean_sample_name(c) for c in cols]
        except Exception as e:
            warn(f"could not read {fc_summary}: {e}")

    counts = first_existing(
        os.path.join(run_dir, "featurecounts", "featurecounts.readcounts.ann.tsv"),
        os.path.join(run_dir, "featurecounts", "featurecounts.readcounts.tsv"),
    )
    if counts:
        try:
            with open(counts, errors="ignore") as fh:
                first = fh.readline()
                if first.startswith("#"):
                    first = fh.readline()
            cols = [c for c in first.rstrip("\n").split("\t") if c not in FC_META_COLS]
            if cols:
                return [clean_sample_name(c) for c in cols]
        except Exception as e:
            warn(f"could not read {counts}: {e}")

    for fs in ("top", "all"):
        pcs = os.path.join(run_dir, "plots", f"PCA_{fs}_PCs.csv")
        if os.path.exists(pcs):
            try:
                return [clean_sample_name(x) for x in pd.read_csv(pcs)["sample"]]
            except Exception:
                pass

    return []


def collect_metrics(run_dir, samples, fastqc_dirs=()):
    """Merge every per-sample source into {sample: {metric: value}} plus coverage curves."""
    metrics = {s: {"sample": s} for s in samples}
    coverage = {}
    sources = {}

    for s in samples:
        rec = metrics[s]

        star_log = os.path.join(run_dir, "alignment", f"{s}.Log.final.out")
        hisat_log = os.path.join(run_dir, "alignment", f"{s}.hisat.summary.txt")
        if os.path.exists(star_log):
            rec.update(parse_star_log(star_log))
            sources.setdefault("aligner", star_log)
        elif os.path.exists(hisat_log):
            rec.update(parse_hisat_summary(hisat_log))
            sources.setdefault("aligner", hisat_log)

        picard_rna = os.path.join(run_dir, "metrics", "picard", f"{s}.picard.rna.metrics.txt")
        if os.path.exists(picard_rna):
            vals, cov = parse_picard_rnaseq(picard_rna)
            rec.update(vals)
            if cov:
                coverage[s] = cov
            sources.setdefault("picard_rnaseq", picard_rna)

        mkdup = os.path.join(run_dir, "markdup", f"{s}.mkdup.log.txt")
        if os.path.exists(mkdup):
            rec.update(parse_picard_markdup(mkdup))
            sources.setdefault("markdup", mkdup)

        trim = os.path.join(run_dir, "trimming", f"{s}.cutadapt.report")
        if os.path.exists(trim):
            rec.update(parse_cutadapt_report(trim))
            sources.setdefault("cutadapt", trim)

        flagstat = os.path.join(run_dir, "alignment", "stats", f"{s}.srt.bam.flagstat")
        if os.path.exists(flagstat):
            rec.update(parse_flagstat(flagstat))
            sources.setdefault("flagstat", flagstat)

    fc_summary = first_existing(
        os.path.join(run_dir, "featurecounts", "featurecounts.readcounts.raw.tsv.summary"),
        *glob.glob(os.path.join(run_dir, "featurecounts", "*.summary")),
    )
    if fc_summary:
        sources["featurecounts"] = fc_summary
        for s, rec in parse_featurecounts_summary(fc_summary).items():
            if s in metrics:
                metrics[s].update(rec)
    else:
        warn("no featureCounts .summary found - assigned-read metrics will be blank")

    ribo = os.path.join(run_dir, "ribodetector", "rrna_norrna_pct_mqc.tsv")
    if os.path.exists(ribo):
        sources["ribodetector"] = ribo
        for s, rec in parse_ribodetector_tsv(ribo).items():
            if s in metrics:
                metrics[s].update(rec)

    fastqc, fq_path = parse_multiqc_fastqc(run_dir, fastqc_dirs)
    if fq_path:
        sources["fastqc"] = fq_path
    for s, rec in fastqc.items():
        if s in metrics:
            metrics[s].update(rec)

    gs, gs_path = parse_multiqc_general_stats(run_dir)
    if gs_path:
        sources["multiqc_general_stats"] = gs_path
    for s, rec in gs.items():
        if s not in metrics:
            continue
        # General stats only fills gaps; parsed tool output always wins.
        for k, v in rec.items():
            metrics[s].setdefault(k, v)

    return metrics, coverage, sources


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COUNT MATRIX, PCA, CORRELATION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FC_META_COLS = ["Ensembl ID", "Gene Name", "Geneid", "Chr", "Start", "End", "Strand", "Length"]


def read_counts_table(path):
    """Read an (annotated or plain) featureCounts matrix -> (meta DataFrame, numeric DataFrame)."""
    with open(path, errors="ignore") as fh:
        skiprows = 1 if fh.readline().startswith("#") else 0
    df = pd.read_csv(path, sep="\t", skiprows=skiprows, low_memory=False)

    meta_cols = [c for c in FC_META_COLS if c in df.columns]
    sample_cols = [c for c in df.columns if c not in meta_cols]

    meta = df[meta_cols].copy()
    mat = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    mat.columns = [clean_sample_name(c) for c in mat.columns]

    idx_col = "Ensembl ID" if "Ensembl ID" in meta.columns else ("Geneid" if "Geneid" in meta.columns else None)
    if idx_col:
        meta.index = df[idx_col].astype(str)
        mat.index = df[idx_col].astype(str)
    return meta, mat


def median_of_ratios(mat):
    """DESeq2-style size-factor normalization, matching scripts/pca_plotting.py."""
    m = mat.astype(float)
    nonzero = m[(m != 0).all(axis=1)]
    if nonzero.shape[0] == 0:
        return m
    geo = np.exp(np.log(nonzero).mean(axis=1))
    size_factors = np.median(nonzero.div(geo, axis=0), axis=0)
    size_factors = np.where(size_factors > 0, size_factors, 1.0)
    return m.div(size_factors, axis=1)


def select_hvgs(logmat, requested=None):
    """
    Pick highly variable genes.

    Mirrors pca_plotting.py: fit an exponential decay to the ranked variance curve
    and cut at the plateau, so the dashboard PCA agrees with the PNGs the pipeline
    already ships. Falls back to a fixed top-N when scipy is unavailable or the fit
    does not converge.
    """
    gene_var = logmat.var(axis=1).sort_values(ascending=False)
    if requested:
        n = min(int(requested), len(gene_var))
        return gene_var.index[:n], gene_var, n, "user-specified"

    y = gene_var.values
    if len(y) < 10:
        return gene_var.index, gene_var, len(y), "all genes (too few to threshold)"

    plateau, how = None, "top 500 (fallback)"
    try:
        from scipy.optimize import curve_fit

        def decay(x, a, b, c):
            return a * np.exp(-b * x) + c

        params, _ = curve_fit(
            decay, np.arange(len(y)), y, p0=[y.max() - y.min(), 0.001, y.min()], maxfev=10000
        )
        c = params[2]
        tol = 0.01 * (y.max() - y.min())
        hits = np.where(y <= c + tol)[0]
        if len(hits):
            plateau, how = int(hits[0]) + 1, "variance-plateau (auto)"
    except Exception:
        pass

    if not plateau or plateau < 10:
        plateau = min(500, len(y))
        how = "top 500 (fallback)"
    return gene_var.index[:plateau], gene_var, plateau, how


def compute_pca(mat, n_comp=10, hvg=None):
    """
    Fallback PCA, for runs with no precomputed pca_plotting.py output.

    Deliberately mirrors that script step for step - including keeping all-zero genes,
    which only ever land below the HVG cut - so a recomputed PCA lands on the same
    numbers the pipeline would have produced. Returns the JSON payload plus the gene
    index that was used, so the distance panel can be built from the same subset.
    """
    if mat.shape[1] < 2 or mat.shape[0] < 2:
        return None, None

    norm = np.log2(median_of_ratios(mat) + 1)
    genes, gene_var, n_hvg, how = select_hvgs(norm, hvg)
    top = norm.loc[genes]

    X = top.T.values.astype(float)
    X = X - X.mean(axis=0)
    sd = X.std(axis=0, ddof=0)
    X = X / np.where(sd > 0, sd, 1.0)

    n_comp = int(min(n_comp, X.shape[0], X.shape[1]))
    # SVD rather than sklearn: keeps the dependency surface to numpy + pandas.
    U, S, _ = np.linalg.svd(X, full_matrices=False)
    scores = U[:, :n_comp] * S[:n_comp]
    var_ratio = (S ** 2) / np.sum(S ** 2) * 100.0

    # Variance-vs-rank curve, thinned for the browser (log-spaced, ~400 points).
    n_genes_total = len(gene_var)
    if n_genes_total > 400:
        idx = np.unique(np.round(np.logspace(0, np.log10(n_genes_total), 400)).astype(int) - 1)
        idx = idx[(idx >= 0) & (idx < n_genes_total)]
    else:
        idx = np.arange(n_genes_total)
    var_curve = [[int(i + 1), round(float(gene_var.values[i]), 5)] for i in idx]

    return {
        "samples": list(top.columns),
        "scores": [[round(float(v), 4) for v in row] for row in scores],
        "var_explained": [round(float(v), 3) for v in var_ratio[:n_comp]],
        "n_comp": n_comp,
        "n_hvg": int(n_hvg),
        "n_genes_total": int(n_genes_total),
        "hvg_method": how,
        "var_curve": var_curve,
        "groups": None,
        "precomputed": False,
    }, top


def _thin_curve(ranks, values, n=400):
    """Log-space thin a rank/value curve so 60k genes stay a sane payload."""
    total = len(values)
    if total > n:
        idx = np.unique(np.round(np.logspace(0, np.log10(total), n)).astype(int) - 1)
        idx = idx[(idx >= 0) & (idx < total)]
    else:
        idx = np.arange(total)
    return [[int(ranks[i]), round(float(values[i]), 5)] for i in idx]


def read_pca_outputs(plots_dir, feature_set="top"):
    """
    Load the PCA that scripts/pca_plotting.py already computed.

    The pipeline's PCA is the authoritative one - it is what the delivered PNGs show -
    so the dashboard reads its scores rather than recomputing and risking a figure that
    disagrees with the PNG next to it. `top` is the high-variance (HVG) feature set;
    `all` is every gene. Returns None when the CSVs are absent, which is the signal to
    fall back to computing in-process.
    """
    prefix = f"PCA_{feature_set}"
    pcs_path = os.path.join(plots_dir, f"{prefix}_PCs.csv")
    var_path = os.path.join(plots_dir, f"{prefix}_variance_explained.csv")
    if not (os.path.exists(pcs_path) and os.path.exists(var_path)):
        return None

    try:
        pcs = pd.read_csv(pcs_path)
        var = pd.read_csv(var_path)
    except Exception as e:
        warn(f"could not read precomputed PCA ({e}); recomputing")
        return None

    if "sample" not in pcs.columns:
        warn(f"{pcs_path} has no 'sample' column; recomputing")
        return None

    pc_cols = [c for c in pcs.columns if re.fullmatch(r"PC\d+", str(c))]
    pc_cols.sort(key=lambda c: int(c[2:]))
    if not pc_cols:
        warn(f"{pcs_path} has no PC columns; recomputing")
        return None

    var_col = next((c for c in ("Percent_Variance", "Variance", "percent_variance") if c in var.columns), None)
    if var_col is None:
        warn(f"{var_path} has no variance column; recomputing")
        return None

    samples = [clean_sample_name(x) for x in pcs["sample"]]
    scores = [[round(float(v), 4) for v in row] for row in pcs[pc_cols].to_numpy()]
    var_explained = [round(float(v), 3) for v in var[var_col].to_numpy()[: len(pc_cols)]]

    # pca_plotting.py writes `group` only when the sample sheet had that column.
    groups = None
    if "group" in pcs.columns:
        groups = [None if pd.isna(g) else str(g) for g in pcs["group"]]
        if not any(groups):
            groups = None

    # HVG counts come from the script's own log so the caption matches its cut exactly.
    n_hvg = n_total = None
    log_path = os.path.join(plots_dir, "pca_hvg_log.txt")
    if os.path.exists(log_path):
        txt = open(log_path, errors="ignore").read()
        m = re.search(r"Total genes:\s*(\d+)", txt)
        if m:
            n_total = int(m.group(1))
        m = re.search(r"Selected genes \(plateau\):\s*(\d+)", txt)
        if m:
            n_hvg = int(m.group(1))

    var_curve = []
    gv_path = os.path.join(plots_dir, "PCA_gene_variance.csv")
    if os.path.exists(gv_path):
        try:
            gv = pd.read_csv(gv_path)
            if {"Rank", "Variance"}.issubset(gv.columns):
                var_curve = _thin_curve(gv["Rank"].to_numpy(), gv["Variance"].to_numpy())
                n_total = n_total or int(len(gv))
        except Exception as e:
            warn(f"could not read {gv_path}: {e}")

    if feature_set == "all":
        n_hvg = n_total
    label = "high-variance gene set" if feature_set == "top" else "all genes"

    return {
        "samples": samples,
        "scores": scores,
        "var_explained": var_explained,
        "n_comp": len(pc_cols),
        "n_hvg": int(n_hvg or len(scores[0] if scores else [])),
        "n_genes_total": int(n_total or 0),
        "hvg_method": f"{label}, from pca_plotting.py",
        "var_curve": var_curve,
        "groups": groups,
        "precomputed": True,
    }


def read_pca_normalized_counts(plots_dir, feature_set="top"):
    """The log2 median-of-ratios matrix pca_plotting.py fed to the PCA, if it wrote one."""
    path = os.path.join(plots_dir, f"PCA_{feature_set}_normalized_counts.tsv")
    if not os.path.exists(path):
        return None, None
    try:
        mat = pd.read_csv(path, sep="\t", index_col=0)
    except Exception as e:
        warn(f"could not read {path}: {e}")
        return None, None
    mat.columns = [clean_sample_name(c) for c in mat.columns]
    return mat.apply(pd.to_numeric, errors="coerce").fillna(0.0), path


def compute_dissimilarity(norm_mat, method="euclidean"):
    """
    Sample-to-sample dissimilarity from the normalised matrix that fed the PCA.

    Using the PCA's own input rather than re-deriving from raw counts keeps the two
    panels telling the same story: samples that sit together in the PCA are the ones
    that read as close here. `euclidean` is the DESeq2 sample-distance convention;
    the correlation-based options are offered as 1 - r for the same orientation
    (0 = identical, larger = more different).
    """
    if norm_mat is None or norm_mat.shape[1] < 2 or norm_mat.shape[0] < 2:
        return None

    X = norm_mat.astype(float)
    samples = list(X.columns)

    if method == "euclidean":
        V = X.to_numpy().T
        # (a-b)^2 = a.a + b.b - 2a.b, clipped because round-off can go slightly negative
        sq = np.einsum("ij,ij->i", V, V)
        d2 = np.maximum(sq[:, None] + sq[None, :] - 2.0 * (V @ V.T), 0.0)
        D = np.sqrt(d2)
        np.fill_diagonal(D, 0.0)
        label = "Euclidean distance"
    else:
        corr = X.corr(method=method)
        D = (1.0 - corr).to_numpy()
        np.fill_diagonal(D, 0.0)
        label = f"1 − {method} correlation"

    D = pd.DataFrame(D, index=samples, columns=samples)

    # Average-linkage leaf order, so blocks of similar samples read as blocks.
    order = samples
    try:
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import squareform

        sym = (D.to_numpy() + D.to_numpy().T) / 2.0
        np.fill_diagonal(sym, 0.0)
        order = [samples[i] for i in leaves_list(linkage(squareform(sym, checks=False), "average"))]
    except Exception:
        pass

    D = D.loc[order, order]
    return {
        "samples": list(D.columns),
        "matrix": [[round(float(v), 4) for v in row] for row in D.values],
        "method": method,
        "label": label,
        "n_genes": int(norm_mat.shape[0]),
    }


MITO_CHR = {"M", "MT", "MITO", "CHRM", "CHRMT"}


def is_mito(chrom, gene_name):
    """A gene is mitochondrial if its chromosome is chrM/MT, or its name is MT-/mt-."""
    if chrom:
        first = str(chrom).split(";")[0].strip()
        if first.upper().replace("CHR", "") in {"M", "MT", "MITO"} or first.upper() in MITO_CHR:
            return True
    if gene_name:
        g = str(gene_name)
        if g.upper().startswith("MT-") or g.startswith("mt-"):
            return True
        if re.match(r"^(MT|mt)-?(CO|ND|ATP|CYB|RNR|TF|TV)", g):
            return True
    return False


def build_top_genes(tpm_meta, tpm_mat, counts_mat=None, top_n=100):
    """Top N genes by mean TPM across samples, with a mitochondrial flag on each."""
    if tpm_mat.shape[0] == 0:
        return None

    mean_tpm = tpm_mat.mean(axis=1)
    order = mean_tpm.sort_values(ascending=False).index[:top_n]
    samples = list(tpm_mat.columns)

    name_col = "Gene Name" if "Gene Name" in tpm_meta.columns else None
    chr_col = "Chr" if "Chr" in tpm_meta.columns else None

    rows = []
    for gid in order:
        gname = str(tpm_meta.loc[gid, name_col]) if name_col else gid
        if isinstance(gname, pd.Series):
            gname = str(gname.iloc[0])
        chrom = str(tpm_meta.loc[gid, chr_col]) if chr_col else ""
        if isinstance(chrom, pd.Series):
            chrom = str(chrom.iloc[0])
        chrom_short = chrom.split(";")[0].strip()

        tpms = [round(float(v), 2) for v in tpm_mat.loc[gid, samples].values.ravel()[: len(samples)]]
        row = {
            "id": str(gid),
            "gene": gname,
            "chr": chrom_short,
            "mito": bool(is_mito(chrom, gname)),
            "mean_tpm": round(float(mean_tpm.loc[gid]), 2),
            "tpm": tpms,
        }
        if counts_mat is not None and gid in counts_mat.index:
            cnt = counts_mat.loc[gid, [s for s in samples if s in counts_mat.columns]]
            row["counts"] = [int(v) for v in np.ravel(cnt.values)]
            row["mean_count"] = round(float(np.mean(np.ravel(cnt.values))), 1)
        rows.append(row)

    mito_total = float(tpm_mat.loc[[g for g in tpm_mat.index if is_mito(
        str(tpm_meta.loc[g, chr_col]) if chr_col else "",
        str(tpm_meta.loc[g, name_col]) if name_col else g,
    )]].sum(axis=0).mean()) if (chr_col or name_col) else None

    return {
        "samples": samples,
        "rows": rows,
        "has_counts": counts_mat is not None,
        "n_mito_in_top": sum(1 for r in rows if r["mito"]),
        "mito_pct_tpm": round(mito_total / 10000.0, 2) if mito_total else None,
    }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# METRIC DEFINITIONS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# `group` drives the table's column sections, `fmt` the cell rendering, and
# `chart` marks the metrics that get their own bar panel in the QC section.

METRIC_DEFS = [
    # --- Library / trimming
    dict(key="raw_reads", label="Raw reads", group="Library", fmt="reads",
         desc="Reads (or pairs) input to Cutadapt", chart=True),
    dict(key="trimmed_reads", label="Reads post-trim", group="Library", fmt="reads",
         desc="Reads (or pairs) written by Cutadapt, i.e. what the aligner receives"),
    dict(key="pct_reads_retained", label="% Reads retained", group="Library", fmt="pct",
         desc="Reads post-trim as a fraction of raw reads"),
    dict(key="pct_gc", label="% GC", group="Library", fmt="pct1",
         desc="Mean GC content from FastQC (raw reads)"),
    # --- Alignment
    dict(key="total_reads", label="Aligner input", group="Alignment", fmt="reads",
         desc="Reads presented to the aligner"),
    dict(key="pct_unique", label="% uniquely mapped", group="Alignment", fmt="pct",
         desc="Uniquely mapped reads as a fraction of input", chart=True),
    dict(key="uniquely_mapped", label="Uniquely mapped", group="Alignment", fmt="reads",
         desc="Absolute uniquely mapped read count"),
    dict(key="pct_multimapped", label="% multi-mapped", group="Alignment", fmt="pct",
         desc="Reads mapping to more than one locus"),
    dict(key="pct_unmapped", label="% unmapped", group="Alignment", fmt="pct",
         desc="Reads the aligner could not place (incl. too-short)"),
    dict(key="pct_duplication", label="% duplicates", group="Alignment", fmt="pct",
         desc="Picard MarkDuplicates PERCENT_DUPLICATION", chart=True),
    # --- RNA composition
    dict(key="pct_mrna", label="% mRNA bases", group="RNA composition", fmt="pct",
         desc="Picard PCT_MRNA_BASES (coding + UTR)", chart=True),
    dict(key="pct_coding", label="% coding", group="RNA composition", fmt="pct",
         desc="Picard PCT_CODING_BASES"),
    dict(key="pct_utr", label="% UTR", group="RNA composition", fmt="pct",
         desc="Picard PCT_UTR_BASES"),
    dict(key="pct_intronic", label="% intronic", group="RNA composition", fmt="pct",
         desc="Picard PCT_INTRONIC_BASES", chart=True),
    dict(key="pct_intergenic", label="% intergenic", group="RNA composition", fmt="pct",
         desc="Picard PCT_INTERGENIC_BASES", chart=True),
    dict(key="pct_ribosomal", label="% rRNA", group="RNA composition", fmt="pct",
         desc="Picard PCT_RIBOSOMAL_BASES", chart=True),
    dict(key="pct_rrna_ribodetector", label="% rRNA (RiboDetector)", group="RNA composition", fmt="pct",
         desc="Reads classified as rRNA before filtering"),
    dict(key="pct_correct_strand", label="% correct strand", group="RNA composition", fmt="pct",
         desc="Picard PCT_CORRECT_STRAND_READS", chart=True),
    dict(key="median_cv_coverage", label="Median CV coverage", group="RNA composition", fmt="num2",
         desc="Picard MEDIAN_CV_COVERAGE - evenness of transcript coverage"),
    dict(key="median_5prime_to_3prime_bias", label="5'/3' bias", group="RNA composition", fmt="num2",
         desc="Picard MEDIAN_5PRIME_TO_3PRIME_BIAS"),
    # --- Quantification
    dict(key="assigned_reads", label="Assigned reads", group="Quantification", fmt="reads",
         desc="featureCounts reads assigned to a gene", chart=True),
    dict(key="pct_assigned", label="% assigned", group="Quantification", fmt="pct",
         desc="Assigned as a fraction of all featureCounts records", chart=True),
    dict(key="pct_unassigned_nofeatures", label="% no feature", group="Quantification", fmt="pct",
         desc="Reads overlapping no annotated gene"),
    dict(key="pct_unassigned_ambiguity", label="% ambiguous", group="Quantification", fmt="pct",
         desc="Reads overlapping more than one gene"),
]


def build_summary_cards(metrics, samples, top_genes, pca):
    """The headline stat tiles - deliberately descriptive, no pass/fail verdicts."""

    def col(key):
        return [metrics[s].get(key) for s in samples if metrics[s].get(key) is not None]

    def fmt_reads(v):
        if v is None:
            return "-"
        return f"{v / 1e6:.1f}M" if v >= 1e6 else f"{v:,.0f}"

    def reads_card(label, values, note=None):
        if not values:
            return None
        return dict(
            label=label, value=fmt_reads(float(np.mean(values))),
            sub=note or f"range {fmt_reads(min(values))} - {fmt_reads(max(values))}",
        )

    def pct_card(label, values, dp=1, note=None):
        if not values:
            return None
        return dict(
            label=label, value=f"{np.median(values):.{dp}f}%",
            sub=note or f"range {min(values):.{dp}f} - {max(values):.{dp}f}%",
        )

    cards = [dict(label="Samples", value=str(len(samples)), sub="in this run")]

    # Read counts as a funnel - obtained, uniquely aligned, assigned to a gene - so the
    # drop-off between stages is readable at a glance. Means, since these are depths.
    # Preference order is Cutadapt input, then aligner input, then the featureCounts
    # summary; the last is flagged in the sub-line because it counts BAM records rather
    # than input reads (fragments, for a paired-end run).
    from_fc = "from featureCounts summary"
    obtained, note = col("raw_reads") or col("total_reads"), None
    if not obtained:
        obtained, note = col("fc_total"), from_fc
    cards.append(reads_card("Mean reads obtained", obtained, note))

    unique_n, note = col("uniquely_mapped"), None
    if not unique_n:
        unique_n, note = col("fc_unique_aligned"), from_fc
    cards.append(reads_card("Mean uniquely aligned", unique_n, note))

    cards.append(reads_card("Mean assigned reads", col("assigned_reads")))

    # Rates stay medians: one collapsed library should not drag the headline figure.
    cards.append(pct_card("Median unique alignment", col("pct_unique")))
    cards.append(pct_card("Median mRNA bases", col("pct_mrna")))
    cards.append(pct_card("Median correct strand", col("pct_correct_strand")))
    rrna = col("pct_ribosomal")
    cards.append(pct_card("Median rRNA bases", rrna, dp=2,
                          note=f"max {max(rrna):.2f}%" if rrna else None))
    cards.append(pct_card("Median duplication", col("pct_duplication")))

    if pca:
        cards.append(dict(
            label="PC1 variance", value=f"{pca['var_explained'][0]:.1f}%",
            sub=f"{pca['n_hvg']:,} variable genes",
        ))

    return [c for c in cards if c]


def _load_config_file(path):
    """Top-level keys of one config file, YAML if available and a scalar reader if not."""
    try:
        import yaml  # noqa

        with open(path) as fh:
            return yaml.safe_load(fh) or {}
    except ImportError:
        pass
    except Exception as e:
        warn(f"could not parse {path}: {e}")
        return {}

    # PyYAML is not guaranteed in every container. This fallback reads top-level
    # scalars only, and duplicate keys resolve last-wins exactly as Snakemake sees them.
    cfg = {}
    with open(path, errors="ignore") as fh:
        for line in fh:
            m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*):\s*(.*?)\s*(?:#.*)?$", line)
            if m and m.group(2) != "":
                cfg[m.group(1)] = m.group(2).strip().strip("'\"")
    return cfg


def read_run_config(run_dir, explicit=None):
    """
    Pull organism/aligner/layout context out of the pipeline config.

    Several paths may be given and are merged later-wins, mirroring how Snakemake
    layers a `--configfile` override on top of the `configfile:` directive - passing
    only the override would blank out any key it does not restate.
    """
    paths = list(explicit or [])
    if not paths:
        auto = first_existing(
            os.path.join(run_dir, "config.yaml"), os.path.join(run_dir, "config.yml")
        )
        if auto:
            paths = [auto]

    cfg, used = {}, []
    for p in paths:
        if not os.path.exists(p):
            warn(f"config not found, skipping: {p}")
            continue
        cfg.update(_load_config_file(p))
        used.append(p)
    return cfg, (used[-1] if used else None)


# Gene ID prefixes are the one organism statement that comes from the data itself
# rather than from a config file that may not be the one this run used.
GENE_ID_ORGANISM = [
    ("ENSMUSG", "Mouse"),
    ("ENSRNOG", "Rat"),
    ("ENSDARG", "Zebrafish"),
    ("ENSG", "Human"),
    ("WBGene", "C. elegans"),
    ("FBgn", "Drosophila"),
]


def organism_from_gene_ids(index):
    """Organism implied by the gene IDs in the count matrix, or None if unrecognised."""
    if index is None or len(index) == 0:
        return None
    sample = [str(g) for g in list(index[:200])]
    for prefix, name in GENE_ID_ORGANISM:
        if sum(1 for g in sample if g.startswith(prefix)) > len(sample) * 0.5:
            return name
    return None


def is_test_reference(cfg):
    """True when the config points at the repo's CI subset rather than a real genome."""
    blob = " ".join(str(cfg.get(k, "")) for k in ("annotation_gtf", "reference_fa", "aligner_index"))
    return bool(re.search(r"sample_ref/|hg38_chr567", blob))


def locate_run_dir(path):
    """
    Allow pointing one level above the results folder.

    Someone handed a project directory rather than the run directory inside it is a
    common slip, and the fix is unambiguous when exactly one child looks like a run.
    """
    markers = ("featurecounts", "alignment", "metrics", "plots")
    if any(os.path.isdir(os.path.join(path, m)) for m in markers):
        return path
    try:
        children = sorted(
            os.path.join(path, d) for d in os.listdir(path)
            if os.path.isdir(os.path.join(path, d))
        )
    except OSError:
        return path
    hits = [c for c in children if any(os.path.isdir(os.path.join(c, m)) for m in markers)]
    if len(hits) == 1:
        log(f"no pipeline outputs directly in {path}; using {hits[0]}")
        return hits[0]
    return path


def guess_organism(cfg):
    blob = " ".join(str(cfg.get(k, "")) for k in ("annotation_gtf", "aligner_index", "reference_fa"))
    if re.search(r"GRCm39|vM\d+|mm39|mouse", blob, re.I):
        return "Mouse (GRCm39)"
    if re.search(r"GRCh38|hg38|gencode\.v\d+|human", blob, re.I):
        return "Human (GRCh38)"
    return None


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RENDER
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ASSET_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(ASSET_DIR)

LOGO_URL = "https://sites.dartmouth.edu/cqb/"


def embed_logo(explicit=None):
    """
    Inline the CQB mark as a data URI.

    A linked file would break the moment the report is moved or shared, which is the
    normal case for a WebShare deliverable, so the bytes travel with the page. The
    pre-scaled PNG is preferred over the full-size JPG: the source art is ~227 KB and
    would otherwise triple the size of a typical report.
    """
    if explicit and not os.path.exists(explicit):
        warn(f"--logo {explicit} not found; falling back to the bundled mark")
        explicit = None

    path = first_existing(
        explicit,
        os.path.join(REPO_DIR, "img", "cqb_logo_small.png"),
        os.path.join(REPO_DIR, "img", "cqb_logo.jpg"),
    )
    if not path:
        warn("no logo found under img/; header logo omitted")
        return ""

    mime = "image/png" if path.lower().endswith(".png") else "image/jpeg"
    with open(path, "rb") as fh:
        data = base64.b64encode(fh.read()).decode("ascii")
    return (
        f'<a class="logo-wrap" href="{LOGO_URL}" target="_blank" rel="noopener">'
        f'<img src="data:{mime};base64,{data}" '
        f'alt="Dartmouth Center for Quantitative Biology"></a>'
    )


def render_html(payload, title, subtitle, footer, logo=""):
    tpl_path = os.path.join(ASSET_DIR, "dashboard_template.html")
    js_path = os.path.join(ASSET_DIR, "dashboard.js")
    for p in (tpl_path, js_path):
        if not os.path.exists(p):
            sys.exit(f"[dashboard] ERROR: missing asset {p}")

    tpl = open(tpl_path, encoding="utf-8").read()
    js = open(js_path, encoding="utf-8").read()

    # </script> inside the JSON payload would close the host <script> tag early.
    data = json.dumps(payload, allow_nan=False, separators=(",", ":")).replace("</", "<\\/")

    out = tpl
    out = out.replace("__TITLE__", html.escape(title))
    out = out.replace("__SUBTITLE__", html.escape(subtitle))
    out = out.replace("__FOOTER__", footer)
    out = out.replace("__LOGO__", logo)
    out = out.replace("__PAYLOAD__", data)
    out = out.replace("__JS__", js)
    return out


def sanitize(obj):
    """json.dumps(allow_nan=False) is strict; scrub NaN/Inf that survived parsing."""
    if isinstance(obj, dict):
        return {k: sanitize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [sanitize(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        f = float(obj)
        return None if (np.isnan(f) or np.isinf(f)) else f
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAIN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Build a self-contained interactive HTML QC dashboard from "
                    "GDSC-RNAseq pipeline outputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("-d", "--run-dir", default=".",
                    help="Pipeline run directory (the one holding alignment/, metrics/, featurecounts/)")
    ap.add_argument("-o", "--output", default="RNAseq_Dashboard.html", help="Output HTML file")
    ap.add_argument("-c", "--config", action="append", default=None,
                    help="Pipeline config YAML. Repeatable; later files win, mirroring "
                         "Snakemake's --configfile override (default: <run-dir>/config.yaml)")
    ap.add_argument("-s", "--sample-csv", default=None, help="Sample sheet CSV (default: from config)")
    ap.add_argument("--fastqc-dir", action="append", default=[],
                    help="Extra directory to search for a raw-read MultiQC (repeatable)")
    ap.add_argument("--plots-dir", default=None,
                    help="Where pca_plotting.py wrote its outputs (default: <run-dir>/plots)")
    ap.add_argument("--feature-set", default="top", choices=["top", "all"],
                    help="Which precomputed PCA to show: 'top' = high-variance genes, "
                         "'all' = every gene")
    ap.add_argument("--recompute-pca", action="store_true",
                    help="Ignore the precomputed PCA outputs and recompute from the counts matrix")
    ap.add_argument("--logo", default=None,
                    help="Logo image to inline in the header (default: img/cqb_logo_small.png, "
                         "falling back to img/cqb_logo.jpg). Pass an empty string for none.")
    ap.add_argument("--title", default="GDSC Bulk RNA-Seq — QC Dashboard")
    ap.add_argument("--subtitle", default=None, help="Header subtitle (default: run directory name)")
    ap.add_argument("--top-n", type=int, default=100, help="Genes in the expression table")
    ap.add_argument("--pcs", type=int, default=10, help="Principal components to compute")
    ap.add_argument("--hvg", type=int, default=None,
                    help="Fixed number of variable genes for PCA (default: variance-plateau auto-detect)")
    ap.add_argument("--dist-method", default="euclidean",
                    choices=["euclidean", "spearman", "pearson"],
                    help="Sample-to-sample dissimilarity: euclidean distance (DESeq2 "
                         "convention) or 1 - correlation")
    ap.add_argument("--version", action="version", version=f"build_dashboard.py {VERSION}")
    args = ap.parse_args(argv)

    run_dir = os.path.abspath(args.run_dir)
    if not os.path.isdir(run_dir):
        sys.exit(f"[dashboard] ERROR: run directory not found: {run_dir}")
    run_dir = locate_run_dir(run_dir)

    cfg, cfg_path = read_run_config(run_dir, args.config)
    sample_csv = args.sample_csv or (
        os.path.join(run_dir, str(cfg["sample_csv"])) if cfg.get("sample_csv") else None
    )

    samples = discover_samples(run_dir, sample_csv)
    if not samples:
        sys.exit("[dashboard] ERROR: no samples found. Run this inside a finished pipeline "
                 "directory, or pass --run-dir / --sample-csv.")
    log(f"{len(samples)} samples: {', '.join(samples[:6])}{' …' if len(samples) > 6 else ''}")

    metrics, coverage, sources = collect_metrics(run_dir, samples, args.fastqc_dir)
    log(f"parsed sources: {', '.join(sorted(sources)) or 'none'}")

    # ---- count matrices
    counts_path = first_existing(
        os.path.join(run_dir, "featurecounts", "featurecounts.readcounts.ann.tsv"),
        os.path.join(run_dir, "featurecounts", "featurecounts.readcounts.tsv"),
    )
    tpm_path = first_existing(
        os.path.join(run_dir, "featurecounts", "featurecounts.readcounts_tpm.ann.tsv"),
        os.path.join(run_dir, "featurecounts", "featurecounts.readcounts_tpm.tsv"),
    )

    pca = dist = top_genes = None
    counts_meta = counts_mat = None
    if counts_path:
        sources["counts"] = counts_path
        counts_meta, counts_mat = read_counts_table(counts_path)
        counts_mat = counts_mat[[c for c in counts_mat.columns if c in samples]]
        log(f"count matrix: {counts_mat.shape[0]:,} genes x {counts_mat.shape[1]} samples")

    # ---- PCA: the pipeline's own, so the dashboard and the delivered PNGs agree
    plots_dir = args.plots_dir or os.path.join(run_dir, "plots")
    norm_mat = norm_path = None
    if not args.recompute_pca:
        pca = read_pca_outputs(plots_dir, args.feature_set)
        if pca:
            sources["pca"] = os.path.join(plots_dir, f"PCA_{args.feature_set}_PCs.csv")
            log(f"PCA from {sources['pca']} "
                f"({pca['n_hvg']:,} of {pca['n_genes_total']:,} genes, {pca['n_comp']} PCs)")
            missing = [x for x in pca["samples"] if x not in samples]
            if missing:
                warn(f"PCA CSV names {len(missing)} sample(s) not in this run: "
                     f"{', '.join(missing[:5])}")
        norm_mat, norm_path = read_pca_normalized_counts(plots_dir, args.feature_set)
        if norm_path:
            sources["pca_normalized_counts"] = norm_path
            log(f"normalised counts from {norm_path} ({norm_mat.shape[0]:,} genes)")

    # A sample sheet with a `group` column groups the fallback PCA as well, so the two
    # code paths look the same to the reader.
    sheet_groups = {}
    if sample_csv and os.path.exists(sample_csv):
        try:
            sheet = pd.read_csv(sample_csv)
            cols = {str(c).strip().lower(): c for c in sheet.columns}
            if "sample_id" in cols and "group" in cols:
                sheet_groups = {
                    str(a).strip(): str(b).strip()
                    for a, b in zip(sheet[cols["sample_id"]], sheet[cols["group"]])
                    if not pd.isna(b) and str(b).strip()
                }
        except Exception as e:
            warn(f"could not read groups from {sample_csv}: {e}")

    if pca is None and counts_mat is not None:
        why = "--recompute-pca" if args.recompute_pca else "no precomputed PCA found"
        log(f"{why}; computing PCA from the counts matrix")
        pca, fallback_norm = compute_pca(counts_mat, n_comp=args.pcs, hvg=args.hvg)
        # Same subset the fallback PCA used, so the distance caption stays true.
        if norm_mat is None and fallback_norm is not None:
            norm_mat = fallback_norm

    if pca and not pca.get("groups") and sheet_groups:
        mapped = [sheet_groups.get(x) for x in pca["samples"]]
        if any(mapped):
            pca["groups"] = mapped

    if norm_mat is None and counts_mat is not None:
        norm_mat = np.log2(median_of_ratios(counts_mat) + 1)
    if norm_mat is not None:
        norm_mat = norm_mat[[c for c in norm_mat.columns if c in samples]]
        dist = compute_dissimilarity(norm_mat, method=args.dist_method)

    if counts_path is None and pca is None:
        warn("no featureCounts matrix and no precomputed PCA - PCA and distances omitted")

    if tpm_path:
        sources["tpm"] = tpm_path
        tpm_meta, tpm_mat = read_counts_table(tpm_path)
        tpm_mat = tpm_mat[[c for c in tpm_mat.columns if c in samples]]
        top_genes = build_top_genes(tpm_meta, tpm_mat, counts_mat, top_n=args.top_n)
        if top_genes:
            log(f"top {len(top_genes['rows'])} genes by mean TPM "
                f"({top_genes['n_mito_in_top']} mitochondrial)")
    else:
        warn("no TPM matrix found - expression table will be omitted")

    # ---- run context shown in the Overview card
    #
    # The config is only as trustworthy as its provenance: a run driven by a prebuilt
    # config still leaves the repo's default config.yaml sitting in the directory, and
    # reading that one would describe the CI test reference as though it were the run's.
    # Gene IDs come from the data, so they arbitrate.
    cfg_organism = guess_organism(cfg)
    data_organism = organism_from_gene_ids(counts_mat.index if counts_mat is not None else None)

    run_info = [("Run directory", os.path.basename(run_dir) or run_dir)]
    conflicts = []
    if data_organism:
        run_info.append(("Organism", f"{data_organism} — from gene IDs"))
        if cfg_organism and not cfg_organism.startswith(data_organism):
            conflicts.append(f"gene IDs are {data_organism.lower()}, but "
                             f"{os.path.basename(cfg_path or 'the config')} says {cfg_organism}")
    elif cfg_organism:
        run_info.append(("Organism", f"{cfg_organism} — from config"))

    if cfg and is_test_reference(cfg):
        conflicts.append(f"{os.path.basename(cfg_path or 'the config')} points at the "
                         f"pipeline's CI test reference (sample_ref/)")
    aligners = {metrics[s].get("aligner") for s in samples if metrics[s].get("aligner")}
    run_info += [
        ("Aligner", ", ".join(sorted(a for a in aligners if a)) or str(cfg.get("aligner_name", "—"))),
        ("Library layout", str(cfg.get("layout", "—"))),
        ("Annotation", os.path.basename(str(cfg.get("annotation_gtf", "—")))),
        ("Reference", os.path.basename(str(cfg.get("reference_fa", "—")))),
        ("featureCounts strand", str(cfg.get("featurecounts_strand", "—"))),
        ("Picard strand", str(cfg.get("picard_strand", "—"))),
        ("rRNA filtering", "RiboDetector" if cfg.get("remove_rRNA") in (True, "true", "True") else "off"),
        ("Samples", str(len(samples))),
        ("Generated", datetime.now().strftime("%Y-%m-%d %H:%M")),
    ]
    if cfg_path:
        run_info.append(("Config", os.path.basename(cfg_path)))
    for c in conflicts:
        warn(f"run context: {c}. Pass --config with the config this run actually used.")
        run_info.append(("⚠ Check config", c))

    payload = {
        "samples": samples,
        "metrics": metrics,
        "metric_defs": METRIC_DEFS,
        "coverage": coverage,
        "pca": pca,
        "dist": dist,
        "top_genes": top_genes,
        "cards": build_summary_cards(metrics, samples, top_genes, pca),
        "run_info": run_info,
        "sources": sources,
    }

    subtitle = args.subtitle or (os.path.basename(run_dir) or "run")
    footer = (
        "Generated by <code>scripts/build_dashboard.py</code> v" + VERSION +
        " · Dartmouth Genomic Data Science Core · "
        '<a href="mailto:gdsc@groups.dartmouth.edu">gdsc@groups.dartmouth.edu</a><br>'
        "PCA and correlation shown here are exploratory; use VST or rlog normalisation "
        "for differential expression."
    )

    out_path = os.path.abspath(args.output)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        logo = "" if args.logo == "" else embed_logo(args.logo)
        fh.write(render_html(sanitize(payload), args.title, subtitle, footer, logo))

    size = os.path.getsize(out_path) / 1024.0
    log(f"wrote {out_path} ({size:,.0f} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
