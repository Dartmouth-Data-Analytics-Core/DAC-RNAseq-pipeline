#!/usr/bin/env python3
"""QC cutoffs from MultiQC output; writes qc_metrics.tsv, qc_summary.md, qc_email_draft.md
and downstream_analysis_samples.txt (only when samples are excluded)."""
import argparse, json, re, sys
from io import StringIO
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from mqc_lz import decompress_from_base64  # noqa: E402
import qc_email_text  # noqa: E402

PIPELINE_DIR = Path(__file__).resolve().parent.parent
CUTOFFS_DIR = PIPELINE_DIR / "additional_rules" / "qc_email" / "cutoffs"
TEMPLATES_DIR = PIPELINE_DIR / "additional_rules" / "qc_email" / "templates"
ORGANISM_PATTERNS = {
    "human":      r"GRCh3[78]|\bhg(19|38)|homo[_ ]sapiens|gencode\.v\d",
    "mouse":      r"GRCm3[89]|\bmm(10|39)|mus[_ ]musculus|gencode\.vM",
    "zebrafish":  r"GRCz1[12]|danRer1[01]|danio[_ ]rerio",
    "c_elegans":  r"WBcel235|\bce11|caenorhabditis",
    "drosophila": r"BDGP6|\bdm6|drosophila",
    "yeast":      r"(?<![A-Za-z0-9])R64(?![0-9])|sacCer3|saccharomyces",
}
ZEBRAFISH_ANNOTATION = {
    "zebrafish_lawson":  r"V4\.3|lawson",
    "zebrafish_ensembl": r"danio_rerio\.|ensembl",
}
LOW_QC_ASSIGNED_FRAC = 0.8   # low QC = failed a cutoff or assigned < 0.8 x median


class GateError(SystemExit):
    """Print the message and exit with code 2."""
    def __init__(self, msg):
        print(f"[qc_email] ERROR: {msg}", file=sys.stderr)
        super().__init__(2)


def detect_organism(cfg):
    idx, gtf = str(cfg.get("aligner_index") or ""), str(cfg.get("annotation_gtf") or "")
    hits = [o for o, p in ORGANISM_PATTERNS.items() if re.search(p, f"{idx} {gtf}", re.I)]
    if len(hits) != 1:
        raise GateError(
            f"could not detect a single organism from aligner_index={idx!r} and "
            f"annotation_gtf={gtf!r} (matches: {hits or 'none'}). "
            "Set qc_organism in config.yaml (e.g. human, mouse, zebrafish_lawson).")
    if hits[0] != "zebrafish":
        return hits[0]
    for name, p in ZEBRAFISH_ANNOTATION.items():
        if re.search(p, gtf, re.I):
            return name
    raise GateError(f"zebrafish detected but the annotation source of {gtf!r} is unknown "
                    "(Lawson V4.3 or Ensembl). Set qc_organism: zebrafish_lawson | zebrafish_ensembl.")


def load_gate(cutoffs_dir, kit, organism):
    """(file, cutoffs, section, unverified keys) for this kit and organism."""
    f = Path(cutoffs_dir) / f"{kit}.yaml"
    if not f.exists():
        have = ", ".join(sorted(p.stem for p in Path(cutoffs_dir).glob("*.yaml"))) or "none"
        raise GateError(f"no cutoff file for kit {kit!r} ({f} not found). "
                        f"Available kits: {have}. Set qc_kit in config.yaml.")
    spec = yaml.safe_load(f.read_text()) or {}
    section = (spec.get("organisms") or {}).get(organism)
    if section is None:
        raise GateError(f"organism {organism!r} is not covered by {f} "
                        f"(covered: {', '.join(spec.get('organisms') or {})}). "
                        "Add a section to that file, or set qc_organism to a covered organism "
                        "and adjust its thresholds with qc_gate in config.yaml.")
    sections = spec.get("sections") or {}

    def resolve(name, seen=()):
        if name not in sections or name in seen:
            raise GateError(f"{f}: section {name!r} missing or inherits in a loop")
        s = dict(sections[name])
        parent = s.pop("inherits", None)
        return {**(resolve(parent, seen + (name,)) if parent else {}), **s}

    gate = resolve(section)
    unverified = gate.pop("unverified", None) or []
    return f, gate, section, list(unverified)


WANT = {   # our column -> MultiQC / Picard column names
    "unique_pct":   ["uniquely_mapped_percent"],
    "unique_M":     ["uniquely_mapped"],
    "multi_pct":    ["multimapped_percent"],
    "tooshort_pct": ["unmapped_tooshort_percent"],
    "other_pct":    ["unmapped_other_percent"],
    "avg_mapped_len": ["avg_mapped_read_length"],
    "assigned_M":   ["Assigned"],
    "assigned_pct": ["percent_assigned"],
    "mrna_pct":     ["PCT_MRNA_BASES"],
    "coding_pct":   ["PCT_CODING_BASES"],
    "utr_pct":      ["PCT_UTR_BASES"],
    "intronic_pct": ["PCT_INTRONIC_BASES"],
    "intergenic_pct": ["PCT_INTERGENIC_BASES"],
    "rrna_pct":     ["PCT_RIBOSOMAL_BASES"],
    "strand_pct":   ["PCT_CORRECT_STRAND_READS"],
    "bias_5to3":    ["MEDIAN_5PRIME_TO_3PRIME_BIAS"],
    "dup_pct":      ["PERCENT_DUPLICATION"],
}
FRACTION_COLS = {"mrna_pct", "coding_pct", "utr_pct", "intronic_pct", "intergenic_pct",
                 "rrna_pct", "strand_pct", "dup_pct"}  # Picard reports 0-1


# ----------------------------------------------------------------- reading
def _read_tsv(p):
    return pd.read_csv(p, sep="\t", dtype={"Sample": str}).set_index("Sample")


def _pick(df, names):
    low = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in low:
            return df[low[n.lower()]]
    return None


def _hisat2(frames, have):
    """HISAT2 has no unique % column: aligned exactly once / total."""
    if "unique_pct" in have:
        return {}
    for df in frames:
        for prefix in ("paired", "unpaired"):
            tot, one, multi = (_pick(df, [f"{prefix}_{x}"]) for x in ("total", "aligned_one", "aligned_multi"))
            if tot is not None and one is not None and pd.to_numeric(tot, errors="coerce").sum() > 0:
                tot, one = pd.to_numeric(tot, errors="coerce"), pd.to_numeric(one, errors="coerce")
                res = {"unique_pct": 100 * one / tot, "unique_M": one}
                if multi is not None:
                    res["multi_pct"] = 100 * pd.to_numeric(multi, errors="coerce") / tot
                rate = _pick(df, ["overall_alignment_rate"])
                if rate is not None:
                    res["overall_aln_pct"] = pd.to_numeric(rate, errors="coerce")
                return res
    return {}


def _picard_rna_files(folder):
    """CollectRnaSeqMetrics values per sample, from the Picard files."""
    rows = {}
    for f in sorted(Path(folder).glob("*.picard.rna.metrics.txt")):
        lines = f.read_text(errors="ignore").splitlines()
        i = next((k for k, l in enumerate(lines) if l.startswith("## METRICS CLASS")), None)
        if i is None or i + 2 >= len(lines):
            continue
        rows[f.name.replace(".picard.rna.metrics.txt", "")] = dict(
            zip(lines[i + 1].split("\t"), lines[i + 2].split("\t")))
    if not rows:
        return None
    return pd.DataFrame(rows).T.apply(pd.to_numeric, errors="coerce")


def _picard_coverage(folder):
    """Gene-body coverage per sample (positions 0-100), from the Picard files."""
    rows = {}
    for f in sorted(Path(folder).glob("*.picard.rna.metrics.txt")):
        lines = f.read_text(errors="ignore").splitlines()
        i = next((k for k, l in enumerate(lines) if l.startswith("## HISTOGRAM")), None)
        if i is None:
            continue
        vals = {}
        for l in lines[i + 2:]:
            parts = l.split("\t")
            if len(parts) < 2 or not parts[0].strip().isdigit():
                break
            vals[int(parts[0])] = pd.to_numeric(parts[1], errors="coerce")
        if vals:
            rows[f.name.replace(".picard.rna.metrics.txt", "")] = vals
    return pd.DataFrame(rows).T if rows else None


def from_multiqc_data(aln_dir, raw_dir):
    aln_dir = Path(aln_dir)
    frames = []
    for f in sorted(aln_dir.glob("multiqc_*.txt")):
        if f.name in ("multiqc_sources.txt", "multiqc_software_versions.txt",
                      "multiqc_general_stats.txt", "multiqc_citations.txt"):
            continue
        try:
            frames.append(_read_tsv(f))
        except Exception:
            pass
    cols = {}
    for key, names in WANT.items():
        for df in frames:
            s = _pick(df, names)
            if s is not None:
                s = pd.to_numeric(s, errors="coerce")
                if key in FRACTION_COLS and s.dropna().le(1.0).all():
                    s = s * 100
                cols[key] = s
                break
    cols.update(_hisat2(frames, cols))
    # RustQC runs leave Picard RNA metrics out of MultiQC: read the Picard files
    picard = _picard_rna_files(aln_dir.parent / "metrics" / "picard")
    if picard is not None:
        added = []
        for key, names in WANT.items():
            s = _pick(picard, names)
            if key not in cols and s is not None and s.notna().any():
                cols[key] = s * 100 if key in FRACTION_COLS and s.dropna().le(1.0).all() else s
                added.append(key)
        if added:
            print(f"Picard RNA metrics not in the MultiQC data; read {len(picard)} files from "
                  f"metrics/picard/ for: {', '.join(added)}")
    out = pd.concat(cols, axis=1)
    if "assigned_M" in out:
        out["assigned_M"] = out["assigned_M"] / 1e6
    if "unique_M" in out:
        out["unique_M"] = out["unique_M"] / 1e6
    out.index = out.index.astype(str)

    strand_f = next(iter(aln_dir.glob("mqc_picard_rnaseqmetrics_strand_plot*.txt")), None)
    if strand_f and ("strand_pct" not in out or out["strand_pct"].isna().all()):
        sd = _read_tsv(strand_f)
        c, i = _pick(sd, ["Correct"]), _pick(sd, ["Incorrect"])
        if c is not None and i is not None:
            out["strand_pct"] = 100 * c / (c + i)

    # gene-body coverage -> 3'/5' ratio (mean of 3' half / mean of 5' half)
    cov_f = next(iter(aln_dir.glob("mqc_picard_rna_coverage*.txt")), None)
    cov = _read_tsv(cov_f).apply(pd.to_numeric, errors="coerce") if cov_f else None
    if cov is None:
        cov = _picard_coverage(aln_dir.parent / "metrics" / "picard")
        if cov is not None:
            print(f"gene-body coverage read from {len(cov)} Picard files in metrics/picard/")
    if cov is not None:
        pos = pd.to_numeric(cov.columns, errors="coerce")
        out["cov_3to5_ratio"] = cov.loc[:, pos >= 50].mean(axis=1) / cov.loc[:, pos < 50].mean(axis=1)

    raw = {}
    if raw_dir:
        fq = Path(raw_dir) / "multiqc_fastqc.txt"
        if fq.exists():
            d = _read_tsv(fq)
            raw = pd.DataFrame({
                "raw_M": pd.to_numeric(_pick(d, ["Total Sequences"]), errors="coerce") / 1e6,
                "gc_pct": pd.to_numeric(_pick(d, ["%GC"]), errors="coerce"),
                "read_len": _pick(d, ["avg_sequence_length", "Sequence length"]),
            })
            raw.index = raw.index.astype(str)
            raw = raw.join(fastqc_curves(raw_dir), how="left")
    return out, raw, None


def _fastqc_plot(raw_dir, tsv_glob, json_id):
    """A FastQC plot as a table, from mqc_<plot>.txt or multiqc_data.json."""
    f = next(iter(sorted(Path(raw_dir).glob(tsv_glob))), None)
    if f:
        t = pd.read_csv(f, sep="\t", dtype={"Sample": str}).set_index("Sample")
        t.columns = pd.to_numeric(t.columns, errors="coerce")
        return t.apply(pd.to_numeric, errors="coerce")
    j = Path(raw_dir) / "multiqc_data.json"
    if not j.exists():
        return None
    p = (json.loads(j.read_text()).get("report_plot_data") or {}).get(json_id)
    if not p:
        return None
    ds = p["datasets"][0]
    if "samples" in p:      # bar plot
        return pd.DataFrame({s["name"]: s["data"] for s in ds}, index=p["samples"][0])
    return pd.DataFrame({s["name"]: dict((float(x), y) for x, y in s["data"]) for s in ds}).T


def fastqc_curves(raw_dir):
    """Per sample: mean Phred, main GC peak, % reads >70% GC, overrepresented %, max adapter %."""
    res = {}
    q = _fastqc_plot(raw_dir, "mqc_fastqc_per_base_sequence_quality_plot*.txt",
                     "fastqc_per_base_sequence_quality_plot")
    if q is not None:
        res["mean_phred"] = q.mean(axis=1)
    gc = _fastqc_plot(raw_dir, "mqc_fastqc_per_sequence_gc_content_plot_Percentages*.txt",
                      "fastqc_per_sequence_gc_content_plot")
    if gc is None:
        gc = _fastqc_plot(raw_dir, "mqc_fastqc_per_sequence_gc_content_plot*.txt",
                          "fastqc_per_sequence_gc_content_plot")
    if gc is not None:
        gc = gc.div(gc.sum(axis=1), axis=0) * 100
        res["gc_main_peak"] = gc.loc[:, gc.columns <= 60].idxmax(axis=1)
        res["high_gc_pct"] = gc.loc[:, gc.columns > 70].sum(axis=1)
    ov = _fastqc_plot(raw_dir, "mqc_fastqc_overrepresented_sequences_plot*.txt",
                      "fastqc_overrepresented_sequences_plot")
    if ov is not None:
        res["overrep_top_pct"] = ov.iloc[:, 0]
        res["overrep_total_pct"] = ov.sum(axis=1)
    ad = _fastqc_plot(raw_dir, "mqc_fastqc_adapter_content_plot*.txt", "fastqc_adapter_content_plot")
    if ad is not None:     # rows: "<sample> - <adapter>"
        ad.index = ad.index.astype(str).str.replace(r" - [^-]+$", "", regex=True)
        res["adapter_max_pct"] = ad.max(axis=1).groupby(level=0).max()
    out = pd.DataFrame(res)
    out.index = out.index.astype(str)
    return out


def _html_tables(path):
    html = Path(path).read_text(encoding="utf-8", errors="ignore")
    m = re.search(r'id="mqc_compressed_plotdata">([^<]+)<', html)
    plots = json.loads(decompress_from_base64(m.group(1))) if m else {}
    tabs = pd.read_html(StringIO(html))
    return html, plots, tabs


def _num(s):
    return pd.to_numeric(s.astype(str).str.replace(r"[%a-zA-Z ]", "", regex=True), errors="coerce")


def from_html(aln_html, raw_html):
    """Fallback: General Statistics table only."""
    _, plots, tabs = _html_tables(aln_html)
    g = tabs[0].set_index("Sample Name")
    out = pd.DataFrame(index=g.index.astype(str))
    m = {"% Aligned": "unique_pct", "M Aligned": "unique_M", "% Assigned": "assigned_pct",
         "M Assigned": "assigned_M", "Duplication": "dup_pct", "rRNA": "rrna_pct",
         "mRNA": "mrna_pct", "M Reads": "bam_reads_M"}
    for col, key in m.items():
        if col in g:
            out[key] = _num(g[col]).values
    for pid, p in plots.items():
        if "strand_plot" in pid:
            ds = p["datasets"][0]
            samples = p.get("samples", [[]])[0]
            c = dict(zip(samples, next(d["data"] for d in ds if d["name"] == "Correct")))
            i = dict(zip(samples, next(d["data"] for d in ds if d["name"] == "Incorrect")))
            out["strand_pct"] = [100 * c[s] / (c[s] + i[s]) if s in c else None for s in out.index]
    raw, overrep = {}, None
    if raw_html:
        _, _, rt = _html_tables(raw_html)
        r = rt[0].set_index("Sample Name")
        raw = pd.DataFrame({"raw_M": _num(r["M Seqs"]), "gc_pct": _num(r["% GC"]),
                            "dup_fastqc_pct": _num(r["% Dups"])})
        for t in rt:
            if "Overrepresented sequence" in t.columns:
                overrep = t.head(10)
    return out, raw, overrep


def collapse_mates(raw, regex):
    """Merge FastQC R1/R2 rows per sample."""
    if not isinstance(raw, pd.DataFrame) or raw.empty:
        return None
    raw = raw.copy()
    raw.index = raw.index.astype(str).str.replace(regex, "", regex=True)
    agg = {c: "mean" for c in raw.columns if c != "read_len"}
    if "raw_M" in agg:
        agg["raw_M"] = "first"   # mates have the same read count
    return raw.groupby(level=0).agg(agg)


# ----------------------------------------------------------------- gate
def resolve_auto(df, g, notes):
    """Replace "auto" cutoffs with values from this run."""
    if g.get("min_assigned_M") == "auto" and "assigned_M" in df:
        frac = g.get("auto_assigned_frac", 0.2)
        g["min_assigned_M"] = round(frac * df["assigned_M"].median(), 2)
        notes.append(f"min_assigned_M=auto -> {g['min_assigned_M']}M ({frac:g} x run median)")
    if g.get("min_mrna") == "auto" and "mrna_pct" in df:
        x = df["mrna_pct"].dropna()
        med = x.median()
        raw = med - 3 * (x - med).abs().median()
        lo, hi = g.get("auto_mrna_floor", 20), g.get("auto_mrna_cap", 40)
        g["min_mrna"] = round(min(max(raw, lo), hi), 1)
        notes.append(f"min_mrna=auto -> {g['min_mrna']}% (median-3*MAD={raw:.1f}, clamped to [{lo}, {hi}])")
        if raw < lo:
            notes.append(f"WARNING: mRNA% is low across the whole run (median {med:.1f}%); "
                         f"the auto threshold hit its floor of {lo}%")
    return g


def apply_gate(df, gate):
    g = dict(gate)
    notes, missing = [], []
    g = resolve_auto(df, g, notes)
    # rRNA 0 everywhere = not measured (empty rRNA interval list)
    if "rrna_pct" in df and df["rrna_pct"].notna().any() and df["rrna_pct"].fillna(0).eq(0).all():
        df = df.drop(columns="rrna_pct")
        missing.append("rrna_pct (0 in every sample: rRNA interval list is probably empty)")
    reasons = pd.Series([[] for _ in range(len(df))], index=df.index)

    def rule(key, col, ok, label):
        if g.get(key) is None:
            return
        if col not in df or df[col].isna().all():
            missing.append(col)
            return
        bad = ~ok(df[col], g[key]) & df[col].notna()
        for s in df.index[bad]:
            reasons[s].append(label.format(g[key]))

    rule("min_strand", "strand_pct", lambda x, t: x >= t, "strand<{}%")
    rule("min_mrna", "mrna_pct", lambda x, t: x >= t, "mRNA<{}%")
    rule("min_assigned_M", "assigned_M", lambda x, t: x >= t, "assigned<{}M")
    rule("min_unique", "unique_pct", lambda x, t: x >= t, "unique<{}%")
    rule("max_rrna", "rrna_pct", lambda x, t: x <= t, "rRNA>{}%")
    rule("max_dup", "dup_pct", lambda x, t: x <= t, "dup>{}%")
    df = df.copy()
    df["gate_include"] = reasons.map(len).eq(0)
    df["gate_reasons"] = reasons.map(lambda r: ";".join(r))
    return df, g, sorted(set(missing)), notes


# ----------------------------------------------------------------- PCA cross-check
def pca_crosscheck(df, pcs_csv, n_pcs=3):
    """Where low-QC samples sit on PC1-3 (robust z, rank, PC extremes). Flag only."""
    df = df.copy()
    low = ~df["gate_include"]
    if "assigned_M" in df:
        low |= df["assigned_M"] < LOW_QC_ASSIGNED_FRAC * df["assigned_M"].median()
    df["low_qc"] = low
    if not pcs_csv or not Path(pcs_csv).exists():
        return df, f"PCA table not found ({pcs_csv}); cross-check not done."
    p = pd.read_csv(pcs_csv)
    p = p.set_index(p.columns[0])
    p.index = p.index.astype(str)
    pcs = [c for c in p.columns if re.fullmatch(r"PC\d+", str(c), re.I)][:n_pcs]
    X = p[pcs].apply(pd.to_numeric, errors="coerce")
    dev = (X - X.median()).abs()
    z = dev / (1.4826 * dev.median())
    rank = dev.rank(ascending=False, method="min").astype(int)
    ext = {s: [] for s in X.index}
    for c in pcs:
        ext[X[c].idxmin()].append(f"{c}-")
        ext[X[c].idxmax()].append(f"{c}+")
    best = z.idxmax(axis=1)
    pc_info = pd.DataFrame({
        "pca_z": z.max(axis=1).round(1),
        "pca_z_pc": best,
        "pca_rank": [rank.at[s, best[s]] for s in X.index],
        "pca_extreme": [",".join(ext[s]) for s in X.index],
    })
    unmatched = sorted(set(df.index) - set(pc_info.index))
    df = df.join(pc_info, how="left")

    hvg = Path(pcs_csv).parent / "pca_hvg_log.txt"
    m = re.search(r"Selected genes.*?:\s*(\d+)", hvg.read_text()) if hvg.exists() else None
    n = len(X)
    lines = [f"PCA table: {pcs_csv} ({n} samples"
             + (f", {m.group(1)} genes" if m else "") + "). "
             "z = |score - median| / (1.4826 x MAD) per PC; rank 1 = farthest from the median. "
             "Descriptive only: with few samples a clean run often has some z > 3, so z alone "
             "does not make a sample an outlier, and PCA position never excludes a sample."]
    lq = df[df["low_qc"]]
    if lq.empty:
        lines.append("No low-QC samples (all pass the QC cutoffs and none has assigned < "
                     f"{LOW_QC_ASSIGNED_FRAC:g} x median).")
    else:
        lines.append(f"Low-QC samples ({len(lq)}; failed a cutoff or assigned < "
                     f"{LOW_QC_ASSIGNED_FRAC:g} x run median) and where they sit on the PCA:")
        for s, r in lq.sort_values("assigned_M").iterrows():
            why = r["gate_reasons"] or f"assigned {r['assigned_M']:.2f}M"
            if pd.isna(r.get("pca_z")):
                lines.append(f"- {s}: {why}; not in the PCA table")
            else:
                e = f", extreme of {r['pca_extreme']}" if r["pca_extreme"] else ""
                lines.append(f"- {s}: {why}; max z {r['pca_z']} on {r['pca_z_pc']} "
                             f"(rank {int(r['pca_rank'])}/{n}){e}")
        others = df[~df["low_qc"] & df["pca_extreme"].fillna("").ne("")]
        if len(others):
            lines.append("Other samples at a PC extreme (QC not low): " + "; ".join(
                f"{s} ({r['pca_extreme']}, z {r['pca_z']})" for s, r in others.iterrows()))
    if unmatched:
        lines.append(f"Samples missing from the PCA table (name mismatch?): {', '.join(unmatched)}")
    return df, "\n".join(lines)


# ----------------------------------------------------------------- summary
RANGE_COLS = {
    "raw_M": "raw reads (M)", "gc_pct": "GC %", "unique_pct": "uniquely mapped %",
    "multi_pct": "multi-mapped %", "tooshort_pct": "too short %", "assigned_M": "assigned (M)",
    "mrna_pct": "mRNA %", "intronic_pct": "intronic %", "intergenic_pct": "intergenic %",
    "rrna_pct": "rRNA %", "strand_pct": "correct strand %", "cov_3to5_ratio": "3'/5' coverage",
    "dup_pct": "duplication %",
}


def _md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(x) for x in r) + " |" for r in rows]
    return "\n".join(out)


def run_context(cfg):
    keys = ["aligner_index", "annotation_gtf", "picard_refflat", "picard_rrna_list",
            "remove_rRNA", "layout", "aligner_name", "picard_strand", "read_length", "umi"]
    return "\n".join(f"- {k}: {cfg.get(k)}" for k in keys if k in cfg) or "(config.yaml not provided)"


def build_summary(df, gate_used, missing, cfg, overrep, pca_var, source, gate_info, crosscheck):
    n, k = len(df), int(df["gate_include"].sum())
    parts = [f"# QC summary\n\nSource: {source}. Samples: {n}. Pass the QC cutoffs: {k}. "
             f"Excluded: {n - k}. Per-sample values: `qc_metrics.tsv`."]

    shown = {key: v for key, v in gate_used.items() if not key.startswith("auto_")}
    parts.append("## QC cutoffs\n" + "\n".join(f"- {l}" for l in gate_info.splitlines()) + "\n\n"
                 + _md_table(["cutoff", "value"],
                             [(key, "no cutoff" if v is None else v) for key, v in shown.items()]))

    if n - k:
        ex = df[~df["gate_include"]]
        parts.append(f"## Excluded samples ({n - k})\n"
                     + _md_table(["sample", "reasons"], list(ex["gate_reasons"].items())))
    else:
        parts.append("## Excluded samples\nNone: every sample passes the QC cutoffs "
                     "(downstream_analysis_samples.txt is not written).")

    rows = []
    for col, label in RANGE_COLS.items():
        if col in df and df[col].notna().any():
            x = df[col].dropna()
            rows.append((label, f"{x.min():.2f}", f"{x.median():.2f}", f"{x.max():.2f}"))
    parts.append("## Metric ranges\n" + _md_table(["metric", "min", "median", "max"], rows))

    notes = [f"Not assessed in this run: {', '.join(missing) or 'none'}."]
    partial = df.drop(columns=["gate_reasons"]).isna().sum()
    partial = partial[(partial > 0) & (partial < len(df))]
    if len(partial):
        notes.append("Missing for some samples (possible sample-name mismatch between tools): "
                     + ", ".join(f"{c} ({m})" for c, m in partial.items()))
    parts.append("## Missing metrics\n" + "\n".join(f"- {x}" for x in notes))

    if crosscheck:
        parts.append("## Low-QC samples vs PCA\n" + crosscheck)
    if pca_var and Path(pca_var).exists():
        v = pd.read_csv(pca_var).head(3)
        parts.append("## PCA variance explained\n"
                     + _md_table(list(v.columns), v.round(1).values.tolist()))
    if overrep is not None:
        parts.append("## Top overrepresented sequences (FastQC)\n"
                     + _md_table(list(overrep.columns), overrep.values.tolist()))
    parts.append("## Run context (effective config.yaml)\n" + run_context(cfg))
    return "\n\n".join(parts) + "\n"


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kit", required=True, help="watchmaker-polya | lexogen-3prime "
                                                 "(name of a file in the cutoffs folder)")
    ap.add_argument("--cutoffs-dir", default=str(CUTOFFS_DIR), help="folder with <kit>.yaml cutoff files")
    ap.add_argument("--aln-data", help="multiqc_data/ folder of the alignment report")
    ap.add_argument("--raw-data", help="multiqc_data/ folder of the raw (FastQC) report")
    ap.add_argument("--aln-html", help="fallback: alignment_QC.html")
    ap.add_argument("--raw-html", help="fallback: raw_QC.html")
    ap.add_argument("--config", help="pipeline config.yaml")
    ap.add_argument("--organism", help="skip detection: human, mouse, zebrafish_lawson, "
                                       "zebrafish_ensembl, c_elegans, drosophila, yeast")
    ap.add_argument("--gate", help="JSON with thresholds that override the gate file")
    ap.add_argument("--mate-regex", default=r"(_Forward|_Reverse|_R[12](_001)?|_[12])$",
                    help="suffix that marks R1/R2 in FastQC sample names")
    ap.add_argument("--pca-pcs", help="optional: PCA coordinates CSV (e.g. plots/PCA_top_PCs.csv)")
    ap.add_argument("--pca-var", help="optional: PCA variance CSV (e.g. plots/PCA_top_variance_explained.csv)")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

    # cutoffs: <kit>.yaml section for the organism, then remove_rRNA, then --gate
    # (yaml keeps the last value of a duplicated key, which is the effective one)
    cfg = (yaml.safe_load(Path(a.config).read_text()) if a.config and Path(a.config).exists()
           else None) or {}
    user_gate = json.loads(a.gate or "{}")
    organism = a.organism or detect_organism(cfg)
    gate_file, gate, section, unverified = load_gate(a.cutoffs_dir, a.kit, organism)
    if cfg.get("remove_rRNA") and "max_rrna" not in user_gate:
        if gate.get("max_rrna_filtered") is not None:
            gate["max_rrna"] = gate["max_rrna_filtered"]
        elif gate.get("max_rrna") is not None:
            unverified.append("max_rrna (remove_rRNA is true but this kit has no post-filter cutoff)")
    gate.pop("max_rrna_filtered", None)
    gate.update(user_gate)
    unverified = [u for u in unverified if u.split()[0] not in user_gate]

    if a.aln_data:
        aln, raw, overrep = from_multiqc_data(a.aln_data, a.raw_data)
        source = "multiqc_data"
    else:
        aln, raw, overrep = from_html(a.aln_html, a.raw_html)
        source = "HTML (fallback; General Statistics only)"
    raw = collapse_mates(raw, a.mate_regex)
    df = aln.join(raw, how="left") if raw is not None else aln
    df.index.name = "sample"

    rrna_reliable = gate.pop("rrna_reliable", True)
    expected_unique = gate.pop("expected_unique_pct", 70)
    auto_keys = [k for k, v in gate.items() if v == "auto"]
    df, gate_used, missing, notes = apply_gate(df, gate)
    df, crosscheck = pca_crosscheck(df, a.pca_pcs)
    df.to_csv(out / "qc_metrics.tsv", sep="\t")

    how = "from --organism/qc_organism" if a.organism else "detected from config.yaml"
    gate_lines = [f"cutoffs: {gate_file} -> organism {organism} ({how}), section '{section}'"]
    if user_gate:
        gate_lines.append(f"overridden by qc_gate: {json.dumps(user_gate)}")
    gate_lines += notes
    gate_lines.append("unverified values (not in SKILL.md; provisional): "
                      + (", ".join(unverified) or "none"))
    if not rrna_reliable:
        gate_lines.append("rRNA % for this organism is unreliable (see the skill); "
                          "a low value does not prove that rRNA is absent.")
    gate_info = "\n".join(gate_lines)
    print(gate_info)
    print(f"{len(df)} samples, {int(df.gate_include.sum())} pass the QC cutoffs; "
          f"low-QC for the PCA cross-check: {int(df.low_qc.sum())}")

    dl = out / "downstream_analysis_samples.txt"
    if (~df.gate_include).any():
        dl.write_text("\n".join(df.index[df.gate_include]) + "\n")
    elif dl.exists():
        dl.unlink()   # nothing excluded: no list

    summary = build_summary(df, gate_used, missing, cfg, overrep, a.pca_var, source,
                            gate_info, crosscheck)
    (out / "qc_summary.md").write_text(summary)
    written = [out / "qc_metrics.tsv", out / "qc_summary.md"] + ([dl] if dl.exists() else [])

    template = TEMPLATES_DIR / "email.md"
    email = (yaml.safe_load(Path(gate_file).read_text()) or {}).get("email") or {}
    if template.exists():
        hvg = Path(a.pca_pcs).parent / "pca_hvg_log.txt" if a.pca_pcs else None
        log = hvg.read_text() if hvg and hvg.exists() else ""
        sel = re.search(r"Selected genes.*?:\s*(\d+)", log)
        tot = re.search(r"Total genes:\s*(\d+)", log)
        draft = qc_email_text.render(template, df, dict(
            kit=a.kit, email=email, organism=organism, cfg=cfg, gate=gate_used, missing=missing,
            rrna_reliable=rrna_reliable, patterns=ORGANISM_PATTERNS, pipeline_dir=PIPELINE_DIR,
            pca_var=a.pca_var,
            pca_genes=(sel.group(1) if sel else None, tot.group(1) if tot else None),
            expected_unique=expected_unique, auto_keys=auto_keys,
            run_dir=Path(a.aln_data).resolve().parent if a.aln_data else Path(".")))
        (out / "qc_email_draft.md").write_text(draft)
        written.append(out / "qc_email_draft.md")
    else:
        print(f"WARNING: {template} not found; qc_email_draft.md not written")
    print("written: " + ", ".join(str(p) for p in written))


if __name__ == "__main__":
    main()
