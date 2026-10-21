#!/usr/bin/env python3
"""
Draft the RNA-seq QC client email from MultiQC output + a skill (SKILL.md).

Flow (a single API call per run):
  1. Extract per-sample metrics from multiqc_data/ (preferred) or the HTML report (fallback).
  2. Apply the inclusion gate DETERMINISTICALLY (thresholds from the kit defaults / config).
  3. Send SKILL.md (system prompt, cached) + a compact metrics table to Claude.
  4. Write: qc_metrics.tsv, downstream_analysis_samples.txt (only if samples are excluded),
     qc_assessment.md, qc_email_draft.md and a token-usage log.

Without ANTHROPIC_API_KEY (or with --dry-run) it only writes the prompt, for inspection.
"""
import argparse, base64, json, os, re, sys
from io import StringIO
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from mqc_lz import decompress_from_base64  # noqa: E402

# Default gate per kit (starting values from the skills; override with config: qc_gate)
GATE_DEFAULTS = {
    "watchmaker-polya": dict(min_strand=80, min_mrna=40, min_assigned_M=5, min_unique=40,
                             max_dup=None, max_rrna=15),
    "lexogen-3prime":   dict(min_strand=80, min_mrna=40, min_assigned_M="auto", min_unique=40,
                             max_dup=None, max_rrna=15),
}
AUTO_ASSIGNED_FRAC = 0.2   # "auto" = 20% of the run's median assigned reads

# Columns of interest, with alternative names across MultiQC versions
WANT = {
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
    out = pd.DataFrame()
    for key, names in WANT.items():
        for df in frames:
            s = _pick(df, names)
            if s is not None:
                s = pd.to_numeric(s, errors="coerce")
                if key in FRACTION_COLS and s.dropna().le(1.0).all():
                    s = s * 100
                out[key] = s
                break
    if "assigned_M" in out:
        out["assigned_M"] = out["assigned_M"] / 1e6
    if "unique_M" in out:
        out["unique_M"] = out["unique_M"] / 1e6

    # Strand: if Picard ran without STRAND, PCT_CORRECT_STRAND_READS is empty -> use plot data
    strand_f = next(iter(aln_dir.glob("mqc_picard_rnaseqmetrics_strand_plot*.txt")), None)
    if strand_f and ("strand_pct" not in out or out["strand_pct"].isna().all()):
        sd = _read_tsv(strand_f)
        c, i = _pick(sd, ["Correct"]), _pick(sd, ["Incorrect"])
        if c is not None and i is not None:
            out["strand_pct"] = 100 * c / (c + i)

    # Gene-body coverage -> 3'/5' ratio (mean of 3' half / mean of 5' half)
    cov_f = next(iter(aln_dir.glob("mqc_picard_rna_coverage*.txt")), None)
    if cov_f:
        cov = _read_tsv(cov_f).apply(pd.to_numeric, errors="coerce")
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
    return out, raw, None


def _html_tables(path):
    html = Path(path).read_text(encoding="utf-8", errors="ignore")
    m = re.search(r'id="mqc_compressed_plotdata">([^<]+)<', html)
    plots = json.loads(decompress_from_base64(m.group(1))) if m else {}
    tabs = pd.read_html(StringIO(html))
    return html, plots, tabs


def _num(s):
    return pd.to_numeric(s.astype(str).str.replace(r"[%a-zA-Z ]", "", regex=True), errors="coerce")


def from_html(aln_html, raw_html):
    """Fallback: General Statistics table only. For large cohorts MultiQC renders
    strand/coverage as static PNGs, so those metrics will be missing."""
    _, plots, tabs = _html_tables(aln_html)
    g = tabs[0].set_index("Sample Name")
    out = pd.DataFrame(index=g.index.astype(str))
    m = {"% Aligned": "unique_pct", "M Aligned": "unique_M", "% Assigned": "assigned_pct",
         "M Assigned": "assigned_M", "Duplication": "dup_pct", "rRNA": "rrna_pct",
         "mRNA": "mrna_pct", "M Reads": "bam_reads_M"}
    for col, key in m.items():
        if col in g:
            out[key] = _num(g[col]).values
    # strand from plot data, when present (small cohorts)
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
    """FastQC lists R1/R2 separately; merge them per sample."""
    if not isinstance(raw, pd.DataFrame) or raw.empty:
        return None
    raw = raw.copy()
    raw.index = raw.index.astype(str).str.replace(regex, "", regex=True)
    agg = {c: "mean" for c in raw.columns if c != "read_len"}
    if "raw_M" in agg:
        agg["raw_M"] = "first"   # mates have the same read count
    return raw.groupby(level=0).agg(agg)


# ----------------------------------------------------------------- gate
def apply_gate(df, gate):
    g = dict(gate)
    if g.get("min_assigned_M") == "auto":
        g["min_assigned_M"] = round(AUTO_ASSIGNED_FRAC * df["assigned_M"].median(), 2)
    reasons = pd.Series([[] for _ in range(len(df))], index=df.index)
    missing = []

    def rule(col, ok, label):
        if col not in df or df[col].isna().all():
            missing.append(col)
            return
        bad = ~ok(df[col]) & df[col].notna()
        for s in df.index[bad]:
            reasons[s].append(label)

    rule("strand_pct", lambda x: x >= g["min_strand"], f"strand<{g['min_strand']}%")
    rule("mrna_pct", lambda x: x >= g["min_mrna"], f"mRNA<{g['min_mrna']}%")
    rule("assigned_M", lambda x: x >= g["min_assigned_M"], f"assigned<{g['min_assigned_M']}M")
    rule("unique_pct", lambda x: x >= g["min_unique"], f"unique<{g['min_unique']}%")
    rule("rrna_pct", lambda x: x <= g["max_rrna"], f"rRNA>{g['max_rrna']}%")
    if g.get("max_dup") is not None:
        rule("dup_pct", lambda x: x <= g["max_dup"], f"dup>{g['max_dup']}%")
    df = df.copy()
    df["gate_include"] = reasons.map(len).eq(0)
    df["gate_reasons"] = reasons.map(lambda r: ";".join(r))
    return df, g, sorted(set(missing))


# ----------------------------------------------------------------- prompt
def run_context(cfg_path):
    if not cfg_path or not Path(cfg_path).exists():
        return "config.yaml not provided."
    # yaml.safe_load keeps the LAST value of duplicated keys (= the effective one)
    cfg = yaml.safe_load(Path(cfg_path).read_text()) or {}
    keys = ["aligner_index", "annotation_gtf", "picard_refflat", "picard_rrna_list",
            "remove_rRNA", "umi", "paired_end", "strandedness"]
    return "\n".join(f"{k}: {cfg.get(k)}" for k in keys if k in cfg) or "(no relevant keys)"


def build_prompt(df, gate_used, missing, ctx, overrep, pca_tsv, group_regex, source):
    cols = [c for c in ["raw_M", "gc_pct", "unique_pct", "unique_M", "multi_pct", "tooshort_pct",
                        "assigned_M", "assigned_pct", "mrna_pct", "intronic_pct",
                        "intergenic_pct", "rrna_pct", "strand_pct", "cov_3to5_ratio",
                        "bias_5to3", "dup_pct", "gate_include", "gate_reasons"] if c in df]
    tsv = df[cols].round(2).to_csv(sep="\t")
    n, k = len(df), int(df["gate_include"].sum())
    parts = [
        "The data below HAVE ALREADY BEEN EXTRACTED from MultiQC by a script; there is no need to "
        "decode plot data or read HTML. The inclusion gate HAS ALREADY BEEN COMPUTED "
        "deterministically: use `gate_include` as is, do not recompute it.",
        f"Source: {source}. Samples: {n}. Pass the gate: {k}. Excluded: {n - k}.",
        f"Gate thresholds used: {json.dumps(gate_used)}",
        f"Metrics NOT available in this run: {', '.join(missing) or 'none'}. "
        "Do not invent values for them; state in the email what could not be assessed.",
        f"## Run context (effective config.yaml)\n{ctx}",
        f"## Per-sample metrics (percentages in %, counts in millions)\n{tsv}",
    ]
    if group_regex:
        g = df.index.to_series().str.extract(group_regex)[0]
        tab = (pd.DataFrame({"total": g.value_counts(),
                             "included": g[df.gate_include.values].value_counts()})
               .fillna(0).astype(int).rename_axis("group"))
        parts.append("## Samples per group (total and included by the gate)\n" + tab.to_csv(sep="\t"))
    if overrep is not None:
        parts.append("## Top overrepresented sequences (FastQC)\n" + overrep.to_csv(sep="\t", index=False))
    if pca_tsv and Path(pca_tsv).exists():
        parts.append("## PCA (coordinates; first comment line gives % variance explained)\n"
                     + Path(pca_tsv).read_text())
    parts.append(
        "## Task\nFollow the skill. Reply with exactly two sections, in this order:\n"
        "`=== ASSESSMENT ===`: per-sample verdict ONLY for FAIL or PASS WITH CAVEATS samples "
        "(the skill's short format), then the cohort summary.\n"
        "`=== EMAIL ===`: the client email draft, ready for review.")
    return "\n\n".join(parts)


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skill", required=True, help="skill folder (contains SKILL.md)")
    ap.add_argument("--aln-data", help="multiqc_data/ folder of the alignment report")
    ap.add_argument("--raw-data", help="multiqc_data/ folder of the raw (FastQC) report")
    ap.add_argument("--aln-html", help="fallback: alignment_QC.html")
    ap.add_argument("--raw-html", help="fallback: raw_QC.html")
    ap.add_argument("--config", help="pipeline config.yaml")
    ap.add_argument("--kit", help="watchmaker-polya | lexogen-3prime (default: skill folder name)")
    ap.add_argument("--gate", help="JSON with thresholds that override the kit defaults")
    ap.add_argument("--group-regex",
                    help="regex with 1 capture group to get the group from the sample name, "
                         "e.g. '(Kidney|Liver|Spleen)$'")
    ap.add_argument("--mate-regex", default=r"(_Forward|_Reverse|_R[12](_001)?|_[12])$",
                    help="suffix that marks R1/R2 in FastQC sample names")
    ap.add_argument("--pca-tsv", help="optional: PCA coordinates")
    ap.add_argument("--pca-png", help="optional: PCA image (sent as an image, ~1-2k tokens)")
    ap.add_argument("--model", default="claude-sonnet-5-5")
    ap.add_argument("--max-tokens", type=int, default=6000)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--dry-run", action="store_true", help="write the prompt without calling the API")
    a = ap.parse_args()

    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    skill_dir = Path(a.skill)
    kit = a.kit or skill_dir.name

    if a.aln_data:
        aln, raw, overrep = from_multiqc_data(a.aln_data, a.raw_data)
        source = "multiqc_data"
    else:
        aln, raw, overrep = from_html(a.aln_html, a.raw_html)
        source = "HTML (fallback; General Statistics only)"
    raw = collapse_mates(raw, a.mate_regex)
    df = aln.join(raw, how="left") if raw is not None else aln
    df.index.name = "sample"

    gate = {**GATE_DEFAULTS.get(kit, GATE_DEFAULTS["watchmaker-polya"]), **json.loads(a.gate or "{}")}
    df, gate_used, missing = apply_gate(df, gate)
    df.to_csv(out / "qc_metrics.tsv", sep="\t")

    dl = out / "downstream_analysis_samples.txt"
    if (~df.gate_include).any():
        keep = df.index[df.gate_include]
        lines = list(keep)
        if a.group_regex:
            g = df.index.to_series().str.extract(a.group_regex)[0]
            cnt = g[df.gate_include.values].value_counts().reindex(g.unique(), fill_value=0)
            lines += ["", "# per group"] + [f"# {grp}: {c}" for grp, c in cnt.items()]
        dl.write_text("\n".join(lines) + "\n")
    elif dl.exists():
        dl.unlink()   # every sample passed: the skill says NOT to write this file

    system = "\n\n".join(p.read_text() for p in sorted(skill_dir.rglob("*.md")))
    user_text = build_prompt(df, gate_used, missing, run_context(a.config), overrep,
                             a.pca_tsv, a.group_regex, source)
    (out / "prompt_user.md").write_text(user_text)

    content = [{"type": "text", "text": user_text}]
    if a.pca_png:
        content.insert(0, {"type": "image", "source": {
            "type": "base64", "media_type": "image/png",
            "data": base64.b64encode(Path(a.pca_png).read_bytes()).decode()}})
        content.insert(1, {"type": "text", "text": "Image above: PCA (PC1 vs PC2) for this run."})

    est = (len(system) + len(user_text)) // 4
    log = out / "qc_email_usage.log"
    if a.dry_run or not os.environ.get("ANTHROPIC_API_KEY"):
        msg = f"dry-run: prompt written; ~{est} input tokens estimated (API not called)"
        log.write_text(msg + "\n"); print(msg); return

    import anthropic
    client = anthropic.Anthropic()
    resp = client.messages.create(
        model=a.model, max_tokens=a.max_tokens,
        system=[{"type": "text", "text": system, "cache_control": {"type": "ephemeral"}}],
        messages=[{"role": "user", "content": content}],
    )
    text = "".join(b.text for b in resp.content if b.type == "text")
    assessment, _, email = text.partition("=== EMAIL ===")
    (out / "qc_assessment.md").write_text(assessment.replace("=== ASSESSMENT ===", "").strip() + "\n")
    (out / "qc_email_draft.md").write_text((email or text).strip() + "\n")
    u = resp.usage
    with log.open("a") as f:
        f.write(f"model={a.model} input={u.input_tokens} output={u.output_tokens} "
                f"cache_write={getattr(u, 'cache_creation_input_tokens', 0)} "
                f"cache_read={getattr(u, 'cache_read_input_tokens', 0)} stop={resp.stop_reason}\n")
    print(log.read_text())


if __name__ == "__main__":
    main()
