"""Fill templates/<kit>.md with sentences built from the QC table; [brackets] are for the analyst."""
import re
from pathlib import Path

import pandas as pd

import run_info

# wording rules (not from SKILL.md)
HIGH_GC_SAMPLE_PCT = 5.5     # % of reads above 70% GC
OVERREP_FAIL_PCT = 1.0       # FastQC fail limits
ADAPTER_FAIL_PCT = 10.0
THREE_PRIME_ENRICHED = 1.5   # 3'/5' coverage ratio
UNIFORM_COVERAGE = (0.8, 1.25)
LOW_QC_GAP = 1.10            # jump in assigned reads needed to name low samples
STRAND_GOOD = 90
PHRED_DROP = 1.0

ORGANISM_LABEL = {"human": "human", "mouse": "mouse", "zebrafish_lawson": "zebrafish",
                  "zebrafish_ensembl": "zebrafish", "c_elegans": "C. elegans",
                  "drosophila": "Drosophila", "yeast": "yeast"}
ASSEMBLY_ALIAS = {"danrer11": "GRCz11 (danRer11)", "danrer10": "GRCz10 (danRer10)",
                  "hg38": "GRCh38 (hg38)", "hg19": "GRCh37 (hg19)",
                  "mm39": "GRCm39 (mm39)", "mm10": "GRCm38 (mm10)",
                  "ce11": "WBcel235 (ce11)", "dm6": "BDGP6 (dm6)", "saccer3": "R64 (sacCer3)"}
REASON_LABEL = {"strand": ("strand_pct", "correct strand", "%"),
                "mRNA": ("mrna_pct", "mRNA", "%"),
                "assigned": ("assigned_M", "assigned reads", "M"),
                "unique": ("unique_pct", "unique alignment", "%"),
                "rRNA": ("rrna_pct", "rRNA", "%"),
                "dup": ("dup_pct", "duplication", "%")}


# ----------------------------------------------------------------- helpers
MAX_LISTED = 10   # longer lists become a count


def _names(idx, cap=MAX_LISTED):
    idx = list(idx)
    if len(idx) > cap:
        return f"{len(idx)} samples (listed in qc_metrics.tsv)"
    return idx[0] if len(idx) == 1 else ", ".join(idx[:-1]) + " and " + idx[-1]


def _rng(s, fmt="{:.0f}", unit="%"):
    s = s.dropna()
    lo, hi = fmt.format(s.min()), fmt.format(s.max())
    return f"{lo}{unit}" if lo == hi else f"{lo}{unit} to {hi}{unit}"


def _has(df, col):
    return col in df and df[col].notna().any()


def _num_word(n):
    return ["no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine",
            "ten"][n] if n <= 10 else str(n)


# ----------------------------------------------------------------- header
def reference(cfg, organism, patterns):
    idx, gtf = str(cfg.get("aligner_index") or ""), str(cfg.get("annotation_gtf") or "")
    base = "zebrafish" if organism.startswith("zebrafish") else organism
    m = re.search(patterns.get(base, r"$^"), f"{idx} {gtf}", re.I)
    if not m:
        return "[organism + assembly]"
    primary = " primary assembly" if "primary" in idx.lower() else ""
    name = ASSEMBLY_ALIAS.get(m.group(0).lower(), m.group(0))
    return f"{ORGANISM_LABEL.get(organism, organism)} {name}{primary}"


def annotation(cfg, organism):
    gtf = Path(str(cfg.get("annotation_gtf") or "")).name
    if not gtf:
        return "[annotation]"
    m = re.search(r"V4\.3[\d.]*\d", gtf)
    if organism == "zebrafish_lawson" and m:
        return f"Lawson lab {m.group(0)}"
    m = re.search(r"gencode\.(vM?\d+)", gtf, re.I)
    if m:
        return f"GENCODE {m.group(1)}"
    m = re.search(r"\.(\d{2,3})\.(chr\.)?gtf", gtf)
    if m:
        return f"Ensembl release {m.group(1)}"
    return gtf


def pipeline_version(pipeline_dir):
    f = Path(pipeline_dir) / "Snakefile"
    m = re.search(r"v(\d+(?:\.\d+)+)", f.read_text()[:500]) if f.exists() else None
    return f" v{m.group(1)}" if m else ""


def steps(cfg, ref, annot, run_dir="."):
    aligner = {"star": "STAR", "hisat2": "HISAT2"}.get(str(cfg.get("aligner_name", "star")).lower(),
                                                      str(cfg.get("aligner_name")))
    cut = run_info.cutadapt_phrases(run_info.cutadapt_params(run_dir, cfg))
    trim = ("Trimming: reads were processed with CutAdapt"
            + (" to " + run_info.join_phrases(cut) if cut else ""))
    s = [trim]
    if cfg.get("remove_rRNA"):
        s.append("Ribodetection: ribosomal RNA reads were identified with RiboDetector and removed")
    s += [
          f"Alignment: reads were aligned to the {ref} using {aligner}",
          f"Quantification: gene expression was quantified with featureCounts (Subread) "
          f"using the {annot} annotation"]
    return "\n".join(f"{i}. {x}" for i, x in enumerate(s, 1))


# ----------------------------------------------------------------- sections
def _named_low(df, inc):
    """Largest bottom group of low-QC samples set apart by a LOW_QC_GAP jump; (named, rest)."""
    if not _has(df, "assigned_M"):
        return df.iloc[0:0], df[inc]
    rec = df[inc].sort_values("assigned_M")
    a = rec["assigned_M"].tolist()
    m = int(rec["low_qc"].sum())
    for k in range(min(m, len(a) - 1), 0, -1):
        if a[k] >= LOW_QC_GAP * a[k - 1]:
            return rec.iloc[:k], rec.iloc[k:]
    return rec.iloc[0:0], rec


def summary(df):
    n = len(df)
    inc = df["gate_include"]
    ex = df[~inc]
    out = []
    if len(ex):
        out.append(f"After evaluating all QC metrics, we recommend limiting downstream analysis to "
                   f"{int(inc.sum())} of the {n} samples, listed in "
                   "`downstream_analysis_samples.txt`.")
        if len(ex) <= 10:
            items = []
            for s, r in ex.iterrows():
                parts = []
                for code in str(r["gate_reasons"]).split(";"):
                    key = re.split(r"[<>]", code)[0]
                    col, label, unit = REASON_LABEL.get(key, (None, code, ""))
                    v = r.get(col) if col else None
                    parts.append(f"{label} {v:.1f}{unit}" if v is not None and pd.notna(v) else label)
                items.append(f"- {s}: {', '.join(parts)}")
            out.append("The following samples did not meet our QC thresholds:\n" + "\n".join(items))
        else:
            out.append(f"{len(ex)} samples did not meet our QC thresholds; see "
                       "`downstream_analysis_samples.txt` for the recommended set.")
    low, rest = _named_low(df, inc)
    if df.loc[inc, "low_qc"].any() and _has(df, "assigned_M") and not len(low):
        out.append(f"Reads assigned to genes varied across the recommended samples "
                   f"(~{_rng(df.loc[inc, 'assigned_M'], '{:.1f}', 'M')}), with some samples notably "
                   "lower; downstream analysis should check whether these metrics drive any "
                   "expression differences. [Interpretation / recommendation.]")
        low = low.iloc[0:0]
    if len(low) and _has(df, "assigned_M"):
        lead = ("No samples were excluded by our QC thresholds, but" if not len(ex)
                else "Among the recommended samples,")
        s = (f"{lead} {_num_word(len(low))} sample{'s' if len(low) > 1 else ''} "
             f"({_names(low.sort_values('assigned_M').index)}) "
             f"{'have' if len(low) > 1 else 'has'} concerningly low numbers of reads assigned to genes "
             f"(~{_rng(low['assigned_M'], '{:.1f}', 'M')}")
        if len(rest):
            s += f", compared with {_rng(rest['assigned_M'], '{:.1f}', 'M')} in the other {len(rest)} samples"
        s += ")"
        extra = []
        if _has(df, "unique_pct") and len(rest) and low["unique_pct"].max() < rest["unique_pct"].min():
            extra.append(f"the lowest alignment rates ({_rng(low['unique_pct'])})")
        if _has(df, "tooshort_pct") and len(rest) and \
                low["tooshort_pct"].min() > rest["tooshort_pct"].max():
            extra.append(f"the highest fraction of reads too short to map "
                         f"({_rng(low['tooshort_pct'])})")
        if extra:
            s += ", with " + " and ".join(extra)
        out.append(s + ". [Interpretation / recommendation, e.g. exclude if a replicate passed, "
                       "sequence deeper, or keep with a caveat.]")
    if not len(ex) and not df["low_qc"].any():
        out.append(f"After evaluating all QC metrics, all {n} samples passed and are recommended "
                   "for downstream analysis.")
    out.append("To run differential expression we would need the sample metadata to define the "
               "experimental groups and contrasts; generally there should be at least 3 replicates "
               "in each group.")
    return "\n\n".join(out)


def _scope(df, rec):
    """Recommended samples and a suffix, when some were excluded."""
    if rec is None or len(rec) == len(df):
        return df, ""
    return rec, " across the recommended samples"


def raw_qc(df, cfg, rec=None):
    if not _has(df, "raw_M"):
        return ("[Raw FASTQ QC not available: set qc_raw_multiqc_data in config.yaml, "
                "e.g. fastqc_results/multiqc_data]")
    unit = "read pairs" if str(cfg.get("layout", "")).lower() == "paired" else "reads"
    d, note = _scope(df, rec)
    out = [f"Library depth averaged ~{d['raw_M'].mean():.1f}M {unit}{note} "
           f"(range: {d['raw_M'].min():.1f}M to {d['raw_M'].max():.1f}M)."]
    if _has(df, "mean_phred"):
        q = df["mean_phred"]
        cand = q[q < q.median() - PHRED_DROP]
        rest = q.drop(cand.index)
        named = cand[cand.round() < round(rest.min())] if len(rest) else cand.iloc[0:0]
        if q.min() >= 30 and 0 < len(named) <= 3:
            ok = q.drop(named.index)
            lo, hi = f"{ok.min():.0f}", f"{ok.max():.0f}"
            out.append(f"Mean Phred scores were ~Q{lo}{'' if lo == hi else '-Q' + hi} across the "
                       f"read body in {len(ok)} of {len(q)} samples; "
                       + ", ".join(f"{s} was slightly lower (~Q{v:.0f})" for s, v in named.items())
                       + " but still acceptable.")
        elif q.min() >= 30:
            lo, hi = f"{q.min():.0f}", f"{q.max():.0f}"
            out.append(f"Mean Phred scores were ~Q{lo}{'' if lo == hi else '-Q' + hi} across the "
                       "read body in all samples.")
        else:
            out.append(f"Mean Phred scores ranged Q{q.min():.0f}-Q{q.max():.0f}; below Q30: "
                       f"{_names(df.index[q < 30])}.")
    if _has(df, "gc_main_peak"):
        out.append(f"The main GC peak sat at ~{_rng(df['gc_main_peak'])} in all samples.")
    if _has(df, "high_gc_pct"):
        hi = df[df["high_gc_pct"] > HIGH_GC_SAMPLE_PCT].sort_values("high_gc_pct", ascending=False)
        if len(hi):
            who = (f"All {len(df)} samples" if len(hi) == len(df) else
                   f"{_num_word(len(hi)).capitalize()} sample{'s' if len(hi) > 1 else ''}")
            s = (f"{who} also showed a secondary peak at high GC (>70% GC; "
                 f"{_rng(hi['high_gc_pct'], '{:.0f}')} of reads)")
            extra = []
            adapters = _has(hi, "adapter_max_pct") and (hi["adapter_max_pct"] > ADAPTER_FAIL_PCT).all()
            if adapters:
                extra.append("adapter-containing")
            if _has(hi, "overrep_top_pct") and (hi["overrep_top_pct"] > OVERREP_FAIL_PCT).all():
                extra.append("overrepresented")
            if extra:
                s += f" together with {' and '.join(extra)} sequences"
            top = hi.index[:5] if len(hi) > MAX_LISTED else hi.index
            s += f", most pronounced in {_names(top)}."
            if adapters:
                s += (" [Suggested, please review: this is consistent with adapter dimers and other "
                      "fragments with very short or no insert being sequenced, which reduces the "
                      "number of usable reads for these samples.]")
            elif cfg.get("remove_rRNA"):
                s += (" [Suggested, please review: GC-rich residual ribosomal sequences in the raw "
                      "reads, which were filtered in silico before alignment.]")
            else:
                s += " [Possible cause: please review.]"
            out.append(s)
    return " ".join(out)


def alignment_qc(df, email, expected_unique, rec=None):
    full = df
    df, note = _scope(df, rec)
    n, out = len(df), []
    kind = " recommended" if note else ""
    if _has(df, "unique_pct"):
        u = df["unique_pct"]
        below = int((u < expected_unique).sum())
        if below:
            where = f"all {n}" if below == n else f"{below} of {n}"
            out.append(f"Unique alignment rates were below the ~{expected_unique:g}% typically "
                       f"expected in {where}{kind} samples, ranging from {_rng(u)}.")
        else:
            out.append(f"Unique alignment rates ranged from {_rng(u)}{note}.")
    df = full
    if _has(df, "tooshort_pct") and _has(df, "other_pct") and \
            df["tooshort_pct"].median() > df["other_pct"].median():
        out.append(f"Most unaligned reads were classified as \"too short\" by STAR, "
                   f"{email.get('tooshort_cause', 'which usually points to short inserts')}.")
    d, note = _scope(df, rec)
    if _has(d, "assigned_M"):
        out.append(f"Reads assigned to genes ranged from ~{_rng(d['assigned_M'], '{:.1f}', 'M')}{note}.")
    return " ".join(out)


def coverage(df, mode):
    """mode: three_prime (3' tag kits) or uniform (full-length kits)."""
    if not _has(df, "cov_3to5_ratio"):
        return "Gene-body coverage was not available in this run."
    r = df["cov_3to5_ratio"]
    if mode == "three_prime":
        flat = df.index[r <= THREE_PRIME_ENRICHED]
        if not len(flat):
            return (f"Reads showed 3' end enrichment in all {len(df)} samples, as expected for the "
                    "QuantSeq 3' tag design, confirming the enrichment chemistry worked.")
        return (f"3' end enrichment was weak or absent in {_names(flat)} (a flat profile can mean "
                "gDNA contamination or failed 3' enrichment).")
    lo, hi = UNIFORM_COVERAGE
    skew = df.index[(r < lo) | (r > hi)]
    if not len(skew):
        return "Gene-body coverage was uniform in all samples, consistent with intact RNA."
    return f"Gene-body coverage showed a 3'/5' imbalance in {_names(skew)} [check for degradation]."


def strand(df):
    if not _has(df, "strand_pct"):
        return "Strand specificity could not be assessed in this run."
    s = df["strand_pct"]
    if s.min() >= STRAND_GOOD:
        return f"Correct-strand reads were >={s.min():.0f}% in all samples, consistent with a stranded library."
    low = s[s < STRAND_GOOD].sort_values(ascending=False)
    good = s[s >= STRAND_GOOD]
    if len(good) and len(low) <= MAX_LISTED:
        return (f"Strand specificity was high (correct-strand >={good.min():.0f}% in {len(good)} of "
                f"{len(s)} samples; " + ", ".join(f"{v:.1f}% for {k}" for k, v in low.items())
                + "), consistent with a stranded library.")
    return (f"Correct-strand reads ranged {_rng(s)}; below {STRAND_GOOD}%: {_names(low.index)}.")


def content(df, email, rec=None):
    df, note = _scope(df, rec)
    out = []
    if _has(df, "mrna_pct"):
        out.append(f"mRNA (exonic) percentages ranged ~{_rng(df['mrna_pct'])}{note}")
    parts = [(lab, _rng(df[c])) for c, lab in (("intronic_pct", "intronic"),
                                                 ("intergenic_pct", "intergenic")) if _has(df, c)]
    if parts:
        s = f"{parts[0][0]} content ranged ~{parts[0][1]}"
        if len(parts) == 2:
            s += f" and {parts[1][0]} ~{parts[1][1]}"
        if _has(df, "unique_pct"):
            other = df[[c for c in ("intronic_pct", "intergenic_pct") if c in df]].sum(axis=1)
            if df["unique_pct"].rank().corr(other.rank()) < -0.5:
                cause = email.get("intronic_cause", "gDNA carryover")
                s += f", highest in the lowest-aligning samples [possible cause: {cause}]"
        out.append(s)
    return ("; ".join(out) + ".") if out else ""


def rrna(df, missing, reliable):
    if any(m.startswith("rrna_pct") for m in missing) or not _has(df, "rrna_pct"):
        return "Ribosomal content could not be assessed with this annotation."
    r = df["rrna_pct"]
    s = ("Ribosomal content was <0.1% in all samples." if r.max() < 0.1
         else f"Ribosomal content ranged {_rng(r, '{:.1f}')}.")
    return s + ("" if reliable else " [rRNA % is less reliable for this organism.]")


def thresholds(g, auto_keys=()):
    run = " (set from this run's distribution)"
    t = []
    if g.get("min_strand") is not None:
        t.append(f"correct strand >={g['min_strand']}%")
    if g.get("min_mrna") is not None:
        cap = g.get("auto_mrna_cap")
        at_cap = cap is not None and float(g["min_mrna"]) >= float(cap)
        t.append(f"mRNA (PCT_MRNA_BASES) >={g['min_mrna']:g}%"
                 + (run if "min_mrna" in auto_keys and not at_cap else ""))
    if g.get("min_assigned_M") is not None:
        t.append(f"assigned reads >={g['min_assigned_M']}M" + (run if "min_assigned_M" in auto_keys else ""))
    if g.get("min_unique") is not None:
        t.append(f"unique alignment >={g['min_unique']}%")
    if g.get("max_rrna") is not None:
        t.append(f"rRNA <={g['max_rrna']}%")
    if g.get("max_dup") is not None:
        t.append(f"duplication <={g['max_dup']}%")
    return ("Taking the distribution of all QC metrics into consideration, we used the following "
            "thresholds to define samples with adequate quality for downstream analysis: "
            + ", ".join(t) + ".")


def pca(df, pca_var, pca_genes):
    out = []
    sel, total = pca_genes if isinstance(pca_genes, tuple) else (pca_genes, None)
    if sel and total and int(sel) >= int(total):
        out.append(f"PCA was performed on all {int(total):,} genes (the automatic selection of "
                   "variable genes kept every gene).")
    elif sel:
        out.append(f"PCA was performed on the {int(sel):,} most variable genes"
                   + (f" (of {int(total):,})." if total else "."))
    if pca_var and Path(pca_var).exists():
        v = pd.read_csv(pca_var).iloc[:3, 1].round(1).tolist()
        if len(v) == 3:
            out.append(f"PC1, PC2 and PC3 explained {v[0]}%, {v[1]}% and {v[2]}% of the variance.")
    if "pca_extreme" in df:
        low = df[df["low_qc"]]
        ext = low[low["pca_extreme"].fillna("").ne("")]
        if len(ext):
            out.append("Of the samples with the lowest QC, "
                       + "; ".join(f"{s} sits at the extreme of {r['pca_extreme']}" for s, r in ext.iterrows())
                       + ".")
        elif len(low):
            out.append("The samples with the lowest QC do not sit at the extremes of PC1-PC3.")
    out.append("[Interpretation: which samples separate, along which PC; whether QC drives it.]")
    return " ".join(out) if out else "[PCA not available]"


# ----------------------------------------------------------------- render
class _Keep(dict):
    def __missing__(self, key):
        return "[" + key + "]"


def render(template, df, ctx):
    cfg = ctx["cfg"]
    ref, annot = reference(cfg, ctx["organism"], ctx["patterns"]), annotation(cfg, ctx["organism"])
    excluded = (~df["gate_include"]).any()
    rec = df[df["gate_include"]] if excluded else None
    run_dir = Path(ctx.get("run_dir", "."))
    pipe = run_dir if (run_dir / "Snakefile").exists() else ctx["pipeline_dir"]
    email = ctx.get("email") or {}
    dup_range = _rng(df["dup_pct"]) if _has(df, "dup_pct") else "[duplication range]"
    f = _Keep(
        n_samples=len(df),
        kit_description=email.get("kit", "[library prep kit]"),
        duplication=email.get("duplication", "").format(dup_range=dup_range),
        tpm_note=email.get("tpm_note", "[Note on the top genes.]"),
        pipeline_version=pipeline_version(pipe),
        steps=steps(cfg, ref, annot, ctx.get("run_dir", ".")),
        summary=summary(df),
        raw_qc=raw_qc(df, cfg, rec),
        alignment_qc=alignment_qc(df, email, ctx["expected_unique"], rec),
        coverage=coverage(df, email.get("coverage", "uniform")),
        strand=strand(df),
        content=content(df, email, rec),
        rrna=rrna(df, ctx["missing"], ctx["rrna_reliable"]),
        thresholds=thresholds(ctx["gate"], ctx.get("auto_keys", ())),
        pca=pca(df, ctx["pca_var"], ctx["pca_genes"]),
        downstream_file_line=("* downstream_analysis_samples.txt (samples we recommend limiting "
                              "downstream analysis to)" if excluded else ""),
    )
    text = Path(template).read_text().format_map(f)
    # spaces left by empty fragments (keeps list indents)
    text = re.sub(r"(?m)^ (?=[^ ])", "", text)
    text = re.sub(r"(?<=\S)  +", " ", text)
    return re.sub(r"\n{3,}", "\n\n", text)
