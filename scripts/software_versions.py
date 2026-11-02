#!/usr/bin/env python3
"""
Write software_versions.txt: the steps that ran, the tool and version used for each, and the
genome / annotation names. Meant to be attached to the client results email.

Versions are read from what actually ran (STAR Log.out, cutadapt reports, BAM @PG header,
RustQC logs, FastQC / MultiQC data, the DESeq2 note written by rlog_transform.R). When a step
leaves no version behind, the version pinned in env_config/ is used.
"""
import argparse, datetime, gzip, json, re, subprocess, sys
from pathlib import Path

import yaml

PIPELINE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).parent))
from qc_email import ORGANISM_PATTERNS, ZEBRAFISH_ANNOTATION  # noqa: E402
import qc_email_text  # noqa: E402
import run_info  # noqa: E402


def _first(pattern, files):
    for f in files:
        try:
            m = re.search(pattern, Path(f).read_text(errors="ignore"))
        except OSError:
            continue
        if m:
            return m.group(1)
    return None


def _env_pin(env_file, package):
    f = PIPELINE_DIR / "env_config" / env_file
    if not f.exists():
        return None
    m = re.search(rf"^\s*-\s*{re.escape(package)}\s*==?\s*([\w.]+)", f.read_text(), re.M)
    return m.group(1) if m else None


def _bam_pg(run):
    """@PG entries (tool id -> version) from the first duplicate-marked BAM header."""
    bam = next(iter(sorted(run.glob("markdup/*.mkdup.bam"))), None)
    if not bam:
        return {}
    try:
        with gzip.open(bam, "rb") as f:
            head = f.read(400000)
    except OSError:
        return {}
    pairs = re.findall(rb"@PG\tID:([^\t\n]+)\t(?:PN:[^\t\n]+\t)?(?:PP:[^\t\n]+\t)?VN:([^\t\n]+)", head)
    return {k.decode().split(".")[0].lower(): v.decode() for k, v in pairs}


def _multiqc_version(folder):
    j = folder / "multiqc_data.json"
    try:
        return json.loads(j.read_text()).get("config_version")
    except (OSError, ValueError):
        return None


def _organism(cfg):
    idx, gtf = str(cfg.get("aligner_index") or ""), str(cfg.get("annotation_gtf") or "")
    hits = [o for o, p in ORGANISM_PATTERNS.items() if re.search(p, f"{idx} {gtf}", re.I)]
    if len(hits) != 1:
        return None
    if hits[0] != "zebrafish":
        return hits[0]
    return next((n for n, p in ZEBRAFISH_ANNOTATION.items() if re.search(p, gtf, re.I)), "zebrafish_lawson")


def steps(cfg, run):
    pg = _bam_pg(run)
    # raw FastQC MultiQC folder: config, else the names seen in pipeline runs
    raw_dir = next((run / d for d in (cfg.get("qc_raw_multiqc_data"), "fastqc_results/multiqc_data",
                                      "fastqc_results/fastqc_multiqc_data")
                    if d and (run / d).is_dir()), run / "fastqc_results/multiqc_data")
    aligner = str(cfg.get("aligner_name", "star")).lower()
    rows = [("Raw read QC", "FastQC",
             _first(r"FastQC\t([\d.]+)", [raw_dir / "multiqc_software_versions.txt"])
             or _env_pin("rna_fastqc.yaml", "fastqc"))]
    if cfg.get("remove_rRNA"):
        rows.append(("rRNA read filtering", "RiboDetector", _env_pin("rna_ribodetector.yaml", "ribodetector")))
    rows.append(("Adapter and quality trimming", "cutadapt",
                 _first(r"This is cutadapt ([\d.]+)", sorted(run.glob("trimming/*.cutadapt.report"))[:3])
                 or _env_pin("rna_cutadapt.yaml", "cutadapt")))
    if aligner == "hisat2":
        rows.append(("Alignment", "HISAT2", pg.get("hisat2") or _env_pin("rna_alignment.yaml", "hisat2")))
    else:
        rows.append(("Alignment", "STAR",
                     _first(r"STAR version=([\w.]+)", sorted(run.glob("alignment/*.Log.out"))[:3])
                     or pg.get("star") or _env_pin("rna_alignment.yaml", "star")))
    picard = pg.get("markduplicates") or _env_pin("rna_picard.yaml", "picard-slim")
    rows += [("Duplicate marking", "Picard MarkDuplicates", picard),
             ("Alignment metrics", "Picard CollectRnaSeqMetrics", picard),
             ("BAM indexing and statistics", "SAMtools", _env_pin("rna_samtools.yaml", "samtools"))]
    if cfg.get("run_rustqc"):
        rows.append(("Gene-body coverage and TIN", "RustQC",
                     _first(r"RustQC v([\d.]+)", sorted(run.glob("logs/rustqc/*.log"))[:3])))
    rows.append(("Gene quantification", "featureCounts (Subread)",
                 _first(r"v(\d+\.\d+\.\d+)", sorted(run.glob("slurm_logs/log_featurecounts_*.out"))[:1])
                 or _env_pin("rna_featurecounts.yaml", "subread")))
    if cfg.get("run_rsem"):
        rows.append(("Transcript quantification", "RSEM", _env_pin("rna_rsem.yaml", "rsem")))
    deseq = run / "plots" / "deseq2_version.txt"
    dv = _first(r"DESeq2 ([\d.]+)", [deseq]) or _env_pin("rna_pcaplot.yaml", "bioconductor-deseq2")
    method = _first(r"method (.+)", [deseq]) or "rlog"
    rows.append(("Normalization for PCA", f"DESeq2 ({method})", dv))
    rows.append(("QC reports", "MultiQC",
                 _multiqc_version(run / "multiqc_report_data") or _multiqc_version(run / "multiqc_data")
                 or _env_pin("rna_multiqc.yaml", "multiqc")))
    return rows


def methods(cfg, run, rows, genome, annot):
    """Methods paragraphs for the client, built from what ran (config + the run's own logs)."""
    ver = {tool: v for _, tool, v in rows}

    def v(tool, prefix="v"):
        return f" ({prefix}{ver[tool]})" if ver.get(tool) else ""

    paired = str(cfg.get("layout", "")).lower() == "paired"
    paras = []

    # raw QC, trimming, rRNA filtering
    s = (f"Raw {'paired-end' if paired else 'single-end'} FASTQ files were assessed for base "
         f"quality and nucleotide composition with FastQC{v('FastQC')}.")
    cut = run_info.cutadapt_phrases(run_info.cutadapt_params(run, cfg))
    s += (f" Reads were processed with Cutadapt{v('cutadapt')}"
          + (f" to {run_info.join_phrases(cut)}." if cut else "."))
    if cfg.get("remove_rRNA"):
        s += (f" Ribosomal RNA reads were identified with RiboDetector{v('RiboDetector')} "
              "and removed before alignment.")
    paras.append(s)

    # alignment and alignment QC
    if str(cfg.get("aligner_name", "star")).lower() == "hisat2":
        s = f"Reads were aligned to the {genome} using HISAT2{v('HISAT2')} with default parameters."
    else:
        st = run_info.star_params(run) or {}
        s = f"Reads were aligned to the {genome} using STAR{v('STAR')}"
        s += " in two-pass mode" if st.get("twopass") else ""
        opts = []
        if st.get("multimap"):
            opts.append(f"allowing up to {st['multimap']} multimapping loci")
        if st.get("mismatch_frac"):
            opts.append(f"mismatches up to {float(st['mismatch_frac']) * 100:g}% of read length")
        s += (", " + " and ".join(opts) if opts else "") + "."
    picard = f" (v{ver['Picard MarkDuplicates']}; " if ver.get("Picard MarkDuplicates") else " ("
    s += (" Post-alignment quality metrics (duplication rate, strand specificity, and "
          "coding/UTR/intronic/intergenic/rRNA read distribution) were computed with Picard Tools"
          f"{picard}MarkDuplicates and CollectRnaSeqMetrics)")
    if cfg.get("run_rustqc"):
        s += f", gene-body coverage and transcript integrity with RustQC{v('RustQC')}"
    else:
        s += f", alignment statistics with SAMtools{v('SAMtools')}"
    s += f", and aggregated across samples with MultiQC{v('MultiQC')}."
    paras.append(s)

    # quantification
    strand = {"0": "unstranded mode -s 0", "1": "forward-stranded mode -s 1",
              "2": "reverse-stranded mode -s 2"}.get(str(cfg.get("featurecounts_strand")),
                                                    f"-s {cfg.get('featurecounts_strand')}")
    fc_ver = f"Subread v{ver['featureCounts (Subread)']}; " if ver.get("featureCounts (Subread)") else ""
    s = ("Gene-level counts were generated from coordinate-sorted BAM files using featureCounts "
         f"({fc_ver}-t exon, -g gene_id, {strand}{', -p, counting read pairs' if paired else ''}) "
         f"with the {annot} annotation, excluding ambiguous and multi-mapped reads.")
    if cfg.get("run_rsem"):
        s += (f" Transcript- and gene-level expression was also estimated with RSEM{v('RSEM')} "
              f"from the STAR alignments to the transcriptome ({cfg.get('rsem_strandedness', 'none')} "
              "strandedness), which assigns multi-mapping reads probabilistically.")
    paras.append(s)

    # PCA
    deseq = next((t for t in ver if t.startswith("DESeq2")), None)
    dv = f" (v{ver[deseq]})" if deseq and ver.get(deseq) else ""
    if deseq and "log2" in deseq:
        s = (f"For principal component analysis, counts were normalized with DESeq2{dv} size "
             "factors and transformed as log2(normalized counts + 1)")
    else:
        s = (f"For principal component analysis, counts were transformed with the DESeq2{dv} "
             "regularized log (rlog, blind to sample groups)")
    s += ", centered without scaling, and restricted to the most variable genes."
    paras.append(s)
    return "\n\n".join(paras) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--run-dir", default=".", help="pipeline run folder (default: current folder)")
    ap.add_argument("--out", default="software_versions.txt")
    ap.add_argument("--methods-out", default="methods.txt",
                    help="methods paragraphs for the client (empty string to skip)")
    a = ap.parse_args()

    run = Path(a.run_dir)
    # yaml.safe_load keeps the LAST value of duplicated keys (= the effective one)
    cfg = yaml.safe_load((run / a.config).read_text()) or {}
    organism = _organism(cfg)
    if organism:
        genome = qc_email_text.reference(cfg, organism, ORGANISM_PATTERNS)
        annot = qc_email_text.annotation(cfg, organism)
    else:
        genome = Path(str(cfg.get("aligner_index") or "unknown")).name
        annot = Path(str(cfg.get("annotation_gtf") or "unknown")).name
    # the run folder is normally the pipeline clone that ran; fall back to this script's pipeline
    src = run if (run / "Snakefile").exists() else PIPELINE_DIR
    try:
        commit = subprocess.run(["git", "-C", str(src), "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, timeout=10).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        commit = ""

    rows = steps(cfg, run)
    w1 = max(len(r[0]) for r in rows) + 3
    w2 = max(len(r[1]) for r in rows) + 3
    version = qc_email_text.pipeline_version(src).strip() or "(version unknown)"
    lines = [f"GDSC-RNAseq Pipeline {version}" + (f" (commit {commit})" if commit else ""),
             f"Generated: {datetime.date.today().isoformat()}", "",
             f"{'Step':<{w1}}{'Tool':<{w2}}Version"]
    lines += [f"{s:<{w1}}{t:<{w2}}{v or '—'}" for s, t, v in rows]
    lines += ["", "References",
              f"{'Genome':<{w1}}{genome}",
              f"{'Annotation':<{w1}}{annot}"]
    out = run / a.out
    out.write_text("\n".join(lines) + "\n")
    print(f"written: {out}")
    if a.methods_out:
        m = run / a.methods_out
        m.write_text(methods(cfg, run, rows, genome, annot))
        print(f"written: {m}")


if __name__ == "__main__":
    main()
