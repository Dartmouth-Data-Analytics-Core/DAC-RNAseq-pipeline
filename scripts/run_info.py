"""
Read what a pipeline run actually did from its own logs, so text written for the client
(methods.txt, the email draft) follows the real parameters instead of assumed ones.
"""
import re
from pathlib import Path


def _first_file(run, pattern):
    return next(iter(sorted(Path(run).glob(pattern))), None)


def _flag(args, name):
    """Value of `--name=v`, `--name v` or `-n v` in an argument list, else None."""
    for i, a in enumerate(args):
        if a.startswith(name + "="):
            return a.split("=", 1)[1]
        if a == name and i + 1 < len(args):
            return args[i + 1]
    return None


# ----------------------------------------------------------------- cutadapt
def _rule_tokens(path, rule=None):
    """Whitespace tokens of a rule file (or of one rule inside it), for flag lookups."""
    f = Path(path)
    if not f.exists():
        return None
    text = f.read_text(errors="ignore")
    if rule:
        m = re.search(rf"^rule {rule}:\n(.*?)(?=^rule |\Z)", text, re.S | re.M)
        if not m:
            return None
        text = m.group(1)
    return text.replace("\\\n", " ").split()


def cutadapt_params(run, cfg=None):
    """Trimming options from the first cutadapt report ("Command line parameters: ...").
    Without reports (trimming/ removed after the run), fall back to the trimming rule of the
    pipeline copy in the run folder plus the configurable NextSeq flag from config.yaml."""
    f = _first_file(run, "trimming/*.cutadapt.report")
    m = re.search(r"Command line parameters: (.+)", f.read_text(errors="ignore")) if f else None
    if m:
        args = m.group(1).split()
    else:
        args = _rule_tokens(Path(run) / "Snakefile", "trimming")
        if args is None:
            return None
        args = [a for a in args if not a.startswith("{params.nextseq")]
        args += str((cfg or {}).get("cutadapt_nextseq_flag") or "").split()
    return {
        "nextseq": _flag(args, "--nextseq-trim"),
        "quality": _flag(args, "-q") or _flag(args, "--quality-cutoff"),
        "adapters": any(a in args for a in ("-a", "-A", "-g", "-G", "-b", "-B", "--adapter")),
        "poly_a": "--poly-a" in args,
        "trim_n": "--trim-n" in args,
        "max_n": _flag(args, "--max-n"),
        "min_len": _flag(args, "-m") or _flag(args, "--minimum-length"),
    }


def cutadapt_phrases(p):
    """Each thing cutadapt did, as a verb phrase: 'trim terminal N bases', ..."""
    if not p:
        return []
    out = []
    if p["nextseq"]:
        out.append(f"quality-trim 3′ ends with the NextSeq two-colour setting "
                   f"(--nextseq-trim={p['nextseq']}, which also removes poly-G tails)")
    elif p["quality"]:
        out.append(f"quality-trim 3′ ends (-q {p['quality']})")
    if p["adapters"]:
        out.append("remove adapter sequences")
    if p["poly_a"]:
        out.append("remove poly(A) tails")
    if p["trim_n"]:
        out.append("trim terminal N bases")
    if p["max_n"]:
        v = float(p["max_n"])
        out.append(f"discard reads with >{v * 100:g}% N content" if v < 1
                   else f"discard reads with more than {v:g} N bases")
    if p["min_len"] and int(p["min_len"]) > 1:
        out.append(f"discard reads shorter than {p['min_len']} bp")
    return out


def join_phrases(items):
    items = list(items)
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


# ----------------------------------------------------------------- STAR
def star_params(run):
    """STAR options from the command line recorded in the first Log.out; without logs, from the
    STAR rule of the pipeline copy in the run folder."""
    f = _first_file(run, "alignment/*.Log.out")
    m = re.search(r"##### Command Line:\n(.+)", f.read_text(errors="ignore")) if f else None
    if m:
        args = m.group(1).split()
    else:
        args = _rule_tokens(Path(run) / "additional_rules" / "alignment" / "star.smk", "alignment")
        if args is None:
            return None
    return {
        "twopass": (_flag(args, "--twopassMode") or "None") != "None",
        "multimap": _flag(args, "--outFilterMultimapNmax"),
        "mismatch_frac": _flag(args, "--outFilterMismatchNoverReadLmax"),
        "transcriptome": "TranscriptomeSAM" in args,
    }
