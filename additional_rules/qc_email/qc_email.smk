#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QC inclusion gate + summary for the client email
#
# Included by the Snakefile only when `run_qc_email: true`. It does not add
# a rule: it runs once at the end, after the whole pipeline (including the
# MultiQC step in rule all) finished successfully.
#
# Reads multiqc_report_data/ + the PCA tables and applies the inclusion gate
# deterministically (gates/<qc_kit>.yaml, section chosen by organism).
# Outputs go to qc_email/: qc_metrics.tsv, qc_summary.md, qc_email_draft.md
# (templates/<qc_kit>.md filled by rules, no model), qc_email.log and
# downstream_analysis_samples.txt (only when some samples are excluded).
#
# Requirements:
#   conda env create -f env_config/rna_qc_email.yaml   (once)
#   qc_email_python: <that env>/bin/python             (config.yaml; find it with
#                                                       `conda run -n rna_qc_email which python`)
#   qc_kit: watchmaker-polya | lexogen-3prime          (config.yaml)
# If qc_email_python is empty, the Python running snakemake is used instead.
# A failure here never fails the pipeline.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import json
import os as _os
import subprocess as _sp
import sys as _sys


_QC_GATE_FILE = f"additional_rules/qc_email/gates/{config.get('qc_kit')}.yaml"


def _qc_email_is_stale():
    """Re-run only when MultiQC, the gate file or the config changed since the last run,
    so a repeated `snakemake` call ("Nothing to be done") does not redo it."""
    marker = "qc_email/qc_metrics.tsv"
    if not _os.path.exists(marker):
        return True
    sources = ["multiqc_report.html", _QC_GATE_FILE, "config.yaml",
               f"additional_rules/qc_email/templates/{config.get('qc_kit')}.md"]
    newest = max(_os.path.getmtime(f) for f in sources if _os.path.exists(f))
    return newest > _os.path.getmtime(marker)


onsuccess:
    if not config.get("qc_kit"):
        print("[qc_email] skipped: set qc_kit in config.yaml (watchmaker-polya | lexogen-3prime)")
    elif not _os.path.isdir("multiqc_report_data"):
        print("[qc_email] skipped: multiqc_report_data/ not found")
    elif not _qc_email_is_stale():
        print("[qc_email] up to date, skipped")
    else:
        _py = config.get("qc_email_python") or _sys.executable
        cmd = [
            _py, config["qc_email_script"],
            "--kit", config["qc_kit"],
            "--aln-data", "multiqc_report_data",
            "--config", "config.yaml",
            "--pca-pcs", "plots/PCA_top_PCs.csv",
            "--pca-var", "plots/PCA_top_variance_explained.csv",
            "--outdir", "qc_email",
        ]
        if config.get("qc_raw_multiqc_data"):
            cmd += ["--raw-data", config["qc_raw_multiqc_data"]]
        if config.get("qc_organism"):   # skip detection from aligner_index / annotation_gtf
            cmd += ["--organism", config["qc_organism"]]
        if config.get("qc_gate"):       # overrides single keys of gates/<qc_kit>.yaml
            _g = config["qc_gate"]
            cmd += ["--gate", _g if isinstance(_g, str) else json.dumps(_g)]
        _os.makedirs("qc_email", exist_ok=True)
        print("[qc_email] applying the QC gate...")
        with open("qc_email/qc_email.log", "w") as _log:
            try:
                rc = _sp.call(cmd, stdout=_log, stderr=_sp.STDOUT)
            except OSError as e:   # e.g. qc_email_python points to a missing file
                _log.write(f"could not start {_py}: {e}\n")
                rc = 1
        if rc == 0:
            print("[qc_email] done: see qc_email/qc_summary.md")
        else:
            print(f"[qc_email] FAILED (exit {rc}); see qc_email/qc_email.log. Pipeline results are unaffected.")
