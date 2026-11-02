#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QC cutoffs and client email draft (run_qc_email: true, qc_kit set).
# Called from the Snakefile's onsuccess block; writes to qc_email/.
# A failure here never fails the pipeline.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import json
import os as _os
import subprocess as _sp
import sys as _sys


_QC_GATE_FILE = f"additional_rules/qc_email/cutoffs/{config.get('qc_kit')}.yaml"


def _qc_email_is_stale():
    """True when MultiQC, the cutoffs, the template or the config changed since the last run."""
    marker = "qc_email/qc_metrics.tsv"
    if not _os.path.exists(marker):
        return True
    sources = ["multiqc_report.html", _QC_GATE_FILE, "config.yaml",
               f"additional_rules/qc_email/templates/{config.get('qc_kit')}.md"]
    newest = max(_os.path.getmtime(f) for f in sources if _os.path.exists(f))
    return newest > _os.path.getmtime(marker)


def run_qc_email_module():
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
        if config.get("qc_organism"):
            cmd += ["--organism", config["qc_organism"]]
        if config.get("qc_gate"):
            _g = config["qc_gate"]
            cmd += ["--gate", _g if isinstance(_g, str) else json.dumps(_g)]
        _os.makedirs("qc_email", exist_ok=True)
        print("[qc_email] applying the QC cutoffs...")
        with open("qc_email/qc_email.log", "w") as _log:
            try:
                rc = _sp.call(cmd, stdout=_log, stderr=_sp.STDOUT)
            except OSError as e:
                _log.write(f"could not start {_py}: {e}\n")
                rc = 1
        if rc == 0:
            print("[qc_email] done: see qc_email/qc_summary.md")
        else:
            print(f"[qc_email] FAILED (exit {rc}); see qc_email/qc_email.log. Pipeline results are unaffected.")
