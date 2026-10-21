# =============================================================================
# rules/qc_email.smk  —  QC email draft with Claude (one API call per run)
#
# In the main Snakefile:
#     include: "rules/qc_email.smk"
#
#     rule all:
#         input:
#             ...,
#             QC_EMAIL_TARGETS,      # empty when qc_email: false
#
# In config.yaml:
#     qc_email: true
#     qc_skill: /path/to/GDSC-QC-Data-Skills/rnaseq/watchmaker-polya   # or lexogen-3prime
#     qc_email_outdir: results/qc_report                       # optional
#     qc_multiqc_aln_data: results/multiqc/alignment_QC_data   # multiqc_data folder (alignment)
#     qc_multiqc_raw_data: results/multiqc/raw_QC_data         # multiqc_data folder (raw)
#     qc_pca_png: results/Plots/PCA_top_PC1_vs_PC2.png         # optional
#     qc_group_regex: "(Kidney|Liver|Spleen)$"                 # optional
#     qc_gate: '{"min_assigned_M": 5}'                         # optional
#     qc_email_model: claude-sonnet-5-5                        # optional
#
# The API key lives in the environment (export ANTHROPIC_API_KEY=...), never in the config.
# Without the key the rule runs in dry-run mode: it writes the prompt and a token estimate.
# =============================================================================

QC_EMAIL_OUT = config.get("qc_email_outdir", "results/qc_report")

QC_EMAIL_TARGETS = (
    [f"{QC_EMAIL_OUT}/qc_email_draft.md"] if config.get("qc_email", False) else []
)

# Run on the login node: cluster compute nodes often have no internet access.
localrules: qc_email


def _qc_email_optional_args(wildcards):
    args = []
    if config.get("qc_group_regex"):
        args.append(f"--group-regex '{config['qc_group_regex']}'")
    if config.get("qc_gate"):
        args.append(f"--gate '{config['qc_gate']}'")
    if config.get("qc_pca_png"):
        args.append(f"--pca-png {config['qc_pca_png']}")
    return " ".join(args)


rule qc_email:
    input:
        aln=config.get("qc_multiqc_aln_data", "results/multiqc/alignment_QC_data"),
        raw=config.get("qc_multiqc_raw_data", "results/multiqc/raw_QC_data"),
        skill=config["qc_skill"] + "/SKILL.md",   # editing the skill re-runs the email
        cfg="config.yaml",
        script="scripts/qc_email.py",
        pca=config.get("qc_pca_png", []),
    output:
        metrics=f"{QC_EMAIL_OUT}/qc_metrics.tsv",
        email=f"{QC_EMAIL_OUT}/qc_email_draft.md",
    params:
        skill=config["qc_skill"],
        outdir=QC_EMAIL_OUT,
        model=config.get("qc_email_model", "claude-sonnet-5-5"),
        extra=_qc_email_optional_args,
    log:
        "logs/qc_email.log",
    shell:
        r"""
        python {input.script} \
            --skill {params.skill} \
            --aln-data {input.aln} \
            --raw-data {input.raw} \
            --config {input.cfg} \
            --model {params.model} \
            --outdir {params.outdir} \
            {params.extra} > {log} 2>&1

        # in dry-run (no API key) no email is generated; write a notice so Snakemake does not fail
        if [ ! -s {output.email} ]; then
            echo "Email not generated: ANTHROPIC_API_KEY is not set. See {params.outdir}/prompt_user.md" > {output.email}
        fi
        """
