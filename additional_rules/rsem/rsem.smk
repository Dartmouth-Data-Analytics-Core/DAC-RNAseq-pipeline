rule rsem:
    """
    Isoform counting
    """
    input:
        "alignment/{sample}.srt.bam",
    output:
        "rsem/{sample}.genes.results",
        "rsem/{sample}.isoforms.results"
    params:
        sample = lambda wildcards:  wildcards.sample,
        rsem_calc_exp_path = config['rsem_calc_exp_path'],
        rsem_ref_path = config["rsem_ref_path"],
        rsem_strandedness = config["rsem_strandedness"],
        rsem_paired_flag = '--paired-end' if config["layout"]=='paired' else '',
    conda:
        "../../env_config/rna_rsem.yaml",
    container: "docker://ghcr.io/dartmouth-data-analytics-core/rna_rsem:2.0"
    resources: cpus="10", maxtime="8:00:00", mem_mb=60000,
    message: "Counting transcript isoforms for {wildcards.sample} reads with RSEM."
    shell: """
        rsem-calculate-expression \
          {params.rsem_paired_flag} \
          --alignments \
          -p {resources.cpus} \
          --strandedness {params.rsem_strandedness} \
          --no-bam-output \
          alignment/{params.sample}.Aligned.toTranscriptome.out.bam \
          {params.rsem_ref_path} \
          rsem/{params.sample}
 """
