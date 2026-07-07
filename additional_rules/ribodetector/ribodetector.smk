rule ribodetector:
    """
    Ribosomal RNA filtering
    """
    input:
        trimmed_read_1 = "trimming/{sample}.R1.trim.fastq.gz",
        trimmed_read_2 = "trimming/{sample}.R2.trim.fastq.gz" if config["layout"]=="paired" else []
    output:
        filtered_read_1 = "ribodetector/{sample}/{sample}.nonrrna.1.fq.gz",
        rrna_reads_1 = "ribodetector/{sample}/{sample}.rrna.1.fq.gz",
        filtered_read_2 = "ribodetector/{sample}/{sample}.nonrrna.2.fq.gz" if config["layout"]=="paired" else [],
        rrna_reads_2 = "ribodetector/{sample}/{sample}.rrna.2.fq.gz" if config["layout"]=="paired" else []
    resources: cpus="10", maxtime="3:00:00", mem_mb=300000
    message: "Removing rRNA sequences for {wildcards.sample} reads with ribodetector."
    conda:
    	"../../env_config/rna_ribodetector.yaml"
    container: "docker://ghcr.io/dartmouth-data-analytics-core/rna_ribodetector:2.0"
    params:
        sample = lambda wildcards: wildcards.sample,
        read_length = config["read_length"],
        layout = config["layout"],
        filtering_method = "norrna" if config["layout"]=="paired" else "none"
    log:
        stdout = "ribodetector/{sample}/stdout.log",
        stderr = "ribodetector/{sample}/stderr.log"
    shell: """
        if [ "{params.layout}" = "paired" ]; then
            ribodetector_cpu -t {resources.cpus} \
                -l {params.read_length} \
                -i {input.trimmed_read_1} {input.trimmed_read_2} \
                -e norrna \
                -o {output.filtered_read_1} {output.filtered_read_2} \
                -r {output.rrna_reads_1} {output.rrna_reads_2} \
                > {log.stdout} 2> {log.stderr}
        else
            ribodetector_cpu -t {resources.cpus} \
                -l {params.read_length} \
                -i {input.trimmed_read_1} \
                -o {output.filtered_read_1} \
                -r {output.rrna_reads_1} \
                > {log.stdout} 2> {log.stderr}
        fi
    """

rule ribodetector_mqc_summary:
    """
    Generate ribodetector multiqc content
    """
    input:
        expand("ribodetector/{sample}/stderr.log", sample=sample_list)
    output:
        "ribodetector/rrna_norrna_pct_mqc.tsv"
    params:
        ribo_dir = "ribodetector"
    resources: cpus="1", maxtime="10:00", mem_mb=2000
    message: "Generating ribodetector report."
    shell: """
        bash scripts/ribo_stats.sh ribodetector
    """
