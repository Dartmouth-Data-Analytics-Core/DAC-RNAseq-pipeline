rule rustqc:
    """
    RNA-seq QC with RustQC (containerized)
    """
    input:
        bam = "markdup/{sample}.mkdup.bam"
    output:
        qc_dir = directory("qc/{sample}")
    container: "docker://ghcr.io/seqeralabs/rustqc:v0.2.1"
    params:
        gtf = config["annotation_gtf"],
        paired_flag = '-p' if config['layout']=='paired' else ''
    threads: 4
    resources: cpus="4", maxtime="8:00:00", mem_mb=40000
    message: "Running comprehensive QC for {wildcards.sample} with RustQC."
    log: "logs/rustqc/{sample}.log"
    shell: """

        rustqc rna \
            {input.bam} \
            --gtf {params.gtf} \
            {params.paired_flag} \
            --skip-dup-check \
            -o {output.qc_dir} 2> {log}


    """
