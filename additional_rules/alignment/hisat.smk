rule alignment:
    """
    HISAT2 alignment
    """
    input:
        R1_FASTQ_INPUT,
        R2_FASTQ_INPUT if config["layout"] == "paired" else [],
    output:
        "alignment/{sample}.srt.bam",
        "alignment/{sample}.srt.bam.bai",
    params:
        layout = config["layout"],
        sample = lambda wildcards: wildcards.sample,
        aligner_name = config["aligner_name"],
        hisat2 = config["aligner_path"],
        aligner_index = config["aligner_index"],
        samtools = config["samtools_path"],
        fastq_1_flag = '-1' if config['layout'] == 'paired' else '-U',
        fastq_2 = '-2 trimming/{sample}.R2.trim.fastq.gz' if config['layout'] == 'paired' else '',
    conda:
        "../../env_config/rna_alignment.yaml",
    container: "docker://ghcr.io/dartmouth-data-analytics-core/rna_alignment:2.0"
    resources: cpus="4", maxtime="8:00:00", mem_mb=40000,
    message: "aligning {wildcards.sample} reads with Hisat2."
    shell: """
        {params.hisat2} \
            -x {params.aligner_index} \
            --rg-id {params.sample} \
            --rg SM:{params.sample} \
            --rg LB:{params.sample} \
            {params.fastq_1_flag} trimming/{params.sample}.R1.trim.fastq.gz \
            {params.fastq_2} \
            -p {resources.cpus} \
            --summary-file alignment/{params.sample}.hisat.summary.txt | \
            samtools view -@ {resources.cpus} -b | \
            samtools sort -T /scratch/samtools_{params.sample} -@ {resources.cpus} -m 128M - 1> alignment/{params.sample}.srt.bam

        # generate BAM index
        samtools index -@ {resources.cpus} alignment/{params.sample}.srt.bam
    """
