rule pre_alignment:
    """
    Generate aligner index and verify
    """
    output: "alignment/index_status.txt",
    params:
        layout = config["layout"],
        aligner_name = config["aligner_name"],
        aligner_path = "STAR",
        aligner_index = config["aligner_index"],
        samtools = config["samtools_path"],
    conda:
        "../../env_config/rna_alignment.yaml",
    container: "docker://ghcr.io/dartmouth-data-analytics-core/rna_alignment:2.0"
    resources: cpus="10", maxtime="8:00:00", mem_mb=120000,
    message: "Checking STAR alignment index"
    shell: """
        align_folder="sample_ref/STAR_index"
        if [ ! -d "{params.aligner_index}" ]
            then
                if [ ! -d "$align_folder" ]
                    then
                        mkdir "$align_folder"
                fi
                {params.aligner_path} --runThreadN 16 \
                    --runMode genomeGenerate \
                    --genomeDir "$align_folder" \
                    --genomeFastaFiles {params.aligner_index}.fa \
                    --sjdbGTFfile {params.aligner_index}.chr.gtf \
                    --genomeSAindexNbases 10
            else
                align_folder={params.aligner_index}
        fi
        echo "$align_folder" > alignment/index_status.txt
    """

rule alignment:
    """
    STAR alignment
    """
    input:
        R1_FASTQ_INPUT,
        R2_FASTQ_INPUT if config["layout"] == "paired" else [],
        "alignment/index_status.txt",
    output:
        "alignment/{sample}.srt.bam",
        "alignment/{sample}.srt.bam.bai",
        "alignment/{sample}.Aligned.toTranscriptome.out.bam",
    params:
        layout = config["layout"],
        sample = lambda wildcards: wildcards.sample,
        aligner_name = config["aligner_name"],
        aligner_path = "STAR",
        aligner_index = config["aligner_index"],
        samtools = config["samtools_path"],
        readFilesIn = R1_FASTQ_INPUT + (f" {R2_FASTQ_INPUT}" if config["layout"] == "paired" else '')
    conda:
        "../../env_config/rna_alignment.yaml",
    container: "docker://ghcr.io/dartmouth-data-analytics-core/rna_alignment:2.0"
    resources: cpus="5", maxtime="8:00:00", mem_mb=100000,
    message: "aligning {wildcards.sample} reads with STAR."
    shell: """
        align_folder=`cat alignment/index_status.txt`
                {params.aligner_path} \
                    --genomeDir "$align_folder" \
                    --runThreadN {resources.cpus} \
                    --outSAMunmapped Within \
                    --outFilterType BySJout \
                    --outSAMattributes NH HI AS NM MD \
                    --outSAMtype BAM SortedByCoordinate \
                    --outFilterMultimapNmax 10 \
                    --outFilterMismatchNmax 999 \
                    --outFilterMismatchNoverReadLmax 0.04 \
                    --alignIntronMin 20 \
                    --alignIntronMax 1000000 \
                    --alignMatesGapMax 1000000 \
                    --alignSJoverhangMin 8 \
                    --alignSJDBoverhangMin 1 \
                    --readFilesIn {params.readFilesIn} \
                    --twopassMode Basic \
                    --quantMode TranscriptomeSAM \
                    --readFilesCommand zcat \
                    --outFileNamePrefix alignment/{params.sample}.

        #----- rename output BAM
        mv alignment/{params.sample}.Aligned.sortedByCoord.out.bam alignment/{params.sample}.srt.bam

        #----- index BAM
        samtools index -@ 4 alignment/{params.sample}.srt.bam
    """
