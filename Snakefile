import pandas as pd

configfile: "config.yaml"




samples_df = pd.read_table(config["sample_tsv"]).set_index("sample_id", drop=False)
sample_list = list(samples_df['sample_id'])


rule all:
    input:
        expand("alignment/{sample}.srt.bam", sample=sample_list),
        expand("alignment/{sample}.srt.bam.bai", sample=sample_list),
        expand("alignment/stats/{sample}.srt.bam.flagstat", sample=sample_list),
        expand("fastqc/{sample}/{sample}_fastqc.zip", sample=sample_list),
        expand("markdup/{sample}.mkdup.bam", sample=sample_list),
        expand("metrics/picard/{sample}.picard.rna.metrics.txt", sample=sample_list),
        "featurecounts/featurecounts.readcounts.tsv"
    output:
        "multiqc_report.html"
    shell: """
        multiqc fastqc alignment markdup metrics featurecounts
"""

rule fastqc:
    output: "fastqc/{sample}/{sample}_fastqc.html",
            "fastqc/{sample}/{sample}_fastqc.zip" ##fix this to not hardcode 1M
    params:
        sample = lambda wildcards:  wildcards.sample,
        fastqc = config["fastqc_path"],
        fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
    shell: """
            mkdir  -p fastqc/{params.sample}
            {params.fastqc} -t 2 -o fastqc/{params.sample} {params.fastq_file_1}
            mv fastqc/{params.sample}/{params.sample}*fastqc.html fastqc/{params.sample}/{params.sample}_fastqc.html
            mv fastqc/{params.sample}/{params.sample}*fastqc.zip fastqc/{params.sample}/{params.sample}_fastqc.zip
"""



rule alignment:
    output: "alignment/{sample}.srt.bam",
            "alignment/{sample}.srt.bam.bai"
    params:
        fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        fastq_file_2 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_2"] if config["layout"]=="paired" else "None",
        layout = config["layout"],
        sample = lambda wildcards:  wildcards.sample,
        aligner_name = config["aligner_name"],
        aligner = config["aligner_path"],
        aligner_index = config["aligner_index"],
        samtools = config["samtools_path"],
    shell: """

           if [ "{params.aligner_name}" = "hisat" ]
            then
           if [ "{params.layout}" = "single" ]
            then
           {params.aligner} -x {params.aligner_index} --rg ID:{params.sample} --rg SM:{params.sample} --rg LB:{params.sample}  -U {params.fastq_file_1} -p 12  --summary-file alignment/{params.sample}.hisat.summary.txt | {params.samtools} view -@ 2 -b | {params.samtools} sort -@ 8 - 1> alignment/{params.sample}.srt.bam 
           {params.samtools} index -@ 4 alignment/{params.sample}.srt.bam
            else
            {params.aligner} -x {params.aligner_index} --rg ID:{params.sample} --rg SM:{params.sample} --rg LB:{params.sample}  -1 {params.fastq_file_1} -2  {params.fastq_file_2}  -p 12  --summary-file alignment/{params.sample}.hisat.summary.txt | {params.samtools} view -@ 2 -b | {params.samtools} sort -@ 8 - 1> alignment/{params.sample}.srt.bam
           {params.samtools} index -@ 4 alignment/{params.sample}.srt.bam

fi


           elif [ "{params.aligner_name}" = "star" ]
            then
    {params.aligner} --genomeDir {params.aligner_index} \
  --runThreadN 18 \
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
  --readFilesIn {params.fastq_file_1} \
  --readFilesCommand zcat \
  --outFileNamePrefix alignment/{params.sample}.

    mv alignment/{params.sample}.Aligned.sortedByCoord.out.bam alignment/{params.sample}.srt.bam
    {params.samtools} index -@ 4 alignment/{params.sample}.srt.bam

           else
                echo "aligner not supported."
           fi
   """ 





        
rule alignment_metrics:
    input: "alignment/{sample}.srt.bam"
    output: "alignment/stats/{sample}.srt.bam.flagstat",
            "alignment/stats/{sample}.srt.bam.idxstats"
    params:
        samtools = config["samtools_path"],
        sample = lambda wildcards:  wildcards.sample,
    shell: """
            {params.samtools} flagstat alignment/{params.sample}.srt.bam > alignment/stats/{params.sample}.srt.bam.flagstat
            {params.samtools} idxstats alignment/{params.sample}.srt.bam > alignment/stats/{params.sample}.srt.bam.idxstats
           """

rule picard_markdup:
    input: "alignment/{sample}.srt.bam"
    output: "markdup/{sample}.mkdup.bam"
    params:
        sample = lambda wildcards:  wildcards.sample,
        picard = config['picard_path']
    shell: """
            java -Xmx4G -Xms4G -jar {params.picard} MarkDuplicates I=alignment/{params.sample}.srt.bam O=markdup/{params.sample}.mkdup.bam M=markdup/{params.sample}.mkdup.log.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 CREATE_INDEX=false  MAX_RECORDS_IN_RAM=100000 ASSUME_SORTED=true MAX_FILE_HANDLES=768
"""

rule picard_collectmetrics:
    input: "markdup/{sample}.mkdup.bam"
    output: "metrics/picard/{sample}.picard.rna.metrics.txt"
    params:
        sample = lambda wildcards:  wildcards.sample,
        picard = config['picard_path'],
        flatref = config['picard_refflat'],
        rrna_list = config['picard_rrna_list'],
        strand = config['picard_strand'],
    shell: """
        java -Xmx4G -Xms4G -jar {params.picard} CollectRnaSeqMetrics I=markdup/{params.sample}.mkdup.bam O=metrics/picard/{params.sample}.picard.rna.metrics.txt REF_FLAT={params.flatref} STRAND={params.strand} RIBOSOMAL_INTERVALS={params.rrna_list} MAX_RECORDS_IN_RAM=1000000
"""


rule featurecounts:
    input:  expand("alignment/{sample}.srt.bam", sample=sample_list)
    output: "featurecounts/featurecounts.readcounts.tsv",
            "featurecounts/featurecounts.readcounts.ann.tsv",
            "featurecounts/featurecounts.readcounts_tpm.ann.tsv",
            "featurecounts/featurecounts.readcounts_rpkm.ann.tsv",
    params:
        featurecounts = config['featurecounts_path'],
        layout = config["layout"],
        pair_flag = "-p" if config["layout"]=="paired" else "",
        strand = config['featurecounts_strand'],
        gtf = config['annotation_gtf'],
        fc_tpm_script = config['featurecounts_rscript'],
        fc_ann_script = config['featurecounts_annscript'],
    shell: """
        {params.featurecounts} -T 32 {params.pair_flag} -s {params.strand}  -a {params.gtf} -o featurecounts/featurecounts.readcounts.raw.tsv {input}
        sed s/"alignment\/"//g featurecounts/featurecounts.readcounts.raw.tsv| sed s/".srt.bam"//g| tail -n +2 > featurecounts/featurecounts.readcounts.tsv
        Rscript {params.fc_tpm_script} featurecounts/featurecounts.readcounts.tsv
        python {params.fc_ann_script} {params.gtf} featurecounts/featurecounts.readcounts.tsv > featurecounts/featurecounts.readcounts.ann.tsv
        python {params.fc_ann_script} {params.gtf} featurecounts/featurecounts.readcounts_tpm.tsv > featurecounts/featurecounts.readcounts_tpm.ann.tsv
        python {params.fc_ann_script} {params.gtf} featurecounts/featurecounts.readcounts_rpkm.tsv > featurecounts/featurecounts.readcounts_rpkm.ann.tsv
 """


