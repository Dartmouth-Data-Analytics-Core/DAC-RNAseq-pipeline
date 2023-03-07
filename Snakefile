import pandas as pd

configfile: "config.yaml"




samples_df = pd.read_table(config["sample_tsv"]).set_index("sample_id", drop=False)
sample_list = list(samples_df['sample_id'])


rule all:
    input:
        expand("umi_tools/{sample}.R1.umi.fastq.gz", sample=sample_list),
        expand("umi_tools/{sample}.R2.umi.fastq.gz", sample=sample_list),
        expand("trimming/{sample}.R1.trim.fastq.gz", sample=sample_list),
        expand("trimming/{sample}.R2.trim.fastq.gz", sample=sample_list),
        expand("trimming/{sample}.cutadapt.report", sample=sample_list),
        expand("alignment/{sample}.srt.bam", sample=sample_list),
        expand("alignment/{sample}.srt.bam.bai", sample=sample_list),
#        expand("alignment/{sample}.Aligned.toTranscriptome.out.bam", sample=sample_list), #commenting out until condition for STAR exists
        expand("alignment/stats/{sample}.srt.bam.flagstat", sample=sample_list),
        expand("dedup/{sample}.dedup.srt.bam", sample=sample_list),
        expand("fastqc/{sample}/{sample}_fastqc.zip", sample=sample_list),
        expand("markdup/{sample}.mkdup.bam", sample=sample_list),
        expand("metrics/picard/{sample}.picard.rna.metrics.txt", sample=sample_list),
#        expand("rsem/{sample}.genes.results", sample=sample_list), #commenting out until condition for STAR/RSEM exists
#        expand("rsem/{sample}.isoforms.results", sample=sample_list), #commenting out until condition for STAR/RSEM exists
        "featurecounts/featurecounts.readcounts.tsv"
    resources: cpus="1", maxtime="1:00:00", memory="4gb",
    output:
        "multiqc_report.html"
    shell: """
        multiqc fastqc alignment markdup metrics featurecounts umi_tools dedup trimming
"""

rule fastqc:
    output: "fastqc/{sample}/{sample}_fastqc.html",
            "fastqc/{sample}/{sample}_fastqc.zip", ##fix this to not hardcode 1M
    params:
        sample = lambda wildcards: wildcards.sample,
        fastqc = config["fastqc_path"],
        fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
    resources: cpus="1", maxtime="1:00:00", memory="4gb",
    shell: """
            mkdir  -p fastqc/{params.sample}
            {params.fastqc} -t 2 -o fastqc/{params.sample} {params.fastq_file_1}
            mv fastqc/{params.sample}/{params.sample}*fastqc.html fastqc/{params.sample}/{params.sample}_fastqc.html
            mv fastqc/{params.sample}/{params.sample}*fastqc.zip fastqc/{params.sample}/{params.sample}_fastqc.zip
"""

rule umi_extract:
    output: "umi_tools/{sample}.R1.umi.fastq.gz",
            "umi_tools/{sample}.R2.umi.fastq.gz",
    params:
        sample = lambda wildcards: wildcards.sample,
        umi_tools = config["umi_tools_path"],
        fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        fastq_file_2 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_2"] if config["layout"]=="paired" else "None",
    resources: cpus="10", maxtime="8:00:00", memory="20gb",
    shell: """
            {params.umi_tools} extract \
                -I {params.fastq_file_2} \
                --bc-pattern=NNNNNNNN \
                --extract-method=string \
                --read2-in={params.fastq_file_1} \
                -S umi_tools/{params.sample}.R2.umi.fastq.gz \
                --read2-out=umi_tools/{params.sample}.R1.umi.fastq.gz
"""

rule umi_dedup:
    input: "alignment/{sample}.srt.bam",
    output: "dedup/{sample}.dedup.srt.bam",
            "dedup/{sample}.dedup.srt.bam.bai",
    params:
        sample = lambda wildcards: wildcards.sample,
        umi_tools = config["umi_tools_path"],
        samtools = config["samtools_path"],
    resources: cpus="10", maxtime="8:00:00", memory="20gb",
    shell: """
            {params.umi_tools} dedup \
                -I alignment/{params.sample}.srt.bam \
                -S dedup/{params.sample}.dedup.srt.bam
            {params.samtools} index -@ 4 dedup/{params.sample}.dedup.srt.bam
"""

rule trimming:
    input: "umi_tools/{sample}.R1.umi.fastq.gz",
           "umi_tools/{sample}.R2.umi.fastq.gz",
    output: "trimming/{sample}.R1.trim.fastq.gz",
            "trimming/{sample}.R2.trim.fastq.gz", ##fix this to not hardcode 1M
            "trimming/{sample}.cutadapt.report",
    params:
        sample = lambda wildcards: wildcards.sample,
        cutadapt = config["cutadapt_path"],
    resources: cpus="1", maxtime="2:00:00", memory="4gb",
    shell: """
            {params.cutadapt} \
                -o trimming/{params.sample}.R1.trim.fastq.gz \
                -p trimming/{params.sample}.R2.trim.fastq.gz \
                umi_tools/{params.sample}.R1.umi.fastq.gz \
                umi_tools/{params.sample}.R2.umi.fastq.gz \
                -m 1 \
                -j 4 \
                -u 6 \
                --max-n 0.8 \
                --trim-n > trimming/{params.sample}.cutadapt.report
"""

# Notes from Owen
# all base qualities in test data 2 it seems, causing cutadapt with -q 20 to  output no reads, might be worth
# changing the test dataset for the paired end reads

# this currently only works with paired end

# its probably wirth us adding an ifelse for instrument, as nextseq trim isnt appropriate
# if sequenced on other instrucments
#--nextseq-trim 20

# it could be worth adding a PolyA trimming option for the 3'-end data, as there is often a lot of
# poolyA in these sequences. for full-length data, I don't really ever notice a difference
# when i trim polyA or not and then align with STAR. I have less epxerience on if thats also the case with
# hisat



rule alignment:
    input:	"trimming/{sample}.R1.trim.fastq.gz",
            "trimming/{sample}.R2.trim.fastq.gz",
    output: "alignment/{sample}.srt.bam",
            "alignment/{sample}.srt.bam.bai",
#            "alignment/{sample}.Aligned.toTranscriptome.out.bam" #Commenting out until condition for STAR aligner exists
    params:
     #   fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
      #  fastq_file_2 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_2"] if config["layout"]=="paired" else "None",
        layout = config["layout"],
        sample = lambda wildcards: wildcards.sample,
        aligner_name = config["aligner_name"],
        aligner = config["aligner_path"],
        aligner_index = config["aligner_index"],
        samtools = config["samtools_path"],
    resources: cpus="10", maxtime="8:00:00", memory="100gb",
    shell: """

           if [ "{params.aligner_name}" = "hisat" ]
            then
           if [ "{params.layout}" = "single" ]
            then
           {params.aligner} -x {params.aligner_index} --rg ID:{params.sample} --rg SM:{params.sample} --rg LB:{params.sample}  -U {params.sample}.R1.trim.fastq.gz -p 12  --summary-file alignment/{params.sample}.hisat.summary.txt | {params.samtools} view -@ 2 -b | {params.samtools} sort -T /scratch/samtools_{params.sample} -@ 4 -m 512M - 1> alignment/{params.sample}.srt.bam
           {params.samtools} index -@ 4 alignment/{params.sample}.srt.bam
            else
            {params.aligner} -x {params.aligner_index} --rg ID:{params.sample} --rg SM:{params.sample} --rg LB:{params.sample}  -1 trimming/{params.sample}.R1.trim.fastq.gz -2  trimming/{params.sample}.R2.trim.fastq.gz  -p 12  --summary-file alignment/{params.sample}.hisat.summary.txt | {params.samtools} view -@ 2 -b | {params.samtools} sort -T /scratch/samtools_{params.sample} -@ 4 -m 512M - 1> alignment/{params.sample}.srt.bam
           {params.samtools} index -@ 4 alignment/{params.sample}.srt.bam

fi


           elif [ "{params.aligner_name}" = "star" ]
            then
    {params.aligner} \
      --genomeDir {params.aligner_index} \
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
      --readFilesIn trimming/{params.sample}.R1.fastq.gz trimming/{params.sample}.R2.fastq.gz \
      --twopassMode Basic \
      --quantMode TranscriptomeSAM \
      --readFilesCommand zcat \
      --outFileNamePrefix alignment/{params.sample}.

    mv alignment/{params.sample}.Aligned.sortedByCoord.out.bam alignment/{params.sample}.srt.bam
    {params.samtools} index -@ 4 alignment/{params.sample}.srt.bam

           else
                echo "aligner not supported."
           fi
   """






rule alignment_metrics:
    input: "alignment/{sample}.srt.bam",
    output: "alignment/stats/{sample}.srt.bam.flagstat",
            "alignment/stats/{sample}.srt.bam.idxstats",
    params:
        samtools = config["samtools_path"],
        sample = lambda wildcards: wildcards.sample,
    resources: cpus="2", maxtime="1:00:00", memory="20gb",
    shell: """
            {params.samtools} flagstat alignment/{params.sample}.srt.bam > alignment/stats/{params.sample}.srt.bam.flagstat
            {params.samtools} idxstats alignment/{params.sample}.srt.bam > alignment/stats/{params.sample}.srt.bam.idxstats
"""

rule picard_markdup:
    input: "dedup/{sample}.dedup.srt.bam",
    output: "markdup/{sample}.mkdup.bam",
    params:
        sample = lambda wildcards: wildcards.sample,
        picard = config['picard_path'],
        java = config['java_path'],
    resources: cpus="4", maxtime="8:00:00", memory="20gb",
    shell: """
            {params.java} -Xmx8G -Xms8G -jar {params.picard} MarkDuplicates I=dedup/{params.sample}.dedup.srt.bam O=markdup/{params.sample}.mkdup.bam M=markdup/{params.sample}.mkdup.log.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 CREATE_INDEX=false  MAX_RECORDS_IN_RAM=4000000 ASSUME_SORTED=true MAX_FILE_HANDLES=768
"""

rule picard_collectmetrics:
    input: "markdup/{sample}.mkdup.bam",
    output: "metrics/picard/{sample}.picard.rna.metrics.txt",
    params:
        sample = lambda wildcards: wildcards.sample,
        picard = config['picard_path'],
        java = config['java_path'],
        flatref = config['picard_refflat'],
        rrna_list = config['picard_rrna_list'],
        strand = config['picard_strand'],
    resources: cpus="2", maxtime="2:00:00", memory="20gb",
    shell: """
        {params.java} -Xmx4G -Xms4G -jar {params.picard} CollectRnaSeqMetrics I=markdup/{params.sample}.mkdup.bam O=metrics/picard/{params.sample}.picard.rna.metrics.txt REF_FLAT={params.flatref} STRAND={params.strand} RIBOSOMAL_INTERVALS={params.rrna_list} MAX_RECORDS_IN_RAM=1000000
"""



rule rsem:
    input: expand("alignment/{sample}.Aligned.toTranscriptome.out.bam", sample=sample_list)
    output: "rsem/{sample}.genes.results",
            "rsem/{sample}.isoforms.results"
    params:
        sample = lambda wildcards: wildcards.sample,
        rsem_path = config['rsem_path'],
        rsem_ref_path = config["rsem_ref_path"]
    resources: sample_load=1
    shell: """
        {params.rsem_path}/rsem-calculate-expression \
        --paired-end \
        --alignments \
        -p 16 \
        --no-bam-output \
        alignment/{params.sample}.Aligned.toTranscriptome.out.bam \
        {params.rsem_ref_path} \
        rsem/{params.sample}
 """

#        --strandedness reverse \
# the test data seems to be unstranded so i removed this option, but we can add a parameter for strandedness in future





rule featurecounts:
    input: expand("dedup/{sample}.dedup.srt.bam", sample=sample_list),
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
    resources: cpus="12", maxtime="2:00:00", memory="8gb",
    shell: """
        {params.featurecounts} -T 32 {params.pair_flag} -s {params.strand}  -a {params.gtf} -o featurecounts/featurecounts.readcounts.raw.tsv {input}
        sed s/"dedup\/"//g featurecounts/featurecounts.readcounts.raw.tsv| sed s/".dedup.srt.bam"//g| tail -n +2 > featurecounts/featurecounts.readcounts.tsv
        Rscript {params.fc_tpm_script} featurecounts/featurecounts.readcounts.tsv
        python {params.fc_ann_script} {params.gtf} featurecounts/featurecounts.readcounts.tsv > featurecounts/featurecounts.readcounts.ann.tsv
        python {params.fc_ann_script} {params.gtf} featurecounts/featurecounts.readcounts_tpm.tsv > featurecounts/featurecounts.readcounts_tpm.ann.tsv
        python {params.fc_ann_script} {params.gtf} featurecounts/featurecounts.readcounts_rpkm.tsv > featurecounts/featurecounts.readcounts_rpkm.ann.tsv
 """
