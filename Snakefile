#To do:
#- polyA trimming for 3' data
#- update dataset to have better base qualities
#- add potential handling for rsem quantification after hisat alignment (would require getting a transcriptome alignment from hisat)
# confirm why picard_rrna list needs to be different for hisat

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# setup environment
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import pandas as pd

# set config file
configfile: "config.yaml"

# read in sample data
samples_df = pd.read_table(config["sample_tsv"]).set_index("sample_id", drop=False)
sample_list = list(samples_df['sample_id'])

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define rules
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rule all:
    input:
        expand("trimming/{sample}.R1.trim.fastq.gz", sample=sample_list),
        expand("trimming/{sample}.R2.trim.fastq.gz", sample=sample_list) if config["layout"] == "paired" else [],
        expand("trimming/{sample}.cutadapt.report", sample=sample_list),
        expand("alignment/{sample}.srt.bam", sample=sample_list),
        expand("alignment/{sample}.srt.bam.bai", sample=sample_list),
        expand("alignment/{sample}.srt.bam", sample=sample_list), #commenting out until condition for STAR exists
        expand("alignment/stats/{sample}.srt.bam.flagstat", sample=sample_list),
        expand("markdup/{sample}.mkdup.bam", sample=sample_list),
        expand("metrics/picard/{sample}.picard.rna.metrics.txt", sample=sample_list),
        expand("rsem/{sample}.genes.results", sample=sample_list),
        expand("rsem/{sample}.isoforms.results", sample=sample_list),
        "featurecounts/featurecounts.readcounts.tsv",

    conda:
        "env_config/multiqc.yaml",
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    params:
        layout=config["layout"],
        multiqc=config["multiqc_path"],
        run_rsem=config["run_rsem"],
        aligner_name=config["aligner_name"],

    output:
        "multiqc_report.html"

    shell: """
        #multiqc fastqc alignment markdup metrics featurecounts
        {params.multiqc}  alignment markdup metrics featurecounts

        #remove dummy R2 file (created to meet input rule requirements for rule all:)
        if [ "{params.layout}" = "single" ]
          then
            rm -f trimming/*R2.fastq.gz
        fi

        #remove dummy rsem files (created to meet input rule requirements for rule all:)
        if [ "{params.run_rsem}" = "no" ]
          then
            rm -rf rsem/
        fi

        #remove dummy alignment files (created to meet input rule requirements for rule all:)
        if [ "{params.aligner_name}" = "hisat" ]
          then
            rm -rf alignment/*.Aligned.toTranscriptome.out.bam
        fi
"""


rule trimming:
    output: 
        "trimming/{sample}.R1.trim.fastq.gz",
        "trimming/{sample}.R2.trim.fastq.gz" if config["layout"]=="paired" else [],
        "trimming/{sample}.cutadapt.report"
    params:
        sample = lambda wildcards:  wildcards.sample,
        cutadapt = config["cutadapt_path"],
        fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        fastq_file_2 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_2"] if config["layout"]=="paired" else "None",
        layout=config["layout"],
        nextseq_flag = config["cutadapt_nextseq_flag"]
    conda:
        "env_config/cutadapt.yaml",
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    shell: """
        if  [ "{params.layout}" == "paired" ] 
        then
            {params.cutadapt} \
                -o trimming/{params.sample}.R1.trim.fastq.gz \
                -p trimming/{params.sample}.R2.trim.fastq.gz \
                {params.fastq_file_1} \
                {params.fastq_file_2} \
                -m 1 \
                {params.nextseq_flag} \
                -j {resources.cpus} \
                --max-n 0.8 \
                --trim-n > trimming/{params.sample}.cutadapt.report
        else
            {params.cutadapt} \
                -o trimming/{params.sample}.R1.trim.fastq.gz \
                {params.fastq_file_1} \
                -m 1 \
                {params.nextseq_flag} \
                -j {resources.cpus} \
                --max-n 0.8 \
                --trim-n > trimming/{params.sample}.cutadapt.report
        fi

    """





if config["aligner_name"]=="star":
  rule pre_alignment:
      output: "alignment/index_status.txt",
      params: 
          layout = config["layout"],
          aligner_name = config["aligner_name"],
          aligner = config["aligner_path"],
          aligner_index = config["aligner_index"],
          samtools = config["samtools_path"],
      conda:
          "env_config/alignment.yaml",

      resources: cpus="10", maxtime="8:00:00", mem_mb="120gb",

      shell: """
        align_folder="sample_ref/STAR_index"
        
        if [ ! -d "{params.aligner_index}" ]
            then
                echo "Need to change aligner index"
                if [ ! -d "$align_folder" ]
                    then
                        echo "Current folder does not exist, STAR index needed"
                        mkdir "$align_folder"
                fi
                echo "Begin STAR index creation"
                {params.aligner} --runThreadN 16 \
                    --runMode genomeGenerate \
                    --genomeDir "$align_folder" \
                    --genomeFastaFiles {params.aligner_index}.fa \
                    --sjdbGTFfile {params.aligner_index}.chr.gtf \
                    --genomeSAindexNbases 10
                echo "End STAR index creation"
            else
                echo "STAR index already exists, skip creation"
                align_folder={params.aligner_index}
        fi
        echo "$align_folder" > alignment/index_status.txt
      """

  rule alignment:
      input: "trimming/{sample}.R1.trim.fastq.gz",
             "trimming/{sample}.R2.trim.fastq.gz" if config["layout"] == "paired" else [],
             "alignment/index_status.txt",

      output: "alignment/{sample}.srt.bam",
              "alignment/{sample}.srt.bam.bai",
              "alignment/{sample}.Aligned.toTranscriptome.out.bam",

      params:
          layout = config["layout"],
          sample = lambda wildcards:  wildcards.sample,
          aligner_name = config["aligner_name"],
          aligner = config["aligner_path"],
          aligner_index = config["aligner_index"],
          samtools = config["samtools_path"],
      conda:
          "env_config/alignment.yaml",

      resources: cpus="5", maxtime="8:00:00", mem_mb="100gb",

      shell: """
        align_folder=`cat alignment/index_status.txt`

        if [ "{params.layout}" == "single" ]
            then
                echo "Begin STAR alignment single"
                {params.aligner} --genomeDir "$align_folder" \
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
                        --readFilesIn trimming/{params.sample}.R1.trim.fastq.gz \
                        --twopassMode Basic \
                        --quantMode TranscriptomeSAM \
                        --readFilesCommand zcat \
                        --outFileNamePrefix alignment/{params.sample}.
            else
                echo "Begin STAR alignment paired"
                {params.aligner} \
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
                    --readFilesIn trimming/{params.sample}.R1.trim.fastq.gz trimming/{params.sample}.R2.trim.fastq.gz \
                    --twopassMode Basic \
                    --quantMode TranscriptomeSAM \
                    --readFilesCommand zcat \
                    --outFileNamePrefix alignment/{params.sample}.
        fi
        echo "End STAR alignment"
        # rename output BAM
        mv alignment/{params.sample}.Aligned.sortedByCoord.out.bam alignment/{params.sample}.srt.bam
        
        # index BAM
        {params.samtools} index -@ 4 alignment/{params.sample}.srt.bam
        echo "End alignment rule"
     """



if config["aligner_name"]=="hisat":
  rule alignment:
      input: "trimming/{sample}.R1.trim.fastq.gz",
             "trimming/{sample}.R2.trim.fastq.gz" if config["layout"]=="paired" else [],

      output: "alignment/{sample}.srt.bam",
              "alignment/{sample}.srt.bam.bai",

      params:
          layout = config["layout"],
          sample = lambda wildcards:  wildcards.sample,
          aligner_name = config["aligner_name"],
          hisat2 = config["aligner_path"],
          aligner_index = config["aligner_index"],
          samtools = config["samtools_path"],
      conda:
          "env_config/alignment.yaml",

      resources: cpus="4", maxtime="8:00:00", mem_mb="40gb",

      shell: """
      if [ "{params.layout}" == "single" ]
        then
          # run hisat in single-end mode
          {params.hisat2} \
            -x {params.aligner_index} \
            --rg ID:{params.sample} \
            --rg SM:{params.sample} \
            --rg LB:{params.sample}  \
            -U trimming/{params.sample}.R1.trim.fastq.gz \
            -p {resources.cpus}  \
            --summary-file alignment/{params.sample}.hisat.summary.txt | \
            {params.samtools} view -@ {resources.cpus} -b | \
            {params.samtools} sort -T /scratch/samtools_{params.sample} -@ {resources.cpus} -m 128M - 1> alignment/{params.sample}.srt.bam
        else
          # run hisat in paired-end mode
          {params.hisat2} \
            -x {params.aligner_index} \
            --rg ID:{params.sample} \
            --rg SM:{params.sample} \
            --rg LB:{params.sample}  \
            -1 trimming/{params.sample}.R1.trim.fastq.gz \
            -2 trimming/{params.sample}.R2.trim.fastq.gz \
            -p {resources.cpus}  \
            --summary-file alignment/{params.sample}.hisat.summary.txt | \
            {params.samtools} view -@ {resources.cpus} -b | \
            {params.samtools} sort -T /scratch/samtools_{params.sample} -@ {resources.cpus} -m 128M - 1> alignment/{params.sample}.srt.bam
        fi

        # generate BAM index
        {params.samtools} index -@ {resources.cpus} alignment/{params.sample}.srt.bam

     """



rule alignment_metrics:
    input: "alignment/{sample}.srt.bam",

    output: "alignment/stats/{sample}.srt.bam.flagstat",
            "alignment/stats/{sample}.srt.bam.idxstats",

    params:
        samtools = config["samtools_path"],
        sample = lambda wildcards:  wildcards.sample,
    conda:
        "env_config/samtools.yaml",

    resources: cpus="2", maxtime="8:00:00", mem_mb="20gb",

    shell: """
            {params.samtools} flagstat alignment/{params.sample}.srt.bam > alignment/stats/{params.sample}.srt.bam.flagstat
            {params.samtools} idxstats alignment/{params.sample}.srt.bam > alignment/stats/{params.sample}.srt.bam.idxstats
           """

rule picard_markdup:
    input: "alignment/{sample}.srt.bam",

    output: "markdup/{sample}.mkdup.bam",

    params:
        sample = lambda wildcards:  wildcards.sample,
        picard = config['picard_path'],
        java = config['java_path'],
    conda:
        "env_config/picard.yaml",

    resources: cpus="2", maxtime="30:00", mem_mb="20gb",

    shell: """
            {params.picard} -Xmx2G -Xms2G  \
                 MarkDuplicates \
                I=alignment/{params.sample}.srt.bam \
                O=markdup/{params.sample}.mkdup.bam \
                M=markdup/{params.sample}.mkdup.log.txt \
                OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
                CREATE_INDEX=false  \
                MAX_RECORDS_IN_RAM=4000000 \
                ASSUME_SORTED=true \
                MAX_FILE_HANDLES=768
"""

rule picard_collectmetrics:
    input: "markdup/{sample}.mkdup.bam",

    output: "metrics/picard/{sample}.picard.rna.metrics.txt",

    params:
        sample = lambda wildcards:  wildcards.sample,
        picard = config['picard_path'],
        java = config['java_path'],
        flatref = config['picard_refflat'],
        rrna_list = config['picard_rrna_list'],
        strand = config['picard_strand'],
    conda:
        "env_config/picard.yaml",

    resources: cpus="2", maxtime="8:00:00", mem_mb="20gb",

    shell: """
        {params.picard} -Xmx2G -Xms2G \
             CollectRnaSeqMetrics \
            I=markdup/{params.sample}.mkdup.bam \
            O=metrics/picard/{params.sample}.picard.rna.metrics.txt \
            REF_FLAT={params.flatref} STRAND={params.strand} \
            RIBOSOMAL_INTERVALS={params.rrna_list} \
            MAX_RECORDS_IN_RAM=1000000
"""

if config["run_rsem"]=="yes":
    rule rsem:
        input:  "alignment/{sample}.srt.bam",

        output: "rsem/{sample}.genes.results",
                "rsem/{sample}.isoforms.results"
        params:
            sample = lambda wildcards:  wildcards.sample,
            rsem_path = config['rsem_path'],
            rsem_ref_path = config["rsem_ref_path"],
            rsem_strandedness = config["rsem_strandedness"],
            layout = config["layout"],
        conda:
            "env_config/rsem.yaml",
        resources: cpus="10", maxtime="8:00:00", mem_mb="60gb",

        shell: """
        echo "The BAM length is: $(wc -l alignment/{params.sample}.srt.bam)"
        
        if [ "{params.layout}" == "single" ]
          then
            if ["{params.rsem_path}" == ""]
              then
                rsem-calculate-expression \
                --alignments \
                -p {resources.cpus} \
                --strandedness {params.rsem_strandedness} \
                --no-bam-output \
                alignment/{params.sample}.srt.bam \
                {params.rsem_ref_path} \
                rsem/{params.sample}
              else
                {params.rsem_path}/rsem-calculate-expression \
                --alignments \
                -p {resources.cpus} \
                --strandedness {params.rsem_strandedness} \
                --no-bam-output \
                alignment/{params.sample}.srt.bam \
                {params.rsem_ref_path} \
                rsem/{params.sample}
            fi
        else
          if ["{params.rsem_path}" == ""]
            then
              rsem-calculate-expression \
                --paired-end \
                --alignments \
                -p {resources.cpus} \
                --strandedness {params.rsem_strandedness} \
                --no-bam-output \
                alignment/{params.sample}.srt.bam \
                {params.rsem_ref_path} \
                rsem/{params.sample}
            else
              {params.rsem_path}/rsem-calculate-expression \
                --paired-end \
                --alignments \
                -p {resources.cpus} \
                --strandedness {params.rsem_strandedness} \
                --no-bam-output \
                alignment/{params.sample}.srt.bam \
                {params.rsem_ref_path} \
                rsem/{params.sample}
            fi
        fi
     """
else:
    rule rsem:
      input: "alignment/{sample}.srt.bam",

      output: "rsem/{sample}.genes.results",
              "rsem/{sample}.isoforms.results",

      params:
          sample = lambda wildcards: wildcards.sample,
          layout = config["layout"],

      resources: cpus="10", maxtime="8:00:00", mem_mb="40gb",

      shell: """
        touch rsem/{params.sample}.genes.results
        touch rsem/{params.sample}.isoforms.results
   """

rule featurecounts:
    input:  expand("alignment/{sample}.srt.bam", sample=sample_list),

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
    conda:
        "env_config/featurecounts.yaml",

    resources: cpus="10", maxtime="8:00:00", mem_mb="100gb",

    shell: """
        {params.featurecounts} -T 32 {params.pair_flag} -s {params.strand}  -a {params.gtf} -o featurecounts/featurecounts.readcounts.raw.tsv {input}
        sed s/"alignment\/"//g featurecounts/featurecounts.readcounts.raw.tsv| sed s/".srt.bam"//g| tail -n +2 > featurecounts/featurecounts.readcounts.tsv
        Rscript {params.fc_tpm_script} featurecounts/featurecounts.readcounts.tsv
        python {params.fc_ann_script} {params.gtf} featurecounts/featurecounts.readcounts.tsv > featurecounts/featurecounts.readcounts.ann.tsv
        python {params.fc_ann_script} {params.gtf} featurecounts/featurecounts.readcounts_tpm.tsv > featurecounts/featurecounts.readcounts_tpm.ann.tsv
        python {params.fc_ann_script} {params.gtf} featurecounts/featurecounts.readcounts_rpkm.tsv > featurecounts/featurecounts.readcounts_rpkm.ann.tsv
 """
