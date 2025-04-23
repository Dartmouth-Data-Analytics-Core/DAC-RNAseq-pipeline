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
samples_df = pd.read_table(config["sample_tsv"], sep = ",").set_index("sample_id", drop=False)
sample_list = list(samples_df['sample_id'])

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define rules
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### ADD OUTPUTS CONDITIONALLY IF 
rule all:
    input:
        expand("umi_tools/{sample}.R1.umi.fastq.gz", sample=sample_list),
        expand("umi_tools/{sample}.R2.umi.fastq.gz", sample=sample_list),
        expand("trimming/{sample}.R1.umi.trim.fastq.gz", sample=sample_list),
        expand("trimming/{sample}.R2.umi.trim.fastq.gz", sample=sample_list),
        expand("trimming/{sample}.cutadapt.report", sample=sample_list),
        expand("alignment/{sample}.srt.bam", sample=sample_list),
        expand("alignment/{sample}.srt.bam.bai", sample=sample_list),
        expand("alignment/stats/{sample}.srt.bam.flagstat", sample=sample_list),
        expand("markdup/{sample}.mkdup.bam", sample=sample_list),
        expand("metrics/picard/{sample}.picard.rna.metrics.txt", sample=sample_list),
        expand("metrics/umi_dedup/{sample}.idxstats", sample=sample_list),
        expand("metrics/umi_dedup/{sample}.flagstat", sample=sample_list),
        expand("rsem/{sample}.genes.results", sample=sample_list) if config["run_rsem"] == "yes" else [],
        expand("rsem/{sample}.isoforms.results", sample=sample_list) if config["run_rsem"] == "yes" else [],
        "featurecounts/featurecounts.readcounts.tsv",
        "plots/PCA_Variance_Bar_Plot.png",
        "featurecounts/featurecounts.readcounts.ann.tsv",
        "featurecounts/featurecounts.readcounts_tpm.tsv",
        "featurecounts/featurecounts.readcounts_tpm.ann.tsv",
        "featurecounts/featurecounts.readcounts_rpkm.tsv",
        "featurecounts/featurecounts.readcounts_rpkm.ann.tsv",
        "featurecounts/featurecounts.readcounts_fpkm.tsv",
        "featurecounts/featurecounts.readcounts_fpkm.ann.tsv",
    conda:
        "env_config/multiqc.yaml",
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    params:
        multiqc=config["multiqc_path"],
        run_rsem=config["run_rsem"],
        aligner_name=config["aligner_name"],

    output:
        "multiqc_report.html"

    shell: """
        #multiqc fastqc alignment markdup metrics featurecounts
        {params.multiqc} -p  alignment markdup metrics featurecounts

        #remove dummy alignment files (created to meet input rule requirements for rule all:)
        if [ "{params.aligner_name}" = "hisat" ]
          then
            rm -rf alignment/*.Aligned.toTranscriptome.out.bam
        fi

"""

rule umi_extract:
    input:
        fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        fastq_file_2 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_2"],
    output: 
        "umi_tools/{sample}.R1.umi.fastq.gz",
        "umi_tools/{sample}.R2.umi.fastq.gz"
    params:
        sample = lambda wildcards: wildcards.sample,
        umi_tools = config["umi_tools_path"],
        fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        fastq_file_2 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_2"]
    conda:
        "env_config/umitools.yaml",
    resources: cpus="10", maxtime="8:00:00", memory="20gb",
    shell: """
        # Extract UMIs
            {params.umi_tools} extract \
                -I {params.fastq_file_2} \
                --bc-pattern=NNNNNNNN \
                --extract-method=string \
                --read2-in={params.fastq_file_1} \
                -S umi_tools/{params.sample}.R2.umi.fastq.gz \
                --read2-out=umi_tools/{params.sample}.R1.umi.fastq.gz
    """

rule trimming:
    input: 
        fastq_file_1 = lambda wildcards: "umi_tools/{sample}.R1.umi.fastq.gz",
        fastq_file_2 = lambda wildcards: "umi_tools/{sample}.R2.umi.fastq.gz",
    output: 
        forwardTrim = "trimming/{sample}.R1.umi.trim.fastq.gz",
        reverseTrim = "trimming/{sample}.R2.umi.trim.fastq.gz",
        trimLog = "trimming/{sample}.cutadapt.report" 
    params:
        sample = lambda wildcards:  wildcards.sample,
        cutadapt = config["cutadapt_path"],
        fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        fastq_file_2 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_2"],
        layout=config["layout"],
        nextseq_flag = config["cutadapt_nextseq_flag"],
    conda:
        "env_config/cutadapt.yaml",
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    shell: """
            {params.cutadapt} \
                -o {output.forwardTrim} \
                -p {output.reverseTrim} \
                {input.fastq_file_1} \
                {input.fastq_file_2} \
                -m 1 \
                {params.nextseq_flag} \
                -j {resources.cpus} \
                --max-n 0.8 \
                --trim-n > {output.trimLog}

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
                if [ ! -d "$align_folder" ]
                    then
                        mkdir "$align_folder"
                fi
                {params.aligner} --runThreadN 16 \
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
      input: 
        forward = lambda wildcards: "trimming/{sample}.R1.umi.trim.fastq.gz",
        read2 = lambda wildcards: "trimming/{sample}.R2.umi.trim.fastq.gz",
        status = "alignment/index_status.txt",
      output: "alignment/{sample}.srt.bam",
              "alignment/{sample}.srt.bam.bai",
              "alignment/{sample}.Aligned.toTranscriptome.out.bam",
      params:
          sample = lambda wildcards:  wildcards.sample,
          aligner_name = config["aligner_name"],
          aligner = config["aligner_path"],
          aligner_index = config["aligner_index"],
          samtools = config["samtools_path"],
          readFilesIn = lambda wildcards: "trimming/{sample}.R1.umi.trim.fastq.gz trimming/{sample}.R2.umi.trim.fastq.gz"
      conda:
          "env_config/alignment.yaml",
      resources: cpus="5", maxtime="8:00:00", mem_mb="100gb",
      shell: """
        align_folder=`cat {input.status}`
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
                    --readFilesIn {input.forward} {input.read2} \
                    --twopassMode Basic \
                    --quantMode TranscriptomeSAM \
                    --readFilesCommand zcat \
                    --outFileNamePrefix alignment/{params.sample}.
        
# rename output BAM
        mv alignment/{params.sample}.Aligned.sortedByCoord.out.bam alignment/{params.sample}.srt.bam
        
        # index BAM
        {params.samtools} index -@ 4 alignment/{params.sample}.srt.bam
     """

if config["aligner_name"]=="hisat":
  rule alignment:
      input: 
        forward = lambda wildcards: "trimming/{sample}.R1.umi.trim.fastq.gz",
        read2 = lambda wildcards: "trimming/{sample}.R2.umi.trim.fastq.gz",
      output: "alignment/{sample}.srt.bam",
              "alignment/{sample}.srt.bam.bai",
      params:
          layout = config["layout"],
          sample = lambda wildcards:  wildcards.sample,
          aligner_name = config["aligner_name"],
          hisat2 = config["aligner_path"],
          aligner_index = config["aligner_index"],
          samtools = config["samtools_path"],
          fastq_1_flag = '-1',
          fastq_2 = "-2 trimming/{sample}.R2.umi.trim.fastq.gz"
      conda:
          "env_config/alignment.yaml",
      resources: cpus="4", maxtime="8:00:00", mem_mb="40gb",
      shell: """
          {params.hisat2} \
            -x {params.aligner_index} \
            --rg ID:{params.sample} \
            --rg SM:{params.sample} \
            --rg LB:{params.sample}  \
            {params.fastq_1_flag} {input.forward} \
            {params.fastq_2}  \
            -p {resources.cpus}  \
            --summary-file alignment/{params.sample}.hisat.summary.txt | \
            {params.samtools} view -@ {resources.cpus} -b | \
            {params.samtools} sort -T /scratch/samtools_{params.sample} -@ {resources.cpus} -m 128M - 1> alignment/{params.sample}.srt.bam

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


#----- Run umi_tools dedup if UMIs are true
rule dedup:
    input: 
        srtBam = "alignment/{sample}.srt.bam"
    output: 
        dedup = "markdup/{sample}.mkdup.bam",
        idxStats = "metrics/umi_dedup/{sample}.idxstats",
        flagStat = "metrics/umi_dedup/{sample}.flagstat",
    params:
        sample = lambda wildcards:  wildcards.sample,
        umi_tools = config["umi_tools_path"],
        samtools = config["samtools_path"]
    conda:
        "env_config/umitools.yaml",
    resources: cpus="2", maxtime="30:00", mem_mb="20gb",
    shell: """
            # Deduplicate and index, collect metrics
            umi_tools dedup --method=unique -I {input.srtBam} -S {output.dedup} && \
            samtools index {output.dedup} && \
            samtools idxstats {output.dedup} > {output.idxStats} && \
            samtools flagstat {output.dedup} > {output.flagStat}
            


    """
rule picard:
    input: 
        mkdupBam = "markdup/{sample}.mkdup.bam"
    output: 
        picardMetrics = "metrics/picard/{sample}.picard.rna.metrics.txt",
    params:
        sample = lambda wildcards:  wildcards.sample,
        picard = config['picard_path'],
        flatref = config['picard_refflat'],
        rrna_list = config['picard_rrna_list'],
        strand = config['picard_strand'],
    conda:
        "env_config/picard.yaml",
    resources: cpus="2", maxtime="30:00", mem_mb="20gb",
    shell: """
            
            {params.picard} -Xmx2G -Xms2G \
                CollectRnaSeqMetrics \
                I={input.mkdupBam} \
                O=metrics/picard/{params.sample}.picard.rna.metrics.txt \
                REF_FLAT={params.flatref} STRAND={params.strand} \
                RIBOSOMAL_INTERVALS={params.rrna_list} \
                MAX_RECORDS_IN_RAM=1000000  


    """
    
rule rsem:
    input:  "alignment/{sample}.srt.bam",

    output: "rsem/{sample}.genes.results",
            "rsem/{sample}.isoforms.results"
    params:
        sample = lambda wildcards:  wildcards.sample,
        rsem_calc_exp_path = config['rsem_calc_exp_path'],
        rsem_ref_path = config["rsem_ref_path"],
        rsem_strandedness = config["rsem_strandedness"],
        rsem_paired_flag = '--paired-end' if config["layout"]=='paired' else '',
    conda:
        "env_config/rsem.yaml",
    resources: cpus="10", maxtime="8:00:00", mem_mb="60gb",

    shell: """   
        {params.rsem_calc_exp_path} \
          {params.rsem_paired_flag} \
          --alignments \
          -p {resources.cpus} \
          --strandedness {params.rsem_strandedness} \
          --no-bam-output \
          alignment/{params.sample}.Aligned.toTranscriptome.out.bam \
          {params.rsem_ref_path} \
          rsem/{params.sample}
 """

rule featurecounts:
    input:  expand("alignment/{sample}.srt.bam", sample=sample_list),

    output: "featurecounts/featurecounts.readcounts.tsv",
            "featurecounts/featurecounts.readcounts.ann.tsv",
            "featurecounts/featurecounts.readcounts_tpm.tsv",
            "featurecounts/featurecounts.readcounts_tpm.ann.tsv",
            "featurecounts/featurecounts.readcounts_rpkm.tsv",
            "featurecounts/featurecounts.readcounts_rpkm.ann.tsv",
            "featurecounts/featurecounts.readcounts_fpkm.tsv",
            "featurecounts/featurecounts.readcounts_fpkm.ann.tsv",
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
        python {params.fc_tpm_script} featurecounts/featurecounts.readcounts.tsv {params.layout}
        python {params.fc_ann_script} {params.gtf} featurecounts/featurecounts.readcounts.tsv > featurecounts/featurecounts.readcounts.ann.tsv
        python {params.fc_ann_script} {params.gtf} featurecounts/featurecounts.readcounts_tpm.tsv > featurecounts/featurecounts.readcounts_tpm.ann.tsv
        if [ "{params.layout}" == "single" ]
          then
            python {params.fc_ann_script} {params.gtf} featurecounts/featurecounts.readcounts_rpkm.tsv > featurecounts/featurecounts.readcounts_rpkm.ann.tsv
            touch featurecounts/featurecounts.readcounts_fpkm.tsv
            touch featurecounts/featurecounts.readcounts_fpkm.ann.tsv
        else
            python {params.fc_ann_script} {params.gtf} featurecounts/featurecounts.readcounts_fpkm.tsv > featurecounts/featurecounts.readcounts_fpkm.ann.tsv
            touch featurecounts/featurecounts.readcounts_rpkm.tsv
            touch featurecounts/featurecounts.readcounts_rpkm.ann.tsv
        fi
    """

# The number of genes compared for PCA, chosen by largest variance
num_genes_compared = 500

rule pca_plots:
    input: "featurecounts/featurecounts.readcounts.tsv",

    output:
        "plots/Heatmap_scaled_"+str(num_genes_compared)+"_features.png",
        # there potentially could be more, but this plot must exist. Make sure -p flag has number at least 2 if specified
        "plots/PCA_1_vs_2.png",
        "plots/PCA_Variance_Bar_Plot.png",
        "plots/Gene_Variance_Plot.png",

    params:
        num_genes = num_genes_compared,
        pca_plot_script = config['pca_plot_script'],
        
    conda:
        # uses a subset of the packages that featurecounts does
        "env_config/pcaplot.yaml",

    resources: cpus="1", maxtime="1:00:00", mem_mb="2gb",

    shell: """
        python {params.pca_plot_script} \
        featurecounts/featurecounts.readcounts.tsv \
        plots \
        --genes_considered {params.num_genes} \
        --color_file sample_ref/sample_colors_hex.tsv
    """



####
# Reference path checking
# This rule is not run by the default Snakemake target.
# To run these checks, run snakemake -s Snakefile check_refs
####
rule check_refs:
    params:
        ref_gtf = config["annotation_gtf"],
        aligner_index = config["aligner_index"],
        aligner_name = config["aligner_name"],
        picard_refflat = config["picard_refflat"],
        picard_rrna_list = config["picard_rrna_list"],
        run_rsem = config["run_rsem"],
        rsem_ref = config["rsem_ref_path"],
    resources: cpus="1", maxtime="1:00:00", mem_mb="2gb",
    shell: """   
        
        echo "\nChecking for reference annotation GTF file..."
        if [ -f {params.ref_gtf} ] 
        then
            echo "PASSED -- "{params.ref_gtf}" exists."
        else
            echo "FAILED -- "{params.ref_gtf}" reference annotation GTF is missing!!"
            exit 1
        fi

        if [ {params.aligner_name} == "star" ]
        then
            echo "\nChecking for STAR aligner index files..."
            if [ -f {params.aligner_index}/SA ]
            then
                echo "PASSED -- "{params.aligner_index}" exists."
            else
            echo "FAILED -- "{params.aligner_index}" STAR index directory is missing!!"
            exit 1
            fi
        fi
        if [ {params.aligner_name} == "hisat" ]
        then
            echo "\nChecking for HISAT aligner index files..."
            if [ -f {params.aligner_index}.1.ht2 ]
            then
                echo "PASSED -- "{params.aligner_index}" exists."
            else
            echo "FAILED -- "{params.aligner_index}" HISAT index files not found!!"
            exit 1
            fi
        fi

        if [ {params.run_rsem} == "yes" ]
        then
            echo "\nChecking for RSEM reference files..."
            if [ -f {params.rsem_ref}.n2g.idx.fa ]
            then
                echo "PASSED -- "{params.rsem_ref}" exists."
            else
            echo "FAILED -- "{params.rsem_ref}" RSEM reference index files not found!!"
            exit 1
            fi
        fi

        echo "\nChecking for Picard RefFlat file"
        if [ -f {params.picard_refflat} ] 
        then
            echo "PASSED -- "{params.picard_refflat}" exists."
        else
            echo "FAILED -- "{params.picard_refflat}" Picard RefFlat reference is missing!!"
            exit 1
        fi
        echo "\nChecking for Picard rRNA interval list file"
        if [ -f {params.picard_rrna_list} ] 
        then
            echo "PASSED -- "{params.picard_rrna_list}" exists."
        else
            echo "FAILED -- "{params.picard_rrna_list}" Picard rRNA interval list is missing!!"
            exit 1
        fi



"""




####
# Automatic Reference building
# This rule is not run by the default Snakemake target.
# To run these build commands, run snakemake -s Snakefile build_refs
####
rule build_refs:
    params:
        ref_fa = config["reference_fa"],
        ref_gtf = config["annotation_gtf"],        
        aligner_name = config["aligner_name"],
        aligner_path = config["aligner_path"],
        picard_build_script = config["picard_build_script"],
        run_rsem = config["run_rsem"],
        rsem_prepare_path = config["rsem_prep_ref_path"],
    conda:
          "env_config/build_refs.yaml",
    resources: cpus="12", maxtime="8:00:00", mem_mb="48gb",
    shell: """
            REF_NAME=`basename {params.ref_fa} .fa`
            mkdir -p ref/pipeline_refs
    #        cd ref/pipeline_refs
    #        ln -s {params.ref_fa}
    #        ln -s {params.ref_gtf}

            echo "Building Picard Flat Reference and rRNA Interval List files..."
            chmod +x scripts/picard_ref_builder.sh
            scripts/picard_ref_builder.sh {params.ref_fa} {params.ref_gtf} ref/pipeline_refs/$REF_NAME 

#star 
#hisat
            genome_size=`tail -n1 {params.ref_fa}.fai | awk '{{print $3}}'`
            star_genomeSA_calculation=`echo $genome_size |awk '{{print 14 <((log($1)/log(2))/2)-1?14:((log($1)/log(2))/2)-1}}'`

            if [ {params.aligner_name} == "star" ]
            then
                {params.aligner_path} --runThreadN 12 \
                    --runMode genomeGenerate \
                    --genomeDir ref/pipeline_refs/star_index/$REF_NAME \
                    --genomeFastaFiles {params.ref_fa} \
                    --sjdbGTFfile {params.ref_gtf} \
                    --genomeSAindexNbases $star_genomeSA_calculation
            else
            mkdir ref/pipeline_refs/hisat_index
            {params.aligner_path}-build {params.ref_fa} ref/pipeline_refs/hisat_index/$REF_NAME -p 12
            fi

            if [ {params.run_rsem} == "yes" ]
            then
            mkdir -p ref/pipeline_refs/RSEM_index
            {params.rsem_prepare_path} -p 12 --gtf {params.ref_gtf}  {params.ref_fa} ref/pipeline_refs/RSEM_index/$REF_NAME
            fi


echo "Reference and index building complete."
echo "Paths to use in snakemake config.yaml file"
echo "picard_refflat: \"ref/pipeline_refs/$REF_NAME.refFlat\"" 
echo "picard_rrna_list: \"ref/pipeline_refs/$REF_NAME.rRNA.interval.list\"" 
echo "aligner_index: \"ref/pipeline_refs/{params.aligner_name}_index/$REF_NAME\""
if [ {params.run_rsem} == "yes" ]
then
echo "rsem_ref_path: \"ref/pipeline_refs/RSEM_index/$REF_NAME\""
fi


echo "picard_refflat: \"ref/pipeline_refs/$REF_NAME.refFlat\"" >> ref/pipeline_refs/$REF_NAME.entries.yaml
echo "picard_rrna_list: \"ref/pipeline_refs/$REF_NAME.rRNA.interval.list\"" >> ref/pipeline_refs/$REF_NAME.entries.yaml
echo "aligner_index: \"ref/pipeline_refs/{params.aligner_name}_index/$REF_NAME\"" >> ref/pipeline_refs/$REF_NAME.entries.yaml

if [ {params.run_rsem} == "yes" ]
then
echo "rsem_ref_path: \"ref/pipeline_refs/RSEM_index/$REF_NAME\"" >> ref/pipeline_refs/$REF_NAME.entries.yaml
fi


"""
