#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GDSC Bulk RNASeq Pipieline v2.0
#
# Pipeline for the preprocessing and QC of Bulk RNASeq data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import pandas as pd
from snakemake.utils import validate

#----- set config file
configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

#----- read in sample data
samples_df = pd.read_csv(config["sample_csv"]).set_index("sample_id", drop=False)
sample_list = list(samples_df['sample_id'])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RULE INCLUSION AND CONDITIONALS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#----- Set additional rules based on params
REMOVE_rRNA = config.get("remove_rRNA", False)
if REMOVE_rRNA:
    include: "additional_rules/ribodetector/ribodetector.smk"
RUN_RUSTQC = config.get("run_rustqc", False)
if RUN_RUSTQC:
    include: "additional_rules/RustQC/rustqc.smk"
RUN_RSEM = config.get("run_rsem", False)
if RUN_RSEM:
    include: "additional_rules/rsem/rsem.smk"

#----- Conditionally set inputs
if REMOVE_rRNA:
    R1_FASTQ_INPUT = "ribodetector/{sample}/{sample}.nonrrna.1.fq.gz" 
    R2_FASTQ_INPUT = "ribodetector/{sample}/{sample}.nonrrna.2.fq.gz" if config["layout"]=="paired" else None
else:
    R1_FASTQ_INPUT = "trimming/{sample}.R1.trim.fastq.gz"
    R2_FASTQ_INPUT = "trimming/{sample}.R2.trim.fastq.gz" if config["layout"]=="paired" else None

#----- Conditionally set aligner
if config["aligner_name"] == "star":
    include: "additional_rules/alignment/star.smk"
elif config["aligner_name"] == "hisat2":
    include: "additional_rules/alignment/hisat.smk"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RULE ALL INPUTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

all_inputs = []

#----- Trimming
all_inputs += expand("trimming/{sample}.R1.trim.fastq.gz", sample=sample_list)
all_inputs += expand("trimming/{sample}.cutadapt.report", sample=sample_list)
if config["layout"] == "paired":
    all_inputs += expand("trimming/{sample}.R2.trim.fastq.gz", sample=sample_list)

#----- rRNA filtering (optional)
if REMOVE_rRNA:
    all_inputs += expand("ribodetector/{sample}/{sample}.nonrrna.1.fq.gz", sample=sample_list)
    all_inputs += ["ribodetector/rrna_norrna_pct_mqc.tsv"]
    if config["layout"] == "paired":
        all_inputs += expand("ribodetector/{sample}/{sample}.nonrrna.2.fq.gz", sample=sample_list)

#----- Alignment
all_inputs += expand("alignment/{sample}.srt.bam", sample=sample_list)
all_inputs += expand("alignment/{sample}.srt.bam.bai", sample=sample_list)

#----- Deduplication
all_inputs += expand("markdup/{sample}.mkdup.bam", sample=sample_list)

#----- Picard metrics
all_inputs += expand("metrics/picard/{sample}.picard.rna.metrics.txt", sample=sample_list)

#----- QC
if RUN_RUSTQC:
    all_inputs += expand("qc/{sample}", sample=sample_list)
else:
    all_inputs += expand("alignment/stats/{sample}.srt.bam.flagstat", sample=sample_list)

#----- RSEM isoform quantification (optional)
if RUN_RSEM:
    all_inputs += expand("rsem/{sample}.genes.results", sample=sample_list)
    all_inputs += expand("rsem/{sample}.isoforms.results", sample=sample_list)

#----- featureCounts
all_inputs += [
    "featurecounts/featurecounts.readcounts.tsv",
    "featurecounts/featurecounts.readcounts.ann.tsv",
    "featurecounts/featurecounts.readcounts_tpm.tsv",
    "featurecounts/featurecounts.readcounts_tpm.ann.tsv",
    "featurecounts/featurecounts.readcounts_rpkm.tsv",
    "featurecounts/featurecounts.readcounts_rpkm.ann.tsv",
    "featurecounts/featurecounts.readcounts_fpkm.tsv",
    "featurecounts/featurecounts.readcounts_fpkm.ann.tsv",
]

#----- PCA
all_inputs += [
    "plots/PCA_top_PC1_vs_PC2.png",
    "plots/PCA_top_PCA_variance_bar.png",
]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PIPELINE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#----- Main pipeline execution
rule all:
    input:
        all_inputs
    conda:
        "env_config/rna_multiqc.yaml",
    container: "docker://ghcr.io/dartmouth-data-analytics-core/rna_multiqc:2.0"
    resources: cpus="10", maxtime="2:00:00", mem_mb=60000,
    params:
        layout=config["layout"],
        multiqc=config["multiqc_path"],
        aligner_name=config["aligner_name"],
        multiqc_dirs=(
            ("qc " if RUN_RUSTQC else "") +
            "alignment markdup metrics featurecounts" +
            (" ribodetector" if REMOVE_rRNA else "")
        )
    output:
        report = "multiqc_report.html",
    shell: """
        #-----multiqc fastqc alignment markdup metrics featurecounts
        multiqc \
            -c multiqc_config.yaml \
            -p \
            {params.multiqc_dirs} \
            -n {output.report}
        
        #-----remove dummy R2 file (created to meet input rule requirements for rule all:)
        # also remove dummy rpkm and fpkm files from featurecounts normalization
        if [ "{params.layout}" = "single" ]
          then
            rm -f trimming/*R2.fastq.gz
            rm -f featurecounts/featurecounts.readcounts_fpkm.tsv
            rm -f featurecounts/featurecounts.readcounts_fpkm.ann.tsv
          else
            rm -f featurecounts/featurecounts.readcounts_rpkm.tsv
            rm -f featurecounts/featurecounts.readcounts_rpkm.ann.tsv
        fi

        #-----remove dummy alignment files (created to meet input rule requirements for rule all:)
        if [ "{params.aligner_name}" = "hisat" ]
          then
            rm -rf alignment/*.Aligned.toTranscriptome.out.bam
        fi

"""

#----- Rule to execute trimming
rule trimming:
    """
    Read trimming
    """
    output: 
        trim_R1 = "trimming/{sample}.R1.trim.fastq.gz",
        trim_R2 = "trimming/{sample}.R2.trim.fastq.gz" if config["layout"]=="paired" else [],
        trim_log = "trimming/{sample}.cutadapt.report"
    params:
        sample = lambda wildcards:  wildcards.sample,
        cutadapt = config["cutadapt_path"],
        fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        fastq_file_2 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_2"] if config["layout"]=="paired" else "None",
        layout=config["layout"],
        nextseq_flag = config["cutadapt_nextseq_flag"]
    conda:
        "env_config/rna_cutadapt.yaml",
    container: "docker://ghcr.io/dartmouth-data-analytics-core/rna_cutadapt:2.0"
    resources: cpus="10", maxtime="2:00:00", mem_mb=60000,
    message: "Trimming {wildcards.sample} reads with cutadapt."
    shell: """
        if  [ "{params.layout}" == "paired" ] 
        then
            cutadapt \
                -o {output.trim_R1} \
                -p {output.trim_R2} \
                {params.fastq_file_1} \
                {params.fastq_file_2} \
                -m 1 \
                {params.nextseq_flag} \
                -j {resources.cpus} \
                --max-n 0.8 \
                --trim-n > {output.trim_log}
        else
            cutadapt \
                -o {output.trim_R1} \
                {params.fastq_file_1} \
                -m 1 \
                {params.nextseq_flag} \
                -j {resources.cpus} \
                --max-n 0.8 \
                --trim-n > {output.trim_log}
        fi

    """

#----- Rule to get alignment metrics
rule alignment_metrics:
    """
    Samtools metrics (when rustQC is false)
    """
    input: "alignment/{sample}.srt.bam",
    output: "alignment/stats/{sample}.srt.bam.flagstat",
            "alignment/stats/{sample}.srt.bam.idxstats",
    params:
        samtools = config["samtools_path"],
        sample = lambda wildcards:  wildcards.sample,
    conda:
        "env_config/rna_samtools.yaml",
    container: "docker://ghcr.io/dartmouth-data-analytics-core/rna_samtools:2.0"
    resources: cpus="2", maxtime="8:00:00", mem_mb=20000,
    message: "Running flagstats and idxstats QC for {wildcards.sample} with Samtools."
    shell: """
        samtools flagstat alignment/{params.sample}.srt.bam > alignment/stats/{params.sample}.srt.bam.flagstat
        samtools idxstats alignment/{params.sample}.srt.bam > alignment/stats/{params.sample}.srt.bam.idxstats
    """

#----- Rule to mark duplicates
rule picard_markdup:
    """
    Deduplication
    """
    input: 
        sorted_bam = "alignment/{sample}.srt.bam",
    output: 
        mkdups = "markdup/{sample}.mkdup.bam",
        picard_log = "markdup/{sample}.mkdup.log.txt"
    params:
        sample = lambda wildcards:  wildcards.sample,
        picard = config['picard_path'],
    conda:
        "env_config/rna_picard.yaml",
    container: "docker://ghcr.io/dartmouth-data-analytics-core/rna_picard:2.0"
    resources: cpus="2", maxtime="4:00:00", mem_mb=20000,
    message: "Deduplicating reads for {wildcards.sample} reads with Picard."
    shell: """
            picard -Xmx2G -Xms2G  \
                MarkDuplicates \
                I={input.sorted_bam} \
                O={output.mkdups} \
                M={output.picard_log} \
                OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
                CREATE_INDEX=false  \
                MAX_RECORDS_IN_RAM=4000000 \
                ASSUME_SORTED=true \
                MAX_FILE_HANDLES=768
"""

#----- Rule to collect RNASeq metrics
rule picard_collectmetrics:
    """
    Picard metrics
    """
    input: 
        mkdup_bam = "markdup/{sample}.mkdup.bam",
    output: 
        picard_metrics = "metrics/picard/{sample}.picard.rna.metrics.txt",
    params:
        sample = lambda wildcards:  wildcards.sample,
        picard = config['picard_path'],
        flatref = config['picard_refflat'],
        rrna_list = config['picard_rrna_list'],
        strand = config['picard_strand'],
    conda:
        "env_config/rna_picard.yaml",
    container: "docker://ghcr.io/dartmouth-data-analytics-core/rna_picard:2.0"
    resources: cpus="2", maxtime="8:00:00", mem_mb=20000,
    message: "Collecting RNASeq metrics for {wildcards.sample} reads with Picard."
    shell: """
        picard -Xmx2G -Xms2G \
            CollectRnaSeqMetrics \
            I={input.mkdup_bam} \
            O={output.picard_metrics} \
            REF_FLAT={params.flatref} \
            STRAND={params.strand} \
            RIBOSOMAL_INTERVALS={params.rrna_list} \
            MAX_RECORDS_IN_RAM=1000000

        awk -F'\\t' 'BEGIN{{OFS="\\t"; state=0}}
        /^## METRICS CLASS/ {{print; state=1; next}}
        state==1 {{print; ncols=NF; state=2; next}}
        state==2 {{
            if (NF < ncols) {{
                pad = ncols - NF
                line = $0
                for (i = 0; i < pad; i++) line = line "\\t"
                print line
            }} else {{
                print
            }}
            state=0; next
        }}
        {{print}}' {output.picard_metrics} > {output.picard_metrics}.tmp && mv {output.picard_metrics}.tmp {output.picard_metrics}
    """


#----- Rule to count features
rule featurecounts:
    """
    Count reads
    """
    input:  
        expand("alignment/{sample}.srt.bam", sample=sample_list),
    output: 
        "featurecounts/featurecounts.readcounts.tsv",
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
        "env_config/rna_featurecounts.yaml",
    container: "docker://ghcr.io/dartmouth-data-analytics-core/rna_featurecounts:2.0"
    resources: cpus="10", maxtime="8:00:00", mem_mb=100000,
    message: "Counting reads with FeatureCounts."
    shell: """
        featureCounts \
            -T 32 \
            {params.pair_flag} \
            -s {params.strand} \
            -a {params.gtf} \
            -o featurecounts/featurecounts.readcounts.raw.tsv \
            {input}

        #----- Clean
        sed 's|alignment/||g' featurecounts/featurecounts.readcounts.raw.tsv| sed 's|.srt.bam||g'| tail -n +2 > featurecounts/featurecounts.readcounts.tsv
        
        #----- Annotate
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

#----- Rule to plot PCA
rule pca_plots:
    """
    PCA
    """
    input: 
        "featurecounts/featurecounts.readcounts.tsv",
    output:
        "plots/PCA_top_PC1_vs_PC2.png",
        "plots/PCA_top_PCA_variance_bar.png",
    params:
        pca_plot_script = config['pca_plot_script'],   
    conda:
        "env_config/rna_pcaplot.yaml",
    container: "docker://ghcr.io/dartmouth-data-analytics-core/rna_pcaplot:2.0"
    resources: cpus="1", maxtime="1:00:00", mem_mb=2000,
    message: "Running PCA"
    shell: """
        python {params.pca_plot_script} \
        featurecounts/featurecounts.readcounts.tsv \
        plots
    """

#----- CI/CD directives
include: "additional_rules/check_refs/check_refs.smk"
include: "additional_rules/build_refs/build_refs.smk"

