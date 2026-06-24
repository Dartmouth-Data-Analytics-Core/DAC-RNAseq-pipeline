####
# Automatic Reference building
# This rule is not run by the default Snakemake target.
# To run these build commands, run snakemake -s Snakefile build_refs
####
rule build_refs:
    """
    Building refs
    """
    params:
        ref_fa = config["reference_fa"],
        ref_gtf = config["annotation_gtf"],
        aligner_name = config["aligner_name"],
        aligner_path = config["aligner_path"],
        picard_build_script = config["picard_build_script"],
        run_rsem = "yes" if RUN_RSEM else "no",
        rsem_prepare_path = config["rsem_prep_ref_path"],
    conda:
        "../../env_config/rna_build_refs.yaml",
    container: "docker://ghcr.io/dartmouth-data-analytics-core/rna_build_refs:2.0"
    resources: cpus="12", maxtime="8:00:00", mem_mb=48000,
    message: "Building references."
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
