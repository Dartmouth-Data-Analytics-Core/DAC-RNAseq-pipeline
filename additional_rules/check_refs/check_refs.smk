####
# Reference path checking
# This rule is not run by the default Snakemake target.
# To run these checks, run snakemake -s Snakefile check_refs
####
rule check_refs:
    """
    Verify ref paths
    """
    params:
        ref_gtf = config["annotation_gtf"],
        aligner_index = config["aligner_index"],
        aligner_name = config["aligner_name"],
        picard_refflat = config["picard_refflat"],
        picard_rrna_list = config["picard_rrna_list"],
        run_rsem = "yes" if RUN_RSEM else "no",
        rsem_ref = config["rsem_ref_path"],
    resources: cpus="1", maxtime="1:00:00", mem_mb=2000,
    message: "Checking reference paths..."
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
