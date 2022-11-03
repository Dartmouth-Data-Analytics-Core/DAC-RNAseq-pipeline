
#Specify inputs and outputs
REF_FASTA_IN=$1
REF_GTF_IN=$2
OUT_PREFIX=$3

#Specify software paths
gtfToGenePred_path=scripts/gtfToGenePred
samtools_path=samtools

#Specify what biotype string to search for within GTF file.  Usually either "gene_type" or "gene_biotype"
gene_type_string="gene_biotype"
#Space-separated list of rRNA types to search
types_to_search="rRNA Mt_rRNA"

echo -e "\nBuilding refFlat file and rRNA interval list file for use with Picard CollectRNASeqMetrics, using these inputs:"
echo -e "Reference fasta input:     ${REF_FASTA_IN}"
echo -e "Gene annotation GTF input: ${REF_GTF_IN}\n"

#Generate refFlat file
date
echo "Generating refFlat file with gtfToGenePred..."
$gtfToGenePred_path -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons  $REF_GTF_IN /dev/stdout| awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > ${OUT_PREFIX}.refFlat
echo -e "Done.\n"


#Generate rRNA interval list

#Get the name and length of each reference contig and put it in sequence header format
date
echo "Generating reference sequence header index..."
$samtools_path faidx $REF_FASTA_IN
echo -e "@HD\tVN:1.6\tSO:coordinate" > .sequence_header_part.tmp
cat ${REF_FASTA_IN}.fai| awk '{print "@SQ\tSN:"$1"\tLN:"$2}' >> .sequence_header_part.tmp
echo -e "Done.\n"


#Get the transcripts labelled as rRNA from the .gtf filei, write them in interval list format
date
echo "Extracting rRNA intervals from GTF..."
for biotype in $types_to_search
do
echo "Searching for genes of biotype: ${biotype}..."
grep "$gene_type_string \"$biotype" $REF_GTF_IN |  awk '$3 == "transcript"' |  cut -f1,4,5,7,9 | perl -lane ' /transcript_id "([^"]+)"/ or die "no transcript_id on $."; print join "\t", (@F[0,1,2,3], $1)' | sort -k1V -k2n -k3n >> .rrna_gene_intervals_parts.tmp
done
sort .rrna_gene_intervals_parts.tmp | uniq | sort -k1V -k2n -k3n > .rrna_gene_intervals_combined.tmp
echo -e "Done.\n"

#Concatenate the sequence headers and rRNA intervals
date
echo "Building rRNA Interval List..."
cat .sequence_header_part.tmp .rrna_gene_intervals_combined.tmp > ${OUT_PREFIX}.rRNA.interval.list
echo -e "Done.\n"

date
echo "File Generation Complete."
echo "The full paths to the generated refFlat and rRNA interval lists are:"
realpath ${OUT_PREFIX}.refFlat
realpath ${OUT_PREFIX}.rRNA.interval.list
