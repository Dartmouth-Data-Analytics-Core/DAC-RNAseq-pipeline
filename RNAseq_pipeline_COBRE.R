#-------------------------------------#
# [1] RNA-seq pipeline for COBRE #
# From fastq to raw counts           #
#-------------------------------------#

rm(list = ls())
library("glue")
whoseData <- "Meisam_PattabiramanLab"
myRaw <- "/dartfs-hpc/rc/lab/G/GSR_Active/Labs/Pattabiraman/Meisam_RNAseq_10-14-19/"
myQC <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/rnaseq_pipeline/frank/fastqc/"
myTri <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/rnaseq_pipeline/frank/trim/"
myAli <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/rnaseq_pipeline/frank/alignment/"
myExp <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/rnaseq_pipeline/frank/rawcounts/"
myTmp <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/rnaseq_pipeline/frank/tmp/"
#--
givenInd <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/genomic_references/mouse/STAR/GRCm38_index"
givenGen <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/genomic_references/mouse/ensembl-annotations/Mus_musculus.GRCm38.97.gtf"
#--
mytar <- sampleID

#-------------------------
for(i in 1:length(mytar)){
	#-------
	myfiles <- dir(myRaw, pattern = mytar[i])
	
	sample_id_1 <- unique(sapply(strsplit(myfiles[1], "_"), "[", 1))
	sample_id_2 <- unique(sapply(strsplit(myfiles[1], "_"), "[", 2))
	sample_id <- paste(sample_id_1, sample_id_2, sep = "_")
	my_R1 <- myfiles[grep("_R1_", myfiles)]
	my_R2 <- myfiles[grep("_R2_", myfiles)]

	myoutf1 <- glue(myTmp, whoseData, "_", mytar[i], ".sp")
	conOut <- file(myoutf1, "w")
	tmp <- "#!/bin/bash"
	writeLines(tmp, conOut)
	writeLines("\n", conOut)
	#-----------------
	# Quality control
	#-----------------
	comm <- "/isi/whitfield/ywang/SofWar/FastQC/fastqc myFASTQ --outdir="
	tmp <- gsub("myFASTQ", glue(myRaw, my_R1), comm)
	tmp <- glue(tmp, myQC)
	writeLines(tmp, conOut)
	writeLines("\n", conOut)
	#--
	tmp <- gsub("myFASTQ", glue(myRaw, my_R2), comm)
	tmp <- glue(tmp, myQC)
	writeLines(tmp, conOut)
	writeLines("\n", conOut)
	#-------
	# Trime
	#-------
	comm_1 <- "cutadapt -A 'A{76}' -g 'T{76}'"
	comm_2 <- "-o R1_fastq_out -p R2_fastq_out R1_fastq_in R2_fastq_in"
	comm_3 <- "-m 5 --nextseq-trim=20 > Sample_ID_L000.report"  #-j 8 need python 3

	my_R1_in <- paste(myRaw, sample_id, "_R1_001.fastq.gz", sep = "")
	my_R2_in <- paste(myRaw, sample_id, "_R2_001.fastq.gz", sep = "")

	my_R1_out <- paste(myTri, sample_id, "_R1_001.fastq.gz", sep = "")
	my_R2_out <- paste(myTri, sample_id, "_R2_001.fastq.gz", sep = "")

	comm_2 <- gsub("R1_fastq_out", my_R1_out, comm_2)
	comm_2 <- gsub("R2_fastq_out", my_R2_out, comm_2)
	comm_2 <- gsub("R1_fastq_in", my_R1_in, comm_2)
	comm_2 <- gsub("R2_fastq_in", my_R2_in, comm_2)

	tag <- paste(myTri, sample_id, sep = "")
	comm_3 <- gsub("Sample_ID", tag, comm_3)

	comm <- paste(comm_1, comm_2, comm_3, sep = " ")
	tmp <- comm
	writeLines(tmp, conOut)
	writeLines("\n", conOut)
	#----------
	# Alignment
	#----------
	myInd <- givenInd
	myGen <- givenGen

	comm_1 <- "/isi/whitfield/ywang/SofWar/STAR-2.7.3a/source/STAR --quantMode TranscriptomeSAM --genomeDir myind --sjdbGTFfile mygene"
	comm_2 <- "--runThreadN 4 --twopassMode Basic --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1"
	comm_3 <- "--readFilesIn trimed_R1_fastq trimed_R2_fastq --readFilesCommand zcat --outFileNamePrefix Sample_ID."

	comm_1 <- gsub("myind", myInd, comm_1)
	comm_1 <- gsub("mygene", myGen, comm_1)

	myTrim_R1 <- paste(myTri, sample_id, "_L000_R1.fastq.gz", sep = "")
	myTrim_R2 <- paste(myTri, sample_id, "_L000_R2.fastq.gz", sep = "")
	tag <- paste(myAli, sample_id, sep = "")

	comm_3 <- gsub("trimed_R1_fastq", myTrim_R1, comm_3)
	comm_3 <- gsub("trimed_R2_fastq", myTrim_R2, comm_3)
	comm_3 <- gsub("Sample_ID", tag, comm_3)

	comm <- paste(comm_1, comm_2, comm_3, sep = " ")
	tmp <- comm
	writeLines(tmp, conOut)
	writeLines("\n", conOut)
	#------------------
	# QC for BAM file
	#------------------
	comm <- "java -Xmx32G -jar /isi/whitfield/ywang/SofWar/picard.jar CollectRnaSeqMetrics I=myInBam O=myOutBam REF_FLAT=myTxt STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=interval_list MAX_RECORDS_IN_RAM=1000000"
	myBam_in <- ""
	myBam_out <- ""
	myFlat <- ""
	comm <- gsub("myInBam", myBam_in, comm)
	comm <- gsub("myOutBam", myBam_out, comm)
	comm <- gsub("myTxt", myFlat, comm)
	#tmp <- comm
	#writeLines(tmp, conOut)
	#writeLines("\n", conOut)
	#--
	comm <- "java -Xmx32G -jar /isi/whitfield/ywang/SofWar/picard.jar MarkDuplicates I=myInBam O=myOutBam M=myTxt OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 CREATE_INDEX=false"
	myBam_in <- ""
	myBam_out <- ""
	myFlat <- ""
	comm <- gsub("myInBam", myBam_in, comm)
	comm <- gsub("myOutBam", myBam_out, comm)
	comm <- gsub("myTxt", myFlat, comm)
	#tmp <- comm
	#writeLines(tmp, conOut)
	#writeLines("\n", conOut)
	#-----------
	# Expression
	#-----------
	myIn <- paste(myAli, sample_id, ".Aligned.toTranscriptome.out.bam", sep = "")
	tag <- paste(myExp, sample_id, sep = "")

	comm <- "rsem-calculate-expression --paired-end --alignments --strandedness reverse -p 4 input_bam mygene Sample_ID"

	comm <- gsub("input_bam", myIn, comm)
	comm <- gsub("mygene", myGen, comm)
	comm <- gsub("Sample_ID", tag, comm)
	
	tmp <- comm
	writeLines(tmp, conOut)
	close(conOut)
}

mysub <- paste(myTmp, "submit.sp", sep="")
conOut <- file(mysub, "w")
curLine <- paste("qsub -A=NCCC -l nodes=1:ppn=10 -l walltime=20:0:0 ", whoseData, "_", mytar, ".sp", sep = "")
writeLines(curLine, conOut)
close(conOut)

comm <- paste("chmod u+x ", mysub, sep="")
system(comm)
"bash submit.sp"


###################
# [2] Multiqc for every step #
###################
"https://github.com/ewels/MultiQC/blob/master/docs/usage.md"
