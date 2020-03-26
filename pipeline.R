DAC_RNAseq_process <- function(Lab, FastqRaw, SamNames, SeqMethod, AlignInd, AlignRef, PicardInt, PicardRef, QuantRef, CondaEnv, OutputFolder){
	myQC <- paste(OutputFolder, "fastqc/", sep = "")
	myTri <- paste(OutputFolder, "trim/", sep = "")
	myAli <- paste(OutputFolder, "alignment/", sep = "")
	myExp <- paste(OutputFolder, "quantification/", sep = "")
	myTmp <- paste(OutputFolder, "tmp/", sep = "")
	writeLines("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
	writeLines("#This RNA-seq pipeline is funded by the COBRE grand (grand number '1P20GM130454')")
	writeLines("#If you are using this pipeline for processing your own data and preparing publications,")
	writeLines("#please CITE us (grand number '1P20GM130454'). Thank you for your cooperation!!!!")
	writeLines("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
	#-------------------------
	mytar <- SamNames
	#-------------------------
	for(i in 1:length(mytar)){
		#-------
		myfiles <- dir(FastqRaw, pattern = mytar[i])
	
		sample_id_1 <- unique(sapply(strsplit(myfiles[1], "_"), "[", 1))
		sample_id_2 <- unique(sapply(strsplit(myfiles[1], "_"), "[", 2))
		sample_id <- paste(sample_id_1, sample_id_2, sep = "_")
		my_R1 <- myfiles[grep("_R1_", myfiles)]
		if(SeqMethod == "fullLength"){
			my_R2 <- myfiles[grep("_R2_", myfiles)]
		}
	
		myoutf1 <- paste(myTmp, Lab, "_", mytar[i], ".sp", sep = "")
		conOut <- file(myoutf1, "w")
		tmp <- "#!/bin/bash"
		writeLines(tmp, conOut)
		writeLines("\n", conOut)
		writeLines(tmp, conOut)
		writeLines(paste("conda activate ", CondaEnv, sep = ""), conOut)
		#-----------------
		# Quality control
		#-----------------
		writeLines("echo 'QC step'", conOut)
		comm <- "fastqc myFASTQ --outdir="
		tmp <- gsub("myFASTQ", paste(FastqRaw, my_R1, sep = ""), comm)
		tmp <- paste(tmp, myQC, sep = "")
		writeLines(tmp, conOut)
		writeLines("\n", conOut)
		#--
		if(SeqMethod == "fullLength"){
			tmp <- gsub("myFASTQ", paste(FastqRaw, my_R2, sep = ""), comm)
			tmp <- paste(tmp, myQC, sep = "")
			writeLines(tmp, conOut)
			writeLines("\n", conOut)
		}
		#------
		# Trim
		#------
		writeLines("echo 'Trim step'", conOut)
		if(SeqMethod == "fullLength"){
			comm_1 <- "cutadapt -A 'A{76}' -g 'T{76}'"
			comm_2 <- "-o R1_fastq_out -p R2_fastq_out R1_fastq_in R2_fastq_in"
			comm_3 <- "-m 5 --nextseq-trim=20 > Sample_ID.report"

			my_R1_in <- paste(FastqRaw, sample_id, "_R1_001.fastq.gz", sep = "")
			my_R2_in <- paste(FastqRaw, sample_id, "_R2_001.fastq.gz", sep = "")

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
		}else if(SeqMethod == "3Prime"){
			comm_1 <- "cutadapt -a 'A{76}' -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -n 2"
			comm_2 <- "-o R1_fastq_out R1_fastq_in"
			comm_3 <- "-m 5 --nextseq-trim=20 > Sample_ID.report"
		
			my_R1_in <- paste(FastqRaw, sample_id, "_R1_001.fastq.gz", sep = "")
			my_R1_out <- paste(myTri, sample_id, "_R1_001.fastq.gz", sep = "")
			
			comm_2 <- gsub("R1_fastq_out", my_R1_out, comm_2)
			comm_2 <- gsub("R1_fastq_in", my_R1_in, comm_2)

			tag <- paste(myTri, sample_id, sep = "")
			comm_3 <- gsub("Sample_ID", tag, comm_3)
			
			comm <- paste(comm_1, comm_2, comm_3, sep = " ")
			tmp <- comm
		}else{
			tmp <- "Wrong Sequencing Method was referred."
		}
		writeLines(tmp, conOut)
		writeLines("\n", conOut)
		#----------
		# Alignment
		#----------
		writeLines("echo 'Alignment step'", conOut)
		writeLines("date", conOut)
		if(SeqMethod == "fullLength"){
			comm_1 <- "STAR --quantMode TranscriptomeSAM --genomeDir myind --sjdbGTFfile mygene"
			comm_2 <- "--runThreadN 4 --twopassMode Basic --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1"
			comm_3 <- "--readFilesIn trimed_R1_fastq trimed_R2_fastq --readFilesCommand zcat --outFileNamePrefix Sample_ID."

			comm_1 <- gsub("myind", AlignInd, comm_1)
			comm_1 <- gsub("mygene", AlignRef, comm_1)

			myTrim_R1 <- paste(myTri, sample_id, "_R1_001.fastq.gz", sep = "")
			myTrim_R2 <- paste(myTri, sample_id, "_R2_001.fastq.gz", sep = "")
			tag <- paste(myAli, sample_id, sep = "")

			comm_3 <- gsub("trimed_R1_fastq", myTrim_R1, comm_3)
			comm_3 <- gsub("trimed_R2_fastq", myTrim_R2, comm_3)
			comm_3 <- gsub("Sample_ID", tag, comm_3)

			comm <- paste(comm_1, comm_2, comm_3, sep = " ")
			tmp <- comm
		}else if(SeqMethod == "3Prime"){
			comm_1 <- "STAR --genomeDir myind --sjdbGTFfile mygene"
			comm_2 <- "--runThreadN 4 --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 10 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1"
			comm_3 <- "--readFilesIn trimed_R1_fastq --readFilesCommand zcat --outFileNamePrefix Sample_ID."
		
			comm_1 <- gsub("myind", AlignInd, comm_1)
			comm_1 <- gsub("mygene", AlignRef, comm_1)

			myTrim_R1 <- paste(myTri, sample_id, "_R1_001.fastq.gz", sep = "")
			tag <- paste(myAli, sample_id, sep = "")

			comm_3 <- gsub("trimed_R1_fastq", myTrim_R1, comm_3)
			comm_3 <- gsub("Sample_ID", tag, comm_3)

			comm <- paste(comm_1, comm_2, comm_3, sep = " ")
			tmp <- comm
		}else{
			tmp <- "Wrong Sequencing Method was referred."
		}
		writeLines(tmp, conOut)
		writeLines("\n", conOut)
		#------------------
		# QC for BAM file
		#------------------
		writeLines("echo 'Alignment QC step'", conOut)
		if(SeqMethod == "fullLength"){
			comm <- "java -Xmx32G -jar ~/.conda/pkgs/picard-2.21.7-0/share/picard-2.21.7-0/picard.jar CollectRnaSeqMetrics I=myInBam O=myOut REF_FLAT=myTxt STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=interval_list MAX_RECORDS_IN_RAM=1000000"
		}else{
			comm <- "java -Xmx32G -jar ~/.conda/pkgs/picard-2.21.7-0/share/picard-2.21.7-0/picard.jar CollectRnaSeqMetrics I=myInBam O=myOut REF_FLAT=myTxt STRAND=FIRST_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=interval_list MAX_RECORDS_IN_RAM=1000000"
		}
		myBam_in <- paste(myAli, sample_id, ".Aligned.sortedByCoord.out.bam", sep = "")
		my_out <- paste(myAli, sample_id, ".RNA_Metrics", sep = "")
		comm <- gsub("myInBam", myBam_in, comm)
		comm <- gsub("myOut", my_out, comm)
		comm <- gsub("myTxt", PicardRef, comm)
		comm <- gsub("interval_list", PicardInt, comm)
		tmp <- comm
		writeLines(tmp, conOut)
		writeLines("\n", conOut)
		#--
		comm <- "java -Xmx32G -jar ~/.conda/pkgs/picard-2.21.7-0/share/picard-2.21.7-0/picard.jar MarkDuplicates I=myInBam O=myOutBam M=myTxt OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 CREATE_INDEX=false"
		myBam_in <- paste(myAli, sample_id, ".Aligned.sortedByCoord.out.bam", sep = "")
		myBam_out <- paste(myAli, sample_id, ".marked_duplicates.bam", sep = "")
		mytxt <- paste(myAli, sample_id, ".marked_dup_metrics.txt", sep = "")
		comm <- gsub("myInBam", myBam_in, comm)
		comm <- gsub("myOutBam", myBam_out, comm)
		comm <- gsub("myTxt", mytxt, comm)
		tmp <- comm
		writeLines(tmp, conOut)
		writeLines("\n", conOut)
		#-----------
		# Expression
		#-----------
		writeLines("echo 'Quantification step'", conOut)
		if(SeqMethod == "fullLength"){
			myIn <- paste(myAli, sample_id, ".Aligned.toTranscriptome.out.bam", sep = "")
			tag <- paste(myExp, sample_id, sep = "")
			comm <- "rsem-calculate-expression --paired-end --alignments --strandedness reverse -p 8 input_bam mygene Sample_ID"

			comm <- gsub("input_bam", myIn, comm)
			comm <- gsub("mygene", QuantRef, comm)
			comm <- gsub("Sample_ID", tag, comm)
			tmp <- comm
		}else if(SeqMethod == "3Prime"){
			myIn <- paste(myAli, sample_id, ".Aligned.sortedByCoord.out.bam", sep = "")
			tag <- paste(myExp, sample_id, sep = "")
		
			comm <- "htseq-count -f bam -s yes input_bam mygene > myout.htseq-counts"

			comm <- gsub("input_bam", myIn, comm)
			comm <- gsub("mygene", QuantRef, comm)
			comm <- gsub("myout", tag, comm)
	
			tmp <- comm
		}else{
			tmp <- "Wrong Sequencing Method was referred."
		}
	
		writeLines(tmp, conOut)
		writeLines("conda deactivate", conOut)
		close(conOut)
	}

	mysub <- paste(myTmp, "submit.sp", sep="")
	conOut <- file(mysub, "w")
	writeLines("echo 'This RNA-seq pipeline is funded by the COBRE grand (grand number '1P20GM130454')'", conOut)
	writeLines("echo 'If you are using this pipeline for processing your own data and preparing publications,' ", conOut)
	writeLines("echo 'please cite us. Thank you for your cooperation!!!!'                                                      ", conOut)
	writeLines("\n", conOut)
	curLine <- paste("mksub -l nodes=1:ppn=10 -l walltime=20:0:0 ", Lab, "_", mytar, ".sp", sep = "")
	writeLines(curLine, conOut)
	close(conOut)

	comm <- paste("chmod u+x ", mysub, sep="")
	system(comm)
}

cleanFolders <- function(myTar){
	myQC <- paste(myTar, "fastqc/", sep = "")
	myTri <- paste(myTar, "trim/", sep = "")
	myAli <- paste(myTar, "alignment/", sep = "")
	myExp <- paste(myTar, "quantification/", sep = "")
	myTmp <- paste(myTar, "tmp/", sep = "")
	
	comm <- paste("rm -rf ", myQC, sep = "")
	system(comm)
	comm <- paste("rm -rf ", myTri, sep = "")
	system(comm)
	comm <- paste("rm -rf ", myAli, sep = "")
	system(comm)
	comm <- paste("rm -rf ", myExp, sep = "")
	system(comm)
	comm <- paste("rm -rf ", myTmp, sep = "")
	system(comm)
	comm <- paste("mkdir ", myQC, sep = "")
	system(comm)
	comm <- paste("mkdir ", myTri, sep = "")
	system(comm)
	comm <- paste("mkdir ", myAli, sep = "")
	system(comm)
	comm <- paste("mkdir ", myExp, sep = "")
	system(comm)
	comm <- paste("mkdir ", myTmp, sep = "")
	system(comm)
}
