#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Title: Testing reading xlsx file into R to extract sample names
# Description: Need to link metadata with sample sheet
#
# Author: Mike Martinez
# Lab: GDSC
# Project: Development
# Date created: 6/24/2025
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD LIBRARIES AND SET PATHS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----- Libraries
library(openxlsx)
library(data.table)
library(stringr)

#----- README
# Run this code from within the directory containing your fastq files and specify argument 1 as "."

#----- Set command args
args <- commandArgs(trailingOnly = TRUE)

#----- Set args
workingDir <- args[1]
slimsMeta <- args[2]
sampleSheet <- args[3]


#----- Set working directory
wd <- workingDir

message("#--------------------------------------------------#")
message("Starting meta data linkage...")
message(paste0("Working directory: ", workingDir))
message(paste0("Slims Metadata: ", slimsMeta))
message(paste0("Sample sheet: ", sampleSheet))
message("#--------------------------------------------------#\n")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN THE EXCEL FILE TO SEE HOW IT'S FORMATTED
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

message("Reading in slims metadata...\n")
meta <- read.xlsx(slimsMeta,
                  startRow = 3,
                  colNames = TRUE,
                  rowNames = FALSE)

#----- Now, read in a mock sample sheet
message("Reading in sample sheet.../n")
samples <- fread(sampleSheet)
nSamples <- nrow(samples)
message(paste0(nSamples, " file names found!"))
filenames <- samples$fastq_1

message("Cleaning filenames...\n")
sampleBase <- str_split_fixed(basename(filenames), "_S", 2)[, 1]
sampleBase <- gsub("./", "", sampleBase)

# Extract the matched group
message("Cleaned sample names: ")
print(sampleBase)
message("\n")

#----- Add sample base to sample sheet
samples$base <- sampleBase

#----- Ensure that they are in the same order
#----- Strip DIL1 from IDs
meta$Id <- gsub("_DIL1", "", meta$Id)
rownames(meta) <- meta$Id
meta <- meta[samples$base,]
message("Checking all cleaned file names are present in slims metadata...")
all(rownames(meta) == samples$base)

#----- Match sample base with metadata
message("Linking external names...\n")
samples$sample_id <- meta$External.Name


message("Writing final sample sheet...")
if ("fastq_2" %in% colnames(samples)) {
  sampleDF <- data.frame(sample_id = samples$sample_id,
                         fastq_1 = samples$fastq_1,
                         fastq_2 = samples$fastq_2)
  write.csv(sampleDF, file = "sample_fastq_list_paired.csv", row.names = FALSE, quote = FALSE)
  message("Sample sheet saved to sample_fastq_list_paired.csv")
} else {
  sampleDF <- data.frame(sample_id = samples$sample_id,
                         fastq_1 = samples$fastq_1)
  write.csv(sampleDF, file = "sample_fastq_list_single.csv", row.names = FALSE, quote = FALSE)
  message("Sample sheet saved to sample_fastq_list_single.csv")
}




