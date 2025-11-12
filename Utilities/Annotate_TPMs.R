# ===== Annotate TPMs ===== #
# Usage: Rscript Annotate_TPMs.R <featurecounts.readcounts_tpm.ann.tsv> <output path>

#----- Libraries
library(openxlsx)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PARSE ARGUMENTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript Annotate_TPMs.R <input_tpm_file> <output_excel_path>")
}

input_file <- args[1]
output_file <- args[2]

if (!grepl("\\.xlsx$", output_file, ignore.case = TRUE)) {
  stop("❌ Output file must have the '.xlsx' extension.")
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ AND PROCESS INPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
input <- read.csv(input_file, sep = "\t")

#----- Get average TPM and sort
input$Avg_TPM <- rowMeans(input[8:ncol(input)])
input <- input[order(input$Avg_TPM, decreasing = TRUE),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CREATE AND ANNOTATE WORKBOOK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Create workbook
wb <- createWorkbook()
addWorksheet(wb, "TPM_Ordered")
writeData(wb, sheet = "TPM_Ordered", x = input)

#----- Define heatmap-style color scale
heat_style <- createStyle(
  fontColour = "#000000", 
  bgFill = "#FFFFFF"
)

#----- Define numeric column range
start_col <- 8
end_col <- ncol(input)

#----- Create conditional formatting (gradient red scale)
conditionalFormatting(
  wb, 
  sheet = "TPM_Ordered",
  cols = start_col:end_col,
  rows = 2:(nrow(input) + 1),  # +1 for header
  type = "colorScale",
  style = c("forestgreen", "yellow", "red") # white → light red → dark red
)

#----- Highlight mitochondrial genes (Gene.Name starts with MT or mt)
chr_col <- which(colnames(input) == "Chr")
chr_col_letter <- openxlsx::int2col(chr_col)  # Convert to Excel-style letter, e.g., "C"

#----- Highlight ChrM genes
conditionalFormatting(
  wb,
  sheet = "TPM_Ordered",
  cols = chr_col,
  rows = 2:(nrow(input) + 1),
  rule = sprintf('OR(LEFT($%s2,2)="MT",LEFT($%s2,2)="mt")', chr_col_letter, chr_col_letter),
  style = createStyle(fontColour = "red", fgFill = "#FF0000")
)

saveWorkbook(wb, file = output_file, overwrite = TRUE)
cat("✅ Annotated TPM Excel file saved to:", output_file, "\n")


