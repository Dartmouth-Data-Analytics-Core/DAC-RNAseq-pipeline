#!/usr/bin/env Rscript
# Usage: Rscript rlog_transform.R <raw_counts.tsv> <rlog_out.tsv>
# Input must be RAW counts (genes x samples); rlog does the median-of-ratios normalisation itself.
suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
counts <- round(as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE)))
coldata <- data.frame(row.names = colnames(counts), dummy = rep(1, ncol(counts)))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ 1)
rld <- tryCatch(rlog(dds, blind = TRUE), error = function(e) {
  message("rlog failed (", conditionMessage(e), "); falling back to log2(normalized counts + 1)")
  normTransform(estimateSizeFactors(dds))
})
write.table(data.frame(Gene = rownames(rld), assay(rld), check.names = FALSE),
            args[2], sep = "\t", quote = FALSE, row.names = FALSE)