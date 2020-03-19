#################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Title: Differential expression analysis annotation & visualziation using DESeq2


# Description:  Example code to perform basic differential expression (DE) analysis, 
# results annotation, and generation of exploratory and diagnostic plots using DE results  
# 
#   - Inputs: 
#         - DESeq2 class object that has been normalized and transformed 
#         - names of sample comparison and indvidual groups 
#    - Outputs: 
#         - volcano plot of DE results 
#         - dispersion estimates plot
#         - independent filtering threshold plot 
#         - Median-average (MA) plot 
#         - .csv files containing annotated DEG results (complete and FDR<5% )
#         - expression heatmap (Z-scores) of hierachical clustering on 
#           significantly differentially expressed genes 
#
#################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tximport)
library(DESeq2)
library(biomaRt)
library(vsn)
library(dplyr)
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(ComplexHeatmap)
library(readr)
library(circlize)
library(EnhancedVolcano)
library(apeglm)

# output directory for input and output files  
dir1 <- ""
# output directory for output figures  
dir2 <- ""

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in dataset
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in DEseq2 object containing raw, filtered, normalized, and transformed counts 
# (generated in explkoratory_data_analysis.R)
load(paste0(dir1, "DESeq.rdata"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# extract DEG analysis results 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# name comparison being assesed (for output file names, e.g. WT_vs_mutant)
comparison <- ""
# name sample groups being compared individually (e.g. WT & mutant)
g1 <- ""
g2 <- ""

# get results for DEG analysis (and order by Pval) by specifying design 
res <- results(dds, alpha = 0.05, 
  contrast = c("design", paste0(g1), paste0(g2)), 
  lfcThreshold = 0, 
  altHypothesis = "greaterAbs")

# remove results w/ NA for adjusted Pval
res <- res[!is.na(res$padj),]

# order by adj Pval 
res_ord <- res[order(res$pvalue),] 

# quick check for how many DEGs with significance @ 5% level in either FC direction 
sum(res$padj < 0.05 & res$log2FoldChange>2, na.rm=TRUE)
sum(res$padj < 0.05 & res$log2FoldChange<-2, na.rm=TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate basic visualizations 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load custom volcano plotting function 
source(paste0("volcano.R"))

# make volcano plot 
volcano_plot(dat = res, 
             out_dir = paste0(dir1), 
             out_name = paste0(comparison, "-volcano.png"), 
             title = paste0(comparison, " - Differential gene expression"),
             alpha=0.05, 
             fc_line = 2)

# plot disperson estimates of DESeq models
ppi=300
png(paste0(dir1, "", comparison, "-dispersion-estimates.png"), width=5*ppi, height=5*ppi, res=ppi)
plotDispEsts(dds, main=paste0(comparison, ""))
dev.off()

# check indepenet filtering threshold 
metadata(res_ord)$alpha
metadata(res_ord)$filterThreshold
png(paste0(dir1, "", comparison, "-independent-filtering-threshold.png"), width=5*ppi, height=5*ppi, res=ppi)
plot(metadata(res_ord)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter", main = paste0(comparison, "Independent filtering threshold"))
lines(metadata(res_ord)$lo.fit, col="red")
abline(v=metadata(res_ord)$filterTheta)
dev.off()

# calculate shrunken fold change estimate
resLFC <- lfcShrink(dds, 
                    coef=paste0(resultsNames(dds)[which(resultsNames(dds)==paste0("design_", g1, "_vs_", g2))]), 
                    type="apeglm")

# plot shrunken estimates against unshrunken using MA plot 
png(paste0(dir1, "", comparison, "-MA-plot.png"), width=6*ppi, height=10*ppi, res=ppi)
par(mfrow=c(2,1))
plotMA(res_ord, ylim=c(-6,6), main = paste0(comparison, "Raw Log2 Fold change"))
plotMA(resLFC, ylim=c(-6,6), main = paste0(comparison, "Raw Log2 Fold change"))
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# annotate gene features in results 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# select mart to use from BioMart databases   
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# store attributes of this Mart for viewing in df 
ensembl_genes_atts <- listAttributes(mart)

# obtain HUGO gene symbol using Ensembl gene IDs 
anno <- getBM(filters = "ensembl_gene_id", 
              attributes = c("ensembl_gene_id", "hgnc_symbol"),
              values = rownames(res_ord),
              mart = mart)

# order annotation results by Ensembl gene ID  
anno_ord <- anno[order(anno$ensembl_gene_id),]

# reorder DEG results by ENSG ID  
res_ord_2 <- res_ord[order(rownames(res_ord)),]

# identify any genes mapping to more than one Ensembl ID 
dup <- anno_ord[which(duplicated(anno_ord$ensembl_gene_id)),]

# drop duplicates
ind1 <- which(duplicated(anno_ord$ensembl_gene_id))
anno_ord <- anno_ord[-ind1,]

# remove entries from DEG results that don't have an annotation identified 
res_ord_3 <- res_ord_2[-which(!rownames(res_ord_2) %in% anno_ord$ensembl_gene_id),]

# check results and annotatioins are in same order, by Ensembl ID, before merging 
all(rownames(res_ord_3) == anno_ord$ensembl_gene_id)

# add HUGO gene ID to results 
res_ord_3$gene <- NA
res_ord_3$gene <- anno_ord$hgnc_symbol
res_ord_3$ENSG_gene <- anno_ord$ensembl_gene_id

# reorder results by pval 
res_ord_4 <- res_ord_3[order(res_ord_3$pvalue),]

# subset @ 5% adjusted pval sig. level 
res_ord_4_FDR_05 <- res_ord_4[res_ord_4$padj<0.05,]

# write out to csv 
write.csv(as.data.frame(res_ord_4), file=paste0(dir1, comparison, "-Diff-exp_results.csv"))
write.csv(as.data.frame(res_ord_4_FDR_05), file=paste0(dir1, comparison, "-Diff-exp_results-FDR_0.05.csv"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# hierachcial clustering using sig. DE genes 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# subset transformed count values for sig. DEGs 
DEGenes <- rownames(assay(rld))[rownames(assay(rld)) %in% rownames(res_ord_4_FDR_05)]
DEGenes <- DEGenes[order(DEGenes)]

# set up gene expression matrix 
mat1 <- assay(rld)[DEGenes,]
rownames(mat1) <- res_ord_4_FDR_05$gene[order(rownames(res_ord_4_FDR_05))]

# scale matrix  
mat_scaled = t(apply(mat1, 1, scale))

# set up colors for heatmap 
colramp = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
cols1 <- brewer.pal(12, "Paired")
cols2 <- brewer.pal(9, "Greens")

# set up variable w/ levels for sample group to be used in annotation bar 
genos <- as.character(dds$design)
genos <- factor(genos, levels = c(paste0(g2), paste0(g1)))
ha1 = HeatmapAnnotation(Genotype = genos, 
                        col = list(Genotype = c("WT" = cols1[7], "mutant" = cols1[8])), show_legend = TRUE)

# set up column annotation labels (samples) for bottom annotation of individual samples 
ha = columnAnnotation(x = anno_text(rownames(dds@colData), 
                                    which="column", rot = 45, 
                                    gp = gpar(fontsize = 10)))

# generate heatmap object 
ht1 = Heatmap(mat_scaled, name = "Expression", col = colramp, 
              bottom_annotation = c(ha),
              top_annotation = c(ha1),
              show_row_names = FALSE)

# plot the heatmap 
ppi=300
png(paste0(dir1, "heatmap_DE_genes_", comparison, ".png"), width=7.5*ppi, height=7.5*ppi, res=ppi)
draw(ht1, row_title = "Genes", column_title = paste0("DE genes, FDR < 0.05, ", comparison))
dev.off()

