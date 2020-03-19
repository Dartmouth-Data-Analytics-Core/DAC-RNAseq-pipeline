#################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Title: exploratory anlysis of RNAseq data with DESeq2 (from HTSeqCount)


# Description: Example code for exploratory data analysis of RNA-seq data quantified with HTSeq-Count, 
# using several approaches to assess sample similairty and groupings 
# 
#   - Inputs: covariate/sample sheet .csv, .htseq-counts quantifications
#   - Construction of DESeq2 data set object, removal of low count features
#   - DESeq2 normalization, rlog transformation 
#   - exploratory analyses: 
#         - sample distance clustering 
#         - principal components analysis 
#         - hierachical clustering 
#    - Outputs: 
#         - pre-processed, normalized, & transformed DESeq2 class object 
#         - heatmap of hierachically clustered sample distances 
#         - feature variance distribution 
#         - barplot of PC variance explained 
#         - PCA plot (PC 1 vs 2)
#         - expression heatmap (Z-scores) of hierachical clustering on top X varaible features 
#
#################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(dplyr)
library(ggplot2)
library(tximport)
library(DESeq2)
library(biomaRt)
library(vsn)
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

# directory of htseq-count & covariate/sample sheet files  
dir1 <- ""
# output directory for files  
dir2 <- ""
# output directory for figures  
dir3 <- ""

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load data 
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in sample info (must include variable called 'design')
samples <- read.csv(paste0(dir1, "covariates.csv"), stringsAsFactors = FALSE)

# order entries by sample name 
samples <- samples[order(samples$sample),]

# assign sample vars as factor 
samples$siRNA <- factor(samples$siRNA, levels = c("WT", "mutant"))

# get file names 
sampleFiles <- list.files(dir1, pattern = "*htseq-counts")

# order files by sample name 
sampleFiles <- sampleFiles[order(sampleFiles)]

# sample file names must match sample names and order of them in covariate/sample sheet file 
all(sampleFiles == samples$sample)

# set df up containing these
sampleTable <- data.frame(sampleName = samples$sample,
                          fileName = sampleFiles, 
                          design = samples$siRNA)

# read into DEseq2 object 
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = dir2, 
                                  design = ~design)

# run the DEseq2 analysis 
dds <- DESeq(dds)

# drop genes with low counts 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
rld <- rlog(dds, blind = FALSE)

# save DESeq object ( to be used in differential expression analysis)
save(dds, rld, file = paste0(dir2, "DESeq2.rdata"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample distance 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###### sample distance (euclidean) plot ###### 
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( colnames(rld))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# produce plot
png(paste0(dir3, "euclidean_distance_matrix.png"), width=7*ppi, height=6*ppi, res=ppi)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, main = "Euclidean distance")
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Principal components analysis (PCA)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate gene expression level variance between samples (to determine how 
# many genes to use for dimensionality reduction and clustering)
var <- rev(rowVars(assay(rld))[order(rowVars(assay(rld)))])

# plot variance for genes accross samples
png(paste0(dir3, "per_gene_variance.png"), width=5*ppi, height=5*ppi, res=ppi)
plot(var, las = 1, main="Sample gene expression variance", xlab = "Gene", ylab = "Variance")
abline(v=500, col="red") ; abline(v=500, col="green")
dev.off()

# modify variable feature number to be used in PCA and hierachical clutering based on no. of most variable features 
var_feature_n <- 500

# perform PCA and order by variance 
rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(var_feature_n, length(rv)))]
pca <- prcomp(t(assay(rld)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
names(percentVar)[1:10] <- c("PC1", "PC2", "PC3", "PC4", "PC5", 
                             "PC6", "PC7", "PC8", "PC19", "PC10")

# plot variance for top 10 PCs 
png(paste0(dir3, "PCA_variance-top_", var_feature_n, "_genes.png"), width=5*ppi, height=5*ppi, res=ppi)
barplot(percentVar[1:10], col = "indianred", las = 1, ylab = "Percent Variance explained", cex.lab = 1.2)
dev.off()

# construct data frame w/ PC loadings and add sample labels 
pca_df <- as.data.frame(pca$x)
pca_df$group <- dds$design
pca_df$sample_ids <- colnames(dds)

# add colors for plotting to df 
pca_df$col <- NA
for(i in 1:length(levels(pca_df$group))){
  ind1 <- which(pca_df$group == levels(pca_df$group)[i])
  pca_df$col[ind1] <- i
}

# plot PC1 vs PC2
png(paste0(dir3, "pca_top_", var_feature_n, "_genes_PC1vsPC2.png"), width=7*ppi, height=7*ppi, res=ppi)
plot(pca_df[, 1], pca_df[, 2], 
     xlab = paste0("PC1 (", (round(percentVar[1], digits=3)*100), "% variance)"), 
     ylab = paste0("PC2 (", (round(percentVar[2], digits=3)*100), "% variance)"),
     main=("PC1 vs PC2 for ", var_feature_n, " most variable genes"),
     pch=16, cex=1.35, cex.lab=1.3, cex.axis = 1.15, las=1, 
     panel.first = grid(),
     col=pca_df$col)
text((pca_df[, 2])~(pca_df[, 1]), labels = pca_df$sample_ids, cex=0.6, font=2, pos=4)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# hierachical clustering  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# select rlog normalized count values for top X no. of variable genes 
topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), var_feature_n)

# set up gene expression matrix 
mat1 <- assay(rld)[topVarGenes,]

# scale matrix 
mat_scaled = t(apply(mat1, 1, scale))

# set up colors for heatmap 
col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
cols1 <- brewer.pal(12, "Paired")

# set up variable w/ levels for sample group to be used in annotation bar 
genos <- as.character(dds$design)
genos <- factor(genos, levels = c("WT", "mutant"))
ha1 = HeatmapAnnotation(Genotype = genos, 
                        col = list(Genotype = c("WT" = cols1[1], "mutant" = cols1[2]), 
                          show_legend = TRUE)

# set up column annotation labels (samples) for bottom annotation of individual samples 
ha = columnAnnotation(x = anno_text(colnames(dds), 
                                    which="column", rot = 45, 
                                    gp = gpar(fontsize = 10)))

# generate heatmap object 
ht1 = Heatmap(mat_scaled,name = "Expression", col = col, 
              bottom_annotation = c(ha),
              top_annotation = c(ha1),
              show_row_names = FALSE)

# plot the heatmap & write to file 
png(paste0(dir3, "heatmap_top_", var_feature_n, "_genes.png"), width=7.5*ppi, height=5*ppi, res=ppi)
draw(ht1, row_title = "Genes", column_title = "Hierachical clustering, ", var_feature_n, " most variable genes")
dev.off()
