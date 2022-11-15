#! /usr/bin/env Rscript

## functions for rpkm and tpm
## from https://gist.github.com/slowkow/c6ab0348747f86e2748b#file-counts_to_tpm-r-L44
## from https://www.biostars.org/p/171766/
rpkm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(counts) * 1e9
}

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

## read table from featureCounts output
args <- commandArgs(T)
tag <- tools::file_path_sans_ext(args[1])
ftr.cnt <- read.table(args[1], sep="\t", stringsAsFactors=FALSE,
  header=TRUE, check.names=FALSE)
library(dplyr)
library(tidyr)

ftr.rpkm <- ftr.cnt %>%
  gather(sample, cnt, 7:ncol(ftr.cnt)) %>%
  group_by(sample) %>%
  mutate(rpkm=rpkm(cnt, Length)) %>%
  select(-cnt) %>%
  spread(sample, rpkm)
write.table(ftr.rpkm, file=paste0(tag, "_rpkm.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

ftr.tpm <- ftr.cnt %>%
  gather(sample, cnt, 7:ncol(ftr.cnt)) %>%
  group_by(sample) %>%
  mutate(tpm=tpm(cnt, Length)) %>%
  select(-cnt) %>%
  spread(sample, tpm)
write.table(ftr.tpm, file=paste0(tag, "_tpm.tsv"), sep="\t", row.names=FALSE, quote=FALSE)


