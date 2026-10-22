Hi [client first name],
The results for these {n_samples} samples are available at the link below:
[webshare link]
If you'd like to chat about any downstream analyses, let us know, we are happy to help.
These data were generated with the Watchmaker mRNA library prep kit (poly(A)-selected, stranded) and processed with the GDSC-RNAseq Pipeline{pipeline_version}, which includes the following:

{steps}

Summary
{summary}

Raw FASTQ QC
{raw_qc}

Alignment QC
{alignment_qc}

{coverage} {strand} {content} {rrna} Duplication ranged {dup_range}; high duplication is expected at low input or high sequencing depth.

{thresholds}

PCA
{pca}

Files in the webshare
* QC_Reports (multiqc and raw fastqc)
* Featurecounts/ (counts files):
   * TPMs and raw
   * Annotated_TPMs.xlsx (TPMs with mitochondrial genes flagged)
* PCA Plots
{downstream_file_line}

Cheers,
[your name]
