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

Annotated TPMs
Transcript-per-million normalized counts were ranked by decreasing mean TPM across samples to identify the top most highly expressed genes. [Note on the top genes, e.g. mitochondrial genes among the most highly expressed.]

Files in the webshare
* QC_Reports (multiqc and raw fastqc)
* Featurecounts/ (counts files):
   * TPMs and raw
   * Annotated_TPMs.xlsx (TPMs with mitochondrial genes flagged)
* PCA Plots
* software_versions.txt (pipeline steps, tools, versions and references)
* methods.txt (methods paragraphs for publications)
{downstream_file_line}

Citing the core
If you use these results in a publication, please cite our grant number in the acknowledgements. Below is some recommended text for citing the core: Research reported in this publication was supported through Geisel School of Medicine at Dartmouth's Center for Quantitative Biology through a grant from the National Institute of General Medical Sciences of the National Institutes of Health under Award Number P20GM130454.

Cheers,
[your name]
