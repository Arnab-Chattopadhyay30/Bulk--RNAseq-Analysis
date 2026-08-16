## Bulk RNA-seq Transcriptomic Analysis of Down Syndrome

A bulk RNA-seq dataset (GSE305135) was obtained from the NCBI GEO database. The dataset was available in SRA format and was converted into FASTQ files using `fasterq-dump` with the `--split-files` option to generate forward and reverse reads. The quality of the raw reads was assessed using FastQC. The reads showed no apparent adapter contamination and had Phred quality scores above 30; therefore, trimming was not required.

The high-quality reads were aligned to the human reference genome (GRCh38) using HISAT2, a splice-aware alignment tool. The resulting alignment files were processed using SAMtools, and gene-level read quantification was performed using the corresponding GTF annotation file to generate a count matrix for downstream analysis.

The count matrix was subsequently analyzed in R. Differential gene expression analysis between the Down syndrome (DS) and control (CON) conditions was performed using DESeq2. Genes with a log2 fold change > 1 and adjusted p-value < 0.05 were considered upregulated, while genes with a log2 fold change < -1 and adjusted p-value < 0.05 were considered downregulated. Using these criteria, 321 genes were identified as upregulated and 140 genes as downregulated. A volcano plot was generated using ggplot2 to visualize the differentially expressed genes.

Principal Component Analysis (PCA) was performed to evaluate variation between the DS and CON samples. The first two principal components captured approximately 94% of the observed variation, showing clear separation between the conditions. A heatmap was also generated to visualize gene expression patterns and assess clustering across samples.

Gene Ontology (GO) enrichment analysis was performed to identify enriched Biological Process, Molecular Function, and Cellular Component categories. Upregulated genes were primarily enriched in extracellular matrix organization, extracellular matrix structural constituent, and extracellular matrix-related categories. Downregulated genes were enriched in cognition, 3’,5’-cyclic-AMP phosphodiesterase activity, and neuronal cell body categories.

KEGG pathway analysis identified pathways associated with cytoskeletal regulation in muscle cells among the significantly enriched pathways in the DS samples. Machine learning approaches, including Random Forest and Support Vector Machine (SVM), were subsequently applied to identify potential biomarker genes. Random Forest analysis identified GBP1, ALPK1, LNCOG, and CPZ based on Mean Decrease Accuracy. Comparison of the Random Forest and SVM results identified common candidate genes, including HOXA5, HOXA7, HOXD3, HOXD8, and TMEM88.

Weighted Gene Co-expression Network Analysis (WGCNA) was also performed to identify gene co-expression patterns and hub genes. CAV1 and STS were identified as hub genes, while CALB1, CAV1, and CYP26A1 were highlighted as potential biomarker candidates.

Overall, this study provides an integrated transcriptomic analysis of Down syndrome by combining differential expression, dimensionality reduction, functional enrichment, pathway analysis, machine learning, and co-expression network analysis. The analysis identified differentially expressed genes, enriched biological pathways, hub genes, and candidate biomarkers that may provide insights into the molecular landscape of Down syndrome. These candidate biomarkers require further validation using independent cohorts and experimental studies to establish their biological functions and potential clinical relevance.

