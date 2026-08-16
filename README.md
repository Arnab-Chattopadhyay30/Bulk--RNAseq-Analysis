## Bulk RNA-seq Analysis of Down Syndrome

A bulk RNA-seq dataset (GSE305135) was obtained from the NCBI Gene Expression Omnibus (GEO) database. The sequencing data were available in SRA format and were converted into FASTQ files using `fasterq-dump` with the `--split-files` option to generate forward and reverse reads.

Raw reads were initially assessed using FastQC. The reads showed good sequencing quality, with Phred scores above 30 and no significant adapter contamination; therefore, additional trimming was not required. The quality-controlled reads were aligned to the human reference genome (GRCh38) using HISAT2, a splice-aware sequence aligner. The resulting alignments were processed and quantified against the corresponding GTF annotation using featureCounts to generate a gene-level count matrix.

The count matrix was subsequently analyzed in R. Differential gene expression analysis between Down syndrome (DS) and control (CON) samples was performed using DESeq2. Genes with `log2FoldChange > 1` and `adjusted p-value < 0.05` were considered significantly upregulated, while genes with `log2FoldChange < -1` and `adjusted p-value < 0.05` were considered significantly downregulated. Using these criteria, 321 genes were identified as upregulated and 140 genes as downregulated.

PCA was performed to evaluate sample-level variation, showing approximately 94% variation across the first two principal components. A heatmap was generated to visualize differential gene-expression patterns and sample clustering. Volcano plots were also generated using ggplot2 to visualize significantly differentially expressed genes.

Gene Ontology (GO) enrichment analysis revealed that upregulated genes were significantly associated with extracellular matrix organization, extracellular matrix structural constituent activity, and extracellular matrix components. Downregulated genes were enriched in cognition, 3′,5′-cyclic-AMP phosphodiesterase activity, and neuronal cell body-related functions.

KEGG pathway analysis identified significant enrichment of pathways associated with cytoskeletal regulation in muscle cells, among other pathways, providing insights into molecular changes associated with Down syndrome.

To explore potential biomarker candidates, machine-learning approaches using Random Forest and Support Vector Machine (SVM) models were applied. Random Forest analysis identified genes including GBP1, ALPK1, LNCOG, and CPZ based on Mean Decrease Accuracy. Several genes were also identified using SVM, with common candidates between the two approaches including HOXA5, HOXA7, HOXD3, HOXD8, and TMEM88.

Weighted Gene Co-expression Network Analysis (WGCNA) was additionally performed to identify co-expression modules and potential hub genes. CAV1 and STS were identified as hub genes, while CALB1, CAV1, and CYP26A1 were considered potential biomarker candidates.

Overall, this study provides an integrated transcriptomic characterization of Down syndrome by combining differential expression, dimensionality reduction, functional enrichment, pathway analysis, machine learning, and co-expression network analysis. The identified candidate biomarkers provide a basis for further investigation; however, independent cohorts and experimental validation are required to establish their biological and clinical significance.
