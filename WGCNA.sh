library(DESeq2)
library(WGCNA)
library(dplyr)
library(dynamicTreeCut)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

###############################
# Read data
###############################

count <- read.csv("final_counts_annotated.csv",
                  header = TRUE,
                  check.names = FALSE)

coldata <- read.csv("metadata_new.csv",
                    header = TRUE)

coldata_modified <- coldata %>%
  dplyr::select(1, 8)

rownames(coldata_modified) <- coldata_modified$Run

rownames(count) <- make.unique(as.character(count$gene_name))
count <- count[, colnames(count) != "gene_name"]

## Ensure sample order matches
coldata_modified <- coldata_modified[colnames(count), ]

stopifnot(all(colnames(count) == rownames(coldata_modified)))

###############################
# DESeq2 object
###############################

dds <- DESeqDataSetFromMatrix(
  countData = count,
  colData = coldata_modified,
  design = ~ Group
)

###############################
# Filter low-count genes
###############################

keep <- rowSums(counts(dds) >= 10) >= (ncol(dds)/2)

dds <- dds[keep, ]

###############################
# Differential expression
###############################

dds <- DESeq(dds)

res <- results(dds)

deg <- as.data.frame(res)
deg$Gene <- rownames(deg)
deg$category <- "NS"

deg$category[
    deg$log2FoldChange > 1 &
    deg$padj < 0.05
] <- "Up"

deg$category[
    deg$log2FoldChange < -1 &
    deg$padj < 0.05
] <- "Down"
deg_sig <- subset(deg,
                  padj < 0.05 &
                   abs(log2FoldChange) > 1)

###############################
# Variance Stabilization
###############################

vsd <- varianceStabilizingTransformation(
  dds,
  blind = FALSE
)

expression.data <- t(assay(vsd))

###############################
# Trait
###############################

trait <- data.frame(
  Group = as.numeric(as.factor(coldata_modified$Group))
)

rownames(trait) <- rownames(coldata_modified)

trait <- trait[rownames(expression.data), , drop = FALSE]

###############################
# Quality control
###############################

gsg <- goodSamplesGenes(expression.data)

if (!gsg$allOK) {
  expression.data <- expression.data[gsg$goodSamples,
                                     gsg$goodGenes]
}

###############################
# Sample clustering
###############################

sampleTree <- hclust(
  dist(expression.data),
  method = "average"
)

###############################
# Soft threshold
###############################

spt <- pickSoftThreshold(
  expression.data,
  networkType = "signed",
  verbose = 5
)
#Plot the results
par(mfrow=c(1,2))
 #Plot 1: Scale-free topology fit
plot(spt$fitIndices[,1],-sign(spt$softIndices[,3])*spt$softIndices[,2],xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit (R²)",
     type = "n", 
     main = "Scale Independence")

text(spt$fitIndices[,1], 
     -sign(spt$fitIndices[,3]) * spt$fitIndices[,2],
     labels = powers, 
     cex = 0.9, 
     col = "red")

# Plot 2: Mean connectivity
plot(spt$fitIndices[,1], 
     spt$fitIndices[,5],
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", 
     type = "n",
     main = "Mean Connectivity")

text(spt$fitIndices[,1], 
     spt$fitIndices[,5], 
     labels = powers, 
     cex = 0.9, 
     col = "red")

softPower <- spt$powerEstimate

if (is.na(softPower))
  softPower <- 6

###############################
# Network construction
###############################

adjacency <- adjacency(
  expression.data,
  power = softPower,
  type = "signed"
)

#  Measures not just direct correlation between two genes, but also how many neighbors they share.
TOM <- TOMsimilarity(adjacency,TOMType = "signed") 

 dissTOM <- 1 - TOM         #To find out the dissimilarity from the TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree)
###############################
# Module detection
###############################

dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = 30
)
table(dynamicMods)
dynamicColors <- labels2colors(dynamicMods)

###############################
# Merge similar modules
###############################

merge <- mergeCloseModules(
  expression.data,
  dynamicColors,
  cutHeight = 0.25,
  verbose = 3
)
print("#######The columns present in the merge data is :########")
print(colnames(merge))
moduleColors <- merge$colors
MEs <- merge$newMEs
head(MEs)
MElist <- moduleEigengenes(expression.data, colors = moduleColors)
print("########List of Module Eigengenes###########")
print(MElist)
ME_eigen <- MElist$eigengenes
head(ME_eigen)
gene_color_table <- data.frame(
  Gene_Name = colnames(expression.data),
  Module_Color = moduleColors
)
write.csv(gene_color_table, "Genes_with_colors.csv", row.names  = FALSE)
ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete")
METree = hclust(as.dist(ME.dissimilarity), method = "average")
par(mar = c(0,4,2,0))
par(cex = 0.6);
plot(METree)
abline(h=.25, col = "red")


###### Dendogram analysis############################
plotDendroAndColors(
  dendro = geneTree,
  colors = moduleColors,
  groupLabels = "merged Modules",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)

mergedColors = merge$colors
mergedMEs = merge$newMEs
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Original Module", "Merged Module"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors for original and merged modules")


###############################
# Module-trait correlation
###############################

moduleTraitCor <- cor(
  MEs,
  trait,
  use = "p"
)

moduleTraitPvalue <- corPvalueStudent(
  moduleTraitCor,
  nSamples = nrow(expression.data)
)

cat("\nModule-Trait Correlations:\n")
print(moduleTraitCor)

cat("\nModule-Trait P-values:\n")
print(moduleTraitPvalue)

write.csv(
  moduleTraitCor,
  "Module_Trait_Correlation.csv"
)

write.csv(
  moduleTraitPvalue,
  "Module_Trait_Pvalues.csv"
)

###############################
# Automatically select module
###############################

moduleCor <- moduleTraitCor[,1]
moduleP <- moduleTraitPvalue[,1]

moduleNames <- gsub("^ME","",rownames(moduleTraitCor))

sigModules <- which(moduleP < 0.05)

if(length(sigModules)==0){

  stop("No module significantly correlated with the trait.")

}

bestModule <- sigModules[which.max(abs(moduleCor[sigModules]))]

module <- moduleNames[bestModule]

cat("\nSelected module:", module,"\n")

###############################
# Module Membership
###############################

moduleGenes <- moduleColors == module

MM <- cor(
  expression.data,
  MEs[,paste0("ME",module)],
  use = "p"
)

MMPvalue <- corPvalueStudent(
  MM,
  nrow(expression.data)
)


###############################
# Gene Significance
###############################

GS <- cor(
  expression.data,
  trait,
  use = "p"
)

GSPvalue <- corPvalueStudent(
  GS,
  nrow(expression.data)
)

###############################
# Hub genes
###############################

hubGenes <- data.frame(
  Gene = colnames(expression.data)[moduleGenes],
  MM = MM[moduleGenes],
  GS = GS[moduleGenes],
  MM_pvalue = MMPvalue[moduleGenes],
  GS_pvalue = GSPvalue[moduleGenes]
)

hubGenes <- hubGenes[
  order(abs(hubGenes$MM), decreasing = TRUE),
]

hubGenes_filtered <- subset(
  hubGenes,
  abs(MM) > 0.5 &
    abs(GS) > 0.2
)

write.csv(
  hubGenes,
  "HubGenes.csv",
  row.names = FALSE
)

write.csv(
  hubGenes_filtered,
  "HubGenes_Filtered.csv",
  row.names = FALSE
)

###############################
# Candidate hub DEGs
###############################

candidateGenes <- merge(
  hubGenes_filtered,
  deg_sig,
  by = "Gene"
)

cat("\nNumber of candidate hub DEGs:",
    nrow(candidateGenes), "\n")

print(candidateGenes)

write.csv(
  candidateGenes,
  "Candidate_hub_DEGs.csv",
  row.names = FALSE
)


hub_gene_indices <- match(hubGenes_filtered$Gene, colnames(expression.data))

