library(dplyr)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(cowplot)

count <- read.csv("final_counts_annotated.csv",header = TRUE,check.names = FALSE)
coldata <- read.csv("metadata_new.csv",header=TRUE)

coldata_modified <- coldata %>%
  dplyr::select(1,8)
 
rownames(coldata_modified) <- coldata_modified$Run

rownames(count) <- make.unique(as.character(count$gene_name))
count <- count[, colnames(count) !="gene_name"]
all(colnames(count) %in% rownames(coldata_modified))

dds <- DESeqDataSetFromMatrix(countData = count,colData = coldata_modified,design = ~ Group)
dds <- DESeq(dds)
res <- results(dds)

res <- as.data.frame(res)
res$category <- "NS"
res$category[res$log2FoldChange > 1 & res$padj < 0.05] <- "Up"
res$category[res$log2FoldChange < -1 & res$padj < 0.05] <- "Down"

ggplot(res,aes(x = log2FoldChange, y = -log10(padj),color = category)) +
  geom_point() +
  scale_color_manual(values = c("Up" = "red","Down" = "blue", "NS" = "grey")) +
  
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed") +
  
  theme_classic() +
  labs(title = "DS vs CON Volcano Plot",x = "log2 Fold Change",y = "-log10 adjusted p-value")

res$Ensembl_ID <- sub("\\..*", "", rownames(res))

symbols <- mapIds(
  org.Hs.eg.db,
  keys = res$Ensembl_ID,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

res$GeneSymbol <- unname(symbols)


top_up <- res %>%
  filter(category == "Up") %>%
  arrange(padj) 

top_down <- res %>%
  filter(category == "Down") %>%
  arrange(padj)

top_genes <- rbind(top_up, top_down)


ggplot(res,aes(x = log2FoldChange, y = -log10(padj),color = category)) +
  geom_point(aes(color=category),alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Up" = "red","Down" = "blue", "NS" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed") +
  geom_text_repel(data=top_genes,aes(label = GeneSymbol),
                  size = 3.5,box.padding = 0.1,       
                  point.padding = 0.1,     
                  force = 1,               
                  max.overlaps = 20,       
                  min.segment.length = 0)+
  guides(color = guide_legend(
    override.aes = list(
      shape = 16,  
      size =3,    
      alpha = 1)))+
  
  theme_minimal() +
  labs(title = "DS vs CON Volcano Plot",x = "log2 Fold Change",y = "-log10 adjusted p-value")

write.csv(res,"res.csv",row.names = FALSE)
vsd <- vst(dds,blind=FALSE)
assay_vsd <- t(assay(vsd))
expr_ml <- as.data.frame(assay_vsd)
expr_ml$condition <- colData(dds)$condition
write.csv(assay_vsd,"Expression data.csv")
pca_data <- plotPCA(vsd,intgroup = "Group",returnData=TRUE)
plotPCA(vsd,intgroup = "Group")
ggplot(pca_data, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = rownames(pca_data))) +
  theme_classic()+
  labs(title="PCA plot with sample name")

distance <- dist(t(assay(vsd,method = "correlation")))
cluster <- hclust(distance)

plot(cluster)


reorder <- res[order(res$padj),]
top_genes <- rownames(reorder)[1:50]
mat <- assay(vsd)[top_genes,]
mat <- t(scale(t(mat)))
 


pheatmap(mat,annotation_col = annotation_df,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = cols2,
         angle_col = 90,fontsize_row = 8,
         fontsize_col = 10,width = 20,height = 70,border_color = "NA")

cols <- brewer.pal(9, "YlOrRd")
pheatmap(mat,labels_row = symbols,annotation_col = annotation_df,
         color = cols2,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_cols = TRUE,border_color = "NA",
         treeheight_row = 0.1,
         fontsize_row = 10)

cols2 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

sig_genes_up <- subset(res, padj < 0.05 & log2FoldChange > 1)
gene_df <- bitr(
  sig_genes_up$Ensembl_ID,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
gene_list <- unique(gene_df$ENTREZID)
ego <- enrichGO(gene = gene_list,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
dotplot(ego)

ego_MF <- enrichGO(gene = gene_list,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = "MF",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
dotplot(ego_MF)


ego_cc <- enrichGO(gene = gene_list,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = "CC",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
dotplot(ego_cc)



sig_genes_down <- subset(res, padj < 0.05 & log2FoldChange < -1)

gene_df_down <- bitr(sig_genes_down$Ensembl_ID,fromType="SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
gene_list_down <- unique(gene_df_down$ENTREZID)
ego_down <- enrichGO(gene = gene_list_down,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
dotplot(ego_down)

ego_down_mf <- enrichGO(gene = gene_list_down,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = "MF",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
dotplot(ego_down_mf)

ego_down_cc <- enrichGO(gene = gene_list_down,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = "CC",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
dotplot(ego_down_cc)



ego_down <- enrichGO(gene = gene_list_down,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
dotplot(ego_down)





plotMA(res)

kegg_up <- enrichKEGG(gene=gene_list,org='hsa',pvalueCutoff=0.05,pAdjustMethod="BH")
dotplot(kegg_up)

kegg_down <- enrichKEGG(gene=gene_list_down,organism='hsa',keyType = "ncbi-geneid",pvalueCutoff = 0.2,pAdjustMethod="BH")
dotplot(kegg_down)
kegg_down_readable <- setReadable(kegg_down,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
write.csv(as.data.frame(kegg_down_readable),"KEGG down regulating.csv",row.names = FALSE)

KEGG__up_readable <- setReadable(kegg_up,OrgDb=org.Hs.eg.db,keyType="ENTREZID")
head(KEGG__up_readable)
head(as.data.frame(KEGG__up_readable))
write.csv(KEGG__up_readable,"KEGG__up associated genes.csv",row.names=FALSE)




p1_up <- dotplot(ego)
p1_up_MF <- dotplot(ego_MF)
p2_up <- dotplot(kegg_up)
plot_grid(p1_up_MF,p2_up,labels = c("GO_UP(MF)","KEGG_up"))





