library(randomForest)
install.packages("glmnet")
library(glmnet)
install.packages("VennDiagram")
library(VennDiagram)

vsd <- vst(dds,blind=FALSE)
assay_vsd <- t(assay(vsd))
expr_ml <- as.data.frame(assay_vsd)
expr_ml$Group <- colData(dds)$Group

res <- results(dds)
deg <- rownames(subset(res,padj < 0.05 & log2FoldChange > 1))
expr_deg <- expr_ml[, colnames(expr_ml) %in% deg]

expr_deg$Group <- colData(dds)$Group
colnames(expr_deg) <- make.names(colnames(expr_deg))
rf <- randomForest(Group ~.,data = expr_deg,importance=TRUE,ntree=2000)
varImpPlot(rf)
imp <- importance(rf)

x <- as.matrix(expr_deg[, colnames(expr_deg) != "Group"])
y <- expr_deg$Group
y=ifelse(expr_deg$Group=="DS",0,1)
set.seed(123)
 cv_fit <- cv.glmnet(x=x,y=y,family="binomial",alpha=1,nfolds=4)
 plot(cv_fit)
best_model <- cv_fit$lambda.min
cv_fit$lambda.1se

coef_min <- coef(cv_fit,s="lambda.min")

imp_df <- data.frame(Gene=rownames(imp),importance=imp[,"MeanDecreaseGini"])
imp_df <- imp_df[order(imp_df$importance,decreasing =TRUE),]
rf_genes <- imp_df$Gene

selected_genes <- rownames(coef_min)[coef_min[,1] != 0]

selected_genes <- selected_genes[selected_genes != "(Intercept)"]

overlap_genes <- intersect(selected_genes,rf_genes)
res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)


lasso_df <- data.frame(
  Gene = rownames(coef_min),
  Coefficient = as.numeric(coef_min)
)
merged <- merge(res_df,imp_df,by="Gene")
merged <- merge(merged,lasso_df,by="Gene")

final <- merged[merged$Gene %in% overlap_genes,]

prob <- predict(
  cv_fit,
  s = "lambda.min",
  newx = x,
  type = "response"
)

write.csv(final,"Biomarker of upregulated genes in DS.CSV",row.names = FALSE)



grid.newpage()
venn_plot <- draw.pairwise.venn(
  area1 = length(selected_genes),
  area2 = length(rf_genes),
  cross.area = length(overlap_genes),
  category = c("Lasso Regression", "Random Forest"),
  fill = c("#3B82F6", "#10B981"),
  alpha = c(0.5, 0.5),
  lty = "blank",
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  fontfamily = "sans",
  cat.fontfamily = "sans"
)


plot(cv_fit)
abline(v = log(cv_fit$lambda.min), col="red")
abline(v = log(cv_fit$lambda.1se), col="blue")
predict(cv_fit, type="coefficients", s="lambda.min")
