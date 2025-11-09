
library(DESeq2)
library(limma)

expr <- read.csv("./BeatAML2/expression_matrix.csv", row.names = 1)
meta <- read.csv("./BeatAML2/metadata.csv")  
summary(as.numeric(as.matrix(expr)))
min_val <- min(as.numeric(as.matrix(expr)), na.rm = TRUE); min_val

rownames(expr)

length(intersect(colnames(expr), meta$sample_id))

common <- intersect(colnames(expr), meta$sample_id)
expr <- expr[, common]
meta <- meta[match(common, meta$sample_id), ] 
identical(colnames(expr), meta$sample_id) 

meta$timepoint <- factor(meta$timepoint, levels=c("pre","post"))
meta$response <- factor(meta$response, levels=c("sensitive","resistant"))
design <- model.matrix(~ timepoint + response, data=meta)
 
dim(expr)     # genes × samples
dim(design)   # samples × (predictors)
ncol(expr) == nrow(design)  



fit <- lmFit(expr, design)
fit <- eBayes(fit, trend = TRUE)

# post vs pre
res_prepost <- topTable(fit, coef = "timepointpost", number = Inf, sort.by = "P")
deg_prepost <- subset(res_prepost, adj.P.Val < 0.05 & abs(logFC) > 1)
write.csv(deg_prepost, "./BeatAML2/DEG_pre_vs_post_limma.csv", row.names = TRUE)

# responder vs nonresponder
if ("responseresistant" %in% colnames(coef(fit))) {
  res_resp <- topTable(fit, coef = "responseresistant", number = Inf, sort.by = "P")
  deg_resp <- subset(res_resp, adj.P.Val < 0.05 & abs(logFC) > 1)
  gene_list <- rownames(deg_resp)
  write.csv(data.frame(gene = gene_list),
            "./BeatAML2/DEG_list.csv", row.names = FALSE, quote = FALSE)
  write.csv(deg_resp, "./BeatAML2/DEG_responder_vs_nonresponder_limma.csv", row.names = TRUE)
}




