# 安装必要包（如未安装）
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ReactomePA"))

library(clusterProfiler)
library(org.Hs.eg.db)      # 人类注释数据库
library(ReactomePA)
library(enrichplot)

# 基因列表
genes <- c("CALN1", "SOX4", "URI1", "NRIP1", "RPS29", "RPS10", "RPS14", 
           "COA5", "MBD2", "NPDC1", "RPL41", "RPS26")

# 转换为 Entrez ID
gene_df <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
entrez_genes <- gene_df$ENTREZID
ego <- enrichGO(gene         = entrez_genes,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",     # 或 "MF", "CC"
                pAdjustMethod= "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)

dotplot(ego, showCategory=10)

ekegg <- enrichKEGG(gene         = entrez_genes,
                    organism     = 'hsa',  # human
                    pvalueCutoff = 0.05)

dotplot(ekegg, showCategory=10)

ereact <- enrichPathway(gene=entrez_genes, pvalueCutoff=0.05, readable=TRUE)
dotplot(ereact, showCategory=10)




# 查看每个基因所参与的通路
# 安装并加载 KEGGREST
if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("KEGGREST")
}
library(KEGGREST)

# 你的基因列表
genes <- c("CALN1", "SOX4", "URI1", "NRIP1", "RPS29", 
           "RPS10", "RPS14", "COA5", "MBD2", "NPDC1", 
           "RPL41", "RPS26")

# 将 gene symbol 转为 Entrez ID
library(org.Hs.eg.db)
library(clusterProfiler)
gene_df <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# 遍历每个 Entrez ID 查询 KEGG 路径
gene_pathways <- lapply(gene_df$ENTREZID, function(entrez_id) {
  tryCatch({
    pathways <- keggLink("pathway", paste0("hsa:", entrez_id))
    data.frame(Gene = gene_df$SYMBOL[gene_df$ENTREZID == entrez_id],
               PathwayID = sub("path:", "", pathways))
  }, error = function(e) NULL)
})

# 合并所有结果
result_df <- do.call(rbind, gene_pathways)

# 添加 KEGG 通路名称
pathway_names <- keggList("pathway", "hsa")
result_df$PathwayName <- pathway_names[result_df$PathwayID]

# 查看结果
head(result_df)
write.csv(result_df, "gene_pathway_mapping.csv", row.names = FALSE)
