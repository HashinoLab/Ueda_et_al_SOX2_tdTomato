library(Seurat)
library(DESeq2)

message("Library loaded.")

day2030405060 <- readRDS("/N/slate/yosueda/20220208/Output/day2030405060_dim9_res0.7.rds")

message("Data loaded.")

day2030405060 <- subset(day2030405060, ident = c("1","13"))
count <- GetAssayData(day2030405060, assay = "RNA", slot = "counts")
count <- as.matrix(count)
group <- data.frame(con = factor(c(day2030405060$seurat_clusters)))

message("Count and group created.")

saveRDS(count, "/N/slate/yosueda/20220208/Output/count.rds")
saveRDS(group, "/N/slate/yosueda/20220208/Output/group.rds")

message("Count and group saved.")

dds <- DESeqDataSetFromMatrix(count, colData = group, design = ~ con)
dds <- DESeq(dds)
res <- results(dds)
head(res)

message("DEseq done.")

DESeq2::plotMA(res, alpha = 0.01)

message("DEseq done.")

write.table(res, file = "/N/slate/yosueda/20220208/Output/result.txt", row.names = T, col.names = T, sep = "\t")

message("Table saved.")

write.csv(res, file = "/N/slate/yosueda/20220208/Output/result.csv")

message("CSV saved.")
