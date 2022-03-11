####  Run on Quartz  ####
library(scran)
library(Seurat, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(zinbwave, quietly = TRUE)
library(SingleCellExperiment, quietly = TRUE)
library(matrixStats, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(biomaRt, quietly = TRUE)
library(edgeR, quietly = TRUE)
library(iDEA, quietly = TRUE)
library(zingeR, quietly = TRUE)
library(DESeq2, quietly = TRUE)
library(BiocParallel, quietly = TRUE)

message("library loaded")

register(MulticoreParam(8))

dds <- readRDS(paste0("/N/slate/yosueda/20220208/Output/dds_SC.rds"))

message("data loaded")

dds <- DESeq(dds,
             test = "Wald", 
             sfType = "poscounts", 
             useT = T, 
             betaPrior = T, 
             modelMatrixType="expanded",
             parallel = T)

message("DESeq done")

res <- results(dds, cooksCutoff = F, parallel = T)
res$log2FoldChange %>% summary()
res$lfcSE %>% summary()

#Save DESeq results
full.results <- results(dds, parallel = T)
full.results.df <- data.frame(gene = row.names(full.results),
                              log2FoldChange = full.results$log2FoldChange,
                              pvalue = full.results$pvalue,
                              padj = full.results$padj,
                              baseMean = full.results$baseMean,
                              logFC_SE = full.results$lfcSE)
write.csv(full.results.df, file = "/N/slate/yosueda/20220208/Output/DESeq2_SC.csv")
saveRDS(dds, "/N/slate/yosueda/20220208/Output/dds_SC_new.rds")
message("data saved")
