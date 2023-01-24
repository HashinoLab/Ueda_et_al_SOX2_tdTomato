library(Seurat)
library(tidyverse)


#### Input path ####
Proj <- "SOX2 project"
File_path <- "C:/Users/yosueda/Documents/RStudio/"


#### Input day ####
day <- "d20"


#### Load scRNA-seq data ####
RNA_dat <- Read10X(paste0(File_path, Proj, "/", day, "/Data/"))


#### Convert scRNA-seq data into Seurat object ####
RNA_Seu <- CreateSeuratObject(RNA_dat, project = Proj, min.cells = 0, min.features = 0)


#### Calculate %mtRNA ####
RNA_QC <- PercentageFeatureSet(RNA_Seu, pattern = "^MT-", col.name = "Pct_mito")


#### Visualize QC metrics ####
tiff(paste0(File_path, Proj, "/Output/QC_plot.tiff"), width = 8.5, height = 8.5, units = "in", res = 600)
  VlnPlot(RNA_QC, features = c("nFeature_RNA", "Pct_mito"), pt.size = 0.5, ncol = 2)
dev.off()


#### Input control limits ####
# d20, 1200/7000/12; d30, 1500/7000/12; d40, 2000/6000/10; d50, 1500/7000/12; d60, 1200/6500/12
nFeature_min <- 1200
nFeature_max <- 7000
Pct_mito_max <- 12


#### Filter out low-quality cells using control limits ####
RNA_QC <- subset(RNA_QC, subset = nFeature_RNA > nFeature_min & nFeature_RNA < nFeature_max & Pct_mito < Pct_mito_max)


#### Perform SCTransform ####
RNA_SCT <- SCTransform(RNA_QC, do.correct.umi = T, variable.features.n = 3000,
                       vars.to.regress = "Pct_mito", do.scale = F, do.center = T, 
                       return.only.var.genes = T, verbose = T)


#### Construct and export variable feature plot ####
Top_genes <- head(VariableFeatures(RNA_SCT), 20)
  tiff(paste0(File_path, Proj, "/Output/Feat_plot.tiff"), width = 11, height = 8.5, units = "in", res = 600)
  LabelPoints(VariableFeaturePlot(RNA_SCT, pt.size = 0.5, selection.method = "sctransform"),
              points = Top_genes, repel = TRUE, xnudge = 0.0, ynudge = 0.0)
dev.off()


#### Perform PCA ####
RNA_PCA <- RunPCA(RNA_SCT, npcs = 50, weight.by.var = TRUE, verbose = FALSE, approx = TRUE)


#### Construct elbow plot to determine dimensionality ####
tiff(paste0(File_path, Proj, "/Output/Elb_plot.tiff"), width = 11, height = 8.5, units = "in", res = 600)
  ElbowPlot(RNA_PCA, ndims = 50, reduction = "pca")
dev.off()


#### Save Seurat object ####
saveRDS(RNA_PCA, file = paste0(File_path, Proj, "/Output/RNA_PCA_20.rds"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#    Perform above for each dataset of d20, d30, d40, d50, and d60    %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#### Read each day data ####
d20 <- readRDS(file = paste0(File_path, Proj, "/Output/RNA_PCA_", toString(20), ".rds"))
d30 <- readRDS(file = paste0(File_path, Proj, "/Output/RNA_PCA_", toString(30), ".rds"))
d40 <- readRDS(file = paste0(File_path, Proj, "/Output/RNA_PCA_", toString(40), ".rds"))
d50 <- readRDS(file = paste0(File_path, Proj, "/Output/RNA_PCA_", toString(50), ".rds"))
d60 <- readRDS(file = paste0(File_path, Proj, "/Output/RNA_PCA_", toString(60), ".rds"))


#### Input dimensionality and resolution ####
nPC <- 9
res <- 0.4


#### Construct SNN graph and cluster cells by maximizing modularity ####
d20 <- FindNeighbors(d20, k.param = 110, nn.method = "rann", dims = 1:nPC) %>% FindClusters(resolution = res)
d30 <- FindNeighbors(d30, k.param = 120, nn.method = "rann", dims = 1:nPC) %>% FindClusters(resolution = res)
d40 <- FindNeighbors(d40, k.param = 90, nn.method = "rann", dims = 1:nPC) %>% FindClusters(resolution = res)
d50 <- FindNeighbors(d50, k.param = 100, nn.method = "rann", dims = 1:nPC) %>% FindClusters(resolution = res)
d60 <- FindNeighbors(d60, k.param = 130, nn.method = "rann", dims = 1:nPC) %>% FindClusters(resolution = res)


#### Construct UMAP ####
d20 <- RunUMAP(d20, umap.method = "uwot", n.neighbors = 10, min.dist = 0.3, dims = 1:nPC)
d30 <- RunUMAP(d30, umap.method = "uwot", n.neighbors = 10, min.dist = 0.3, dims = 1:nPC)
d40 <- RunUMAP(d40, umap.method = "uwot", n.neighbors = 10, min.dist = 0.3, dims = 1:nPC)
d50 <- RunUMAP(d50, umap.method = "uwot", n.neighbors = 10, min.dist = 0.3, dims = 1:nPC)
d60 <- RunUMAP(d60, umap.method = "uwot", n.neighbors = 10, min.dist = 0.3, dims = 1:nPC)


#### UMAPs ####
color <- NULL
color[[1]] <- c("#D1BBD7","#F6C141","#F7F056","#4EB265","#882E72","#4EB265","#CAE0AB","#90C987")
color[[2]] <- c("#F6C141","#D1BBD7","#4EB265","#F7F056","#90C987","#CAE0AB","#882E72","#1965B0")
color[[3]] <- c("#D1BBD7","#E8601C","#F1932D","#1965B0","#CAE0AB","#AE76A3","#90C987","#4EB265")
color[[4]] <- c("#E8601C","#F1932D","#AE76A3","#1965B0","#CAE0AB","#4EB265","#90C987","#7BAFDE","#5289C7")
color[[5]] <- c("#F1932D","#1965B0","#E8601C","#7BAFDE","#5289C7","#DC050C","#90C987","#4EB265","#CAE0AB")

for (i in c(1:5)){
  day <- i*10 + 10
tiff(paste0(File_path, Proj, "/d", day, "/Output/UMAP_d", day, "_", nPC, "_", res, ".tiff"),
     width = 8, height = 8, units = "in", res = 600)
print(
  DimPlot(RNA_PCA[[i]], reduction = "umap", pt.size = 1.5, label = F, label.size = 5, cols = color[[i]]) +
  theme(plot.title = element_text(size = 0), axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.background = element_blank(),legend.position = "none")
)
dev.off()
}


#### Perform DEA using Wilcoxon test ####
Wilcox_tab <- NULL; for (i in c(1:5)){
  day <- i*10 + 10
  Wilcox_tab[[i]] <- FindAllMarkers(RNA_PCA[[i]], only.pos = T) %>%
    select(c("cluster", "gene", "pct.1", "pct.2", "p_val_adj", "avg_log2FC")) %>%
    group_by(cluster) %>% top_n(50) %>% arrange(cluster, desc(avg_log2FC))
  write.csv(Wilcox_tab[[i]], paste0(File_path, Proj, "/d", day, "/Output/Wilcox_test_", nPC, "_", res, "_d", day, ".csv"))
  saveRDS(RNA_PCA[[i]], paste0(File_path, Proj, "/d", day, "/Output/RNA_PCA_", nPC, "_", res, "_d", day, ".rds"))
}


#### Feature Plots ####
markers <- c("PAX2", "FBXO2", "NEUROD1", "TOP2A", "S100B")
for (i in c(1:5)){
  day <- i*10 + 10
  for (j in markers){
    tiff(paste0(File_path, Proj, "/d", day, "/Output/d", day, "_", j, ".tiff"), width = 8, height = 8, units = "in", res = 600)
    print(
      FeaturePlot(RNA_PCA[[i]], features = j, pt.size = 2, label = F, max.cutoff = 2) + 
        theme(plot.title = element_text(size = 0),
              axis.text.x = element_blank(), axis.text.y = element_blank(),
              axis.line.x = element_blank(), axis.line.y = element_blank(),
              axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
              axis.title.x = element_blank(), axis.title.y = element_blank(),
              panel.background = element_blank(),
              legend.position = "none")
    )
    dev.off()
  }
}

markers <- c("POU4F3", "PCP4", "ATOH1", "LMO7", "OTOF")
for (i in 5){
  day <- i*10 + 10
  for (j in markers){
    tiff(paste0(File_path, Proj, "/d", day, "/Output/d", day, "_", j, ".tiff"), width = 8, height = 8, units = "in", res = 600)
    print(
      FeaturePlot(RNA_PCA[[i]], features = j, pt.size = 2, label = F, max.cutoff = 2) + 
        theme(plot.title = element_text(size = 0),
              axis.text.x = element_blank(), axis.text.y = element_blank(),
              axis.line.x = element_blank(), axis.line.y = element_blank(),
              axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
              axis.title.x = element_blank(), axis.title.y = element_blank(),
              panel.background = element_blank(),
              legend.position = "none")
    )
    dev.off()
  }
}


#### Dot Plots ####
dotlist <- c("TOP2A","UBE2C","MKI67",
             "PCNA","NASP","TYMS",
             "PAX2","JAG1","TBX1",
             "EMX2","OTX2","OC90",
             "FBXO2","BRICD5","SPARCL1",
             "PCP4","MYO6","POU4F3",
             "DACH1","WNT2B",
	     "PTGDS","SPP1","CITED1",
             "NEUROD1","ELAVL4","POU4F1",
             "S100B","PLP1","PMP22",
             "DSTN","TFAP2B","TFAP2A","DCN","KRT19",
             "KRT18","KRT8","TAGLN",
             "COL1A2","COL1A1","POSTN","FN1","S100A6",
             "TPM1","TRPM3","ACTA2",
             "CDH1","CDH2","EPCAM")

# d20
tiff(
  paste0(File_path, Proj, "/d20/Output/d20_dot.tiff"),
  width = 3, height = 9, units = "in", res = 600)
DotPlot(d20,
        features = rev(dotlist),
  col.min = -2, col.max = 2, dot.min = 0, dot.scale = 8) +
  coord_flip() +
  scale_y_discrete(position = "right",
                   limits = c("6","2","1","0","3","5","7","4")) +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
dev.off()

# d30
tiff(
  paste0(File_path, Proj, "/d30/Output/d30_dot.tiff"),
  width = 3, height = 9, units = "in", res = 600)
DotPlot(d30,
        features = rev(dotlist),
        col.min = -2, col.max = 2, dot.min = 0, dot.scale = 8) +
  coord_flip() +
  scale_y_discrete(position = "right",
                   limits = c("5","3","0","1","2","4","6","7")) +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
dev.off()

# d40
tiff(
  paste0(File_path, Proj, "/d40/Output/d40_dot.tiff"),
  width = 3, height = 9, units = "in", res = 600)
DotPlot(d40,
        features = rev(dotlist),
        col.min = -2, col.max = 2, dot.min = 0, dot.scale = 8) +
  coord_flip() +
  scale_y_discrete(position = "right",
                   limits = c("4","2","1","0","5","7","6","3")) +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
dev.off()

# d50
tiff(
  paste0(File_path, Proj, "/d50/Output/d50_dot.tiff"),
  width = 3.2, height = 9, units = "in", res = 600)
DotPlot(d50,
        features = rev(dotlist),
        col.min = -2, col.max = 2, dot.min = 0, dot.scale = 8) +
  coord_flip() +
  scale_y_discrete(position = "right",
                   limits = c("4","1","0","2","5","6","7","3","8")) +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
dev.off()

# d60
tiff(
  paste0(File_path, Proj, "/d60/Output/d60_dot.tiff"),
  width = 3.2, height = 9, units = "in", res = 600)
DotPlot(d60,
        features = rev(dotlist),
        col.min = -2, col.max = 2, dot.min = 0, dot.scale = 8) +
  coord_flip() +
  scale_y_discrete(position = "right",
                   limits = c("8","0","2","5","7","6","3","1","4")) +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
dev.off()


#### Feature plots for d50 hair cells ####
tiff(paste0(File_path, Proj, "/d50/Output/d50_ATOH1_red.tiff"),
     width = 8, height = 8, units = "in", res = 600)
FeaturePlot(d50, features = "ATOH1", max.cutoff = 2, pt.size = 2, cols = c("#99999980", "red")) +
  theme(plot.title = element_text(size = 0),
        axis.line = element_line(size = 1),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "#333333"),
  )
dev.off()
tiff(paste0(File_path, Proj, "/d50/Output/d50_POU4F3_green.tiff"),
     width = 8, height = 8, units = "in", res = 600)
FeaturePlot(d50, features = "POU4F3", max.cutoff = 2, pt.size = 2, cols = c("#99999980", "green")) +
  theme(plot.title = element_text(size = 0),
        axis.line = element_line(size = 1),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "#333333"),
  )
dev.off()


# End of Code
