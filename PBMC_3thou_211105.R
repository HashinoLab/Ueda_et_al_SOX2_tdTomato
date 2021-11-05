# INPUT project name and folder path
Proj <- "SOX2 project"
File_path <- "C:/Users/yosueda/Documents/RStudio/"

# Load prerequisite R libraries
library(ggrepel)
library(scales)
library(Seurat)
library(tidyverse)

# Load scRNA-seq data
RNA_dat <- Read10X(paste0(File_path, Proj, "/D60-4/Data/"))

# Convert scRNA-seq data into Seurat object
RNA_Seu <- CreateSeuratObject(RNA_dat, project = Proj, min.cells = 0, min.features = 0)

# Calculate percentage of mtRNA
RNA_QC <- PercentageFeatureSet(RNA_Seu, pattern = "^MT-", col.name = "Pct_mito")

# Visualize QC metrics using violin plots
tiff(paste0(File_path, Proj, "/Output/QC_plot.tiff"), width = 8.5, height = 8.5, units = "in", res = 600)
  VlnPlot(RNA_QC, features = c("nFeature_RNA", "Pct_mito"), pt.size = 0.5, ncol = 2)
dev.off()

# STOP here and INPUT control limits
  ## D20: 1200/7000/12    ## D30: 1500/7000/12    ## D40: 2000/6000/10
  ## D50: 1500/7000/12    ## D60: 1200/6500/12
nFeature_min <- 1200
nFeature_max <- 6500
Pct_mito_max <- 12

# Filter out low-quality cells using control limits
RNA_QC <- subset(RNA_QC, subset = nFeature_RNA > nFeature_min & nFeature_RNA < nFeature_max & Pct_mito < Pct_mito_max)

# Perform negative binomial regression, limiting subsequent analysis to 3,000 most variably expressed genes
RNA_SCT <- SCTransform(RNA_QC, do.correct.umi = TRUE, variable.features.n = 3000,
                       vars.to.regress = "Pct_mito", do.scale = FALSE, do.center = TRUE, 
                       return.only.var.genes = TRUE, verbose = TRUE)

# Construct and export variable feature plot
Top_genes <- head(VariableFeatures(RNA_SCT), 20)
  tiff(paste0(File_path, Proj, "/Output/Feat_plot.tiff"), width = 11, height = 8.5, units = "in", res = 600)
  LabelPoints(VariableFeaturePlot(RNA_SCT, pt.size = 0.5, selection.method = "sctransform"),
              points = Top_genes, repel = TRUE, xnudge = 0.0, ynudge = 0.0)
dev.off()

# Perform PCA
RNA_PCA <- RunPCA(RNA_SCT, npcs = 50, weight.by.var = TRUE, verbose = FALSE, approx = TRUE)

# Determine 'dimensionality' of scRNA-seq data
  # Construct and export elbow plot
  tiff(paste0(File_path, Proj, "/Output/Elb_plot.tiff"), width = 11, height = 8.5, units = "in", res = 600)
    ElbowPlot(RNA_PCA, ndims = 50, reduction = "pca")
  dev.off()
  
# Save Seurat object
saveRDS(RNA_PCA, file = paste0(File_path, Proj, "/Output/RNA_PCA_60.rds"))


########## Above has done ################
########## Below should be done ##########

# Input day
day <- 20

RNA_PCA <- readRDS(file = paste0(File_path, Proj, "/Output/RNA_PCA_", toString(day), ".rds"))

# STOP here and INPUT dimensionality and resolution
nPC <- 9
res <- 0.4

# Construct SNN graph
RNA_PCA <- FindNeighbors(RNA_PCA, k.param = 150, compute.SNN = TRUE, prune.SNN = 1/15, 
                         nn.method = "rann", nn.eps = 0.0, verbose = TRUE, reduction = "pca", dims = 1:nPC) %>%
  # Cluster cells by maximizing modularity
  FindClusters(modularity.fxn = 1, resolution = res, algorithm = 1, group.singletons = TRUE, verbose = TRUE)

# Visualize PCA results using UMAP
RNA_PCA <- RunUMAP(RNA_PCA, umap.method = "uwot", n.neighbors = 200, n.components = 2, min.dist = 0.5,
                   spread = 1, set.op.mix.ratio = 1.0, local.connectivity = 1L, 
                   verbose = TRUE, dims = 1:nPC, reduction = "pca")

# Construct and export UMAP plot
tiff(paste0(File_path, Proj, "/Output/UMAP", toString(nPC), "_", toString(res), "_d", toString(day), ".tiff"),
     width = 8.5, height = 8.5, units = "in", res = 600)
  DimPlot(RNA_PCA, pt.size = 1.0, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE)
dev.off()

# Perform DEA using Wilcoxon test
Wilcox_tab <- FindAllMarkers(RNA_PCA, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, verbose = TRUE, only.pos = TRUE, return.thresh = 0.01) %>%
  select(c("cluster", "gene", "pct.1", "pct.2", "p_val_adj", "avg_log2FC")) %>%
    group_by(cluster) %>%
      top_n(30) %>%
        arrange(cluster, desc(avg_log2FC))

# Export Wilcoxon test results
write.csv(Wilcox_tab, file = paste0(File_path, Proj, "/Output/Wilcox_test_", toString(nPC), "_", toString(res), "_d", toString(day), ".csv"))

# Save Seurat object
saveRDS(RNA_PCA, file = paste0(File_path, Proj, "/Output/RNA_PCA_", toString(nPC), "_", toString(res), "_d", toString(day),".rds"))

#########################################################################
## PLZ upload the file "UMAP", ""Wilcox_test", and processed "RNA_PCA" ##
#########################################################################
