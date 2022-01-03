#%%%  EXPERIMENTAL PIPELINE, SEURAT  %%%%
day2030405060 <- readRDS("C:/Users/yosueda/Documents/RStudio/SOX2 project/Combine/Output/day2030405060_9_0.5.rds");gc();gc()

# Based on AutoClustR, mathematically reasonable value is dim=10, res=0.8, k.param=160
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
dim = 10; res = 0.7

ElbowPlot(day2030405060)

day2030405060 <- readRDS("C:/Users/yosueda/Documents/RStudio/SOX2 project/Combine/20211227/Output/day2030405060_10_0.7_20211227.rds")

#### Construct SNN graph ####
day2030405060 <- FindNeighbors(day2030405060, reduction = "pca", dims = 1:dim, k.param = 200,
                               compute.SNN = T, prune.SNN = 1/15, nn.method = "rann", n.trees = 50, nn.eps = 0.0,
                               verbose = T) %>%
  # Cluster cells by maximizing modularity
  FindClusters(modularity.fxn = 1, resolution = res, algorithm = 1, group.singletons = T, verbose = T)

n.n <- 50; dist = 0.25
day2030405060 <- RunUMAP(day2030405060, dims = 1:dim, n.neighbors = n.n, min.dist = dist, n.epochs = 500,
                         metric = "cosine", umap.method = "uwot", n.components = 2, return.model = F,
                         verbose = T)

DimPlot(day2030405060, reduction = "umap", pt.size = 0.3, label = T, label.size = 5, repel = T)
DimPlot(day2030405060, group.by = "Day", pt.size = 0.3, label = T, label.size = 5)

####  Changing cluster names  ####
cluster.annot <- c("0 otic prog", # PAX2 JAG1
                   "1 otic epi", # HES1 OC90 WFDC2 FBXO2
                   "2 G1/S", # CENPF TOP2A MKI67 UBE2C HIST1H4C HIST1H1A HIST1H1D HIST1H1B
                   "3 n.prog", # STMN2 NEUROD1 NEFM SST HES6 GAP43 DCX STMN4 ELAVL4 ELAVL3
                   "4 fibrobl", #MGP S100B FN1 COL1A1, COL1A2, COL3A1, 
                   "5 dark?", # PTGDS+ but NR3B2(ESRRB)- KCNQ1- KCNE1-, SLC12A2-(nkcc1)
                   "6 SC", # OTOL1 AGR2 BRICD5 WFDC2 FBXO2 SPARCL1 USH1C
                   "7 Schwann", # S100B PLP1 PMP22 MPZ
                   "8 KRT+", # S100A2 S100A6 S100A10 S100A11 KRT17 KRT18 KRT8 KRT7 KRT19 AGR2
                   "9 end-lymph.du", # WNT2B DACH1  
                   "10 SC-HC", # CXCL14 PCP4 SPARCL1 CALB1 TAGLN
                   "11 n.crest", # TFAP2B SLC2A3(neuronal glucose transporter 3 (GLUT3))
                                 # ARHGAP29(craniofacial expression)
                                 # S100A11
                   "12 stressed", # HERPUD1 DDIT3 NR4A1 GADD45B GADD45A
                   "13 G2/M", # CENPF TOP2A UBE2C MKI67 CDC20 UBE2S CDK1
                   "14 HC")
names(cluster.annot) <- levels(day2030405060)
day2030405060 <- RenameIdents(day2030405060, cluster.annot)


####  Dotplot  ####
tiff(
  paste0("C:/Users/yosueda/Documents/RStudio/SOX2 project/Combine/20211229/Output/dot.tiff"),
  width = 6, height = 13, units = "in", res = 600)
DotPlot(day2030405060,
        features = rev(c(
          "CENPF","TOP2A","UBE2C","MKI67","CDK1","CDC20","UBE2S", #13
          "HIST1H4C","HIST1H1A","HIST1H1B","HIST1H1D", #2
          "PAX2","JAG1", #0
          "WFDC2","FBXO2","HES1","OC90", #1
          "OTOL1","AGR2","BRICD5","USH1C", #6
          "SPARCL1","CXCL14","CALB1", #10
          "PCP4","HES6","CALM2","MYO6","MYO15A","POU4F3","LMO7", #14
          "S100A2","S100A6","S100A10","S100A11","KRT8","KRT18","KRT19", #8
          "WNT2B","DACH1", #9
          "PTGDS", #5
          "HERPUD1","DDIT3","NR4A1","GADD45B","GADD45A", #12
          "NEUROD1","DCX","GAP43","ELAVL4", #3
          "S100B","PLP1","PMP22","MPZ", #7
          "MGP","FN1","COL1A1","COL1A2","COL3A1", #4
          "TFAP2B","SLC2A3","ARHGAP29" #11
          )),
        col.min = -2, col.max = 2, dot.min = 0, dot.scale = 6
        ) +
  scale_y_discrete(limits = c(
    "13 G2/M",
    "2 G1/S",
    "0 otic prog",
    "1 otic epi",
    "6 SC",
    "10 SC-HC",
    "14 HC",
    "8 KRT+",
    "9 end-lymph.du",  
    "5 dark?",
    "12 stressed",
    "3 n.prog",
    "7 Schwann",
    "4 fibrobl", 
    "11 n.crest"
    ), position = "right"
  ) +
  coord_flip() +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()
####

Wilcox_tab <- FindAllMarkers(day2030405060, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1,
                             verbose = T, only.pos = T, return.thresh = 0.01) %>%
  select(c("cluster", "gene", "pct.1", "pct.2", "p_val_adj", "avg_log2FC")) %>%
  group_by(cluster) %>%
  top_n(50) %>%
  arrange(cluster, desc(avg_log2FC))

