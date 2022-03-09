library(Seurat)
library(tidyverse)


# Path
Path <- "C:/Users/yosueda/Documents/RStudio/"
Proj <- "SOX2 project/Combine"


# Load 10X data in Seurat objects
d20 <- Read10X(paste0(Path, Proj, "/Data/d20/")) %>% CreateSeuratObject(); d20$Day <- "20"
d30 <- Read10X(paste0(Path, Proj, "/Data/d30/")) %>% CreateSeuratObject(); d30$Day <- "30"
d40 <- Read10X(paste0(Path, Proj, "/Data/d40/")) %>% CreateSeuratObject(); d40$Day <- "40"
d50 <- Read10X(paste0(Path, Proj, "/Data/d50/")) %>% CreateSeuratObject(); d50$Day <- "50"
d60 <- Read10X(paste0(Path, Proj, "/Data/d60/")) %>% CreateSeuratObject(); d60$Day <- "60"


#### Filtering functions ####
# Filters for (log-transformed) mitochondrial percentage genes
filterPlot <- function(obj) {
  obj[["mt.per"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  mt.per.vec <- obj$mt.per
  obj <- subset(obj, mt.per < 15)
  return(obj)
}

# Filtering for housekeeping gene expression
# Appends metadata column "hk" that stores the log normalized expression of a selected housekeeping gene
# if subset = TRUE, will subset the seurat object discarding cells ~std.dev~
# standard deviations away from the mean 
housekeeping.filter <- function(object, housekeeping.gene = NULL, 
                                subset = F, std.dev = 2) {
  #Sanitize input
  housekeeping.gene <- as.character(housekeeping.gene)
  # Creates a vector containing the reads found in each cell and log normalizing
  expression.vector <- object@assays$RNA@counts[housekeeping.gene, ] %>% log1p()
  #Add metadata column to the Seurat object to facilitate subsetting
  object <- AddMetaData(object, metadata = expression.vector,
                        col.name = "hk")
  # Determines whether or not the object should be subset such that
  # cells greater or less than ~std.dev~ standard deviations are discarded
  if (subset){
    object <- subset(object, hk > mean(object$hk) - std.dev*sd(object$hk) &
                       hk < mean(object$hk) + std.dev*sd(object$hk))
  }
  return(object)
}


#### Filtering ####
# GroupDatasets into a list #
list.all <- list(d20, d30, d40, d50, d60)

# Filters each dataset in list.all according to mitochondrial percentage,
# gives number of cells removed for each dataset
counts.list.pre <- sapply(list.all, ncol)
list.all <- lapply(list.all, filterPlot)
counts.list.post <- sapply(list.all, ncol)
print(counts.list.pre - counts.list.post)

# Filtering based on housekeeping gene expression; RPL27 was used here
list.all <- lapply(list.all, housekeeping.filter, housekeeping.gene = "RPL27", subset = T)
counts.list.post.2 <- sapply(list.all, ncol)
print(counts.list.post - counts.list.post.2)


#### Merge datasets ####
day2030405060 <- merge(list.all[[1]], c(list.all[[2]], list.all[[3]], list.all[[4]], list.all[[5]]))
saveRDS(day2030405060, file = paste0(Path, Proj, "/Output/day2030405060.rds"), compress = F)


#### SCTransform, PCA ####
day2030405060 <- SCTransform(day2030405060, ncells = 30000, vars.to.regress = "mt.per")
day2030405060 <- RunPCA(day2030405060)

DimPlot(day2030405060)
ElbowPlot(day2030405060, ndims = 50, reduction = "pca")


#### Construct SNN graph and UMAP ####
dim = 9
res = 0.7
day2030405060 <- FindNeighbors(day2030405060, k.param = 250, nn.method = "rann", dims = 1:dim) %>% FindClusters(resolution = res)
day2030405060 <- RunUMAP(day2030405060, dims = 1:dim, n.neighbors = 10, min.dist = 0.25)

saveRDS(day2030405060, paste0(Path, Proj, "/Output/day2030405060_dim9_res0.7.rds")


#### Find Markers ####
Wilcox_tab <- FindAllMarkers(day2030405060, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1,
                             verbose = T, only.pos = T, return.thresh = 0.01) %>%
  select(c("cluster", "gene", "pct.1", "pct.2", "p_val_adj", "avg_log2FC")) %>%
  group_by(cluster) %>% top_n(50) %>% arrange(cluster, desc(avg_log2FC))
write.csv(Wilcox_tab, paste0(Path, Proj, "/Output/Wilcox_test_dim9_res0.7.csv")


#### UMAP ####
# Whole
col <- c("#F6C141", "#E8601C", "#F1932D", "#7BAFDE", "#AE76A3", "#4EB265", "#F7F056",
         "#CAE0AB", "#90C987", "#D1BBD7", "#882E72", "#1965B0", "#5289C7", "#DC050C")
tiff(paste0(Path, Proj, "/Output/merged_dimplot.tiff"), width = 8, height = 8, units = "in", res = 600)
DimPlot(day2030405060, reduction = "umap", pt.size = 0.5, label = F, label.size = 5, repel = T, cols = col) +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        legend.position = "none"
  )
dev.off()

# Grouped by day
col <- c("#F76BF3", "#00B0F6", "#00BF7D", "#A3A500", "#F8766D")
tiff(paste0(Path, Proj, "/Output/merged_dimplot_day.tiff"), width = 8, height = 8, units = "in", res = 600)
DimPlot(day2030405060, reduction = "umap", pt.size = 0.5, label = F, label.size = 5, repel = T, group.by = "Day", cols = col) +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        plot.title = element_blank()
  )
dev.off()

# Split by day
col <- c("#F76BF3", "#00B0F6", "#00BF7D", "#A3A500", "#F8766D")
tiff(paste0(Path, Proj, "/Output/merged_dimplot_day_split.tiff"), width = 40, height = 8.2, units = "in", res = 600)
DimPlot(day2030405060,
        reduction = "umap", pt.size = 0.5, label = F, label.size = 5, group.by = "Day", split.by = "Day", cols = col) +
  theme(plot.title = element_text(size = 0),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.background = element_blank(),
        legend.position = "none"
  )
dev.off()


#### FeaturePlot ####
marker <- c(
            c("TOP2A", "PAX2", "FBXO2", "BRICD5", "POU4F3", "DACH1", "SPP1", "TFAP2B", "NEUROD1", "S100B"),
            c("OC90","TOP2A","UBE2C","MKI67","TAGLN","TWIST1","SNAI2","ZEB1","ZEB2","WNT2B")
            )
for (i in marker) {
  tiff(paste0(Path, Proj, "/Output/merged_",i,".tiff"), width = 4, height = 4, units = "in", res = 600)
  (FeaturePlot(day2030405060, features = i, pt.size = 0.4, max.cutoff = 2, cols = c("#cccccc33","#0000ff99")) +
    theme(#plot.title = element_blank(),
      plot.title = element_blank(),
      axis.text.x = element_blank(), axis.text.y = element_blank(),
      axis.line.x = element_blank(), axis.line.y = element_blank(),
      axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      panel.background = element_blank(), legend.position = "none"
    )) %>% print()
  dev.off()
}


#### DotPlot ####
tiff(
  paste0(Path, Proj, "/Output/merged_dotplot.tiff"), width = 6, height = 10, units = "in", res = 600)
DotPlot(
  day2030405060,
  features = rev(c("TOP2A","UBE2C","MKI67",
                   "PCNA","NASP","TYMS",
                   "PAX2","JAG1","TBX1",
                   "EMX2","OTX2","OC90",
                   "FBXO2","BRICD5","SPARCL1",
                   "PCP4","MYO6","POU4F3",
                   "DACH1","WNT2B",
                   "PTGDS","SPP1","CITED1",
                   "NEUROD1","ELAVL4","POU4F1",
                   "S100B","PLP1","PMP22",
                   "DSTN","TFAP2B", "TFAP2A","DCN","KRT19",
                   "KRT18","KRT8","TAGLN",
                   "COL1A2","COL1A1","POSTN","FN1","S100A6",
                   "TPM1","TRPM3","ACTA2",
                   "CDH1","CDH2","EPCAM")),
  col.min = -2, col.max = 2, dot.min = 0, dot.scale = 8) +
  coord_flip() +
  scale_y_discrete(position = "right",
                   limits = c("7","6","0","2","1","13","9","4","5","8","10","3","11","12")) +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
      axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()


# End of Code
