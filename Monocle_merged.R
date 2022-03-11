library(Seurat)
library(SeuratObject)
library(monocle3)
library(tidyverse)


#### Path ####
Path <- "C:/Users/yosueda/Documents/RStudio/"
Proj <- "SOX2 project/Combine"


#### Load megrged Seurat object ####
day2030405060 <- readRDS(paste0(Path, Proj, "/Output/day2030405060_dim9_res0.7.rds"))
# Select otic-related clusters 
day2030405060 <- subset(day2030405060, ident = c("7","6","0","2","9","4","1","13"))


#### Seurat object --> Monocle CellDataSet
#### (Referred to: "https://github.com/cole-trapnell-lab/monocle-release/issues/388")
# Compose CellDataSet
expression_matrix <- day2030405060@assays$RNA@counts
cell_metadata <- day2030405060@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(day2030405060@assays$RNA), row.names = rownames(day2030405060@assays$RNA))
cds_day2030405060 <- new_cell_data_set(expression_matrix,
                                       cell_metadata = cell_metadata,
                                       gene_metadata = gene_annotation)

# Transfer Seurat meta data
reducedDim(cds_day2030405060, type = "PCA") <- day2030405060@reductions$pca@cell.embeddings
cds_day2030405060@preprocess_aux$prop_var_expl <- day2030405060@reductions$pca@stdev
cds_day2030405060@int_colData@listData$reducedDims$UMAP <- day2030405060@reductions$umap@cell.embeddings
cds_day2030405060@clusters$UMAP_so$clusters <- day2030405060@meta.data$seurat_clusters


#### Monocle pipeline ####
# preprocessing and clustering
cds_day2030405060 <- preprocess_cds(cds_day2030405060, num_dim = 9)
cds_day2030405060 <- cluster_cells(cds_day2030405060, reduction_method = "UMAP", 
                                   k = 250, resolution = 0.004, verbose = T)

# calculate trajectory
cds_day2030405060 <- learn_graph(cds_day2030405060, use_partition = F, close_loop = F, verbose = T)

# Check trajectory
plot_cells(cds_day2030405060, label_groups_by_cluster = F,
           show_trajectory_graph = T, color_cells_by = "cluster", cell_size = 1)

# Order cells after checking cycling clusters
plot_cells(cds_day2030405060, genes=c("TOP2A"),
           label_cell_groups = T, show_trajectory_graph = T, cell_size = 1)
plot_cells(cds_day2030405060, genes=c("UBE2C"),
           label_cell_groups = T, show_trajectory_graph = T, cell_size = 1)
plot_cells(cds_day2030405060, genes=c("KI67"),
           label_cell_groups = T, show_trajectory_graph = T, cell_size = 1)
cds_day2030405060 <- order_cells(cds_day2030405060, reduction_method = "UMAP", verbose = T)  

# Save dataset
saveRDS(cds_day2030405060, paste0(Path, Proj, "/Output/cds_day2030405060_9_0.7_sub.rds"))

#### Plot trajectory and pseudotime ####
# Trajectory colored with pseudotime
tiff(paste0(Path, Proj, "/Output/pseudotime.tiff"), width = 8, height = 8, units = "in", res = 600)
plot_cells(cds_day2030405060, show_trajectory_graph = T, color_cells_by = "pseudotime",
           cell_size = 1, group_label_size = 3,
           trajectory_graph_segment_size = 1.5) +
  scale_color_gradient(low = "#00ff00", high = "#ff0000") +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x = element_line(size = 1),
        axis.line.y = element_line(size = 1),
        legend.position = "none",
        panel.background = element_blank()
  )
dev.off()

# Trajectory color with days
tiff(paste0(Path, Proj, "/Output/day.tiff"), width = 8, height = 8, units = "in", res = 600)
plot_cells(cds_day2030405060, show_trajectory_graph = T, color_cells_by = "Day",
           cell_size = 1, group_label_size = 2,
           trajectory_graph_segment_size = 1.5,
           label_cell_groups = F) +
  scale_color_manual(values = c("#F76BF3", "#00B0F6", "#00BF7D", "#A3A500", "#F8766D")) +
  theme(axis.text.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line.x = element_line(size = 1),
      axis.line.y = element_line(size = 1),
      legend.position = "none",
      panel.background = element_blank()
)
dev.off()

# Trajectory color with Seurat clusters
tiff(paste0(Path, Proj, "/Output/SeuratCluster.tiff"), width = 8, height = 8, units = "in", res = 600)
plot_cells(cds_day2030405060, show_trajectory_graph = T, cell_size = 1, group_label_size = 5,
           color_cells_by = "seurat_clusters",
           trajectory_graph_segment_size = 1.5,
           label_cell_groups = F) + 
  scale_color_manual(values = c("#F6C141","#E8601C","#F1932D","#AE76A3",
                                "#F7F056","#CAE0AB","#D1BBD7","#DC050C")) +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x = element_line(size = 1),
        axis.line.y = element_line(size = 1),
        legend.position = "none",
        panel.background = element_blank()
  )
dev.off()

# Pseudotime plot with otic markers
ot_genes <- c("DLX5", "LMX1A", "TBX2", "PAX2", "JAG1",
              "OC90", "SPP1", "PRSS8", "FBXO2", "WFDC2", 
              "COL9A2", "SCG5", "BRICD5", "OTOL1", "SPARCL1",
              "OTOGL", "USH1C", "MYO6", "PCP4", "LMO7", "HES6", "CALM2",
              "ESPN", "MYO15A", "OTOF", "POU4F3", "ATOH1", "MYO7A", "BDNF")
for (i in ot_genes) {
  ot_lineage_cds <- cds_day2030405060[rowData(cds_day2030405060)$gene_short_name %in% i,
                                      colData(cds_day2030405060)$seurat_clusters %in% c("7","6","0","2","9","4","1","13")]
  tiff(paste0(Path, Proj, "/Output/pt_", i, ".tiff"), width = 2.5, height = 1.5, units = "in", res = 600)
  print(
    plot_genes_in_pseudotime(ot_lineage_cds,
                             color_cells_by = "pseudotime",
                             min_expr = 0.5,
                             cell_size = 0.2) + 
      scale_color_gradient(low = "#00ff00", high = "#ff0000") +
      theme(axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line(size = 5),
            legend.position = "none",
            panel.background = element_blank()
            )) 
dev.off()


#### Subsetting the last branch to supporting cells and hair cells ####
temp <- choose_graph_segments(cds_day2030405060)   # Result of choose_graph_segments won't be used; just an index for choose_cells
cds_subset <- choose_cells(cds_day2030405060)      # Choose cells based on choose_graph_segments result

saveRDS(paste0(Path, Proj, "/Output/cds_subset.rds")

# check subset colored by pseudotime
tiff(paste0(Path, Proj, "/Output/subset", ".tiff"), w0idth = 8, height = 8, units = "in", res = 600)
plot_cells(cds_subset,
           color_cells_by = "pseudotime",
           label_cell_groups=F,
           label_leaves=F,
           label_branch_points=F,
           show_trajectory_graph = F,
           cell_size = 1.5) +
  scale_color_gradient(low = "#00ff00", high = "#ff0000") +
  theme(axis.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x = element_line(size = 1),
        axis.line.y = element_line(size = 1),
        legend.position = "none",
        panel.background = element_blank()
  )
dev.off()


#### Gene expression module calculation ####
subset_pr_test_res <- graph_test(cds_subset, neighbor_graph = "principal_graph", k = 70, verbose = T)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], weight = T, cores = 1,
                                    umap.min_dist = 0.25, umap.n_neighbors = 10, k = 70,
                                    partition_qval = 0.01, resolution = 0.005, # max modality = 0.01
                                    random_seed = 42, verbose = T)


# plot all module
tiff(paste0(Path, Proj, "/Output/module.tiff"), width = 8, height = 8, units = "in", res = 600)
plot_cells(cds_subset, genes = gene_module_df,
           show_trajectory_graph = F, label_cell_groups=F, min_expr = 50,
           label_principal_points = T, cell_size = 1) +
  scale_color_gradientn(colors=c("#0000ff", "#ff0000")) +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5),
        panel.background = element_blank()
  )
dev.off()

# plot only module #4 (hair cell enriched) and #5 (supporting cell enriched)
tiff(paste0(Path, Proj, "/Output/module_4_5.tiff"), width = 8, height = 8, units = "in", res = 600)
  print(
  plot_cells(cds_subset, genes = filter(gene_module_df, module == c(4,5)),
           show_trajectory_graph = F, label_cell_groups=F, min_expr = 1,
           label_principal_points = T, cell_size = 1) +
  scale_color_gradientn(colors=c("#0000ff", "#ff0000")) +
  theme(axis.text.y = element_text(size = 32),
        axis.text.x = element_text(size = 32),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x = element_line(size = 1),
        axis.line.y = element_line(size = 1),
        legend.position = "none",
        panel.background = element_blank()
  )
  )
dev.off()

# Out put gene lists
gene_sub_4 <- subset(gene_module_df, module == 4)
gene_sub_5 <- subset(gene_module_df, module == 5)
gene_sub_9 <- subset(gene_module_df, module == 9)
write.csv(gene_sub_4, paste0(Path, Proj, "/Output/gene_sub_4.csv")
write.csv(gene_sub_5, paste0(Path, Proj, "/Output/gene_sub_5.csv")
write.csv(gene_sub_9, paste0(Path, Proj, "/Output/gene_sub_9.csv")

# save dataset
saveRDS(day2030405060, paste0(Path, Proj, "/Output/merged_1_13.rds")


# remove genes expressed in small number of cells
filtering_rate <- 0.05   # genes expressed in more than this rate will be selected

# module #4
cnt <- NULL; for(j in 1:nrow(gene_sub_4)) {
  genename <- gene_sub_4$id[j]
  cnt <- c(cnt,                                               # make the vector of the number of non-zero cells for genes in gene_sub_n (n = 4, 5)
           sum(cds_subset@assays@data$counts[genename,] > 0)) # <- count the number of non-zero cells 
}
gene_sub_4_1 <- mutate(gene_sub_4, cnt) %>% filter(cnt > length(colnames(cds_subset))*filtering_rate)
write.csv(gene_sub_4_1, paste0(Path, Proj, "/Output/gene_sub_4_", filtering_rate, "-filtered.csv"))

# module #5
cnt <- NULL; for(j in 1:nrow(gene_sub_5)) {
  genename <- gene_sub_5$id[j]
  cnt <- c(cnt,                                               # make the vector of the number of non-zero cells for genes in gene_sub_n (n = 4, 5)
           sum(cds_subset@assays@data$counts[genename,] > 0)) # <- count the number of non-zero cells 
}
gene_sub_5_1 <- mutate(gene_sub_5, cnt) %>% filter(cnt > length(colnames(cds_subset))*filtering_rate)
write.csv(gene_sub_5_1, paste0(Path, Proj, "/Output/gene_sub_5_", filtering_rate, "-filtered.csv"))

# module #9
cnt <- NULL; for(j in 1:nrow(gene_sub_9)) {
  genename <- gene_sub_9$id[j]
  cnt <- c(cnt,                                               # make the vector of the number of non-zero cells for genes in gene_sub_n (n = 4, 5)
           sum(cds_subset@assays@data$counts[genename,] > 0)) # <- count the number of non-zero cells 
}
gene_sub_9_1 <- mutate(gene_sub_9, cnt) %>% filter(cnt > length(colnames(cds_subset))*filtering_rate)
write.csv(gene_sub_9_1, paste0(Path, Proj, "/Output/gene_sub_9_", filtering_rate, "-filtered.csv"))


#### Plot all genes stored in filtered "gene_sub_4_1", "gene_sub_5_1", "gene_sub_9_1" ####
# filtered #4
cnt <-0
for(i in gene_sub_4_1$id){
plot_cds <- cds_day2030405060[rowData(cds_day2030405060)$gene_short_name %in% i,
                              colData(cds_day2030405060)$seurat_clusters %in% c("7","6","0","2","9","4","1","13")]
tiff(paste0(Path, Proj, "/Output/pt_set_", i, ".tiff"), width = 2.5, height = 1.5, units = "in", res = 300)
print(
  plot_genes_in_pseudotime(plot_cds,
                           color_cells_by = "pseudotime",
                           min_expr = 0.5,
                           cell_size = 0.2) + 
    scale_color_gradient(low = "#00ff00", high = "#ff0000") +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_line(size = 5),
          legend.position = "none",
          panel.background = element_blank()
    )
)
dev.off()
cnt <- cnt+1
message(paste0(cnt,"/",length(gene_sub_4_1$id)," done (",i,")"))
}

# filtered #5
cnt <-0
for(i in gene_sub_5_1$id){
  plot_cds <- cds_day2030405060[rowData(cds_day2030405060)$gene_short_name %in% i,
                                colData(cds_day2030405060)$seurat_clusters %in% c("7","6","0","2","9","4","1","13")]
  tiff(paste0(Path, Proj, "/Output/pt_set_", i, ".tiff"), width = 2.5, height = 1.5, units = "in", res = 300)
  print(
    plot_genes_in_pseudotime(plot_cds,
                             color_cells_by = "pseudotime",
                             min_expr = 0.5,
                             cell_size = 0.2) + 
      scale_color_gradient(low = "#00ff00", high = "#ff0000") +
      theme(axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line(size = 5),
            legend.position = "none",
            panel.background = element_blank()
      )
  )
  dev.off()
  cnt <- cnt+1
  message(paste0(cnt,"/",length(gene_sub_5_1$id)," done (",i,")"))
}

# filtered #9
cnt <-0
for(i in gene_sub_9_1$id){
  plot_cds <- cds_day2030405060[rowData(cds_day2030405060)$gene_short_name %in% i,
                                colData(cds_day2030405060)$seurat_clusters %in% c("7","6","0","2","9","4","1","13")]
  tiff(paste0(Path, Proj, "/Output/pt_set_", i, ".tiff"), width = 2.5, height = 1.5, units = "in", res = 300)
  print(
    plot_genes_in_pseudotime(plot_cds,
                             color_cells_by = "pseudotime",
                             min_expr = 0.5,
                             cell_size = 0.2) + 
      scale_color_gradient(low = "#00ff00", high = "#ff0000") +
      theme(axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line(size = 5),
            legend.position = "none",
            panel.background = element_blank()
      )
  )
  dev.off()
  cnt <- cnt+1
  message(paste0(cnt,"/",length(gene_sub_9_1$id)," done (",i,")"))
}


# End of code
