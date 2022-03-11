##-------------------------------------------------------------
## Simulation Plot: Bubble plot
##-------------------------------------------------------------
library(data.table)
library(iDEA)
library(tidyverse)
library(stringr)
library(RColorBrewer)


#### Input path ####
Proj <- "SOX2 project"
File_path <- "C:/Users/yosueda/Documents/RStudio/"


customGeneSetsInfo <- readRDS(file = "C:/Users/yosueda/Documents/RStudio/SOX2 project/Combine/20220112/Output/customGeneSetsInfo.rds")

plotdata<- fread(paste0(Proj, file_path, "/Output/iDEA_dn.csv")
bp.data <- data.table::merge.data.table(plotdata, customGeneSetsInfo, by.x = "annot_id", by.y = "gset")
write.csv(bp.data, paste0(Proj, file_path, "/Output/iDEA_dn_fl.csv")

bp.data$Category <- droplevels(bp.data$gsetBioName)


#### Bubble plot ####

includedCats <- c("GO BIOLOGICAL PROCESS",
                  "GO MOLECULAR FUNCTION",
                  "GO CELLULAR COMPONENT",
                  "REACTOME",
                  "KEGG",
                  "PID",
                  "POSITIONAL",
                  "TRANSCRIPTION FACTORS")

caseCats <- c("GO biological process",
              "GO molecular function",
              "GO cellular component",
              "Reactome",
              "KEGG",
              "PID")

bp.data <- bp.data[toupper(bp.data$Category) %in% includedCats, ] #subset data, removing all geneset categories not in includedCats
bp.data$Category <- droplevels(bp.data$Category) %>% factor(caseCats) # Drop unused levels, reorder factor so that plot legend matches x position
bp.data <- bp.data[order(match(toupper(bp.data$Category), includedCats)), ] #Organizes the data so that categories are grouped on the x axis
bp.data$IDNum <- row.names(bp.data) %>% as.integer() #IDNum will define x position on bubble plot
bp.data$Log10_Pvalue_Louis <- -1*log10(bp.data$pvalue_louis)

##-------------------------------------------------------------
## Configurable Variables
##-------------------------------------------------------------

bp.colors <- c("#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69")

labeled.genesets <-
  c("GO_EPITHELIAL_CELL_DIFFERENTIATION","GO_POSITIVE_REGULATION_OF_CELL_CYCLE",
    "GO_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION","GO_STEM_CELL_DIFFERENTIATION",
    "GO_ACTIN_FILAMENT_BUNDLE","GO_ACTOMYOSIN")

included.categories <- c("GO biological process",
                         "GO molecular function",
                         "GO cellular component",
                         "Reactome",
                         "KEGG",
                         "PID")

data.to.plot <- bp.data

##-------------------------------------------------------------
## Create data frame for for plot labels
##-------------------------------------------------------------

Sig <- data.to.plot[which(data.to.plot$annot_id %in% labeled.genesets),]

Sig$Term <- tolower(Sig$annot_id)

##-------------------------------------------------------------
## ggplot function for bubble plot
##-------------------------------------------------------------
library(ggplot2)
library(ggrepel)

bp <- ggplot(filter(data.to.plot, Category %in% included.categories), aes(x = IDNum, y = Log10_Pvalue_Louis, color = Category)) +
  geom_point(shape = 19, alpha=1, size = 3) + 
  labs(x = "",
       y = expression(paste(-log[10], "(p-value)"))) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24),
        axis.line = element_line(colour = 'black'),
        axis.line.x = element_line(size = 2),
        axis.line.y = element_line(size = 2),
        axis.ticks = element_line(colour = 'black'),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 24, face = 'bold'),
        legend.title = element_text(size = 18,face = 'bold'),
        legend.text = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.key = element_rect(color = "transparent", fill = "transparent")) +
  geom_hline(yintercept = 2, col = 'black', linetype = 2, size=2) +
  scale_color_manual(values = bp.colors) +
  theme(legend.direction = "vertical") +
  theme(legend.position = c(0.1, 0.8)) +
  theme(legend.box = "horizontal") +
  theme(legend.title.align = 0) +
  geom_label_repel(
    data = Sig,
    aes(label = Term),
    col = 'black',
    alpha = 1,
    size = 6,
    nudge_y = 10)


#### Output bubble plot ####

tiff("C:/Users/yosueda/Documents/RStudio/SOX2 project/Combine/20220228/Output/bubble_SC.tiff",
     width = 32, height = 8, units = "in", #pointsize = 6,
     res = 600)
bp
dev.off()

setEPS(width = 32, height = 8)
postscript("C:/Users/yosueda/Documents/RStudio/SOX2 project/Combine/20220228/Output/bubble_SC.eps")
bp
dev.off()
