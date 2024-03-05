# Heatmap E17 - Allie's Script
# Created: 3/1/2024
# Last edited: 3/1/2024

# Packages
require(tidyverse)
require(RColorBrewer)
require(ggplot2)
require(ComplexHeatmap)

# Load object and gene list
load("/Volumes/Rachel Rivero/Manuscript Worthy/E17/E17_labeled.Robj") 
load("/Volumes/Rachel Rivero/Manuscript Worthy/E17/FRL17.Seurat_GeneList.Robj")

object <- FRL17.Seurat

# Heatmap 
Idents(object) <- object$Final35
object2 <- subset(object, downsample = 100)
object2$Final35 <- factor(object2$Final35,levels=c('Alveolar','Cycling','Interstitial'))
object2$Sample <- factor(object2$Sample,levels=c('FB13','FB14','BAL','BEF12','BEFM1',
                                                 'BEFM2','BEFM4','BEFM5','BEFM6'))
object2$Dataset3 <- factor(object2$Dataset3,levels=c('Start','Eng'))
object2 <- ScaleData(object2,features = rownames(object2))
Idents(object2) <- object2$Final35
genes2 <- c('Bst2','Cd9','Xrcc5','Cdc42ep3','Pla2g2d','Slc39a2','Lpl','Rgcc','Ly6al',
            'Crip1','S100a4','Scgb1a1',  # Start_Alveolar
            'Stmn1','Tuba1b','Lgals1','Tmem176b','Tmem176a','Limd2','Ccl2','Cd14','Ccl7',
            'Ccl12','Spp1','Ms4a7',  # Start_Interstitial
            'Slc11a1','Fabp4','Mt1','Mt2A','Cxcl2','Cxcl1','Fos','Clec10a','Vsig4','Retn',
            'Clca4l','S100a9','Lcn2','Mgp',  # Eng_Alveolar
            'Ctsl','Ifitm1','Hopx','Pla2g7','Tf','Slpi','Gsta1','F13a1','Gpx3','Sparc',
            'Col1a1')  # Eng_Interstitial
gene_split2 <- c(rep('Start Alveolar',12),rep('Start Interstitial',12),
                 rep('Engineered Alveolar',14),rep('Engineered Interstitial',11))
genes.use2 <- intersect(rownames(GetAssayData(object2,slot = 'scale.data')),genes2)
object.cells2 <- GetAssayData(object2, slot = "scale.data")[genes.use2, ]
object.cells2 <- t(scale(t(object.cells2)))
unityNormalize <- function(x){(x-min(x))/(max(x)-min(x))}
object.cells2 <- t(apply(as.matrix(object.cells2), 1, unityNormalize))
cellorder2 <- data.frame(Dataset=object2$Dataset3,Final35=object2$Final35,Sample=object2$Sample)
cellorder2 <- cellorder2[order(cellorder2[,1],cellorder2[,2],cellorder2[,3]),]
object.cells2 <- object.cells2[,rownames(cellorder2)]
object.cells2 <- object.cells2[genes2,]
dim(object.cells2)

Dataset <- as.matrix(object2$Dataset3)
Dataset <- as.matrix(Dataset[rownames(cellorder2),])
dataset_split <- c(rep('Start',301),rep('Eng',99))
dataset_split2 <- c(rep('Start_Alveolar',150),rep('Start_Interstitial',76),
                    rep('Eng_Alveolar',50),rep('Eng_Interstitial',24))
CellType <- as.matrix(object2$Final35)
CellType <- as.matrix(CellType[rownames(cellorder2),])
Sample <- as.matrix(object2$Sample)
Sample <- as.matrix(Sample[rownames(cellorder2),])
colors.inferno2 <- colorRamp2(breaks = c(seq(min(object.cells2),max(object.cells2),length.out=60)), inferno(n=60), space = "RGB")
ha2 = HeatmapAnnotation(df = list(data.frame(Dataset),
                                  data.frame(CellType),
                                  data.frame(Sample)),
                        col = list(Dataset=starteng_colors,
                                   CellType=type_colors,
                                   Sample=sample_colors),
                        show_annotation_name = c(T),show_legend = c(F),
                        annotation_name_gp = gpar(fontsize = 16))
lgd1 = Legend(col_fun = colors.inferno2,
              title = "Gene\nexpression", at = c(0,1),
              title_gp = gpar(fontface = "bold",fontsize = 18),
              labels = c("Low","High"),
              labels_gp = gpar(fontsize = 16),
              grid_width = unit(1.5, "cm"),
              legend_height = unit(4, "cm"))
lgd3 = Legend(at = levels(object2$Sample),
              legend_gp = gpar(fill=sample_colors),
              title = "Sample",
              title_gp = gpar(fontface = "bold",fontsize = 18),
              labels_gp = gpar(fontsize = 16),
              grid_width = unit(1, "cm"),
              grid_height = unit(1, "cm"))
lgd4 = Legend(at = levels(object2$Dataset3),
              legend_gp = gpar(fill=starteng_colors),
              title = "Dataset",
              title_gp = gpar(fontface = "bold",fontsize = 18),
              labels_gp = gpar(fontsize = 16),
              grid_width = unit(1, "cm"),
              grid_height = unit(1, "cm"))
lgd5 = Legend(at = levels(object2$Final35),
              legend_gp = gpar(fill=type_colors[levels(object2$Final35)]),
              title = "Cell Type",
              title_gp = gpar(fontface = "bold",fontsize = 18),
              labels_gp = gpar(fontsize = 16),
              grid_width = unit(1, "cm"),
              grid_height = unit(1, "cm"))
pd2 = packLegend(lgd4,lgd5,lgd3,lgd1,column_gap = unit(1, "cm"),max_height = unit(30, "cm"))

object.heatmap2 <- Heatmap(as.matrix(object.cells2),
                           col = colors.inferno2,
                           top_annotation = ha2,
                           column_split = factor(dataset_split2,levels=levels(object2$Final4)),
                           row_split = factor(gene_split2,levels=unique(gene_split2)),
                           row_title_gp = gpar(fontsize = 25),
                           column_title_gp = gpar(fontsize = 0),
                           row_names_gp = gpar(fontsize = 20,fontface='italic'),
                           column_title_rot = 0,
                           cluster_rows=F,
                           cluster_columns=F,
                           cluster_row_slices=F,
                           show_row_names=T,
                           show_column_names = F,
                           show_heatmap_legend = F)
p5 = grid.grabExpr(draw(object.heatmap2,annotation_legend_list = pd2,padding = unit(c(2, 2, 10, 2), "mm")))

res <- 400
png(paste(dir_name,'_heatmap_6.png',sep=''), width=800*res/72, height=1100*res/72,res=res)
draw(object.heatmap2,annotation_legend_list = pd2)
dev.off()
