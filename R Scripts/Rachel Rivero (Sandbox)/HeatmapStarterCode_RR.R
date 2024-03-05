# This demonstrates how to make a heatmap on phenotype data

# Set WD
setwd("/Volumes/Rachel Rivero/Manuscript Worthy/E17/figures")

# Load packages
require(Seurat)
require(dplyr)

# Load local functions (change filepath to wherever you have saved 'CustomHeatmapOnly3.R'on your computer)
source('/Users/rr792/Documents/GitHub/NICHESMethods/R Scripts/Rachel Rivero (Sandbox)/CustomHeatmapOnly3_Phenotype.R')

# Load data
load("/Volumes/Rachel Rivero/Manuscript Worthy/E17/E17_labeled.Robj")

# Inspect data to get a handle on the metadata slot names
names(FRL17.Seurat@meta.data)

# Look at the ordering of different metadata slots
table(FRL17.Seurat$CellType,
      FRL17.Seurat$CellClass)

# Re-order the metadata for plotting, if desired (optional)
FRL17.Seurat$CellClass <- factor(FRL17.Seurat$CellClass,
                                 levels=c('Mesenchymal','Epithelial','Immune','Endothelial')) # feel free to modify the order here
# FRL17.Seurat$CellType <- factor(FRL17.Seurat$CellType,
#                                 levels=c('Epithelial','Endothelial','LEC',
#                                         'Mesothelial','Pericytes','Mesenchyme','Proliferating Mesenchyme','Myofibroblasts/SMCs','Mesenchymal Progenitors','NCC',
#                                          'Immune')) # feel free to modify the order here

# Look at the effect of re-ordering
table(FRL17.Seurat$CellType,
      FRL17.Seurat$CellClass)

# Set up color palette for all your plotting (customize for a given paper/project)
col.pal <- list()
col.pal$CellType <- c('firebrick2','mediumseagreen','#B0DB43','#9C528B', 'slategray1','#2F6690',
                                  '#F1C40F','#0F0326','#E65F5C','#14591D','#726DA8',
                                  'cornsilk3', 'blue', 'deepskyblue',
                                  'darkorange1', 'cyan', 'purple', 'yellow', 'deeppink', '#A40606','pink')
col.pal$Condition <- c('#F1C40F', '#2F6690') 
col.pal$Sample <- c('orange','purple','yellow','black') 
names(col.pal$CellType) <- c('Club cells','Early AT2','Il17+ SMCs/Myofibroblasts',
                             'Mesothelium', 'Pericytes 2', 'Proliferating Wnt2+ Mes',
                             'Sox9- Mesenchymal Progenitors', 'Sox9+ Mesenchymal Progenitors',
                             'Transitional Epithelium', 'Wnt2+ Mes', 'Arterial', 'APCs',
                             'Pericytes 1', 'Immune', 'LEC', 'NCCs', 'NECs', 'Non-ciliated Secretory',
                             'SMCs/Myofibroblasts', 'Transitional Endothelial')
names(col.pal$Condition) <- c('CDH','Normal') #must be named correctly
names(col.pal$Sample) <- c('cdhfrl1','cdhfrl2','frl1','frl2') #must be named correctly

save(col.pal, file = 'colorpalette.Robj')

# Create a downsampled object (optional, can be useful both for speed and also to make plot clearer)
table(FRL17.Seurat$Condition)
Idents(FRL17.Seurat) <- FRL17.Seurat$Condition # set idents to condition [initializing like this will preserve relative cell type ratios in the next line. alternatively, you could choose to standardize number of columns per cell type (here, we standardize number of columns per Condition category)]
downsampled <- subset(FRL17.Seurat, downsample = 9213) # same number of cells from each condition, but allows celltype distributions to stay real
table(downsampled$Condition) # same number of cells from each condition now...good for our purposes here

# Scale the object (important becasue we downsampled, and also important in general for heatmaps and many visualizations)
downsampled <- ScaleData(downsampled,features = rownames(downsampled))

# Create a marker list to pick genes from
Idents(downsampled) <- downsampled$CellType
mark <- FindAllMarkers(downsampled,
                       only.pos=T, # only report positive (upregulated) markers
                       min.pct=0.5, # only report markers with at least 50 percent of each test population expressing the gene
                       logfc.threshold = 0.5) # only report markers with logFC > 0.5
mark$ratio <- mark$pct.1/mark$pct.2
mark$power <- mark$ratio*mark$avg_log2FC

# Define the features to plot. This list will be the exact order of the rows in the heatmap. You should customize this.
top.marker.list <- mark %>% group_by(cluster) %>% top_n(10,power) # This is just a quick and dirty way to demonstrate, you can modify OR ideally should manually curate a full list
GOI <- top.marker.list$gene 
GOI

# Define which rows you would like to emphasize / label by name. The order here does not matter. You should customize this.
row.labels <- c('Sox9','Sox10','Gucy1a2','Sox7','Mmrn1','Sox17','Wt1','Ptprc','Prrx1',
                'Rspo2','Cftr', 'Lum', 'Ebf2', 'Ntrk2', 'Twist1', 'Dcn', 'Il17b', 'Fgf18') # etc etc...add as you see fit...pull from marker list...experiment
  
# Make a PNG output plot at 300 dpi for publication, width and height in real-world inches (feel free to change width and height as desired)

png(file = 'Heatmap_Demo_Output.png', width=8, height=6,units = 'in',res=300)
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'CellType' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'Condition' , # sets the second metadata slot to use (secondary column grouping)
                   tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$CellType, # Needs to be a named list of colors
                   secondary.cols = col.pal$Condition, # Needs to be a named list of colors
                   tertiary.cols = col.pal$Sample, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = GOI, # defined above, customize
                   labels = c('CellType','Condition','Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=row.labels,
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.5) # This can be tuned to change the plotting color range, which will change the look of the plot. Experiment.
dev.off()

