# Mes Heatmap
# Created: 3/6/24
# Modified: 3/6/24

# Set WD
setwd("/Volumes/Rachel Rivero/Manuscript Worthy/E17/figures")

# Load packages
require(Seurat)
require(dplyr)

# Load local functions (change filepath to wherever you have saved 'CustomHeatmapOnly3.R'on your computer)
source('/Users/rr792/Documents/GitHub/NICHESMethods/R Scripts/Rachel Rivero (Sandbox)/CustomHeatmapOnly3_RRedits.R')

# Load data
load("/Volumes/Rachel Rivero/Manuscript Worthy/E17/mes_labeled.Robj")

# Inspect data to get a handle on the metadata slot names
names(mes@meta.data)

# Look at the ordering of different metadata slots
table(mes$CellType,
      mes$CellClass)

# Re-order the metadata for plotting, if desired (optional)
mes$CellClass <- factor(mes$CellClass,
                                 levels=c('Mesenchymal')) 
mes$CellType <- factor(mes$CellType,
                                levels=c('Sox9+ Mesenchymal Progenitors', 'Sox9- Mesenchymal Progenitors',
                                         'Wnt2+ Mes', 'Proliferating Wnt2+ Mes','Prex2+ Mes', 
                                         'SMCs/Myofibroblasts','Il17+ SMCs/Myofibroblasts',  
                                         'Pericytes 1','Pericytes 2', 
                                         'Mesothelium','NCCs', 'NECs')) # feel free to modify the order here

# Look at the effect of re-ordering
table(mes$CellType,
      mes$CellClass)

# Load color palette
# load("/Volumes/Rachel Rivero/Manuscript Worthy/E17/figures/E17colorpalette.Robj")

# # Set up color palette for all your plotting (customize for a given paper/project)
col.pal.mes <- list()
col.pal.mes$CellClass <- c('#3C5488FF')
col.pal.mes$CellType <- c('cornflowerblue','red','mediumseagreen','#B0DB43','orange', 
                      'slategray1','deepskyblue', 'yellow','purple',
                      '#14591D','black','darkgoldenrod')
col.pal.mes$Condition <- c('#F1C40F', '#2F6690')
col.pal.mes$Sample <- c('#BAB700','#EDE3E4','#FF5E5B','#00CECB')

names(col.pal.mes$CellClass) <- c('Mesenchymal') 
names(col.pal.mes$CellType) <- c('Sox9+ Mesenchymal Progenitors', 'Sox9- Mesenchymal Progenitors',
                                 'Wnt2+ Mes', 'Proliferating Wnt2+ Mes','Prex2+ Mes', 
                             'SMCs/Myofibroblasts','Il17+ SMCs/Myofibroblasts',  
                             'Pericytes 1','Pericytes 2', 
                             'Mesothelium','NCCs', 'NECs')
names(col.pal.mes$Condition) <- c('CDH','Normal') 
names(col.pal.mes$Sample) <- c('cdhfrl1','cdhfrl2','frl1','frl2')
 
save(col.pal.mes, file = 'mescolorpalette.Robj')

# Create a downsampled object (optional, can be useful both for speed and also to make plot clearer)
table(mes$Condition)
Idents(mes) <- mes$Condition # set idents to condition [initializing like this will preserve relative cell type ratios in the next line. alternatively, you could choose to standardize number of columns per cell type (here, we standardize number of columns per Condition category)]
downsampled <- subset(mes, downsample = 7522) # same number of cells from each condition, but allows celltype distributions to stay real
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
row.labels <- c('Frem1','Wnt2','Cdc20','Bnc2','Prex2','Myh11','Fgf18','Il17b','Lgr6',
                'Gucy1a2','Cxcr4','Pdgfrb','Capn6','Matn4','Sox9','Col9a2','Prrx1',
                'Lum','Ebf2','Dcn','Ntrk2','Twist1','Wt1','Bnc1') 

# Make a PNG output plot at 300 dpi for publication, width and height in real-world inches (feel free to change width and height as desired)

png(file = 'Heatmap_Demo_Mes.png', width=8, height=6,units = 'in',res=300)
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'CellClass' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'CellType' , # sets the second metadata slot to use (secondary column grouping)
                   tertiary = 'Condition' , # sets the third metadata slot to use (tertiary column grouping)
                   quarternary = 'Sample' , 
                   primary.cols = col.pal.mes$CellClass, # Needs to be a named list of colors, defaulted
                   secondary.cols = col.pal.mes$CellType, # Needs to be a named list of colors
                   tertiary.cols = col.pal.mes$Condition, # Needs to be a named list of colors
                   quarternary.cols = col.pal.mes$Sample, # not used in this specific function
                   features = GOI, # defined above, customize
                   labels = c('CellClass','CellType','Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=row.labels,
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.5) # This can be tuned to change the plotting color range, which will change the look of the plot. Experiment.
dev.off()

