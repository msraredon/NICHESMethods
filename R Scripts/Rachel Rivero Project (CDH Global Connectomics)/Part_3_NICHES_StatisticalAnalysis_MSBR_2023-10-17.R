# PART 3 This script performs statistical analysis of the NICHES data

# Set WD
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Rachel/E17.E19.Merged.Explorations.2023-10-09")

# Packages
require(Seurat)
require(SeuratDisk)
require(SeuratWrappers)
require(NICHES)
require(ggplot2)
require(cowplot)
require(dplyr)
require(Matrix)

# Load from previous save point
# CellToCell data
load('epi.epi.Robj')
load('epi.mes.Robj')
load('mes.epi.Robj')
load('mes.mes.Robj')
# SystemToCell data
load('epi.niche.Robj')
load('mes.niche.Robj')
# CellToSystem data
load('epi.influence.Robj')
load('mes.influence.Robj')

##### Let's look just at mesenchymal to epithelial signaling, at first ####

# Define colors
col.pal <- list()
col.pal$cell_type <- c('#A40606','#9CFFFA','#B0DB43','#9C528B','#2F6690',
                                '#946846','#F1C40F','green','#0F0326','#E65F5C','#14591D','#726DA8',
                                'yellow','purple')
names(col.pal$cell_type) <- c(unique(mes.epi$RNA$ReceivingType),unique(mes.epi$RNA$SendingType))
                                
# Check out mesenchymal to epithelium signaling in RNA slot
mes.epi$RNA <- FindVariableFeatures(mes.epi$RNA)
mes.epi$RNA <- RunPCA(mes.epi$RNA)
mes.epi$RNA <- RunUMAP(mes.epi$RNA,dims = 1:20) # arbitrary, for now

pdf(file='mes.epi.RNA.UMAPS.pdf',width=10,height=8)
DimPlot(mes.epi$RNA,group.by = 'Condition',raster=F,shuffle = T)
DimPlot(mes.epi$RNA,group.by = 'Sample',raster=F,shuffle = T)
DimPlot(mes.epi$RNA,group.by = 'Class.Sending',raster=F,shuffle = T)
DimPlot(mes.epi$RNA,group.by = 'Class.Receiving',raster=F,shuffle = T)
DimPlot(mes.epi$RNA,group.by = 'SendingType',raster=F,cols = col.pal$cell_type,shuffle = T)
DimPlot(mes.epi$RNA,group.by = 'ReceivingType',raster=F,cols = col.pal$cell_type,shuffle = T)
dev.off()

# Check out epi to mesenchymal signaling in Imputed slot
mes.epi$Imputed <- FindVariableFeatures(mes.epi$Imputed)
mes.epi$Imputed <- RunPCA(mes.epi$Imputed)
mes.epi$Imputed <- RunUMAP(mes.epi$Imputed,dims = 1:20) # arbitrary, for now

pdf(file='mes.epi.Imputed.UMAPS.pdf',width=10,height=8)
DimPlot(mes.epi$Imputed,group.by = 'Condition',raster=F,shuffle = T)
DimPlot(mes.epi$Imputed,group.by = 'Sample',raster=F,shuffle = T)
DimPlot(mes.epi$Imputed,group.by = 'Class.Sending',raster=F,shuffle = T)
DimPlot(mes.epi$Imputed,group.by = 'Class.Receiving',raster=F,shuffle = T)
DimPlot(mes.epi$Imputed,group.by = 'SendingType',raster=F,cols = col.pal$cell_type,shuffle = T)
DimPlot(mes.epi$Imputed,group.by = 'ReceivingType',raster=F,cols = col.pal$cell_type,shuffle = T)
dev.off()

# Comments 2023-10-17
# Not really loving any of this. The biggest issue is that there is tremendous batch effect between the two timepoints.
# While that *might* be biological, I recognize it pretty clearly as an effect of information depth:
FeaturePlot(mes.epi$Imputed,'nFeature_CellToCell',max.cutoff = 400,split.by = 'Timepoint.Joint')
VlnPlot(mes.epi$Imputed,'nFeature_CellToCell',pt.size = 0)
VlnPlot(mes.epi$Imputed,'nFeature_CellToCell',pt.size = 0,group.by = 'Condition')
VlnPlot(mes.epi$Imputed,'nFeature_CellToCell',pt.size = 0,group.by = 'Timepoint.Joint')

# I think I would like to try integrating by timepoint
split <- SplitObject(mes.epi$Imputed,split.by = 'Timepoint.Joint')
split <- lapply(X = split, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = split)
mes.epi.anchors <- FindIntegrationAnchors(object.list = split, anchor.features = features)
mes.epi.combined <- IntegrateData(anchorset = mes.epi.anchors)

# temporary save in case R crashes
save(mes.epi.combined,file = 'mes.epi.integrated.Robj')

DefaultAssay(mes.epi.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
mes.epi.combined <- ScaleData(mes.epi.combined, verbose = FALSE)
mes.epi.combined <- RunPCA(mes.epi.combined, npcs = 100, verbose = FALSE)

# look at INTEGRATED pcs
pdf(file='mes.epi.integrated.PCs.pdf',width=10,height=8)
ElbowPlot(mes.epi.combined,ndims = 100)
PCHeatmap(mes.epi.combined,cells=200,balanced=T,dims=1:9)
PCHeatmap(mes.epi.combined,cells=200,balanced=T,dims=10:18)
PCHeatmap(mes.epi.combined,cells=200,balanced=T,dims=19:27)
PCHeatmap(mes.epi.combined,cells=200,balanced=T,dims=28:36)
PCHeatmap(mes.epi.combined,cells=200,balanced=T,dims=37:45)
PCHeatmap(mes.epi.combined,cells=200,balanced=T,dims=46:54)
PCHeatmap(mes.epi.combined,cells=200,balanced=T,dims=55:63)
PCHeatmap(mes.epi.combined,cells=200,balanced=T,dims=64:72)
PCHeatmap(mes.epi.combined,cells=200,balanced=T,dims=73:81)
PCHeatmap(mes.epi.combined,cells=200,balanced=T,dims=82:90)
PCHeatmap(mes.epi.combined,cells=200,balanced=T,dims=91:99)
dev.off()

# check out a first embedding in 2D
mes.epi.combined <- RunUMAP(mes.epi.combined, reduction = "pca", dims = 1:20)
DimPlot(mes.epi.combined,group.by = 'Timepoint.Joint',shuffle = T)
DimPlot(mes.epi.combined,group.by = 'SendingType',shuffle = T)
DimPlot(mes.epi.combined,group.by = 'Condition',shuffle = T)

# Looks ok, let's cluster and see
mes.epi.combined <- FindNeighbors(mes.epi.combined, reduction = "pca", dims = 1:20)
mes.epi.combined <- FindClusters(mes.epi.combined, resolution = 0.2)

# plot umaps
pdf(file='mes.epi.imputed.integrated.UMAPS.PCtuned.pdf',width=10,height=8)
DimPlot(mes.epi.combined,group.by = 'Condition',raster=F,shuffle = T)
DimPlot(mes.epi.combined,group.by = 'Sample',raster=F,shuffle = T)
DimPlot(mes.epi.combined,group.by = 'Class.Sending',raster=F,shuffle = T)
DimPlot(mes.epi.combined,group.by = 'Class.Receiving',raster=F,shuffle = T)
DimPlot(mes.epi.combined,group.by = 'SendingType',raster=F,cols = col.pal$cell_type,shuffle = T)
DimPlot(mes.epi.combined,group.by = 'ReceivingType',raster=F,cols = col.pal$cell_type,shuffle = T)
dev.off()


# Not bad, good start.

# What if we are missing something? What if we need to look at how other cell types are talking to mesenchyme? And.or a mesenchymal receptor deficit??
##### Epithelial to mesenchymal signaling ####

# Define colors
col.pal <- list()
col.pal$cell_type <- c('#A40606','#9CFFFA','#B0DB43','#9C528B','#2F6690',
                                '#946846','#F1C40F','green','#0F0326','#E65F5C','#14591D','#726DA8',
                                'yellow','purple')
names(col.pal$cell_type) <- c(unique(epi.mes$RNA$ReceivingType),unique(epi.mes$RNA$SendingType))

# Check out mesenchymal to epithelium signaling in RNA slot
epi.mes$RNA <- FindVariableFeatures(epi.mes$RNA)
epi.mes$RNA <- RunPCA(epi.mes$RNA)
epi.mes$RNA <- RunUMAP(epi.mes$RNA,dims = 1:20) # arbitrary, for now

pdf(file='epi.mes.RNA.UMAPS.pdf',width=10,height=8)
DimPlot(epi.mes$RNA,group.by = 'Condition',raster=F,shuffle = T)
DimPlot(epi.mes$RNA,group.by = 'Sample',raster=F,shuffle = T)
DimPlot(epi.mes$RNA,group.by = 'Class.Sending',raster=F,shuffle = T)
DimPlot(epi.mes$RNA,group.by = 'Class.Receiving',raster=F,shuffle = T)
DimPlot(epi.mes$RNA,group.by = 'SendingType',raster=F,cols = col.pal$cell_type,shuffle = T)
DimPlot(epi.mes$RNA,group.by = 'ReceivingType',raster=F,cols = col.pal$cell_type,shuffle = T)
dev.off()

# Check out epi to mesenchymal signaling in Imputed slot
epi.mes$Imputed <- FindVariableFeatures(epi.mes$Imputed)
epi.mes$Imputed <- RunPCA(epi.mes$Imputed)
epi.mes$Imputed <- RunUMAP(epi.mes$Imputed,dims = 1:20) # arbitrary, for now

pdf(file='epi.mes.Imputed.UMAPS.pdf',width=10,height=8)
DimPlot(epi.mes$Imputed,group.by = 'Condition',raster=F,shuffle = T)
DimPlot(epi.mes$Imputed,group.by = 'Sample',raster=F,shuffle = T)
DimPlot(epi.mes$Imputed,group.by = 'Class.Sending',raster=F,shuffle = T)
DimPlot(epi.mes$Imputed,group.by = 'Class.Receiving',raster=F,shuffle = T)
DimPlot(epi.mes$Imputed,group.by = 'SendingType',raster=F,cols = col.pal$cell_type,shuffle = T)
DimPlot(epi.mes$Imputed,group.by = 'ReceivingType',raster=F,cols = col.pal$cell_type,shuffle = T)
dev.off()

# Comments 2023-10-17
# Not really loving any of this. The biggest issue is that there is tremendous batch effect between the two timepoints.
# While that *might* be biological, I recognize it pretty clearly as an effect of information depth:
FeaturePlot(epi.mes$Imputed,'nFeature_CellToCell',max.cutoff = 400,split.by = 'Timepoint.Joint')
VlnPlot(epi.mes$Imputed,'nFeature_CellToCell',pt.size = 0)
VlnPlot(epi.mes$Imputed,'nFeature_CellToCell',pt.size = 0,group.by = 'Condition')
VlnPlot(epi.mes$Imputed,'nFeature_CellToCell',pt.size = 0,group.by = 'Timepoint.Joint')

# I think I would like to try integrating by timepoint
split <- SplitObject(epi.mes$Imputed,split.by = 'Timepoint.Joint')
split <- lapply(X = split, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = split)
epi.mes.anchors <- FindIntegrationAnchors(object.list = split, anchor.features = features)
epi.mes.combined <- IntegrateData(anchorset = epi.mes.anchors)

# temporary save in case R crashes
save(epi.mes.combined,file = 'epi.mes.integrated.Robj')

DefaultAssay(epi.mes.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
epi.mes.combined <- ScaleData(epi.mes.combined, verbose = FALSE)
epi.mes.combined <- RunPCA(epi.mes.combined, npcs = 100, verbose = FALSE)

# look at INTEGRATED pcs
pdf(file='epi.mes.integrated.PCs.pdf',width=10,height=8)
ElbowPlot(epi.mes.combined,ndims = 100)
PCHeatmap(epi.mes.combined,cells=200,balanced=T,dims=1:9)
PCHeatmap(epi.mes.combined,cells=200,balanced=T,dims=10:18)
PCHeatmap(epi.mes.combined,cells=200,balanced=T,dims=19:27)
PCHeatmap(epi.mes.combined,cells=200,balanced=T,dims=28:36)
PCHeatmap(epi.mes.combined,cells=200,balanced=T,dims=37:45)
PCHeatmap(epi.mes.combined,cells=200,balanced=T,dims=46:54)
PCHeatmap(epi.mes.combined,cells=200,balanced=T,dims=55:63)
PCHeatmap(epi.mes.combined,cells=200,balanced=T,dims=64:72)
PCHeatmap(epi.mes.combined,cells=200,balanced=T,dims=73:81)
PCHeatmap(epi.mes.combined,cells=200,balanced=T,dims=82:90)
PCHeatmap(epi.mes.combined,cells=200,balanced=T,dims=91:99)
dev.off()

# check out a first embedding in 2D
epi.mes.combined <- RunUMAP(epi.mes.combined, reduction = "pca", dims = 1:20)
DimPlot(epi.mes.combined,group.by = 'Timepoint.Joint',shuffle = T)
DimPlot(epi.mes.combined,group.by = 'SendingType',shuffle = T)
DimPlot(epi.mes.combined,group.by = 'Condition',shuffle = T)

# Looks ok, let's cluster and see
epi.mes.combined <- FindNeighbors(epi.mes.combined, reduction = "pca", dims = 1:20)
epi.mes.combined <- FindClusters(epi.mes.combined, resolution = 0.2)

# plot umaps
pdf(file='epi.mes.imputed.integrated.UMAPS.PCtuned.pdf',width=10,height=8)
DimPlot(epi.mes.combined,group.by = 'Condition',raster=F,shuffle = T)
DimPlot(epi.mes.combined,group.by = 'Sample',raster=F,shuffle = T)
DimPlot(epi.mes.combined,group.by = 'Class.Sending',raster=F,shuffle = T)
DimPlot(epi.mes.combined,group.by = 'Class.Receiving',raster=F,shuffle = T)
DimPlot(epi.mes.combined,group.by = 'SendingType',raster=F,cols = col.pal$cell_type,shuffle = T)
DimPlot(epi.mes.combined,group.by = 'ReceivingType',raster=F,cols = col.pal$cell_type,shuffle = T)
dev.off()

# Really, not a lot of signal here...
##### What about the system level measurements? These are designed to take in variable cell population distributions ####
