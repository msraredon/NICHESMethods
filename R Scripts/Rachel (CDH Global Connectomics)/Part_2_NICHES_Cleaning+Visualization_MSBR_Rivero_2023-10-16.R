# PART 2 This script cleans the NICHES data and performs initial global visualizations of Dr. Rivero's data #

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

# Load Data
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Rachel/E17.E19.Merged.Explorations.2023-10-09/system.to.cell.Robj")
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Rachel/E17.E19.Merged.Explorations.2023-10-09/cell.to.system.Robj")
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Rachel/E17.E19.Merged.Explorations.2023-10-09/cell.to.cell.Robj")
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Rachel/E17.E19.Merged.Explorations.2023-10-09/system.to.cell.imputed.Robj")
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Rachel/E17.E19.Merged.Explorations.2023-10-09/cell.to.system.imputed.Robj")
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Rachel/E17.E19.Merged.Explorations.2023-10-09/cell.to.cell.imputed.Robj")

# Clean each data set
# CTC
p1 <- VlnPlot(cell.to.cell,group.by='Sample',c('nFeature_CellToCell'),raster=F,pt.size = 0.1)+ggtitle('Not Imputed')+ylab('nFeature_CellToCell')
p2 <- VlnPlot(cell.to.cell.imputed,group.by='Sample',c('nFeature_CellToCell'),raster=F,pt.size = 0.1)+ggtitle('Imputed')+ylab('nFeature_CellToCell')
png(filename = 'CellToCell_BeforeFiltration.png',width = 12,height = 7,units = 'in',res=300)
plot_grid(p1,p2,nrow=1)
dev.off()
# STC
p1 <- VlnPlot(system.to.cell,group.by='Sample',c('nFeature_SystemToCell'),raster=F,pt.size = 0.1)+ggtitle('Not Imputed')+ylab('nFeature_SystemToCell')
p2 <- VlnPlot(system.to.cell.imputed,group.by='Sample',c('nFeature_SystemToCell'),raster=F,pt.size = 0.1)+ggtitle('Imputed')+ylab('nFeature_SystemToCell')
png(filename = 'SystemToCell_BeforeFiltration.png',width = 12,height = 7,units = 'in',res=300)
plot_grid(p1,p2,nrow=1)
dev.off()
# CTS
p1 <- VlnPlot(cell.to.system,group.by='Sample',c('nFeature_CellToSystem'),raster=F,pt.size = 0.1)+ggtitle('Not Imputed')+ylab('nFeature_CellToSystem')
p2 <- VlnPlot(cell.to.system.imputed,group.by='Sample',c('nFeature_CellToSystem'),raster=F,pt.size = 0.1)+ggtitle('Imputed')+ylab('nFeature_CellToSystem')
png(filename = 'CellToSystem_BeforeFiltration.png',width = 12,height = 7,units = 'in',res=300)
plot_grid(p1,p2,nrow=1)
dev.off()

# Apply reasonable thresholds as a starting point
cell.to.cell <- subset(cell.to.cell,nFeature_CellToCell>50)
cell.to.cell.imputed <- subset(cell.to.cell.imputed,nFeature_CellToCell>100)
system.to.cell <- subset(system.to.cell,nFeature_SystemToCell>100)
system.to.cell.imputed <- subset(system.to.cell.imputed,nFeature_SystemToCell>100)
cell.to.system <- subset(cell.to.system,nFeature_CellToSystem>100)
cell.to.system.imputed <- subset(cell.to.system.imputed,nFeature_CellToSystem>100)

# CTC
p1 <- VlnPlot(cell.to.cell,group.by='Sample',c('nFeature_CellToCell'),raster=F,pt.size = 0.1)+ggtitle('Not Imputed')+ylab('nFeature_CellToCell')
p2 <- VlnPlot(cell.to.cell.imputed,group.by='Sample',c('nFeature_CellToCell'),raster=F,pt.size = 0.1)+ggtitle('Imputed')+ylab('nFeature_CellToCell')
png(filename = 'CellToCell_AfterFiltration.png',width = 12,height = 7,units = 'in',res=300)
plot_grid(p1,p2,nrow=1)
dev.off()
# STC
p1 <- VlnPlot(system.to.cell,group.by='Sample',c('nFeature_SystemToCell'),raster=F,pt.size = 0.1)+ggtitle('Not Imputed')+ylab('nFeature_SystemToCell')
p2 <- VlnPlot(system.to.cell.imputed,group.by='Sample',c('nFeature_SystemToCell'),raster=F,pt.size = 0.1)+ggtitle('Imputed')+ylab('nFeature_SystemToCell')
png(filename = 'SystemToCell_AfterFiltration.png',width = 12,height = 7,units = 'in',res=300)
plot_grid(p1,p2,nrow=1)
dev.off()
# CTS
p1 <- VlnPlot(cell.to.system,group.by='Sample',c('nFeature_CellToSystem'),raster=F,pt.size = 0.1)+ggtitle('Not Imputed')+ylab('nFeature_CellToSystem')
p2 <- VlnPlot(cell.to.system.imputed,group.by='Sample',c('nFeature_CellToSystem'),raster=F,pt.size = 0.1)+ggtitle('Imputed')+ylab('nFeature_CellToSystem')
png(filename = 'CellToSystem_AfterFiltration.png',width = 12,height = 7,units = 'in',res=300)
plot_grid(p1,p2,nrow=1)
dev.off()

# Visualize each data set

# Scale
cell.to.cell <- ScaleData(cell.to.cell)
cell.to.cell.imputed <- ScaleData(cell.to.cell.imputed)
system.to.cell <- ScaleData(system.to.cell)
system.to.cell.imputed <- ScaleData(system.to.cell.imputed)
cell.to.system <- ScaleData(cell.to.system)
cell.to.system.imputed <- ScaleData(cell.to.system.imputed)

# Find Variable Features
cell.to.cell <- FindVariableFeatures(cell.to.cell)
cell.to.cell.imputed <- FindVariableFeatures(cell.to.cell.imputed)
system.to.cell <- FindVariableFeatures(system.to.cell)
system.to.cell.imputed <- FindVariableFeatures(system.to.cell.imputed)
cell.to.system <- FindVariableFeatures(cell.to.system)
cell.to.system.imputed <- FindVariableFeatures(cell.to.system.imputed)

# Run PCA
cell.to.cell <- RunPCA(cell.to.cell,npcs = 100)
cell.to.cell.imputed <- RunPCA(cell.to.cell.imputed,npcs = 100)
system.to.cell <- RunPCA(system.to.cell,npcs = 100)
system.to.cell.imputed <- RunPCA(system.to.cell.imputed,npcs = 100)
cell.to.system <- RunPCA(cell.to.system,npcs = 100)
cell.to.system.imputed <- RunPCA(cell.to.system.imputed,npcs = 100)

# Check out PCs
pdf(file='cell.to.cell.PCs.pdf',width=10,height=8)
ElbowPlot(cell.to.cell,ndims = 100)
PCHeatmap(cell.to.cell,cells=200,balanced=T,dims=1:9)
PCHeatmap(cell.to.cell,cells=200,balanced=T,dims=10:18)
PCHeatmap(cell.to.cell,cells=200,balanced=T,dims=19:27)
PCHeatmap(cell.to.cell,cells=200,balanced=T,dims=28:36)
PCHeatmap(cell.to.cell,cells=200,balanced=T,dims=37:45)
PCHeatmap(cell.to.cell,cells=200,balanced=T,dims=46:54)
PCHeatmap(cell.to.cell,cells=200,balanced=T,dims=55:63)
PCHeatmap(cell.to.cell,cells=200,balanced=T,dims=64:72)
PCHeatmap(cell.to.cell,cells=200,balanced=T,dims=73:81)
PCHeatmap(cell.to.cell,cells=200,balanced=T,dims=82:90)
PCHeatmap(cell.to.cell,cells=200,balanced=T,dims=91:99)
dev.off()
pdf(file='cell.to.cell.imputed.PCs.pdf',width=10,height=8)
ElbowPlot(cell.to.cell.imputed,ndims = 100)
PCHeatmap(cell.to.cell.imputed,cells=200,balanced=T,dims=1:9)
PCHeatmap(cell.to.cell.imputed,cells=200,balanced=T,dims=10:18)
PCHeatmap(cell.to.cell.imputed,cells=200,balanced=T,dims=19:27)
PCHeatmap(cell.to.cell.imputed,cells=200,balanced=T,dims=28:36)
PCHeatmap(cell.to.cell.imputed,cells=200,balanced=T,dims=37:45)
PCHeatmap(cell.to.cell.imputed,cells=200,balanced=T,dims=46:54)
PCHeatmap(cell.to.cell.imputed,cells=200,balanced=T,dims=55:63)
PCHeatmap(cell.to.cell.imputed,cells=200,balanced=T,dims=64:72)
PCHeatmap(cell.to.cell.imputed,cells=200,balanced=T,dims=73:81)
PCHeatmap(cell.to.cell.imputed,cells=200,balanced=T,dims=82:90)
PCHeatmap(cell.to.cell.imputed,cells=200,balanced=T,dims=91:99)
dev.off()
pdf(file='system.to.cell.PCs.pdf',width=10,height=8)
ElbowPlot(system.to.cell,ndims = 100)
PCHeatmap(system.to.cell,cells=200,balanced=T,dims=1:9)
PCHeatmap(system.to.cell,cells=200,balanced=T,dims=10:18)
PCHeatmap(system.to.cell,cells=200,balanced=T,dims=19:27)
PCHeatmap(system.to.cell,cells=200,balanced=T,dims=28:36)
PCHeatmap(system.to.cell,cells=200,balanced=T,dims=37:45)
PCHeatmap(system.to.cell,cells=200,balanced=T,dims=46:54)
PCHeatmap(system.to.cell,cells=200,balanced=T,dims=55:63)
PCHeatmap(system.to.cell,cells=200,balanced=T,dims=64:72)
PCHeatmap(system.to.cell,cells=200,balanced=T,dims=73:81)
PCHeatmap(system.to.cell,cells=200,balanced=T,dims=82:90)
PCHeatmap(system.to.cell,cells=200,balanced=T,dims=91:99)
dev.off()
pdf(file='system.to.cell.imputed.PCs.pdf',width=10,height=8)
ElbowPlot(system.to.cell.imputed,ndims = 100)
PCHeatmap(system.to.cell.imputed,cells=200,balanced=T,dims=1:9)
PCHeatmap(system.to.cell.imputed,cells=200,balanced=T,dims=10:18)
PCHeatmap(system.to.cell.imputed,cells=200,balanced=T,dims=19:27)
PCHeatmap(system.to.cell.imputed,cells=200,balanced=T,dims=28:36)
PCHeatmap(system.to.cell.imputed,cells=200,balanced=T,dims=37:45)
PCHeatmap(system.to.cell.imputed,cells=200,balanced=T,dims=46:54)
PCHeatmap(system.to.cell.imputed,cells=200,balanced=T,dims=55:63)
PCHeatmap(system.to.cell.imputed,cells=200,balanced=T,dims=64:72)
PCHeatmap(system.to.cell.imputed,cells=200,balanced=T,dims=73:81)
PCHeatmap(system.to.cell.imputed,cells=200,balanced=T,dims=82:90)
PCHeatmap(system.to.cell.imputed,cells=200,balanced=T,dims=91:99)
dev.off()
pdf(file='cell.to.system.PCs.pdf',width=10,height=8)
ElbowPlot(cell.to.system,ndims = 100)
PCHeatmap(cell.to.system,cells=200,balanced=T,dims=1:9)
PCHeatmap(cell.to.system,cells=200,balanced=T,dims=10:18)
PCHeatmap(cell.to.system,cells=200,balanced=T,dims=19:27)
PCHeatmap(cell.to.system,cells=200,balanced=T,dims=28:36)
PCHeatmap(cell.to.system,cells=200,balanced=T,dims=37:45)
PCHeatmap(cell.to.system,cells=200,balanced=T,dims=46:54)
PCHeatmap(cell.to.system,cells=200,balanced=T,dims=55:63)
PCHeatmap(cell.to.system,cells=200,balanced=T,dims=64:72)
PCHeatmap(cell.to.system,cells=200,balanced=T,dims=73:81)
PCHeatmap(cell.to.system,cells=200,balanced=T,dims=82:90)
PCHeatmap(cell.to.system,cells=200,balanced=T,dims=91:99)
dev.off()
pdf(file='cell.to.system.imputed.PCs.pdf',width=10,height=8)
ElbowPlot(cell.to.system.imputed,ndims = 100)
PCHeatmap(cell.to.system.imputed,cells=200,balanced=T,dims=1:9)
PCHeatmap(cell.to.system.imputed,cells=200,balanced=T,dims=10:18)
PCHeatmap(cell.to.system.imputed,cells=200,balanced=T,dims=19:27)
PCHeatmap(cell.to.system.imputed,cells=200,balanced=T,dims=28:36)
PCHeatmap(cell.to.system.imputed,cells=200,balanced=T,dims=37:45)
PCHeatmap(cell.to.system.imputed,cells=200,balanced=T,dims=46:54)
PCHeatmap(cell.to.system.imputed,cells=200,balanced=T,dims=55:63)
PCHeatmap(cell.to.system.imputed,cells=200,balanced=T,dims=64:72)
PCHeatmap(cell.to.system.imputed,cells=200,balanced=T,dims=73:81)
PCHeatmap(cell.to.system.imputed,cells=200,balanced=T,dims=82:90)
PCHeatmap(cell.to.system.imputed,cells=200,balanced=T,dims=91:99)
dev.off()

# Run UMAP
cell.to.cell <- RunUMAP(cell.to.cell,dims = 1:75)
cell.to.cell.imputed <- RunUMAP(cell.to.cell.imputed,dims = 1:75)
system.to.cell <- RunUMAP(system.to.cell,dims = 1:75)
system.to.cell.imputed <- RunUMAP(system.to.cell.imputed,dims = 1:75)
cell.to.system <- RunUMAP(cell.to.system,dims = 1:75)
cell.to.system.imputed <- RunUMAP(cell.to.system.imputed,dims = 1:75)

# Define colors
col.pal <- list()
col.pal$cell_type <- c('#A40606','#9CFFFA','#B0DB43','#9C528B','#2F6690',
                                '#946846','#F1C40F','green','#0F0326','#E65F5C','#14591D','#726DA8',
                                'yellow','purple')
names(col.pal$cell_type) <- unique(system.to.cell.imputed$ReceivingType)

# Plot UMAPS with relevant metadata

# CTC
pdf(file='cell.to.cell.UMAPS.pdf',width=10,height=8)
DimPlot(cell.to.cell,group.by = 'Condition',raster=F,shuffle = T)
DimPlot(cell.to.cell,group.by = 'Sample',raster=F,shuffle = T)
DimPlot(cell.to.cell,group.by = 'Class.Sending',raster=F,shuffle = T)
DimPlot(cell.to.cell,group.by = 'Class.Receiving',raster=F,shuffle = T)
DimPlot(cell.to.cell,group.by = 'SendingType',raster=F,cols = col.pal$cell_type,shuffle = T)
DimPlot(cell.to.cell,group.by = 'ReceivingType',raster=F,cols = col.pal$cell_type,shuffle = T)
dev.off()

pdf(file='cell.to.cell.imputed.UMAPS.pdf',width=10,height=8)
DimPlot(cell.to.cell.imputed,group.by = 'Condition',raster=F,shuffle = T)
DimPlot(cell.to.cell.imputed,group.by = 'Sample',raster=F,shuffle = T)
DimPlot(cell.to.cell.imputed,group.by = 'Class.Sending',raster=F,shuffle = T)
DimPlot(cell.to.cell.imputed,group.by = 'Class.Receiving',raster=F,shuffle = T)
DimPlot(cell.to.cell.imputed,group.by = 'SendingType',raster=F,cols = col.pal$cell_type,shuffle = T)
DimPlot(cell.to.cell.imputed,group.by = 'ReceivingType',raster=F,cols = col.pal$cell_type,shuffle = T)
dev.off()

# STC
pdf(file='system.to.cell.UMAPS.pdf',width=10,height=8)
DimPlot(system.to.cell,group.by = 'Condition',raster=F,shuffle = T)
DimPlot(system.to.cell,group.by = 'Sample',raster=F,shuffle = T)
DimPlot(system.to.cell,group.by = 'Class',raster=F,shuffle = T)
DimPlot(system.to.cell,group.by = 'ReceivingType',raster=F,cols = col.pal$cell_type,shuffle = T)
dev.off()
pdf(file='system.to.cell.imputed.UMAPS.pdf',width=10,height=8)
DimPlot(system.to.cell.imputed,group.by = 'Condition',raster=F,shuffle = T)
DimPlot(system.to.cell.imputed,group.by = 'Sample',raster=F,shuffle = T)
DimPlot(system.to.cell.imputed,group.by = 'Class',raster=F,shuffle = T)
DimPlot(system.to.cell.imputed,group.by = 'ReceivingType',raster=F,cols = col.pal$cell_type,shuffle = T)
dev.off()

# CTS
pdf(file='cell.to.system.UMAPS.pdf',width=10,height=8)
DimPlot(cell.to.system,group.by = 'Condition',raster=F,shuffle = T)
DimPlot(cell.to.system,group.by = 'Sample',raster=F,shuffle = T)
DimPlot(cell.to.system,group.by = 'Class',raster=F,shuffle = T)
DimPlot(cell.to.system,group.by = 'SendingType',raster=F,cols = col.pal$cell_type,shuffle = T)
dev.off()
pdf(file='cell.to.system.imputed.UMAPS.pdf',width=10,height=8)
DimPlot(cell.to.system.imputed,group.by = 'Condition',raster=F,shuffle = T)
DimPlot(cell.to.system.imputed,group.by = 'Sample',raster=F,shuffle = T)
DimPlot(cell.to.system.imputed,group.by = 'Class',raster=F,shuffle = T)
DimPlot(cell.to.system.imputed,group.by = 'SendingType',raster=F,cols = col.pal$cell_type,shuffle = T)
dev.off()

# Save in a convenient format for later analysis
# Initialize
epi.epi <- list(cell.to.cell,cell.to.cell.imputed)
epi.mes <- list(cell.to.cell,cell.to.cell.imputed)
mes.epi <- list(cell.to.cell,cell.to.cell.imputed)
mes.mes <- list(cell.to.cell,cell.to.cell.imputed)
mes.niche <- list(system.to.cell,system.to.cell.imputed)
epi.niche <- list(system.to.cell,system.to.cell.imputed)
mes.influence <- list(cell.to.system,cell.to.system.imputed)
epi.influence <- list(cell.to.system,cell.to.system.imputed)

# Name
names(epi.epi) <- c('RNA','Imputed')
names(epi.mes) <- c('RNA','Imputed')
names(mes.epi) <- c('RNA','Imputed')
names(mes.mes) <- c('RNA','Imputed')
names(mes.niche) <- c('RNA','Imputed')
names(epi.niche) <- c('RNA','Imputed')
names(mes.influence) <- c('RNA','Imputed')
names(mes.influence) <- c('RNA','Imputed')

# Subset
for (i in 1:length(epi.epi)){
  epi.epi[[i]] <- subset(epi.epi[[i]],subset = Class.Sending == 'Epithelium' & Class.Receiving == 'Epithelium')
  epi.epi[[i]] <- ScaleData(epi.epi[[i]])
}
for (i in 1:length(epi.mes)){
  epi.mes[[i]] <- subset(epi.mes[[i]],subset = Class.Sending == 'Epithelium' & Class.Receiving == 'Mesenchyme')
  epi.mes[[i]] <- ScaleData(epi.mes[[i]])
}
for (i in 1:length(mes.mes)){
  mes.mes[[i]] <- subset(mes.mes[[i]],subset = Class.Sending == 'Mesenchyme' & Class.Receiving == 'Mesenchyme')
  mes.mes[[i]] <- ScaleData(mes.mes[[i]])
}
for (i in 1:length(mes.epi)){
  mes.epi[[i]] <- subset(mes.epi[[i]],subset = Class.Sending == 'Mesenchyme' & Class.Receiving == 'Epithelium')
  mes.epi[[i]] <- ScaleData(mes.epi[[i]])
}
for (i in 1:length(mes.niche)){
  mes.niche[[i]] <- subset(mes.niche[[i]],subset = Class == 'Mesenchyme')
  mes.niche[[i]] <- ScaleData(mes.niche[[i]])
}
for (i in 1:length(epi.niche)){
  epi.niche[[i]] <- subset(epi.niche[[i]],subset = Class == 'Epithelium')
  epi.niche[[i]] <- ScaleData(epi.niche[[i]])
}
for (i in 1:length(mes.influence)){
  mes.influence[[i]] <- subset(mes.influence[[i]],subset = Class == 'Mesenchyme')
  mes.influence[[i]] <- ScaleData(mes.influence[[i]])
}
for (i in 1:length(epi.influence)){
  epi.influence[[i]] <- subset(epi.influence[[i]],subset = Class == 'Epithelium')
  epi.influence[[i]] <- ScaleData(epi.influence[[i]])
}

# Save for later
save(epi.epi,file = 'epi.epi.Robj')
save(epi.mes,file = 'epi.mes.Robj')
save(mes.epi,file = 'mes.epi.Robj')
save(mes.mes,file = 'mes.mes.Robj')

save(epi.niche,file = 'epi.niche.Robj')
save(mes.niche,file = 'mes.niche.Robj')

save(epi.influence,file = 'epi.influence.Robj')
save(mes.influence,file = 'mes.influence.Robj')
