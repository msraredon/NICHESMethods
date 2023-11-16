# Set wd
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Rachel/E17.E19.Merged.Explorations.2023-10-09")

# Load packages
require(Seurat)
require(ggplot2)
require(cowplot)
require(dplyr)
library(future) # enables multithreadhing # see https://satijalab.org/seurat/archive/v3.0/future_vignette.html

# change the current plan to access parallelization
plan("multisession", workers = 8)  #adjust based on how many cores you have. I have 10 total so i am reserving 2 for other tasks
plan()
options(future.globals.maxSize= 10000*1024^2) # 10 Gb futures export limit, bug fix for below

# load data
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Rachel/E17.E19.Merged.Explorations.2023-10-09/FRL17.19_merged.2023-10-09.Robj")

# inspect data
names(FRL17.19@meta.data)
table(Idents(FRL17.19))
table(Idents(FRL17.19),FRL17.19$Sample)
table(Idents(FRL17.19),FRL17.19$Condition)
table(Idents(FRL17.19),FRL17.19$CellClass)
table(Idents(FRL17.19),FRL17.19$CellType)
table(FRL17.19$CellType,FRL17.19$Sample)
table(FRL17.19$CellClass,FRL17.19$Sample)

# Create a new object with a shorted handle (quicker to type)
frl <- FRL17.19
rm(FRL17.19)
gc() # clears up space

# Clean up metadata
names(frl@meta.data)
# 1. Condition
table(frl$Sample,frl$Condition)
Idents(frl) <- frl$Sample
frl <- RenameIdents(frl,
                    'CDH19' = 'CDH',
                    'NRL19' = 'Normal',
                    'cdhfrl1'='CDH',
                    'cdhfrl2'='CDH',
                    'frl1'='Normal',
                    'frl2'='Normal')
frl$Condition <- Idents(frl)
table(frl$Sample,frl$Condition)
# 2. Class
table(frl$CellClass,frl$Sample) # Looks good already, no change
# 3. CellType
table(frl$CellType,frl$Sample) # This will work for now, we will update later
# 4. Timepoint (new!)
frl$Timepoint <- NA
Idents(frl) <- frl$Timepoint
cell.e17 <- names(frl$Sample[frl$Sample %in% c('cdhfrl1','cdhfrl2','frl1','frl2')])
frl <- SetIdent(frl, cells = cell.e17, value = 'E17')
cell.e19 <- names(frl$Sample[frl$Sample %in% c('CDH19','NRL19')])
frl <- SetIdent(frl, cells = cell.e19, value = 'E19')
table(Idents(frl))
frl$Timepoint <- Idents(frl)

# clean up unnecessary metadata (saves space)
frl$RNA_snn_res.0.2 <- NULL
frl$RNA_snn_res.0.3 <- NULL
frl$RNA_snn_res.0.4 <- NULL
frl$CellType_Preliminary <- NULL
frl$seurat_clusters <- NULL
gc()

# Now, let's see if we can add the multiplexing information and make anything of it
# Not extracted from Fastqs yet, let's revisit with guilin, and re-visit later

# First let's compare data quality
Idents(frl) <- frl$Sample
pdf('data QC1.pdf', width = 7,height = 5)
VlnPlot(frl,'nFeature_RNA') # not nice ##### we need to re-sequence and we need to make sure we (Dave) are not charged twice.
VlnPlot(frl,'nCount_RNA',log=T)
VlnPlot(frl,'percent.mt')
dev.off()

# In the meantime, let's keep on pluggin here
# First let's scale the data and run PCA
frl <- ScaleData(frl)
frl <- FindVariableFeatures(frl)
frl <- RunPCA(frl,npcs = 100)
# look at pcs
pdf(file='frl.PCs.pdf',width=10,height=8)
ElbowPlot(frl,ndims = 100)
PCHeatmap(frl,cells=200,balanced=T,dims=1:9)
PCHeatmap(frl,cells=200,balanced=T,dims=10:18)
PCHeatmap(frl,cells=200,balanced=T,dims=19:27)
PCHeatmap(frl,cells=200,balanced=T,dims=28:36)
PCHeatmap(frl,cells=200,balanced=T,dims=37:45)
PCHeatmap(frl,cells=200,balanced=T,dims=46:54)
PCHeatmap(frl,cells=200,balanced=T,dims=55:63)
PCHeatmap(frl,cells=200,balanced=T,dims=64:72)
PCHeatmap(frl,cells=200,balanced=T,dims=73:81)
PCHeatmap(frl,cells=200,balanced=T,dims=82:90)
PCHeatmap(frl,cells=200,balanced=T,dims=91:99)
dev.off()

# Try a PC selection and see how it looks in 2D space
frl <- RunUMAP(frl,dims = 1:60)
p1 <- DimPlot(frl,group.by = 'Condition',shuffle = T)
p2 <- DimPlot(frl,group.by = 'Sample',shuffle = T)
png('first embedding.png',width = 10,height = 7,units = 'in',res=300)
plot_grid(p1,p2) # BIG TIME BATCH EFFECT DUE TO READ DEPTH DIFFERENCE, NO SURPIRSE BUT ANNOYING
dev.off()

# Given the read-depth artifact, I would integrate at this point and see if we can stil get useful information out of this combined dataset
# See https://satijalab.org/seurat/articles/integration_introduction

frl.list <- SplitObject(frl, split.by = "Sample") # Try integrating by Sample first, could also do E17 vs. E19...

# normalize and identify variable features for each dataset independently
frl.list <- lapply(X = frl.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = frl.list)
gc()
frl.anchors <- FindIntegrationAnchors(object.list = frl.list, anchor.features = features)
gc()
rm(frl)
rm(frl.list)
gc()
# this command creates an 'integrated' data assay
plan('default') # switch to normal non parallel processing to this line to work
frl.combined <- IntegrateData(anchorset = frl.anchors)
gc()
rm(frl.anchors)
# temporary save in case R crashes
save(frl.combined,file = 'frl.combined.Robj')

# Perform an integrated analysis
# Now we can run a single integrated analysis on all cells!
#   
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(frl.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
frl.combined <- ScaleData(frl.combined, verbose = FALSE)
frl.combined <- RunPCA(frl.combined, npcs = 100, verbose = FALSE)

# look at INTEGRATED pcs
pdf(file='frl.PCs.integrated.pdf',width=10,height=8)
ElbowPlot(frl.combined,ndims = 100)
PCHeatmap(frl.combined,cells=200,balanced=T,dims=1:9)
PCHeatmap(frl.combined,cells=200,balanced=T,dims=10:18)
PCHeatmap(frl.combined,cells=200,balanced=T,dims=19:27)
PCHeatmap(frl.combined,cells=200,balanced=T,dims=28:36)
PCHeatmap(frl.combined,cells=200,balanced=T,dims=37:45)
PCHeatmap(frl.combined,cells=200,balanced=T,dims=46:54)
PCHeatmap(frl.combined,cells=200,balanced=T,dims=55:63)
PCHeatmap(frl.combined,cells=200,balanced=T,dims=64:72)
PCHeatmap(frl.combined,cells=200,balanced=T,dims=73:81)
PCHeatmap(frl.combined,cells=200,balanced=T,dims=82:90)
PCHeatmap(frl.combined,cells=200,balanced=T,dims=91:99)
dev.off()

# check out a first embedding in 2D
frl.combined <- RunUMAP(frl.combined, reduction = "pca", dims = 1:60)

# Looks ok, let's cluster and see
frl.combined <- FindNeighbors(frl.combined, reduction = "pca", dims = 1:60)
frl.combined <- FindClusters(frl.combined, resolution = 0.2)

# Visualization
p1 <- DimPlot(frl.combined,group.by = "Sample",shuffle=T)
p2 <- DimPlot(frl.combined,group.by = "Timepoint",shuffle=T)
p3 <- DimPlot(frl.combined,group.by = "Condition",shuffle=T)
p4 <- DimPlot(frl.combined,group.by = "seurat_clusters",shuffle=T)
png('integrated embedding first pass.png',width = 18,height = 10,units = 'in',res=300)
plot_grid(p1,p2,p3,p4,nrow = 2)
dev.off()

# get some markers
DefaultAssay(frl.combined) <- 'RNA' # all differential testing and visualization of gene expression levels must be done on the raw RNA assay and not on the integrated assay which makes up pseudo-values
mark <- FindAllMarkers(frl.combined,min.pct = 0.5,logfc.threshold = 0.5) # adjust threshold to suit user preferences
mark$ratio <- mark$pct.1/mark$pct.2
mark$power <-  mark$ratio*mark$avg_log2FC
View(mark)

# prelim. labels
DimPlot(frl.combined,group.by = "seurat_clusters",shuffle=T,label=T)
frl.combined <- RenameIdents(frl.combined,
                             '0'='Fibroblasts',
                             '1'='Fibroblasts_Proliferating',
                             '2'='Myofibroblasts',
                             '3'='Epithelium',
                             '4'='Endothelium',
                             '5'='Mesenchymal_Progenitors', # This cluster is very interesting and needs to be subclustered to resolve internal patterns
                             '6'='Myeloid',
                             '7'='Pericytes',
                             '8'='RBC',
                             '9'='Mesothelium',
                             '10'='Lymphoid',
                             '11'='NCC',
                             '12'='Cardiac Artifact',
                             '13'='Fibroblasts',
                             '14'='Cardiac Artifact',
                             '15'='Lymphatic')

FeaturePlot(frl.combined,'Cthrc1',label=T) # just to look at markers to determine above
# stash ident labels in a new metadata slot
frl.combined$prelim.celltypes <- Idents(frl.combined)
# add cell class data based on the above
cell.mes <- names(frl.combined$prelim.celltypes[frl.combined$prelim.celltypes %in% c('Fibroblasts','Fibroblasts_Proliferating',
                                                                   'Myofibroblasts','Mesenchymal_Progenitors',
                                                                   'Pericytes','Mesothelium')])
cell.epi <- names(frl.combined$prelim.celltypes[frl.combined$prelim.celltypes %in% c('Epithelium')])
cell.end <- names(frl.combined$prelim.celltypes[frl.combined$prelim.celltypes %in% c('Endothelium','Lymphatic')])
cell.imm <- names(frl.combined$prelim.celltypes[frl.combined$prelim.celltypes %in% c('Myeloid','Lymphoid')])
cell.misc <- names(frl.combined$prelim.celltypes[frl.combined$prelim.celltypes %in% c('NCC','RBC','Cardiac Artifact')])
frl.combined <- SetIdent(frl.combined, cells = cell.mes, value = 'Mesenchyme')
frl.combined <- SetIdent(frl.combined, cells = cell.epi, value = 'Epithelium')
frl.combined <- SetIdent(frl.combined, cells = cell.end, value = 'Endothelium')
frl.combined <- SetIdent(frl.combined, cells = cell.imm, value = 'Immune')
frl.combined <- SetIdent(frl.combined, cells = cell.misc, value = 'Miscellany')
table(Idents(frl.combined))
frl.combined$Class <- Idents(frl.combined)
sum(is.na(frl.combined$Class)) # needs to equal 0
gc()

# plot for posterity
p1 <- DimPlot(frl.combined,group.by = "Sample",shuffle=T)
p2 <- DimPlot(frl.combined,group.by = "Timepoint",shuffle=T)
p3 <- DimPlot(frl.combined,group.by = "Condition",shuffle=T)
p4 <- DimPlot(frl.combined,group.by = "Class",shuffle=T)
p5 <- DimPlot(frl.combined,group.by = "prelim.celltypes",shuffle=T,label=T)
p6 <- FeaturePlot(frl.combined,"Sox9",label=T)

png('integrated UMAPs annotated.png',width = 24,height = 16,units = 'in',res=300)
plot_grid(p1,p2,p3,p4,p5,p6,nrow = 2)
dev.off()

# subset out mesenchyme
mes <- subset(frl.combined,idents = 'Mesenchyme')
DefaultAssay(mes) <- "integrated" # we switch back here to create an embedding

# Run the standard workflow for visualization and clustering
mes <- ScaleData(mes, verbose = FALSE)
mes <- RunPCA(mes, npcs = 100, verbose = FALSE)

# look at INTEGRATED pcs
pdf(file='mes.PCs.integrated.pdf',width=10,height=8)
ElbowPlot(mes,ndims = 100)
PCHeatmap(mes,cells=200,balanced=T,dims=1:9)
PCHeatmap(mes,cells=200,balanced=T,dims=10:18)
PCHeatmap(mes,cells=200,balanced=T,dims=19:27)
PCHeatmap(mes,cells=200,balanced=T,dims=28:36)
PCHeatmap(mes,cells=200,balanced=T,dims=37:45)
PCHeatmap(mes,cells=200,balanced=T,dims=46:54)
PCHeatmap(mes,cells=200,balanced=T,dims=55:63)
PCHeatmap(mes,cells=200,balanced=T,dims=64:72)
PCHeatmap(mes,cells=200,balanced=T,dims=73:81)
PCHeatmap(mes,cells=200,balanced=T,dims=82:90)
PCHeatmap(mes,cells=200,balanced=T,dims=91:99)
dev.off()

# check out a first embedding in 2D
mes <- RunUMAP(mes, reduction = "pca", dims = 1:22)

# Looks ok, let's cluster and see
mes <- FindNeighbors(mes, reduction = "pca", dims = 1:22)
mes <- FindClusters(mes, resolution = 0.2)

# Visualization
p1 <- DimPlot(mes,group.by = "Sample",shuffle=T)
p2 <- DimPlot(mes,group.by = "Timepoint",shuffle=T)
p3 <- DimPlot(mes,group.by = "Condition",shuffle=T)
p4 <- DimPlot(mes,group.by = "seurat_clusters",shuffle=T,label=T)
png('mes.integrated.umaps.first.pass.png',width = 18,height = 10,units = 'in',res=300)
plot_grid(p1,p2,p3,p4,nrow = 2)
dev.off()

# get some markers
DefaultAssay(mes) <- 'RNA' # all differential testing and visualization of gene expression levels must be done on the raw RNA assay and not on the integrated assay which makes up pseudo-values
mark.mes <- FindAllMarkers(mes,min.pct = 0.5,logfc.threshold = 0.5) # adjust threshold to suit user preferences
mark.mes$ratio <- mark.mes$pct.1/mark.mes$pct.2
mark.mes$power <-  mark.mes$ratio*mark.mes$avg_log2FC
View(mark.mes)

# prelim. labels
DimPlot(mes,group.by = "seurat_clusters",shuffle=T,label=T)
mes <- RenameIdents(mes,
                             '0'='Fibroblasts',
                             '1'='Fibroblasts_Proliferating',
                             '2'='Myofibroblasts',
                             '3'='Mesenchymal_Progenitors',
                             '4'='Pericytes',
                             '5'='Fibroblasts_Proliferating',
                              '6'='Mesothelium',
                              '7'='Fat3+_Fibroblasts')

FeaturePlot(mes,'Myh11',label=T) # just to look at markers to determine above

# stash ident labels in a new metadata slot
mes$subcluster.celltypes <- Idents(mes)

# plot for posterity
p1 <- DimPlot(mes,group.by = "Sample",shuffle=T)
p2 <- DimPlot(mes,group.by = "Timepoint",shuffle=T)
p3 <- DimPlot(mes,group.by = "Condition",shuffle=T)
p4 <- DimPlot(mes,group.by = "subcluster.celltypes",shuffle=T,label=T)
p5 <- FeaturePlot(mes,"Sox9",label=T)
p6 <- VlnPlot(mes,'Sox9',split.by = 'Condition')
png('mesenchyme integrated UMAPs annotated.png',width = 24,height = 16,units = 'in',res=300)
plot_grid(p1,p2,p3,p4,p5,p6,nrow = 2)
dev.off()

# subset out JUST THE MESENCHYMAL PROGENITORS
Idents(mes) <- mes$subcluster.celltypes
pro <- subset(mes,idents = 'Mesenchymal_Progenitors')
DefaultAssay(pro) <- "integrated" # we switch back here to create an embedding

# Run the standard workflow for visualization and clustering
pro <- ScaleData(pro, verbose = FALSE)
pro <- RunPCA(pro, npcs = 100, verbose = FALSE)

# look at INTEGRATED pcs
pdf(file='pro.PCs.integrated.pdf',width=10,height=8)
ElbowPlot(pro,ndims = 100)
PCHeatmap(pro,cells=200,balanced=T,dims=1:9)
PCHeatmap(pro,cells=200,balanced=T,dims=10:18)
PCHeatmap(pro,cells=200,balanced=T,dims=19:27)
PCHeatmap(pro,cells=200,balanced=T,dims=28:36)
PCHeatmap(pro,cells=200,balanced=T,dims=37:45)
PCHeatmap(pro,cells=200,balanced=T,dims=46:54)
PCHeatmap(pro,cells=200,balanced=T,dims=55:63)
PCHeatmap(pro,cells=200,balanced=T,dims=64:72)
PCHeatmap(pro,cells=200,balanced=T,dims=73:81)
PCHeatmap(pro,cells=200,balanced=T,dims=82:90)
PCHeatmap(pro,cells=200,balanced=T,dims=91:99)
dev.off()

# check out a first embedding in 2D
pro <- RunUMAP(pro, reduction = "pca", dims = 1:20)

# Looks ok, let's cluster and see
pro <- FindNeighbors(pro, reduction = "pca", dims = 1:20)
pro <- FindClusters(pro, resolution = 0.2)

# Visualization
p1 <- DimPlot(pro,group.by = "Sample",shuffle=T)
p2 <- DimPlot(pro,group.by = "Timepoint",shuffle=T)
p3 <- DimPlot(pro,group.by = "Condition",shuffle=T)
p4 <- DimPlot(pro,group.by = "seurat_clusters",shuffle=T,label=T)
png('pro.integrated.umaps.first.pass.png',width = 14,height = 7,units = 'in',res=300)
plot_grid(p1,p2,p3,p4,nrow = 2)
dev.off()

# get some markers
DefaultAssay(pro) <- 'RNA' # all differential testing and visualization of gene expression levels must be done on the raw RNA assay and not on the integrated assay which makes up pseudo-values
mark.pro <- FindAllMarkers(pro,min.pct = 0.25,logfc.threshold = 0.25) # adjust threshold to suit user preferences
mark.pro$ratio <- mark.pro$pct.1/mark.pro$pct.2
mark.pro$power <-  mark.pro$ratio*mark.pro$avg_log2FC
View(mark.pro)

# prelim. labels
DimPlot(pro,group.by = "seurat_clusters",shuffle=T,label=T)
pro <- RenameIdents(pro,
                    '0'='Fibroblasts',
                    '1'='Fibroblasts',
                    '2'='Sox9+',
                    '3'='Dcn+',
                    '4'='Cardiac Artifact')

FeaturePlot(pro,'Sox9',label=T) # just to look at markers to determine above

# stash ident labels in a new metadata slot
pro$subcluster.celltypes <- Idents(pro) # oops this overwrites, oh well


# plot cell numbers per each cluster
cell.frac <- table(Idents(pro),pro$Condition)
sums <- colSums(cell.frac)
cell.frac[,1] <- cell.frac[,1]/sums[1]
cell.frac[,2] <- cell.frac[,2]/sums[2]
colSums(cell.frac) # should equal 1
cell.frac <- data.frame(cell.frac)

# Demo for Hahram
distributions2 <- table(Idents(pro),pro@meta.data$Condition)
for(i in 1:nrow(distributions2)){
  distributions2[i,] <- distributions2[i,]/rowSums(distributions2)[i]
}
distributions2 <- data.frame(distributions2)


# plot for posterity
p1 <- DimPlot(pro,group.by = "Sample",shuffle=T)
p2 <- DimPlot(pro,group.by = "Timepoint",shuffle=T)
p3 <- DimPlot(pro,group.by = "Condition",shuffle=T)
p4 <- DimPlot(pro,group.by = "subcluster.celltypes",shuffle=T,label=T)
p5 <- ggplot(cell.frac,aes(x=Var1,y=Freq,group=Var2,fill=Var2))+
              geom_bar(stat = 'identity',position = 'dodge')+
  xlab('SubCluster')+
  ylab('Fraction of SubClustered Dataset')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme_classic()+ggtitle('Distribution across subclusters')
p6 <- ggplot(distributions2,aes(x=Var1,y=Freq,group=Var2,fill=Var2))+
  geom_bar(stat = 'identity',position = 'dodge')+
  xlab('SubCluster')+
  ylab('Fraction of Each Cluster')+
  ggtitle('Condition Contribution to Each Archetype')+theme_classic()
png('mes progenitors integrated UMAPs annotated.png',width = 16,height = 7,units = 'in',res=300)
plot_grid(p1,p2,p3,p4,p5,p6,nrow = 2)
dev.off()


# save progenitor object for later
save(pro,file = 'pro.combined.annotated.2023-10-09.Robj')

# save mesenchyme object for later
save(mes,file = 'mes.combined.annotated.2023-10-09.Robj')

# save combined object for later
save(frl.combined,file = 'frl.combined.annotated.classed.2023-10-09.Robj')

