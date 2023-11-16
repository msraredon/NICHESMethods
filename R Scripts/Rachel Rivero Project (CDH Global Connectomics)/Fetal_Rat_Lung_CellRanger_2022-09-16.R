# Set WD
setwd("/Volumes/Samsung_T5/Fetal_Rat_Lung")

# Load required packages
require(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(NICHES)

# See this vignette for background on this script
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

# Load the data
frl1 <- Read10X(data.dir = "/Volumes/Samsung_T5/Fetal_Rat_Lung/FRL1/raw_feature_bc_matrix")
frl2 <- Read10X(data.dir = "/Volumes/Samsung_T5/Fetal_Rat_Lung/FRL2/raw_feature_bc_matrix")

# Initialize the Seurat objects
frl1.seurat <- CreateSeuratObject(counts = frl1, project = "FetalRatLung", min.cells = 3, min.features = 200)
frl2.seurat <- CreateSeuratObject(counts = frl2, project = "FetalRatLung", min.cells = 3, min.features = 200)

# Tag with sample identity
frl1.seurat$Sample <- 'FRL1'
frl2.seurat$Sample <- 'FRL2'

# Merge into single object
frl <- merge(frl1.seurat,frl2.seurat)

# Remove previous objects (no longer needed) and clear garbage ('gc()') to free up space
rm(frl1)
rm(frl2)
rm(frl1.seurat)
rm(frl2.seurat)
gc()

# Calculate fraction per cell mitochondrial reads
frl[["percent.mt"]] <- PercentageFeatureSet(frl, pattern = "^Mt-")

# First QC metrics
png(filename = 'FRL_QC_1_Unfiltered.png',width = 15,height = 10,res = 200,units = 'in')
VlnPlot(frl,features = c('nCount_RNA','nFeature_RNA','percent.mt'),pt.size = 0.1,log=T,split.by = 'Sample')
dev.off()

# First look temporary filtration
frl <- subset(frl, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 25)

# Normalize data (turn counts into "fraction of transcriptome per cell")
frl <- NormalizeData(frl, normalization.method = "LogNormalize", scale.factor = 10000)

# Scale Data (scaling each ROW) (can do regression steps here to minimize artifacts if they exist)
frl <- ScaleData(frl)

# Find Variable Features
frl <- FindVariableFeatures(frl, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes (unnecessary but for show)
top10 <- head(VariableFeatures(frl), 50)

# Plot variable features with and without labels to get a sense of the data structure
plot1 <- VariableFeaturePlot(frl)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Principle Component Analysis
frl <- RunPCA(frl, features = VariableFeatures(object = frl),npcs = 100)

# Visualize all PCs
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

# Cluster
frl <- FindNeighbors(frl, dims = 1:89)
frl <- FindClusters(frl, resolution = 0.2)

# Embed
frl <- RunUMAP(frl, dims = 1:89)

# Plot for first look
pdf(file='FRL_UMAP_1.pdf')
DimPlot(frl, reduction = "umap")+ggtitle('Fetal Rat Lung')
dev.off()

# Plot more comprehensively
png(file='FRL_UMAP_2.png',width=18,height=6,units = 'in',res=300)
p1 <- DimPlot(frl, reduction = "umap", group.by = "Sample",shuffle = T)
p2 <- DimPlot(frl, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- FeaturePlot(frl,'nFeature_RNA')
p1 + p2 + p3
dev.off()

# Find Markers
mark <- FindAllMarkers(frl,only.pos = T,min.pct = 0.5,logfc.threshold = 0.5)
mark$ratio <- mark$pct.1/mark$pct.2

# Identify cell types
frl <- RenameIdents(frl,
                    '0'='Fibroblasts',
                    '1'='Epithelium',
                    '2'='Ribosomal?',
                    '3'='Endothelial',
                    '4'='Myofibroblasts/SMCs',
                    '5'='Pericytes!',
                    '6'='Rhd_RBC',
                    '7'='Epithelium',
                    '8'='Immune',
                    '9'='Sox9_Mesenchyme',
                    '10'='Akr7a3_RBC',
                    '11'='Rspo1_Mesenchyme',
                    '12'='Sox10_Mesenchyme')
frl$CellTypes_Preliminary <- Idents(frl)


# Plot
png(file='FRL_CellType_Preliminary.png',width=18,height=8,units = 'in',res=300)
p1 <- DimPlot(frl, reduction = "umap", group.by = "Sample",shuffle = T,repel = T)
p2 <- DimPlot(frl,label = T)
p1 + p2
dev.off()

# Save objects
save(frl,file = 'frl.2022-09-16.Robj')
save(mark,file='frl.marker.genes.2022-09-16.Robj')
