# Merging E17 and E19
# Oct 6, 2023

setwd("/Volumes/Rachel Rivero/E19 NL CDH")
set.seed(123)

# Download required packages
require(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(NICHES)

load(file='FRL19.Seurat.unlabeled.2023-09-27.Robj')
DimPlot(FRL19.Seurat, reduction = 'umap')

# Define cell type 
FRL19.Seurat <- RenameIdents(FRL19.Seurat,
                             '0'='Proliferating Mesenchyme',
                             '1'='Proliferating Mesenchyme',
                             '2'='Myofibroblasts/SMCs',
                             '3'='Distal Epithelium',
                             '4'='Endothelial',
                             '5'='Fibroblasts?',
                             '6'='Pericytes',
                             '7'='Macrophages',
                             '8'='Dendritic/APCs',
                             '9'='Prox Epithelium',
                             '10'='B and T cells',
                             '11'='Foxp2+ Epithelium',
                             '12'='RBCs',
                             '13'='Mesothelial',
                             '14'='Cardiac muscle cells',
                             '15'='Cardiac muscle cells',
                             '16'='NCCs',
                             '17'='Lymphatic Endothelial Cells',
                             '18'='Neutrophils',
                             '19'='NECs'
)

# stash cell type
FRL19.Seurat$CellType <- Idents(FRL19.Seurat)
head(FRL19.Seurat) # confirm

# Cell class 
FRL19.Seurat <- RenameIdents(FRL19.Seurat,
                             'Proliferating Mesenchyme'='Mesenchymal',
                             'Proliferating Mesenchyme'='Mesenchymal',
                             'Myofibroblasts/SMCs'='Mesenchymal',
                             'Distal Epithelium'='Epithelial',
                             'Endothelial'='Endothelial',
                             'Fibroblasts?'='Mesenchymal',
                             'Pericytes'='Mesenchymal',
                             'Macrophages'='Immune',
                             'Dendritic/APCs'='Immune',
                             'Prox Epithelium'='Epithelial',
                             'B and T cells'='Immune',
                             'Foxp2+ Epithelium'='Epithelial',
                             'RBCs'='RBCs',
                             'Mesothelial'='Mesenchymal',
                             'Cardiac muscle cells'='Mesenchymal',
                             'Cardiac muscle cells'='Mesenchymal',
                             'NCCs'='Mesenchymal',
                             'Lymphatic Endothelial Cells'='Endothelial',
                             'Neutrophils'='Immune',
                             'NECs'='Mesenchymal')
FRL19.Seurat$CellClass <- Idents(FRL19.Seurat)
head(FRL19.Seurat)

DimPlot(FRL19.Seurat, reduction = 'umap')

# Save object with cell type and cell class
save(FRL19.Seurat, file = 'FRL19.CellClass.2023-09-29.Robj')

# Same as above for E17
# Load unlabeled object
load("/Volumes/Rachel Rivero/E17 NL_CDH/CDH_FRL.merge.Robj")
DimPlot(merge, reduction = 'umap')

# Cell type
merge <- RenameIdents(merge,
                      '0'='Mes_Stem/Fibs',
                      '1'='Prolif_Mes',
                      '2'='Cardiac/SMCs',
                      '3'='Epithelial',
                      '4'='Pericytes',
                      '5'='Immune',
                      '6'='Mes_Prog',
                      '7'='Endothelial',
                      '8'='RBCs',
                      '9'='Mesothelial',
                      '10'='NCC',
                      '11'='Lymphatic_Endo')
merge$CellType <- Idents(merge)
head(merge)

# Cell class
merge <- RenameIdents(merge,
                      'Mes_Stem/Fibs'='Mesenchymal',
                      'Prolif_Mes'='Mesenchymal',
                      'Cardiac/SMCs'='Mesenchymal',
                      'Epithelial'='Epithelial',
                      'Pericytes'='Mesenchymal',
                      'Immune'='Immune',
                      'Mes_Prog'='Mesenchymal',
                      'Endothelial'='Endothelial',
                      'RBCs'='RBCs',
                      'Mesothelial'='Mesenchymal',
                      'NCC'='Mesenchymal',
                      'Lymphatic_Endo'='Endothelial')
merge$CellClass <- Idents(merge)
head(merge)

# Merge labeled E17 and E19
FRL17.19 <- merge(merge, FRL19.Seurat)
table(FRL17.19$Sample)
table(FRL17.19$CellClass)


