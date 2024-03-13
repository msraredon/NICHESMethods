# Set WD
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Collaborations/Fan Lab")

# Load packages
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(dplyr)
require(cowplot)
require(NICHES)
require(ggrepel)

# Load data
gbm1 <- readRDS("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Collaborations/Fan Lab/R_GBM1027N_res0.6.rds")

# Load Cell2Loc data
ad1 <- anndata::read_h5ad("Cell2Loc/GBM1027N/GBM1027N_alpha20_N_cells10.h5ad")

# Take a look at metadata
names(gbm1@meta.data)
table(gbm1$seurat_clusters)

# Set default assay
DefaultAssay(gbm1) <- 'Spatial'

# Look at QC
plot1 <- VlnPlot(gbm1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(gbm1, features = "nCount_Spatial",pt.size.factor = 5,stroke = 0) + theme(legend.position = "right")
wrap_plots(plot1, plot2)

# Add tissue coordinates as metadata (workaround for now)
gbm1.coords <- GetTissueCoordinates(gbm1)
gbm1$x <- gbm1.coords$x
gbm1$y <- gbm1.coords$y

# Pick which genes to impute
num.cells.per.feature <- 25 # this reduces ALRA false-positives, but it is not a perfect approach
genes.to.impute <- rownames(gbm1)[rowSums(gbm1@assays$Spatial@counts>0)>num.cells.per.feature]

# Impute using ALRA
options(warn = 1)
gc()
#gbm1 <- SeuratWrappers::RunALRA(gbm1, genes.use = genes.to.impute)
# Save imputed data for later
gc()
#save(gbm1,file='gbm1.imputed.Robj')
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Collaborations/Fan Lab/gbm1.imputed.Robj")
DefaultAssay(gbm1) <- 'alra'

# Stash cluster and embedding information
DimPlot(gbm1)
SpatialDimPlot(gbm1,pt.size.factor = 5)
gbm1@reductions$umap.ph <- gbm1@reductions$umap
gbm1$phenotype <- Idents(gbm1)

# Run NICHES 
test <- RunNICHES(gbm1,
                  assay = 'alra', # note this change
                  LR.database = 'fantom5',
                  species = 'human',
                  position.x = 'x',
                  position.y = 'y',
                  k = 10,
                  rad.set = NULL,
                  CellToCellSpatial = T,
                  CellToNeighborhood = T,
                  NeighborhoodToCell = T,
                  CellToCell = F,
                  SystemToCell = F,
                  CellToSystem = F)

# Add Niches output as an assay
me.data <- GetAssayData(object =  test$NeighborhoodToCell, slot = 'data')
colnames(me.data) <- test$NeighborhoodToCell$ReceivingCell
me.assay <- CreateAssayObject(data = me.data)
gbm1[['NeighborhoodToCell']] <- me.assay
DefaultAssay(gbm1) <- 'NeighborhoodToCell'
gbm1 <- ScaleData(gbm1)
gbm1 <- FindVariableFeatures(gbm1)
gbm1 <- RunPCA(gbm1,npcs = 100)
ElbowPlot(gbm1,ndims = 100)
gbm1 <- RunUMAP(gbm1,dims = 1:10)
DimPlot(gbm1)
gbm1 <- FindNeighbors(gbm1,dims = 1:10)
gbm1 <- FindClusters(gbm1,res=0.2)
DimPlot(gbm1)
SpatialDimPlot(gbm1,pt.size.factor = 5)
gbm1$me.clusters <- Idents(gbm1)
gbm1@reductions$umap.me <- gbm1@reductions$umap

# Make a clear plot of the data
cols.use <- list()
cols.use[['me']] <- RColorBrewer::brewer.pal(n = 12, name = "Paired")[c(7:12)]
names(cols.use[['me']]) <- unique(gbm1$me.clusters)
cols.use[['ph']] <- RColorBrewer::brewer.pal(n = 12, name = "Paired")[c(1:6)]
names(cols.use[['ph']]) <- unique(gbm1$phenotype)
# Phenotype Plots
png('Phenotype.png',width = 4,height = 4,units = 'in',res=300)
SpatialDimPlot(gbm1,pt.size.factor = 6,stroke = 0,
               cols = cols.use[['ph']],group.by = 'phenotype')+
  ggtitle('Gene Expression')+ 
  guides(fill=guide_legend(title="Cluster"))+
  theme(plot.title = element_text(hjust = 0.5,face='bold',size = 15))+NoLegend()
dev.off()
png('Phenotype UMAP.png',width = 4,height = 4,units = 'in',res=300)
DimPlot(gbm1,cols = alpha(cols.use[['ph']],1),
        pt.size = 1,shuffle = T,
        reduction = 'umap.ph',group.by = 'phenotype',label = T)+
  ggtitle('Transcriptomic UMAP')+ 
  guides(color=guide_legend(title="Cluster"))+NoLegend()
dev.off()

# Microenvironment Plots
png('Microenvironment.png',width = 4,height = 4,units = 'in',res=300)
SpatialDimPlot(gbm1,pt.size.factor = 6,stroke = 0,cols = cols.use[['me']],group.by = 'me.clusters')+
  ggtitle('Microenvironment')+ 
  guides(fill=guide_legend(title="Microenvironment"))+
  theme(plot.title = element_text(hjust = 0.5,face='bold',size = 15))+NoLegend()
dev.off()
png('Microenvironment UMAP.png',width = 4,height = 4,units = 'in',res=300)
DimPlot(gbm1,cols = alpha(cols.use[['me']],0.75),
        pt.size = 1,shuffle = T,
        reduction = 'umap.me',group.by = 'me.clusters',label = T)+
  ggtitle('NICHES UMAP')+ 
  guides(color=guide_legend(title="Microenvironment"))+NoLegend()
dev.off()

# Find microenviroment markers
DefaultAssay(gbm1) <- 'NeighborhoodToCell'
Idents(gbm1) <- gbm1$me.clusters
mark.me <- FindAllMarkers(gbm1,min.pct = 0.1,logfc.threshold = 0.01,only.pos = F)
mark.me$ratio <- mark.me$pct.1/mark.me$pct.2
mark.me$power <- mark.me$ratio*mark.me$avg_log2FC

# Find phenotype me markers
DefaultAssay(gbm1) <- 'NeighborhoodToCell'
Idents(gbm1) <- gbm1$phenotype
mark.ph <- FindAllMarkers(gbm1,min.pct = 0.1,logfc.threshold = 0.01,only.pos = F)
mark.ph$ratio <- mark.ph$pct.1/mark.ph$pct.2
mark.ph$power <- mark.ph$ratio*mark.ph$avg_log2FC

# Scale data, just in case
DefaultAssay(gbm1) <- 'NeighborhoodToCell'
gbm1 <- ScaleData(gbm1)
DefaultAssay(gbm1) <- 'alra'
gbm1 <- ScaleData(gbm1)
DefaultAssay(gbm1) <- 'NeighborhoodToCell'

LRMPlot <- function(object,mechanism){
  # Define two side of the mechanism
  mech.L <- strsplit(mechanism,split = '—')[[1]][[1]]
  mech.R <- strsplit(mechanism,split = '—')[[1]][[2]]
  # First make transcriptomic plots
  DefaultAssay(object) <- 'alra'
  p1 <- SpatialFeaturePlot(object,
                           features = mech.L,
                           slot = 'scale.data',pt.size.factor = 5,stroke = 0)+ggtitle('Ligand')
  p2 <- SpatialFeaturePlot(object,
                           features = mech.R,
                           slot = 'scale.data',pt.size.factor = 5,stroke = 0)+ggtitle('Receptor')
  # Second make microenvironment plots
  DefaultAssay(object) <- 'NeighborhoodToCell'
  p3 <- SpatialFeaturePlot(object,
                           features = mechanism,
                           slot = 'scale.data',pt.size.factor = 5,stroke = 0)+ggtitle('Mechanism')
  # Plot Total
  plot.total <- plot_grid(p1,p2,p3,nrow=1)
  # Output
  print(plot.total)
  return(plot.total)
  plot.total <<- plot.total
}

# Good visual demonstration markers
mech <- c('WNT2—LRP6','GDF7—BMPR2','TNFSF15—TNFRSF25','CCL3—CCR1','HGF—MET')
for(i in 1:length(mech)){
png(filename = paste(mech[i],'.png'),width = 8,height = 4,units = 'in',res = 300)
LRMPlot(object = gbm1,mechanism = mech[i])
dev.off()
}

# Find phenotypic markers
DefaultAssay(gbm1) <- 'alra'
Idents(gbm1) <- gbm1$phenotype
mark.ph.ph <- FindAllMarkers(gbm1,min.pct = 0.5,logfc.threshold = 0.5,only.pos = T)
mark.ph.ph$ratio <- mark.ph.ph$pct.1/mark.ph.ph$pct.2
mark.ph.ph$power <- mark.ph.ph$ratio*mark.ph.ph$avg_log2FC

# Correlate microenvironment and phenotype
DefaultAssay(gbm1) <- 'NeighborhoodToCell'
data.me <- data.frame(feat.me = gbm1@assays$NeighborhoodToCell@data['GDF7—BMPR2',])
data.ph <- data.frame(feat.ph = gbm1@assays$alra@data['MET',])
data <- cbind(data.me,data.ph)
ggplot(data,aes(x=feat.me,y=feat.ph))+geom_point()+theme_classic()

LRMPlot(object = gbm1,mechanism = 'EFNA5—EPHA4')
LRMPlot(object = gbm1,mechanism = 'BMP2—BMPR2')
LRMPlot(object = gbm1,mechanism = 'GDF7—BMPR2')
LRMPlot(object = gbm1,mechanism = 'SHH—PTCH2')
LRMPlot(object = gbm1,mechanism = 'GDF7—BMPR2')
LRMPlot(object = gbm1,mechanism = 'WNT2—LRP6')
LRMPlot(object = gbm1,mechanism = 'FGF2—CD44')
LRMPlot(object = gbm1,mechanism = 'TNC—EGFR')
LRMPlot(object = gbm1,mechanism = 'VEGFA—KDR')

SpatialDimPlot(gbm1,pt.size.factor = 5,stroke = 0.1)
LRMPlot(object = gbm1,mechanism = 'WNT5A—FZD4') #0
LRMPlot(object = gbm1,mechanism = 'TNFSF15—TNFRSF25') #1
LRMPlot(object = gbm1,mechanism = 'BDNF—NTRK2') #2
LRMPlot(object = gbm1,mechanism = 'BMP7—BMPR2') #2

