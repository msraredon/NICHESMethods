# Set WD
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Collaborations/Kyriakides Collab/Yaqing")

# Set seed
set.seed(2)

# Packages
library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(scales)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(cowplot)

#### SetUp ####
# Functions (from Allie) # Some are not fully generalized #
LoadSeuratData <- function(file.names,sample.names){
  data <- list()
  for (i in 1:length(file.names)){
    message(sample.names[i])
    load(file = paste("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Collaborations/Kyriakides Collab/Yaqing/",file.names[i],sep='')) ## Script specific workaround
    data[[i]] <- CreateSeuratObject(counts = output[['Gene']])
    data[[i]][['GeneFull']] <- CreateAssayObject(output[['GeneFull']])
    data[[i]][['spliced']] <- CreateAssayObject(output[['spliced']])
    data[[i]][['unspliced']] <- CreateAssayObject(output[['unspliced']])
    data[[i]][["percent.mt"]] <- PercentageFeatureSet(data[[i]], pattern = "^Mt-")
    data[[i]]$perc.spliced <- 100*(colSums(data[[i]]@assays$spliced)/(colSums(data[[i]]@assays$spliced)+colSums(data[[i]]@assays$unspliced)))
    data[[i]]$Sample <- sample.names[i]
    gc()
  }
  return(data)
}
LoadMetaData <- function(seurat.data){
  data <- data.frame()
  for (i in 1:length(seurat.data)){
    data <- rbind(data,seurat.data[[i]]@meta.data)
    gc()
  }
  return(data)
}
QCViolinPlotter <- function(sample.name,min.counts,min.features,max.mt,tag){
  x <- which(sample.name == names(seurat.data))
  object <- seurat.data[[x]]
  filtered <- subset(object,nFeature_RNA > min.features & nCount_RNA > min.counts & percent.mt < max.mt)
  p <- VlnPlot(object, features = c("nCount_RNA"),pt.size=0.1,log=T) +
    geom_hline(yintercept=min.counts, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Min.Counts = ",min.counts)) + theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
  q <- VlnPlot(object, features = c("nFeature_RNA"),pt.size=0.1,log=T) +
    geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Min.Features = ",min.features)) + theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
  r <- VlnPlot(object, features = c("percent.mt"),pt.size=0.1,y.max=50) +
    geom_hline(yintercept=max.mt, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Max.Mt = ",max.mt)) +
    labs(caption = paste("Barcodes retained with these filters from starting total: ",ncol(filtered),"/",ncol(object))) +
    theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'),plot.caption = element_text(color='red'))
  png(paste(sample.name,'_filterselect_1',tag,'.png',sep=""),width = 1500,height=500)
  print(plot_grid(p,q,r,ncol=3))
  dev.off()
}
TriplePercentMtPlotter <- function(sample.name,min.counts,min.features,max.mt,tag){
  object <- subset(pneum.metadata,Sample %in% sample.name)
  lowpass <- subset(object,percent.mt < max.mt)
  highmt <- subset(object,percent.mt > max.mt)
  filtered <- subset(object,nFeature_RNA > min.features & nCount_RNA > min.counts & percent.mt < max.mt)
  p <- ggplot() +
    geom_point(data=highmt, aes(x=nCount_RNA, y=nFeature_RNA), color = 'gray') +
    geom_point(data=lowpass, aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
    scale_color_viridis(limits = c(0, 50), oob = scales::squish) +
    geom_vline(xintercept=min.counts, linetype='dashed', color='red', size=1) +
    geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
    labs(title = sample.name) +
    labs(subtitle = paste("Min.Counts = ",min.counts," / Min.Features = ",min.features," / Max.Mt = ",max.mt)) +
    labs(caption = paste("Barcodes retained with these filters from starting total: ",nrow(filtered),"/",nrow(object))) +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          plot.subtitle = element_text(color='red'),
          plot.caption = element_text(color='red'),
          axis.title.x = element_text(size = 20,color='black'),
          axis.title.y = element_text(size = 20,color='black'),
          axis.text.x = element_text(size=12,color='black'),
          axis.text.y = element_text(size = 12,color='black'),
          axis.line = element_line(size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  png(paste(sample.name,'_filterselect_2',tag,'.png',sep=""),width = 600,height=500)
  print(p)
  dev.off()
}
TriplePercentSplPlotter <- function(sample.name,min.counts,min.features,max.mt,tag){
  object <- subset(pneum.metadata,Sample %in% sample.name)
  lowpass <- subset(object,percent.mt < max.mt)
  highmt <- subset(object,percent.mt > max.mt)
  filtered <- subset(object,nFeature_RNA > min.features & nCount_RNA > min.counts & percent.mt < max.mt)
  p <- ggplot() +
    geom_point(data=highmt, aes(x=nCount_RNA, y=nFeature_RNA), color = 'gray') +
    geom_point(data=lowpass, aes(x=nCount_RNA, y=nFeature_RNA, color=perc.spliced)) +
    scale_color_viridis(option="magma",direction=-1,limits = c(50, 100), oob = scales::squish) +
    geom_vline(xintercept=min.counts, linetype='dashed', color='red', size=1) +
    geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
    labs(title = sample.name) +
    labs(subtitle = paste("Min.Counts = ",min.counts," / Min.Features = ",min.features," / Max.Mt = ",max.mt)) +
    labs(caption = paste("Barcodes retained with these filters from starting total: ",nrow(filtered),"/",nrow(object))) +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          plot.subtitle = element_text(color='red'),
          plot.caption = element_text(color='red'),
          axis.title.x = element_text(size = 20,color='black'),
          axis.title.y = element_text(size = 20,color='black'),
          axis.text.x = element_text(size=12,color='black'),
          axis.text.y = element_text(size = 12,color='black'),
          axis.line = element_line(size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  png(paste(sample.name,'_filterselect_3',tag,'.png',sep=""),width = 600,height=500)
  print(p)
  dev.off()
}

# More Functions (from Sam)  # Some may not be fully generalized #
theProcess <- function(sample.name,min.counts,min.features,max.mt,tag){
  # Plot QC metrics
  QCViolinPlotter(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)
  TriplePercentMtPlotter(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)
  TriplePercentSplPlotter(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)
  
  # Pull object, apply subset, & dimensional reduction
  x <- which(sample.name == names(seurat.data))
  object <- seurat.data[[x]]
  object <- subset(object, nCount_RNA > min.counts & nFeature_RNA > min.features & percent.mt < max.mt)
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object)
  object <- ScaleData(object) # don't regress on anything for now, see if it matters
  object <- RunPCA(object, npcs = 50, verbose = F)
  png(paste(sample.name,'_elbow.png',tag,sep=""), width = 800, height = 500)
  print(ElbowPlot(object,ndims = 50)+labs(title = sample.name))
  dev.off()
  pdf(paste(sample.name,'_pcheatmap_1',tag,'.pdf',sep=""), width = 12, height = 16)
  print(DimHeatmap(object, dims = 1:15, cells = 200, balanced = T))
  dev.off()
  pdf(paste(sample.name,'_pcheatmap_2',tag,'.pdf',sep=""), width = 12, height = 16)
  print(DimHeatmap(object, dims = 16:30, cells = 200, balanced = T))
  dev.off()
  
  # Cluster, plot, & check QC
  pcs <- 30
  object <- FindNeighbors(object, dims = 1:pcs)
  object <- FindClusters(object, resolution = 0.5)
  object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
  png(paste(sample.name,'_umap_1',tag,'.png',sep=""), width = 500, height = 500)
  print(UMAPPlot(object,label = T) + labs(title = sample.name,subtitle = paste("PC's = ",pcs)) + NoLegend() + NoAxes())
  dev.off()
  png(paste(sample.name,'_fp_lin',tag,'.png',sep=""), width = 800, height = 800)
  print(FeaturePlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'), label = T,repel=T))
  dev.off()
  png(paste(sample.name,'_fp_qc',tag,'.png',sep=""), width = 800, height = 1200)
  print(FeaturePlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),ncol=2,label = T,repel=T))
  dev.off()
  png(paste(sample.name,'_vln_lin',tag,'.png',sep=""), width = 800, height = 800)
  print(VlnPlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'),ncol=2,pt.size=0.1))
  dev.off()
  png(paste(sample.name,'_vln_qc',tag,'.png',sep=""), width = 800, height = 800)
  print(VlnPlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),ncol=2,pt.size=0.1))
  dev.off()
  
  ## Backtrack to original Filtration plots
  nclusters <- length(levels(object$seurat_clusters))
  # ViolinPlots broken apart by cluster
  p <- VlnPlot(object, features = c("nCount_RNA"),pt.size=0.1,log=T) +
    geom_hline(yintercept=min.counts, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Min.Counts = ",min.counts)) + theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
  q <- VlnPlot(object, features = c("nFeature_RNA"),pt.size=0.1,log=T) +
    geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Min.Features = ",min.features)) + theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
  r <- VlnPlot(object, features = c("percent.mt"),pt.size=0.1,y.max=(max.mt+0.5)) +
    geom_hline(yintercept=max.mt, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Max.Mt = ",max.mt)) +
    theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'),plot.caption = element_text(color='red'))
  png(paste(sample.name,'_filterback_1',tag,'.png',sep=""),width = 1500,height=500)
  print(plot_grid(p,q,r,ncol=3))
  dev.off()
  
  # ViolinPlots not broken apart, jitter points colored by cluster
  p <- VlnPlot(object, features = c("nCount_RNA"),pt.size=0,log=T,group.by='orig.ident',cols=1) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    geom_hline(yintercept=min.counts, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Min.Counts = ",min.counts)) + theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
  q <- VlnPlot(object, features = c("nFeature_RNA"),pt.size=0,log=T,group.by='orig.ident',cols=1) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Min.Features = ",min.features)) + theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
  r <- VlnPlot(object, features = c("percent.mt"),pt.size=0,group.by='orig.ident',cols=1,y.max=(max.mt+0.5)) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    geom_hline(yintercept=max.mt, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Max.Mt = ",max.mt)) +
    theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
  png(paste(sample.name,'_filterback_2',tag,'.png',sep=""),width = 1500,height=500)
  print(plot_grid(p,q,r,ncol=3))
  dev.off()
  
  # Scatter Plot colored by cluster
  list <- subset(pneum.metadata,Sample %in% sample.name)
  filtered <- subset(list,nCount_RNA > min.counts & nFeature_RNA > min.features & percent.mt < max.mt)
  
  png(paste(sample.name,'_filterback_3',tag,'.png',sep=""),width = 600,height=500)
  print(
    ggplot() +
      geom_point(data=filtered, aes(x=nCount_RNA, y=nFeature_RNA, color=factor(object$seurat_clusters)), size = 0.5) +
      scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
      geom_vline(xintercept=min.counts, linetype='dashed', color='red', size=1) +
      geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
      labs(title = sample.name) +
      labs(subtitle = paste("Min.Counts = ",min.counts," / Min.Features = ",min.features," / Max.Mt = ",max.mt)) +
      theme(plot.title = element_text(size = 24, hjust = 0.5),
            plot.subtitle = element_text(color='red'),
            plot.caption = element_text(color='red'),
            axis.title.x = element_text(size = 20,color='black'),
            axis.title.y = element_text(size = 20,color='black'),
            axis.text.x = element_text(size=12,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.line = element_line(size = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())
  )
  dev.off()
  
  # Ribosomal Characterization
  object[['percent.ribo']] <- PercentageFeatureSet(object, pattern = "^Rp")
  png(paste(sample.name,'_fp_ribo',tag,'.png',sep=""), width = 500, height = 500)
  print(FeaturePlot(object, c('percent.ribo'), label = T))
  dev.off()
  png(paste(sample.name,'_vln_ribo',tag,'.png',sep=""), width = 500, height = 500)
  print(VlnPlot(object, c('percent.ribo'),pt.size=0.1))
  dev.off()
  
  # Epithelial Characterization
  png(paste(sample.name,'_fp_c5',tag,'.png',sep=""), width = 800, height = 800)
  print(FeaturePlot(object, c('Sftpc','Sftpb','Defb4','Lyz2'), label = T, repel = T))
  dev.off()
  png(paste(sample.name,'_fp_epi',tag,'.png',sep=""), width = 800, height = 800)
  print(FeaturePlot(object, c('Krt5','Sftpc','Hopx','Aqp5'), label = T, repel = T))
  dev.off()
  
  # Pass output to outside for further exploration
  object <<- object # Passes this outside of the function
}
titrationUMAP <- function(object,min.counts,min.features,max.mt){
  object@meta.data$Passing <- 'NoPass'
  cells.Pass <- WhichCells(object,expression = nCount_RNA > min.counts & nFeature_RNA > min.features & percent.mt < max.mt)
  meta.data.pass <- data.frame(Barcode = cells.Pass,Passing = 'Pass')
  rownames(meta.data.pass) <- cells.Pass
  if(length(cells.Pass) < ncol(object)){
    cells.NoPass <- WhichCells(object,cells = cells.Pass,invert = T)
    meta.data.nopass <- data.frame(Barcode = cells.NoPass,Passing = 'NoPass')
    rownames(meta.data.nopass) <- cells.NoPass
  }else{meta.data.nopass <- data.frame()}
  meta.data <- rbind(meta.data.pass,meta.data.nopass)
  object <- AddMetaData(object,metadata = meta.data)
  tag <- paste('min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
  png(paste(sample.name,'_TitrationDimPlot.',tag,'.png',sep=""), width = 500, height = 500)
  print(DimPlot(object,group.by = 'Passing',cols = c('#0529F4','#F42F05'),shuffle = T)
        + labs(title = sample.name,subtitle = paste('nCount_RNA >',min.counts,'& nFeature_RNA >',min.features,'& percent.mt <',max.mt)))
  dev.off()
}

####### Start here, mainly ########

# Load data from 10X raw DGE output folder
tk1 <- Read10X("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Collaborations/Kyriakides Collab/Yaqing/TK1/raw_feature_bc_matrix")
wt1 <- Read10X("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Collaborations/Kyriakides Collab/Yaqing/WT1/raw_feature_bc_matrix")

# stash the gene expression data (will use later as main data framework)
tk1.gene <- tk1$`Gene Expression`
wt1.gene <- wt1$`Gene Expression`

# stash the cell hashing data to manipulate and generate labels
tk1.hash <- tk1$`Antibody Capture`
wt1.hash <- wt1$`Antibody Capture`

# take a look
colSums(tk1.hash)

# convert to df
tk1.hash.df <- data.frame(t(as.matrix(tk1.hash)))
wt1.hash.df <- data.frame(t(as.matrix(wt1.hash)))

# establish label null column
tk1.hash.df$hash.label <- NA
wt1.hash.df$hash.label <- NA
# ratio column
tk1.hash.df$ratio1 <- tk1.hash.df$B0301_anti.mouse_Hashtag_1_Antibody/tk1.hash.df$B0302_anti.mouse_Hashtag_2_Antibody 
tk1.hash.df$ratio2 <- tk1.hash.df$B0302_anti.mouse_Hashtag_2_Antibody/tk1.hash.df$B0301_anti.mouse_Hashtag_1_Antibody
wt1.hash.df$ratio1 <- wt1.hash.df$B0301_anti.mouse_Hashtag_1_Antibody/wt1.hash.df$B0302_anti.mouse_Hashtag_2_Antibody 
wt1.hash.df$ratio2 <- wt1.hash.df$B0302_anti.mouse_Hashtag_2_Antibody/wt1.hash.df$B0301_anti.mouse_Hashtag_1_Antibody

# total column
tk1.hash.df$total <- rowSums(tk1.hash.df[,c(1:2)])
wt1.hash.df$total <- rowSums(wt1.hash.df[,c(1:2)])

# if the total is 0, then label as unknown
tk1.hash.df[tk1.hash.df$total==0,]$hash.label <- 'Unknown'
wt1.hash.df[wt1.hash.df$total==0,]$hash.label <- 'Unknown'

# ratio between the two values > threshold in one direction, then label as that sample
#threshold = 10
tk1.hash.temp <- tk1.hash.df
wt1.hash.temp <- wt1.hash.df

threshold.list <- seq(1, 10, by = 0.1)
output.list <- data.frame()
for (i in 1:length(threshold.list)){
  tk1.hash.temp <- tk1.hash.df
  wt1.hash.temp <- wt1.hash.df
  tk1.hash.temp[!is.na(tk1.hash.temp$ratio1) & tk1.hash.temp$ratio1>threshold.list[i],]$hash.label <- 'TK1.Sample1'
  wt1.hash.temp[!is.na(wt1.hash.temp$ratio1) & wt1.hash.temp$ratio1>threshold.list[i],]$hash.label <- 'WT1.Sample1'
  tk1.hash.temp[!is.na(tk1.hash.temp$ratio2) & tk1.hash.temp$ratio2>threshold.list[i],]$hash.label <- 'TK1.Sample2'
  wt1.hash.temp[!is.na(wt1.hash.temp$ratio2) & wt1.hash.temp$ratio2>threshold.list[i],]$hash.label <- 'WT1.Sample2'
  # if the coSum is >0, and the ratio between the two values is <= threshold.list[i], then label as suspected multiplet
  tk1.hash.temp[is.na(tk1.hash.temp$hash.label),]$hash.label <- 'Multiplet'
  wt1.hash.temp[is.na(wt1.hash.temp$hash.label),]$hash.label <- 'Multiplet'
  # Look at distributions
  df1 <- data.frame(table(tk1.hash.temp$hash.label))
  df1$sample <- 'TK1'
  df2 <- data.frame(table(wt1.hash.temp$hash.label))
  df2$sample <- 'WT1'
  df.out <- rbind(df1,df2)
  df.out$threshold <- threshold.list[i]
  output.list <- rbind(output.list,df.out)
}

View(output.list) # formatted for ggplot2

# take a look at the trend
ggplot(data = output.list[output.list$Var1=='Multiplet',],
       aes(x=threshold,y=Freq,color = sample))+geom_point()

ggplot(data = output.list[output.list$sample=='TK1',],
       aes(x=threshold,y=Freq,color = Var1))+geom_point()+scale_y_log10()
ggplot(data = output.list[output.list$sample=='WT1',],
       aes(x=threshold,y=Freq,color = Var1))+geom_point()+scale_y_log10()

# I think a threshold of 1 actually might make sense here
threshold <- 1
tk1.hash.df[!is.na(tk1.hash.df$ratio1) & tk1.hash.df$ratio1>threshold.list[i],]$hash.label <- 'TK1.Sample1'
wt1.hash.df[!is.na(wt1.hash.df$ratio1) & wt1.hash.df$ratio1>threshold.list[i],]$hash.label <- 'WT1.Sample1'
tk1.hash.df[!is.na(tk1.hash.df$ratio2) & tk1.hash.df$ratio2>threshold.list[i],]$hash.label <- 'TK1.Sample2'
wt1.hash.df[!is.na(wt1.hash.df$ratio2) & wt1.hash.df$ratio2>threshold.list[i],]$hash.label <- 'WT1.Sample2'
# if the coSum is >0, and the ratio between the two values is <= threshold.list[i], then label as suspected multiplet
tk1.hash.df[is.na(tk1.hash.df$hash.label),]$hash.label <- 'Multiplet'
wt1.hash.df[is.na(wt1.hash.df$hash.label),]$hash.label <- 'Multiplet'

# Create Seurat object
tk1.seurat <- CreateSeuratObject(tk1.gene,min.cells = 3,min.features = 200) # Be careful with these thresholds. I know that min.features = 200 is OK here, because this 10X v3 data which has very high UMI
wt1.seurat <- CreateSeuratObject(wt1.gene,min.cells = 3,min.features = 200)

# Add metdata
tk1.seurat$Condition <- 'Tsp2.Knockout'
wt1.seurat$Condition <- 'Wild.Type'

# Add hashtag data
tk1.hash.metadata <- tk1.hash.df
tk1.hash.metadata <- tk1.hash.metadata[colnames(tk1.seurat),]
wt1.hash.metadata <- wt1.hash.df
wt1.hash.metadata <- wt1.hash.metadata[colnames(wt1.seurat),]
?AddMetaData
tk1.seurat <- AddMetaData(tk1.seurat,metadata = tk1.hash.metadata)
wt1.seurat <- AddMetaData(wt1.seurat,metadata = wt1.hash.metadata)

# Calculate percentage mitochondrial reads
tk1.seurat[["percent.mt"]] <- PercentageFeatureSet(tk1.seurat, pattern = "^mt-") # Mouse: "^mt-" | Rat: "^Mt-" | Human: "^MT-"
wt1.seurat[["percent.mt"]] <- PercentageFeatureSet(wt1.seurat, pattern = "^mt-") 

# Merge these two samples into one R object
merge <- merge(tk1.seurat,wt1.seurat)
merge

# check hashing information
table(merge$hash.label)

# Look at QC
png(filename = 'First_Look_QC_Unfiltered.png',width = 7,height = 5,units = 'in',res=300)
VlnPlot(merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)
dev.off()
VlnPlot(merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,log = T)

#### # QC Plots from Luecken 2019, Fig. 2 - Replicas built in R ####
# Panel D - nFeature_RNA vs nUMI, colored by percent.mt
ggplot(merge@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(aes(color=percent.mt)) +
    scale_color_viridis(limits = c(0, 50), oob = scales::squish) +
    labs(title = "Yaqing's Data") +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.title.x = element_text(size = 20,color='black'),
          axis.title.y = element_text(size = 20,color='black'),
          axis.text.x = element_text(size=12,color='black'),
          axis.text.y = element_text(size = 12,color='black'),
          axis.line = element_line(size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())

# Filter the data first pass, be careful
sub <- subset(merge,subset = 
                nFeature_RNA > 1000 &
                percent.mt < 10 &
                percent.mt > 0.1)

VlnPlot(sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,log = T)
VlnPlot(sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,log = T,group.by = 'Condition')


######################################
###### UNADAPTED LEGACY CODE BELOW ######
#### First Look All Project Data ####
sub <- NormalizeData(sub)
sub <- ScaleData(sub)
sub <- FindVariableFeatures(sub)
sub <- RunPCA(sub, npcs = 100)
pdf(file='sub.PCs.pdf',width=10,height=8)
ElbowPlot(sub,ndims = 100)
PCHeatmap(sub,cells=200,balanced=T,dims=1:9)
PCHeatmap(sub,cells=200,balanced=T,dims=10:18)
PCHeatmap(sub,cells=200,balanced=T,dims=19:27)
PCHeatmap(sub,cells=200,balanced=T,dims=28:36)
PCHeatmap(sub,cells=200,balanced=T,dims=37:45)
PCHeatmap(sub,cells=200,balanced=T,dims=46:54)
PCHeatmap(sub,cells=200,balanced=T,dims=55:63)
PCHeatmap(sub,cells=200,balanced=T,dims=64:72)
PCHeatmap(sub,cells=200,balanced=T,dims=73:81)
PCHeatmap(sub,cells=200,balanced=T,dims=82:90)
PCHeatmap(sub,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
sub <- RunUMAP(sub, reduction = "pca", dims = 1:20)
sub <- FindNeighbors(sub, reduction = "pca", dims = 1:20)
sub <- FindClusters(sub, resolution = 0.2)
DimPlot(sub, reduction = "umap", label = TRUE, repel = TRUE)
FeaturePlot(sub,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T) # Note, no obvious multiplet clusters
FeaturePlot(sub,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T) # Cluster 6, 4 , and 2 look low-quality / low-information.
VlnPlot(sub,c('nFeature_RNA','nCount_RNA','percent.mt'))

# Remove garbage clusters, re-scale, and redo!
sub.2 <- subset(sub,idents = c('2','4','6'),invert=T) # from here, go back up to line 333

sub.2 <- ScaleData(sub.2)
sub.2 <- FindVariableFeatures(sub.2)
sub.2 <- RunPCA(sub.2, npcs = 100)
pdf(file='sub.2.PCs.pdf',width=10,height=8)
ElbowPlot(sub.2,ndims = 100)
PCHeatmap(sub.2,cells=200,balanced=T,dims=1:9)
PCHeatmap(sub.2,cells=200,balanced=T,dims=10:18)
PCHeatmap(sub.2,cells=200,balanced=T,dims=19:27)
PCHeatmap(sub.2,cells=200,balanced=T,dims=28:36)
PCHeatmap(sub.2,cells=200,balanced=T,dims=37:45)
PCHeatmap(sub.2,cells=200,balanced=T,dims=46:54)
PCHeatmap(sub.2,cells=200,balanced=T,dims=55:63)
PCHeatmap(sub.2,cells=200,balanced=T,dims=64:72)
PCHeatmap(sub.2,cells=200,balanced=T,dims=73:81)
PCHeatmap(sub.2,cells=200,balanced=T,dims=82:90)
PCHeatmap(sub.2,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster the now clean(er) object
sub.2 <- RunUMAP(sub.2, reduction = "pca", dims = 1:20)
sub.2 <- FindNeighbors(sub.2, reduction = "pca", dims = 1:20)
sub.2 <- FindClusters(sub.2, resolution = 0.2)
pdf(file = 'sub.2.first.look.pdf',width = 7,height = 5)
DimPlot(sub.2, reduction = "umap", label = TRUE, repel = TRUE)
FeaturePlot(sub.2,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T) # Note, no obvious multiplet clusters
FeaturePlot(sub.2,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T) # Cluster 6, 4 , and 2 look low-quality / low-information.
VlnPlot(sub.2,c('nFeature_RNA','nCount_RNA','percent.mt'))
dev.off()

# Remove garbage clusters again, re-scale, and redo!
mark.3 <- FindMarkers(sub.2,ident.1 = '3')
mark.3$ratio <- mark.3$pct.1/mark.3$pct.2
View(mark.3)
sub.3 <- subset(sub.2,idents = c('3'),invert=T)
sub.3
sub.3 <- ScaleData(sub.3)
sub.3 <- FindVariableFeatures(sub.3)
sub.3 <- RunPCA(sub.3, npcs = 100)
pdf(file='sub.3.PCs.pdf',width=10,height=8)
ElbowPlot(sub.3,ndims = 100)
PCHeatmap(sub.3,cells=200,balanced=T,dims=1:9)
PCHeatmap(sub.3,cells=200,balanced=T,dims=10:18)
PCHeatmap(sub.3,cells=200,balanced=T,dims=19:27)
PCHeatmap(sub.3,cells=200,balanced=T,dims=28:36)
PCHeatmap(sub.3,cells=200,balanced=T,dims=37:45)
PCHeatmap(sub.3,cells=200,balanced=T,dims=46:54)
PCHeatmap(sub.3,cells=200,balanced=T,dims=55:63)
PCHeatmap(sub.3,cells=200,balanced=T,dims=64:72)
PCHeatmap(sub.3,cells=200,balanced=T,dims=73:81)
PCHeatmap(sub.3,cells=200,balanced=T,dims=82:90)
PCHeatmap(sub.3,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster the now clean(er) object
sub.3 <- RunUMAP(sub.3, reduction = "pca", dims = 1:30)
sub.3 <- FindNeighbors(sub.3, reduction = "pca", dims = 1:30)
sub.3 <- FindClusters(sub.3, resolution = 0.3)
pdf(file = 'sub.3.first.look.pdf',width = 7,height = 5)
DimPlot(sub.3, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(sub.3, reduction = "umap",group.by = 'Condition',shuffle=T)
FeaturePlot(sub.3,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T) # Note, no obvious multiplet clusters
FeaturePlot(sub.3,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T) # Cluster 6, 4 , and 2 look low-quality / low-information.
VlnPlot(sub.3,c('nFeature_RNA','nCount_RNA','percent.mt'))
dev.off()

# Find markers
mark <- FindAllMarkers(sub.3,only.pos=T)
mark$ratio <- mark$pct.1/mark$pct.2
mark$power <- mark$ratio*mark$avg_log2FC
View(mark)

write.table(mark,file = 'marker.list.clusters.2023-08-08.csv',sep = ',',row.names = T,col.names = NA)

# rename idents preliminary
sub.3 <- RenameIdents(sub.3,
                      '0'='Frzb+_Fibroblasts',
                      '1'='Lgr5+_Fibroblasts',
                      '2'='Fibroblasts_cycling',
                      '3'='Dkk2+_Fibroblasts',
                      '4'='Schwann',
                      '5'='Sfrp2+_Fibroblasts', # https://www.sciencedirect.com/science/article/pii/S0022202X17330695
                      '6'='Schwann_cycling',
                      '7'='Cthrc1+_Fibroblasts',
                      '8'='Epithelium',
                      '9'='Endothelium',
                      '10'='Pericytes',
                      '11'='Immune',
                      '12'='Melanocyte_Progenitors') # https://www.cell.com/cell-reports/pdf/S2211-1247(22)00143-7.pdf

sub.3$first_pass <- Idents(sub.3)

DimPlot(sub.3, reduction = "umap", label = TRUE, repel = TRUE,split.by = 'Condition')
DimPlot(sub.3, reduction = "umap", label = TRUE, repel = TRUE,split.by = 'Condition',group.by = 'first_pass')
DimPlot(sub.3, reduction = "umap", label = TRUE, repel = TRUE,split.by='hash.label')

# Find markers with labeled clusters
mark <- FindAllMarkers(sub.3,only.pos=T)
mark$ratio <- mark$pct.1/mark$pct.2
mark$power <- mark$ratio*mark$avg_log2FC
View(mark)

write.table(mark,file = 'marker.list.clusters.named.2023-08-08.csv',sep = ',',row.names = T,col.names = NA)

# plot cell numbers per each cluster
cell.frac <- table(Idents(sub.3),sub.3$Condition)
sums <- colSums(cell.frac)
cell.frac[,1] <- cell.frac[,1]/sums[1]
cell.frac[,2] <- cell.frac[,2]/sums[2]
colSums(cell.frac) # should equal 1
cell.frac <- data.frame(cell.frac)
pdf(file = 'sub.3.cell.type.distributions.pdf',width = 8,height = 6)
ggplot(cell.frac,aes(x=Var1,y=Freq,group=Var2,fill=Var2))+
  geom_bar(stat = 'identity',position = 'dodge')+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

# some exploratory plots
FeaturePlot(sub.3,'Lgr5',order=T)

# two primary forms of differential analysis in a two-condition single cell project

# 1. For a given cluster which has representation in both conditions, what is the difference in character of the population in one condition vs. the other
lgr5.obj <- subset(sub.3,idents = 'Lgr5+_Fibroblasts')
Idents(lgr5.obj) <- lgr5.obj$Condition
lgr5.mark <- FindAllMarkers(lgr5.obj,only.pos = T)
lgr5.mark$ratio <- lgr5.mark$pct.1/lgr5.mark$pct.2
lgr5.mark$power <- lgr5.mark$ratio*lgr5.mark$avg_log2FC
View(lgr5.mark)

# you can write a for loop that just does this for every cluster and stores the output in a single data table
# useful to experiment with making heatmaps showing the top differenital genes for each cluster across conditions



# 2. For a given cluster which has representation in both conditions,
# what is the difference between the global marker genes in one condition for that population vs. 
# the global marker genes in the other condition for that population?
# HOLD on this untill we might need it.








