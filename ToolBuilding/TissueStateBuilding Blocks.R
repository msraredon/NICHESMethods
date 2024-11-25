# Building blocks for 'tissue embedding'
# MSBR 2024-11-25
# For Sophie & Nuoya

# Multiple potential approaches:

# 1. Relative cell fraction
# a. choice of grouping variable (cell class, cell type, cell state, cycle phase, etc.)
# b. unexplored

# 2. Signal Graph Alignment
# a. choice of graph-graph distance/correlation metric
# b. choice of how to handle zeros in the data
# c. working prototype

# 3. Gene Graph Alignment
# a. same as above
# b. same as above
# c. unexplored

# packages
require(Seurat)

# Setup
load("~/Desktop/unsorted scripts/pplr.in.process.2024-11-19.Robj")

# Inspect
data
table(data$orig.ident)

# Stash for generalization
object <- data

# Split into individual tissues
split <- SplitObject(object,split.by = 'orig.ident')
split
names(split)

# Break each individual tissue into smaller subsamples (bootstrapping)
n.cells.per.chunk <- 500
boot.list <- list()
for(i in 1:length(split)){
  randomly.ordered.barcodes <- sample(size = length(colnames(split[[i]])),colnames(split[[i]]),replace = F)
  barcode.boot.list <- split(randomly.ordered.barcodes, ceiling(seq_along(randomly.ordered.barcodes)/n.cells.per.chunk))
  temp.list <- list()
  for(j in 1:length(barcode.boot.list)){
    temp.list[[j]] <- subset(split[[i]],cells = barcode.boot.list[[j]])
  }
  names(temp.list) <- paste('Sample',c(1:length(temp.list)),sep = '_')
  boot.list[[i]] <- temp.list
}
names(boot.list) <- names(split)
boot.list

# Flatten this list
boot.list.2 <- purrr::list_flatten(boot.list)
length(boot.list.2)


# Initialize all possible features
# Rownames here will be celltypes
grouping.variable <- 'CellType'
row.names <- names(table(object[[grouping.variable]]))
row.names

# Quantify tissue fraction by grouping variable
tissue.sample.prop <- list()
for(i in 1:length(boot.list.2)){ 
    tissue.sample.prop[[i]] <- table(boot.list.2[[i]][[grouping.variable]])/sum(table(boot.list.2[[i]][[grouping.variable]]))
  }
names(tissue.sample.prop) <- names(boot.list.2)

# Convert each fraxtional table to a dataframe
for(i in 1:length(tissue.sample.prop)){
  tissue.sample.prop[[i]] <- data.frame(tissue.sample.prop[[i]])
}

# Make a dummy
blank.table <- data.frame(TissueFraction = rep(0,length(row.names)))
rownames(blank.table) <- row.names
blank.table

# Fill the dummy for each tissue sample
tissue.sample.prop.final <- list()
for(i in 1:length(tissue.sample.prop)){
  temp <- blank.table
  temp$TissueFraction[tissue.sample.prop[[i]]$CellType] <- tissue.sample.prop[[i]]$Freq
  tissue.sample.prop.final[[i]] <- temp
}

# Concatenate into a single dataframe
output <- do.call(cbind, tissue.sample.prop.final)

# Give column names
colnames(output) <- names(tissue.sample.prop)

# Check that it adds to 1 for each tissue
colSums(output)

# Make into feature-barcode matrix
seurat.assay.obj <- CreateAssay5Object(data = output)
tissue.state.obj <- CreateSeuratObject(seurat.assay.obj)
tissue.state.obj

# Make UMAP
# NO NORMALIZATION
tissue.state.obj <- ScaleData(tissue.state.obj)
RunPCA(tissue.state.obj,features = rownames(tissue.state.obj))
#PCHeatmap(tissue.state.obj)
tissue.state.obj <- RunUMAP(tissue.state.obj,features = rownames(tissue.state.obj))
DimPlot(tissue.state.obj)

FeaturePlot(tissue.state.obj,'Mac-Alv')
