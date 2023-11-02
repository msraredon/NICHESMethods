# PART 1 This script performs the initial NICHES connectomics workflow for Dr. Rivero's data #

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
library(future) # enables multithreadhing # see https://satijalab.org/seurat/archive/v3.0/future_vignette.html

# change the current plan to access parallelization
plan("multisession", workers = 8)  #adjust based on how many cores you have. I have 10 total so i am reserving 2 for other tasks
plan()
options(future.globals.maxSize= 10000*1024^2) # 10 Gb futures export limit, bug fix for below

# Load Data
load('frl.combined.annotated.classed.2023-10-09.Robj')
data <- frl.combined
rm(frl.combined)

# Inspect Data
names(data@meta.data)
VlnPlot(data,c('nFeature_RNA','nCount_RNA','percent.mt'))
Idents(data) <- data$prelim.celltypes
table(Idents(data),data$Sample)
table(Idents(data),data$Condition)
table(Idents(data),data$Timepoint)

# Impute the dataset and save the output for later
# Pick which genes to impute
num.cells.per.feature <- 50 # this reduces ALRA false-positives, but it is not a perfect approach
genes.to.impute <- rownames(data)[rowSums(data@assays$RNA@counts>0)>num.cells.per.feature]
# impute using ALRA
options(warn = 1)
gc()
data <- SeuratWrappers::RunALRA(data, genes.use = genes.to.impute)

# Save imputed data for later
gc()
save(data,file='data.imputed.Robj')

####  Run NICHES [IMPUTED] #### 
# Need to split by 'Sample' first so that we are only crossing cells with other cells in the same tissue
gc()
split <- SplitObject(data,split.by='Sample')
names(split)

# Cell To Cell (CTC)
scc.list <- list()
for(i in 1:length(split)){
  print(i)
  scc.list[[i]] <- RunNICHES(split[[i]],
                             LR.database="fantom5", # we could update this ground truth if we like, but this is easy for now
                             species="mouse",
                             assay="alra",
                             cell_types = "prelim.celltypes",
                             meta.data.to.map = names(split[[i]]@meta.data), # this makes sure all of the original metdata maps into the NICHES output(s)
                             SystemToCell = F,
                             CellToCell = T,
                             CellToSystem=F) # Need to standarize the number of cells in the System measurements
}
names(scc.list) <- names(split)

# Merge outputs
temp.list <- list()
for(i in 1:length(scc.list)){
  temp.list[[i]] <- scc.list[[i]]$CellToCell # Isolate CellToCell Signaling
  gc()
}
cell.to.cell.imputed <- merge(temp.list[[1]],temp.list[2:length(temp.list)])
cell.to.cell.imputed

# Format metadata for plotting
cell.to.cell.imputed$Sample <- factor(cell.to.cell.imputed$Sample.Sending)
cell.to.cell.imputed$Condition <- factor(cell.to.cell.imputed$Condition.Sending)

save(cell.to.cell.imputed,file='cell.to.cell.imputed.Robj')

# STC and CTS
# Standardizing cell number for each sample
split <- SplitObject(data,split.by = 'Sample') # redundant but just in case
names(split)
table(data$Sample)
max.cells.per.sample <- min(table(data$Sample)) #**** is the maximum number of cells that can be pulled from each sample evenly

scc.list <- list()
for(i in 1:length(split)){
  print(i)
  Idents(split[[i]]) <- split[[i]]$orig.ident
  split[[i]] <- subset(split[[i]],cells = WhichCells(split[[i]],downsample = max.cells.per.sample))
  scc.list[[i]] <- RunNICHES(split[[i]],
                             LR.database="fantom5",
                             species="rat",
                             assay="alra",
                             cell_types = "prelim.celltypes",
                             meta.data.to.map = names(split[[i]]@meta.data),
                             SystemToCell = T,
                             CellToCell = F,
                             CellToSystem=T) 
}
names(scc.list) <- names(split)

temp.list <- list()
for(i in 1:length(scc.list)){
  temp.list[[i]] <- scc.list[[i]]$SystemToCell # Isolate SystemToCell Signaling
}
system.to.cell.imputed <- merge(temp.list[[1]],temp.list[2:length(temp.list)])
system.to.cell.imputed
system.to.cell.imputed$Sample <- factor(system.to.cell.imputed$Sample)
system.to.cell.imputed$Condition <- factor(system.to.cell.imputed$Condition)
save(system.to.cell.imputed,file='system.to.cell.imputed.Robj')

temp.list <- list()
for(i in 1:length(scc.list)){
  temp.list[[i]] <- scc.list[[i]]$CellToSystem # Isolate CellToSystem Signaling
}
cell.to.system.imputed <- merge(temp.list[[1]],temp.list[2:length(temp.list)])
cell.to.system.imputed
cell.to.system.imputed$Sample <- factor(cell.to.system.imputed$Sample)
cell.to.system.imputed$Condition <- factor(cell.to.system.imputed$Condition)
save(cell.to.system.imputed,file='cell.to.system.imputed.Robj')

####  Run NICHES [NOT IMPUTED] #### 
# Need to split by 'Sample' first so that we are only crossing cells with other cells in the same tissue
gc()
split <- SplitObject(data,split.by='Sample')
names(split)

# Cell To Cell (CTC)
scc.list <- list()
for(i in 1:length(split)){
  print(i)
  scc.list[[i]] <- RunNICHES(split[[i]],
                             LR.database="fantom5", # we could update this ground truth if we like, but this is easy for now
                             species="mouse",
                             assay="RNA",
                             cell_types = "prelim.celltypes",
                             meta.data.to.map = names(split[[i]]@meta.data), # this makes sure all of the original metdata maps into the NICHES output(s)
                             SystemToCell = F,
                             CellToCell = T,
                             CellToSystem=F) # Need to standarize the number of cells in the System measurements
}
names(scc.list) <- names(split)

# Merge outputs
temp.list <- list()
for(i in 1:length(scc.list)){
  temp.list[[i]] <- scc.list[[i]]$CellToCell # Isolate CellToCell Signaling
  gc()
}
cell.to.cell <- merge(temp.list[[1]],temp.list[2:length(temp.list)])
cell.to.cell

# Format metadata for plotting
cell.to.cell$Sample <- factor(cell.to.cell$Sample.Sending)
cell.to.cell$Condition <- factor(cell.to.cell$Condition.Sending)

save(cell.to.cell,file='cell.to.cell.Robj')

# STC and CTS
# Standardizing cell number for each sample
split <- SplitObject(data,split.by = 'Sample') # redundant but just in case
names(split)
table(data$Sample)
max.cells.per.sample <- min(table(data$Sample)) #4255 is the maximum number of cells that can be pulled from each sample evenly

scc.list <- list()
for(i in 1:length(split)){
  print(i)
  Idents(split[[i]]) <- split[[i]]$orig.ident
  split[[i]] <- subset(split[[i]],cells = WhichCells(split[[i]],downsample = max.cells.per.sample))
  scc.list[[i]] <- RunNICHES(split[[i]],
                             LR.database="fantom5",
                             species="rat",
                             assay="RNA",
                             cell_types = "prelim.celltypes",
                             meta.data.to.map = names(split[[i]]@meta.data),
                             SystemToCell = T,
                             CellToCell = F,
                             CellToSystem=T) 
}
names(scc.list) <- names(split)

temp.list <- list()
for(i in 1:length(scc.list)){
  temp.list[[i]] <- scc.list[[i]]$SystemToCell # Isolate SystemToCell Signaling
}
system.to.cell <- merge(temp.list[[1]],temp.list[2:length(temp.list)])
system.to.cell
system.to.cell$Sample <- factor(system.to.cell$Sample)
system.to.cell$Condition <- factor(system.to.cell$Condition)
save(system.to.cell,file='system.to.cell.Robj')

temp.list <- list()
for(i in 1:length(scc.list)){
  temp.list[[i]] <- scc.list[[i]]$CellToSystem # Isolate CellToSystem Signaling
}
cell.to.system <- merge(temp.list[[1]],temp.list[2:length(temp.list)])
cell.to.system
cell.to.system$Sample <- factor(cell.to.system$Sample)
cell.to.system$Condition <- factor(cell.to.system$Condition)
save(cell.to.system,file='cell.to.system.Robj')


