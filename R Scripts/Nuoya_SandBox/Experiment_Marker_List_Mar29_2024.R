# Set WD
setwd("/Volumes/T7/Chris")

# Load data
load("/Volumes/T7/Chris/all_merged.NICHES.Robj")

# Load marker list
load("/Volumes/T7/Chris/me_Marker_List.Robj")

# require packages
require(Seurat)
require(ggplot2)
require(cowplot)

# Inspect object
p1 <- DimPlot(NICHES.integrated)
p2 <- DimPlot(NICHES.integrated,group.by = 'Sample')
p3 <- SpatialDimPlot(NICHES.integrated,ncol = 4)
p.right <- plot_grid(p1,p2,ncol = 1)
p.total <- plot_grid(p.right,p3,rel_widths = c(1,1.5))
p.total

# Set Default assay
DefaultAssay(NICHES.integrated) <- 'integrated'

# Define some Global parameters

Sample_Name_List <- names(me_Marker_List) # the micro-environment marker list, names of each sample
sample_num <- length(Sample_Name_List) # the number of samples

# This is the function to find self-defined marker list.
# least_sample_num: at least this number of samples contain the marker
# cluster_ID: a character variable defining which cluster we are looking into

Experiment_Marker_List <- function(least_sample_num,cluster_ID){
  
  Sample_Comb <- t(combn(Sample_Name_List,least_sample_num)) # create all possible sample combinations
  comb_num <- length(Sample_Comb)/least_sample_num # get number of combinations

  cur_list <- c() # initializing......
  
  # For loop to iterate all combinations to create a list
  for(i in 1:comb_num){
    first_Sample_Name <- Sample_Comb[i,1]
    first_list_flag <- me_Marker_List[[first_Sample_Name]]['cluster']==cluster_ID&me_Marker_List[[first_Sample_Name]]['p_val_adj']<0.001
    first_Sample <- me_Marker_List[[first_Sample_Name]][first_list_flag,]
    select_Sample_Gene <- first_Sample$gene # first list for intersection initializing...
    for(j in 2:least_sample_num){
      cur_Sample_Name <- Sample_Comb[i,j]
      list_flag <- me_Marker_List[[cur_Sample_Name]]['cluster']==cluster_ID&me_Marker_List[[cur_Sample_Name]]['p_val_adj']<0.001
      cur_Sample <- me_Marker_List[[cur_Sample_Name]][list_flag,]
      cur_Sample_Gene <- cur_Sample$gene
      select_Sample_Gene <- intersect(select_Sample_Gene,cur_Sample_Gene)
    }
    cur_list <- union(cur_list,select_Sample_Gene)  # use union to exclude repetition
  }
  return(cur_list)
}

Experiment_Marker_List(least_sample_num=6,cluster_ID='4')


View(me_Marker_List$EWS4550[list2,])

table(Idents(NICHES.integrated),NICHES.integrated$Sample)
p1 <- SpatialDimPlot(NICHES.integrated,group.by = 'me.clusters',ncol=4)+NoLegend()
p2 <- SpatialFeaturePlot(NICHES.integrated,'PGF—NRP2',ncol=4)
plot_grid(p1,p2,ncol = 2)
