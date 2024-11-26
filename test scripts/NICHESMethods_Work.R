# Create NICHESMethods 2024-11-25
# MSBR
require(ggplot2)
require(Seurat)
library(devtools)
packageVersion("devtools")
create_package("/Users/msbr/GitHub/NICHESMethods")

getwd()
load_all()
check()

# sophie's data
load("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Personal_Folders/Sophie/Single Cell/Single Cell Samples - Nuoya Collaboration/BAMC/BAMC.CTC.integrated.Robj")
load("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Personal_Folders/Sophie/Single Cell/Single Cell Samples - Nuoya Collaboration/BAMC/BAMC.integrated_NAs.Robj")

#inspect
color_pal <- c('#A40606','#9CFFFA','#B0DB43','#9C528B','#2F6690',
               '#946846','#F1C40F','green','#0F0326','#E65F5C','#14591D','#726DA8',
               'yellow','purple','blue','red','orange','darkgrey','magenta')
color_pal <- color_pal[1:length(unique(node.object$New_Annotations))]

names(color_pal) <- unique(node.object$New_Annotations)
Seurat::DimPlot(BAMC.CTC.integrated,label = T)
Seurat::DimPlot(BAMC.CTC.integrated,group.by = 'SendingType',cols = color_pal,shuffle=T)
Seurat::DimPlot(BAMC.CTC.integrated,group.by = 'ReceivingType',cols = color_pal,shuffle = T)
FeaturePlot(BAMC.CTC.integrated,'Igf1—Igf1r')
FeaturePlot(BAMC.CTC.integrated,'Areg—Egfr')

mark <- FindAllMarkers(BAMC.CTC.integrated,min.pct = 0.5)
mark$ratio <- mark$pct.1/mark$pct.2

# test
node.object <- DefineNodeObject(transcr.obj = BAMC.integrated,
                 feature = 'Areg')

edge.object <- DefineEdgeObject(connect.obj = BAMC.CTC.integrated,
                 feature = 'Areg—Egfr',
                 assay = 'CellToCell')

node.aggregate <- AggregateNodeData(node.object = node.object,
                    group.by = 'New_Annotations',
                    global.node.list = unique(node.object$New_Annotations))

edge.aggregate <- AggregateEdgeData(edge.object = edge.object,
                                    group.by = 'New_Annotations.Joint')
# test ggCircuit

global.node.list = unique(node.object$New_Annotations)
ggCircuit(node.aggregate = node.aggregate,
          edge.aggregate = edge.aggregate,
          offset = 0.05, # distance between the opposing arrows
          edge.fixed.size = T, # true means edge.scale.factor does nothing
          edge.scale.factor = 1, # bigger means thinner arrows
          graph.angle = 25, # this should be a function of number of nodes
          h = 0.05, # hypotenuse length for trig that spaces the arrows from the circles
          autocrine.offset = 0.01, # size of autocrine cricle diameter
          arrow.head.angle = 10, # fatness of arrows
          arrow.head.length = 0.02, # obvious, size of arrow tip threads
          autocrine.arrow.curvature = 20,# size of autocrine cricle diameter ???
          cols.use = color_pal)

BAMC.CTC.integrated$poor.mans.entropy <- BAMC.CTC.integrated$nFeature_CellToCell/BAMC.CTC.integrated$nCount_CellToCell
FeaturePlot(BAMC.CTC.integrated,features = 'poor.mans.entropy')
