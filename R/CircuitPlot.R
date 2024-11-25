CircuitPlot <- function(transcr.obj,
                        connect.obj,
                        feature,
                        plot.function = 'ggCircuit',
                        group.by = 'CellClass',
                        graph.angle = 45,
                        h = 0.2,
                        offset = 0.05,
                        autocrine.offset = 0.03,
                        edge.scale.factor = 20,
                        arrow.head.angle = 15,
                        arrow.head.length = 0.03,
                        autocrine.arrow.curvature = 8,
                        cols.use = RColorBrewer::brewer.pal(4,'Set2'),
                        edge.fixed.size = 1,
                        global.node.list = c('Endothelium','Epithelium','Mesenchyme','Immune'),
                        min.edge.value = 0,
                        max.edge.value = NULL,
                        ...){
  # DefineObjects
  ligand <- strsplit(feature,split = "—")[[1]][1]
  receptor <- strsplit(feature,split = "—")[[1]][2]
  node.object <- DefineNodeObject(transcr.obj = transcr.obj,
                                  feature = ligand)
  edge.object <- DefineEdgeObject(connect.obj = connect.obj,
                                  feature = feature)
  
  # AggregateObjects
  group.by.edge <- paste(group.by,'Joint',sep='.')
  #options(warn = -1)
  node.aggregate <- AggregateNodeData(node.object = node.object,
                                      group.by = group.by,
                                      global.node.list = global.node.list)
  edge.aggregate <- AggregateEdgeData(edge.object = edge.object,
                                      group.by = group.by.edge)
  #options(warn = 0)
  
  # Plot with desired plot.function
  if(plot.function == 'ggCircuit'){
    circuit.plot <- ggCircuit(node.aggregate = node.aggregate,
                              edge.aggregate = edge.aggregate,
                              graph.angle = graph.angle,
                              h = h,
                              offset = offset,
                              autocrine.offset = autocrine.offset,
                              edge.scale.factor = edge.scale.factor,
                              arrow.head.angle = arrow.head.angle,
                              arrow.head.length = arrow.head.length,
                              autocrine.arrow.curvature = autocrine.arrow.curvature,
                              cols.use = cols.use,
                              edge.fixed.size = edge.fixed.size,
                              min.edge.value = min.edge.value,
                              max.edge.value = max.edge.value)
  }
  
  circuit.plot
  return(circuit.plot)
}
