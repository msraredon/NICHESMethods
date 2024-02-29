ggCircuit <- function(edge.aggregate,
                        node.aggregate,
                        graph.angle,
                        h,
                        offset,
                        autocrine.offset,
                        edge.scale.factor,
                        arrow.head.angle,
                        arrow.head.length,
                        autocrine.arrow.curvature,
                        cols.use,
                        edge.fixed.size,
                        split.by = FALSE,
                        min.edge.value = NULL,
                        max.edge.value = NULL,
                        fixed.node.list = NULL#c('Epithelium','Endothelium','Mesenchyme','Immune')
                        ){
  
  # Define node info and add coordinates for nodes centered around origin
  if(!is.null(fixed.node.list)){
    node.info <- cbind(node.aggregate,
                       netCoin::layoutCircle(data.frame(fixed.node.list),
                                             deg=graph.angle))
  }else{
  node.info <- cbind(node.aggregate,
                 netCoin::layoutCircle(data.frame(node.aggregate$node.label),
                                       deg=graph.angle))
  }
  # Define edge info and add start and end coordinates for each edge
  edge.info <- edge.aggregate
  edge.info$x.start <- node.info[edge.info$sending.label,]$x
  edge.info$y.start <- node.info[edge.info$sending.label,]$y
  edge.info$x.end <- node.info[edge.info$receiving.label,]$x
  edge.info$y.end <- node.info[edge.info$receiving.label,]$y
  
  # Define autocrine and paracrine
  edge.info$category.AP <- NA
  if(length(which(paste(edge.info$x.start,edge.info$y.start)==paste(edge.info$x.end,edge.info$y.end)))>0){
  edge.info[which(paste(edge.info$x.start,edge.info$y.start)==paste(edge.info$x.end,edge.info$y.end)),]$category.AP <- 'Autocrine'
  }
  if(length(which(paste(edge.info$x.start,edge.info$y.start)!=paste(edge.info$x.end,edge.info$y.end)))>0){
  edge.info[which(paste(edge.info$x.start,edge.info$y.start)!=paste(edge.info$x.end,edge.info$y.end)),]$category.AP <- 'Paracrine'
  }
  # Get set up for geometry
  edge.info$dx <- edge.info$x.end-edge.info$x.start
  edge.info$dy <- edge.info$y.end-edge.info$y.start
  edge.info$angle.in.rads <- atan(edge.info$dy/edge.info$dx)
  
  # Define edge angle 'quadrant' based on rotation in space
  edge.info$quadrant <- NA
  if(length(edge.info[which(edge.info$dx>0 & edge.info$dy>0 ),]$quadrant)>0){
    edge.info[which(edge.info$dx>0 & edge.info$dy>0 ),]$quadrant <- 'quad.1'
  } 
  if(length(edge.info[which(edge.info$dx<0 & edge.info$dy>0 ),]$quadrant)>0){
    edge.info[which(edge.info$dx<0 & edge.info$dy>0 ),]$quadrant <- 'quad.2'
  }
  if(length(edge.info[which(edge.info$dx<0 & edge.info$dy<0 ),]$quadrant)>0){
    edge.info[which(edge.info$dx<0 & edge.info$dy<0 ),]$quadrant <- 'quad.3'
  }
  if(length(edge.info[which(edge.info$dx>0 & edge.info$dy<0 ),]$quadrant)>0){
    edge.info[which(edge.info$dx>0 & edge.info$dy<0 ),]$quadrant <- 'quad.4'
  }
  if(length(edge.info[which(edge.info$dx>0 & edge.info$dy==0 ),]$quadrant)>0){
    edge.info[which(edge.info$dx>0 & edge.info$dy==0 ),]$quadrant <- 'quad.1.4'
  }
  if(length(edge.info[which(edge.info$dx<0 & edge.info$dy==0 ),]$quadrant)>0){
    edge.info[which(edge.info$dx<0 & edge.info$dy==0 ),]$quadrant <- 'quad.2.3'
  }
  if(length(edge.info[which(edge.info$dx==0 & edge.info$dy<0 ),]$quadrant)>0){
    edge.info[which(edge.info$dx==0 & edge.info$dy<0 ),]$quadrant <- 'quad.3.4'
  }
  if(length(edge.info[which(edge.info$dx==0 & edge.info$dy>0 ),]$quadrant)>0){
    edge.info[which(edge.info$dx==0 & edge.info$dy>0 ),]$quadrant <- 'quad.1.2'
  }
  
  ## Adjust paracrine edge length
  # Define hypotenuse length
  h = h
  # Define delta (sign determiend below)
  edge.info$x.delta <- h*cos(edge.info$angle.in.rads)
  edge.info$y.delta <- h*sin(edge.info$angle.in.rads)
  
  # Initialize adjusted paracrine edge length
  edge.info$x.start.adj <- NA
  edge.info$y.start.adj <- NA
  edge.info$x.end.adj <- NA
  edge.info$y.end.adj <- NA
  
  # Adjust paracrine edge length
  for(i in 1:nrow(edge.info)){
    if(!is.na(edge.info[i,]$quadrant) & edge.info[i,]$quadrant=='quad.1'){ 
      edge.info[i,]$x.start.adj <- edge.info[i,]$x.start + abs(edge.info[i,]$x.delta)
      edge.info[i,]$y.start.adj <- edge.info[i,]$y.start + abs(edge.info[i,]$y.delta)
      edge.info[i,]$x.end.adj <- edge.info[i,]$x.end - abs(edge.info[i,]$x.delta)
      edge.info[i,]$y.end.adj <- edge.info[i,]$y.end - abs(edge.info[i,]$y.delta)
    }
    if(!is.na(edge.info[i,]$quadrant) & edge.info[i,]$quadrant =='quad.2'){ 
      edge.info[i,]$x.start.adj <- edge.info[i,]$x.start - abs(edge.info[i,]$x.delta)
      edge.info[i,]$y.start.adj <- edge.info[i,]$y.start + abs(edge.info[i,]$y.delta)
      edge.info[i,]$x.end.adj <- edge.info[i,]$x.end + abs(edge.info[i,]$x.delta)
      edge.info[i,]$y.end.adj <- edge.info[i,]$y.end - abs(edge.info[i,]$y.delta)
    }
    if(!is.na(edge.info[i,]$quadrant) & edge.info[i,]$quadrant =='quad.3'){ 
      edge.info[i,]$x.start.adj <- edge.info[i,]$x.start - abs(edge.info[i,]$x.delta)
      edge.info[i,]$y.start.adj <- edge.info[i,]$y.start - abs(edge.info[i,]$y.delta)
      edge.info[i,]$x.end.adj <- edge.info[i,]$x.end + abs(edge.info[i,]$x.delta)
      edge.info[i,]$y.end.adj <- edge.info[i,]$y.end + abs(edge.info[i,]$y.delta)
    }
    if(!is.na(edge.info[i,]$quadrant) & edge.info[i,]$quadrant =='quad.4'){ 
      edge.info[i,]$x.start.adj <- edge.info[i,]$x.start + abs(edge.info[i,]$x.delta)
      edge.info[i,]$y.start.adj <- edge.info[i,]$y.start - abs(edge.info[i,]$y.delta)
      edge.info[i,]$x.end.adj <- edge.info[i,]$x.end - abs(edge.info[i,]$x.delta)
      edge.info[i,]$y.end.adj <- edge.info[i,]$y.end + abs(edge.info[i,]$y.delta)
    }
    if(!is.na(edge.info[i,]$quadrant) & edge.info[i,]$quadrant=='quad.1.2'){
      edge.info[i,]$x.start.adj <- edge.info[i,]$x.start + 0
      edge.info[i,]$y.start.adj <- edge.info[i,]$y.start + abs(edge.info[i,]$y.delta)
      edge.info[i,]$x.end.adj <- edge.info[i,]$x.end - 0
      edge.info[i,]$y.end.adj <- edge.info[i,]$y.end - abs(edge.info[i,]$y.delta)
    }
    if(!is.na(edge.info[i,]$quadrant) & edge.info[i,]$quadrant=='quad.2.3'){
      edge.info[i,]$x.start.adj <- edge.info[i,]$x.start - abs(edge.info[i,]$x.delta)
      edge.info[i,]$y.start.adj <- edge.info[i,]$y.start - 0
      edge.info[i,]$x.end.adj <- edge.info[i,]$x.end + abs(edge.info[i,]$x.delta)
      edge.info[i,]$y.end.adj <- edge.info[i,]$y.end + 0
    }
    if(!is.na(edge.info[i,]$quadrant) & edge.info[i,]$quadrant=='quad.3.4'){
      edge.info[i,]$x.start.adj <- edge.info[i,]$x.start + 0
      edge.info[i,]$y.start.adj <- edge.info[i,]$y.start - abs(edge.info[i,]$y.delta)
      edge.info[i,]$x.end.adj <- edge.info[i,]$x.end - 0
      edge.info[i,]$y.end.adj <- edge.info[i,]$y.end + abs(edge.info[i,]$y.delta)
    }
    if(!is.na(edge.info[i,]$quadrant) & edge.info[i,]$quadrant=='quad.1.4'){
      edge.info[i,]$x.start.adj <- edge.info[i,]$x.start + abs(edge.info[i,]$x.delta)
      edge.info[i,]$y.start.adj <- edge.info[i,]$y.start - 0
      edge.info[i,]$x.end.adj <- edge.info[i,]$x.end - abs(edge.info[i,]$x.delta)
      edge.info[i,]$y.end.adj <- edge.info[i,]$y.end + 0
    }

  }
  
  ## Offset the paracrine edges
  # Define offset length
  offset = offset
  # Define delta (sign determiend below)
  edge.info$x.offset <- offset*cos(edge.info$angle.in.rads)
  edge.info$y.offset <- offset*sin(edge.info$angle.in.rads)
  
  # Initialize offset paracrine edge coords
  edge.info$x.start.offset <- NA
  edge.info$y.start.offset <- NA
  edge.info$x.end.offset <- NA
  edge.info$y.end.offset <- NA
  
  # Offset
  for(i in 1:nrow(edge.info)){
    if(!is.na(edge.info[i,]$quadrant) & edge.info[i,]$quadrant=='quad.1'){ 
      edge.info[i,]$x.start.offset <- edge.info[i,]$x.start.adj + abs(edge.info[i,]$x.offset)
      edge.info[i,]$y.start.offset <- edge.info[i,]$y.start.adj - abs(edge.info[i,]$y.offset)
      edge.info[i,]$x.end.offset <- edge.info[i,]$x.end.adj + abs(edge.info[i,]$x.offset)
      edge.info[i,]$y.end.offset <- edge.info[i,]$y.end.adj - abs(edge.info[i,]$y.offset)
    }
    if(!is.na(edge.info[i,]$quadrant) & edge.info[i,]$quadrant =='quad.2'){ 
      edge.info[i,]$x.start.offset <- edge.info[i,]$x.start.adj + abs(edge.info[i,]$x.offset)
      edge.info[i,]$y.start.offset <- edge.info[i,]$y.start.adj + abs(edge.info[i,]$y.offset)
      edge.info[i,]$x.end.offset <- edge.info[i,]$x.end.adj + abs(edge.info[i,]$x.offset)
      edge.info[i,]$y.end.offset <- edge.info[i,]$y.end.adj + abs(edge.info[i,]$y.offset)
    }
    if(!is.na(edge.info[i,]$quadrant) & edge.info[i,]$quadrant =='quad.3'){ 
      edge.info[i,]$x.start.offset <- edge.info[i,]$x.start.adj -  abs(edge.info[i,]$x.offset)
      edge.info[i,]$y.start.offset <- edge.info[i,]$y.start.adj + abs(edge.info[i,]$y.offset)
      edge.info[i,]$x.end.offset <- edge.info[i,]$x.end.adj - abs(edge.info[i,]$x.offset)
      edge.info[i,]$y.end.offset <- edge.info[i,]$y.end.adj + abs(edge.info[i,]$y.offset)
    }
    if(!is.na(edge.info[i,]$quadrant) & edge.info[i,]$quadrant =='quad.4'){ 
      edge.info[i,]$x.start.offset <- edge.info[i,]$x.start.adj - abs(edge.info[i,]$x.offset)
      edge.info[i,]$y.start.offset <- edge.info[i,]$y.start.adj - abs(edge.info[i,]$y.offset)
      edge.info[i,]$x.end.offset <- edge.info[i,]$x.end.adj - abs(edge.info[i,]$x.offset)
      edge.info[i,]$y.end.offset <- edge.info[i,]$y.end.adj - abs(edge.info[i,]$y.offset)
    }
    if(!is.na(edge.info[i,]$quadrant) & edge.info[i,]$quadrant=='quad.1.2'){
      edge.info[i,]$x.start.offset <- edge.info[i,]$x.start.adj + offset
      edge.info[i,]$y.start.offset <- edge.info[i,]$y.start.adj + 0
      edge.info[i,]$x.end.offset <- edge.info[i,]$x.end.adj + offset
      edge.info[i,]$y.end.offset <- edge.info[i,]$y.end.adj + 0
    }
    if(!is.na(edge.info[i,]$quadrant) & edge.info[i,]$quadrant=='quad.2.3'){
      edge.info[i,]$x.start.offset <- edge.info[i,]$x.start.adj + 0
      edge.info[i,]$y.start.offset <- edge.info[i,]$y.start.adj + offset
      edge.info[i,]$x.end.offset <- edge.info[i,]$x.end.adj - 0
      edge.info[i,]$y.end.offset <- edge.info[i,]$y.end.adj + offset
    }
    if(!is.na(edge.info[i,]$quadrant) & edge.info[i,]$quadrant=='quad.3.4'){
      edge.info[i,]$x.start.offset <- edge.info[i,]$x.start.adj - offset
      edge.info[i,]$y.start.offset <- edge.info[i,]$y.start.adj + 0
      edge.info[i,]$x.end.offset <- edge.info[i,]$x.end.adj - offset
      edge.info[i,]$y.end.offset <- edge.info[i,]$y.end.adj + 0
    }
    if(!is.na(edge.info[i,]$quadrant) & edge.info[i,]$quadrant=='quad.1.4'){
      edge.info[i,]$x.start.offset <- edge.info[i,]$x.start.adj + 0
      edge.info[i,]$y.start.offset <- edge.info[i,]$y.start.adj - offset
      edge.info[i,]$x.end.offset <- edge.info[i,]$x.end.adj + 0
      edge.info[i,]$y.end.offset <- edge.info[i,]$y.end.adj - offset
    }
    
  }
  
  # Define global plot parameters, if not input by user
  if(is.null(min.edge.value)){
  min.edge.value <- min(edge.info$feature.value)
  }else{
    min.edge.value <- min.edge.value
  }
  if(is.null(max.edge.value)){
  max.edge.value <- max(edge.info$feature.value)
  }else{
    max.edge.value <- max.edge.value
  }
  
  # Plot
  b <- ggplot(data = node.info, 
              aes(x = x,y = y))
  
  # PLOT FIXED EDGES
  if(edge.fixed.size){
  circuit.plot <- b + 
    
    geom_segment(data = edge.info[edge.info$category.AP=='Paracrine',],
                 aes(x = x.start.offset, 
                     y = y.start.offset, 
                     xend = x.end.offset, 
                     yend = y.end.offset,
                     # size = feature.value/edge.scale.factor
                     alpha = feature.value),
                 size = edge.fixed.size,
                 arrow = grid::arrow(angle = arrow.head.angle, # how 'fat' the arrowheads are
                                     length = unit(arrow.head.length, "npc"), # size of the arrowheads, including 0
                                     ends = 'last',
                                     type = 'open'))+
    
      geom_curve(data = edge.info[edge.info$category.AP=='Autocrine',],
                 aes(x = x.start, 
                    y = y.start, 
                    xend = x.end+autocrine.offset, 
                    yend = y.end+autocrine.offset,
                    alpha = feature.value),
                 size = edge.fixed.size,
               curvature = autocrine.arrow.curvature,
               ncp = 10,
                 arrow = grid::arrow(angle = arrow.head.angle, # how 'fat' the arrowheads are
                                     length = unit(arrow.head.length, "npc"), # size of the arrowheads, including 0
                                     ends = 'last',
                                     type = 'open'))
  }else{
    # PLOT UNFIXED EDGES
    circuit.plot <- b + 
      
      #scale_size_continuous(limits=c(min.edge.value,max.edge.value))+
      
      geom_segment(data = edge.info[edge.info$category.AP=='Paracrine',],
                   aes(x = x.start.offset, 
                       y = y.start.offset, 
                       xend = x.end.offset, 
                       yend = y.end.offset,
                       size = feature.value/edge.scale.factor,
                       alpha = feature.value),
                   #size = edge.fixed.size,
                   arrow = grid::arrow(angle = arrow.head.angle, # how 'fat' the arrowheads are
                                       length = unit(arrow.head.length, "npc"), # size of the arrowheads, including 0
                                       ends = 'last',
                                       type = 'open'))+
      
      geom_curve(data = edge.info[edge.info$category.AP=='Autocrine',],
                 aes(x = x.start, 
                     y = y.start, 
                     xend = x.end+autocrine.offset, 
                     yend = y.end+autocrine.offset,
                     size = feature.value/edge.scale.factor,
                     alpha = feature.value),
                 #size = edge.fixed.size,
                 curvature = autocrine.arrow.curvature,
                 ncp = 10,
                 arrow = grid::arrow(angle = arrow.head.angle, # how 'fat' the arrowheads are
                                     length = unit(arrow.head.length, "npc"), # size of the arrowheads, including 0
                                     ends = 'last',
                                     type = 'open'))
  }
  
  # Scale edges (allows global scaling with other plots)
  circuit.plot <- circuit.plot +
                  scale_alpha_continuous(range=c(0.01,1),
                                          limits=c(min.edge.value,max.edge.value),
                                            name='Connectivity',trans = 'sqrt')
  # Add nodes on top
  circuit.plot <- circuit.plot +
      geom_point(data = node.info,
             aes(size = ifelse(system.fraction==0, NA, system.fraction),
                 color = node.label))+
      #geom_text(aes(label = node.label),hjust=1, vjust=1,size=2)+
      scale_color_manual(values = cols.use)+
      scale_size_continuous(range = c(0,6),limits = c(0,1),name = 'System Fraction',trans = 'sqrt')+
    
      theme_classic()+Seurat::NoAxes()+Seurat::NoLegend() +xlim(-1,1)+ylim(-1,1)

  circuit.plot
  
  # return output
  return(circuit.plot)
}
