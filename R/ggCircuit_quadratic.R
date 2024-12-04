ggCircuit_quadratic <- function(edge.aggregate,
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
                      title = NULL) {

  #### Step 1: Define node info and add coordinates for nodes centered around origin ####
  node.info <- cbind(node.aggregate,
                     netCoin::layoutCircle(data.frame(node.aggregate$node.label),
                                           deg = graph.angle))

  #### Step 2: Define edge info and add start and end coordinates for each edge ####
  edge.info <- edge.aggregate
  edge.info$x.start <- node.info[edge.info$sending.label, ]$x
  edge.info$y.start <- node.info[edge.info$sending.label, ]$y
  edge.info$x.end <- node.info[edge.info$receiving.label, ]$x
  edge.info$y.end <- node.info[edge.info$receiving.label, ]$y

  #### Step 3: Categorize edges as autocrine or paracrine ####
  edge.info$category.AP <- ifelse(paste(edge.info$x.start, edge.info$y.start) ==
                                    paste(edge.info$x.end, edge.info$y.end),
                                  "Autocrine", "Paracrine")

  #### Step 4: Custom Alpha Scaling ####
  # Define scaling function
  custom_scale <- function(values) {
    scaled <- (values - min(values)) / (max(values) - min(values)) # Normalize to [0, 1]
    scaled <- scaled^2 # SE applying quadratic scaling to emphasize strong signals
    return(scaled * (1 - 0.01) + 0.01) # Scale to [0.01, 1]
  }

  # Apply scaling to feature values
  edge.info$alpha.scaled <- custom_scale(edge.info$feature.value)

  #### Step 5: Define global plot parameters ####
  if (is.null(min.edge.value)) {
    min.edge.value <- min(edge.info$feature.value)
  }
  if (is.null(max.edge.value)) {
    max.edge.value <- max(edge.info$feature.value)
  }
  if (is.null(cols.use)) {
    cols.use <- gg_color_hue(nrow(node.info))
  }

  #### Step 6: Create the ggplot base ####
  b <- ggplot(data = node.info, aes(x = x, y = y))

  #### Step 7: Add edges and nodes ####
  if (edge.fixed.size) {
    circuit.plot <- b +
      # Paracrine edges
      geom_segment(data = edge.info[edge.info$category.AP == 'Paracrine', ],
                   aes(x = x.start, y = y.start,
                       xend = x.end, yend = y.end,
                       alpha = alpha.scaled),
                   size = edge.fixed.size,
                   arrow = grid::arrow(angle = arrow.head.angle,
                                       length = unit(arrow.head.length, "npc"),
                                       ends = 'last', type = 'open')) +
      # Autocrine edges
      geom_curve(data = edge.info[edge.info$category.AP == 'Autocrine', ],
                 aes(x = x.start, y = y.start,
                     xend = x.end + autocrine.offset,
                     yend = y.end + autocrine.offset,
                     alpha = alpha.scaled),
                 size = edge.fixed.size,
                 curvature = autocrine.arrow.curvature,
                 arrow = grid::arrow(angle = arrow.head.angle,
                                     length = unit(arrow.head.length, "npc"),
                                     ends = 'last', type = 'open'))
  } else {
    circuit.plot <- b +
      # Paracrine edges
      geom_segment(data = edge.info[edge.info$category.AP == 'Paracrine', ],
                   aes(x = x.start, y = y.start,
                       xend = x.end, yend = y.end,
                       size = feature.value / edge.scale.factor,
                       alpha = alpha.scaled),
                   arrow = grid::arrow(angle = arrow.head.angle,
                                       length = unit(arrow.head.length, "npc"),
                                       ends = 'last', type = 'open')) +
      # Autocrine edges
      geom_curve(data = edge.info[edge.info$category.AP == 'Autocrine', ],
                 aes(x = x.start, y = y.start,
                     xend = x.end + autocrine.offset,
                     yend = y.end + autocrine.offset,
                     size = feature.value / edge.scale.factor,
                     alpha = alpha.scaled),
                 curvature = autocrine.arrow.curvature,
                 arrow = grid::arrow(angle = arrow.head.angle,
                                     length = unit(arrow.head.length, "npc"),
                                     ends = 'last', type = 'open'))
  }

  #### Step 8: Finalize the plot with nodes and legend ####
  circuit.plot <- circuit.plot +
    scale_alpha_continuous(name = "Connectivity", range = c(0.01, 1)) +
    geom_point(data = node.info,
               aes(size = ifelse(system.fraction == 0, NA, system.fraction),
                   color = node.label)) +
    scale_color_manual(values = cols.use) +
    scale_size_continuous(range = c(0, 8), limits = c(0, 1), name = "System Fraction") +
    theme_classic() +
    Seurat::NoAxes() +
    xlim(-1.1, 1.1) +
    ylim(-1.1, 1.1) +
    ggtitle(title) +
    theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))

  return(circuit.plot)
}
