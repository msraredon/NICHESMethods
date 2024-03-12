# Signaling Mechanism Categorical Tree - Draft

## Measure
# Transcriptomics

## Focus
# Signaling Mechanisms

## Range
# Queried Signaling Mechanisms

## Type
# Matrix / Growth Factors / Cytokines / Metabolic / Direct Cell-Cell Contact / Mechanics / Bioelectrical Field

## Family
# Collagen / Laminin / Fibronectin
# FGF / BMP / TGFB / VEGFA / PGF / NGF
# CC

## Mechanism

DefineCategories <- function(query.set, # A vector of mechanisms, consisting of a sending part and receiving part
                             rule.set.function = rule.set.function, # A function containing categorization rules
                             check.for.consistency = TRUE, # Whether or not to check if rule.set is self-consistent
){
  
  if(check.for.consistency == TRUE){
    CheckConsistency(query.set = query.set,
                     rule.set = rule.set)
  }
  if(check.for.completeness == TRUE){
    CheckCompleteness(query.set = query.set,
                      rule.set = rule.set)
  }
}

CheckConsistency <- function(query.set,
                             rule.set){
  
}

CheckCompleteness <- function(query.set,
                             rule.set){
  
}

