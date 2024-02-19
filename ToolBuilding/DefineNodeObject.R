DefineNodeObject <- function(transcr.obj,
                             feature,
                             assay = 'RNA'){
  # check & initialize
  if(!(feature %in% rownames(transcr.obj))){
    stop('feature not present in transcr.obj')
  }else{}
  
  # extract & store information
  node.object <- data.frame(feature.value = transcr.obj@assays[[assay]]@data[feature,])
  node.object$barcode <- rownames(node.object)
  
  # extract metadata
  node.object <- cbind(node.object,transcr.obj@meta.data)
  return(node.object)
}
