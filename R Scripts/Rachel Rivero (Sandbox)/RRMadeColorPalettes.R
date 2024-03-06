# Color Palettes

# All cells E17
col.pal <- list()
col.pal$CellClass <- c('#3C5488FF', '#4DBBD5FF', '#00A087FF', '#E64B35FF') #Nature color palette
names(col.pal$CellClass) <- c('Mesenchymal','Epithelial', 'Immune', 'Endothelial')
col.pal$CellType <- c('cornflowerblue', 'red','orange',
                                  '#B0DB43','mediumseagreen', 'slategray1','deepskyblue',
                                  '#14591D','yellow', 'purple','black', 'darkgoldenrod',
                                  'slateblue1','thistle', 'blue', 'deeppink',
                                  'navy', '#AF5D63',
                                  'pink','cyan', '#A40606')
                                  names(col.pal$CellType) <- c('Sox9+ Mesenchymal Progenitors', 'Sox9- Mesenchymal Progenitors','Prex2+ Mes',
                                                               'Proliferating Wnt2+ Mes', 'Wnt2+ Mes', 'SMCs/Myofibroblasts', 'Il17+ SMCs/Myofibroblasts',
                                                               'Mesothelium', 'Pericytes 1','Pericytes 2', 'NCCs', 'NECs', #mes
                                                               'Early AT2', 'Transitional Epithelium', 'Club cells', 'Non-ciliated Secretory', #epi
                                                               'APCs','Immune', #immune
                                                               'Arterial','Transitional Endothelial', 'LEC') #endo
col.pal$Condition <- c('#F1C40F', '#2F6690')
names(col.pal$Condition) <- c('CDH','Normal')
col.pal$Sample <- c('#BAB700','#EDE3E4','#FF5E5B','#00CECB')
names(col.pal$Sample) <- c('frl1', 'frl2', 'cdhfrl1', 'cdhfrl2')

save(col.pal, file = 'E17colorpalette.Robj')                                  

# Mes
col.pal.mes <- list()
col.pal.mes$CellClass <- c('#3C5488FF')
col.pal.mes$CellType <- c('cornflowerblue','red','mediumseagreen','#B0DB43','orange', 
                                          'slategray1','deepskyblue', 'yellow','purple',
                                          '#14591D','black','darkgoldenrod')
                                          col.pal.mes$Condition <- c('#F1C40F', '#2F6690')
                                          col.pal.mes$Sample <- c('#BAB700','#EDE3E4','#FF5E5B','#00CECB')
                                          
                                          names(col.pal.mes$CellClass) <- c('Mesenchymal') 
                                          names(col.pal.mes$CellType) <- c('Sox9+ Mesenchymal Progenitors', 'Sox9- Mesenchymal Progenitors',
                                                                           'Wnt2+ Mes', 'Proliferating Wnt2+ Mes','Prex2+ Mes', 
                                                                           'SMCs/Myofibroblasts','Il17+ SMCs/Myofibroblasts',  
                                                                           'Pericytes 1','Pericytes 2', 
                                                                           'Mesothelium','NCCs', 'NECs')
                                          names(col.pal.mes$Condition) <- c('CDH','Normal') 
                                          names(col.pal.mes$Sample) <- c('cdhfrl1','cdhfrl2','frl1','frl2')
                                          
                                          save(col.pal.mes, file = 'mescolorpalette.Robj')