# Rule Set (Creates Feature Metdata)

# RuleSetFunction <- function(query.set = NICHES::ncomms8866$Pair.Name, # A vector of mechanisms, consisting of a sending part and receiving part, separated by a special character
#                               query.set.name = 'FANTOM5',
#                               special.separating.character = '_')
#   
    # Initialize new metadata
    feature.metadata <- data.frame(QUERY.SET = query.set.name,
                                   DEFINITION.LEVEL = NA,
                              LIGAND.CATEGORY = NA,
                              RECEPTOR.CATEGORY = NA,
                              LIGAND.FAMILY = NA,
                              RECEPTOR.FAMILY = NA,
                              MECHANISM = query.set,
                              LIGAND = NA,
                              RECEPTOR = NA)

    # Rule 1: MECHANISM is the feature label itself
    ## Performed above, upon feature.metadata initialization
    
    ## here is a change

    # Rule 2: LIGAND and RECEPTOR definition via special separating character
    feature.metadata$LIGAND <- data.frame(stringr::str_split_fixed(feature.metadata$MECHANISM, 
                                               pattern = special.separating.character,
                                               n = 2))[,1]
    feature.metadata$RECEPTOR <- data.frame(stringr::str_split_fixed(feature.metadata$MECHANISM, 
                                                                   pattern = special.separating.character,
                                                                   n = 2))[,2]
  
    ###### LIGAND.FAMILY  ######
    
    ## Amyloid Beta Precursor
    indices.temp <- grep(pattern = '^APP',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Amyloid Beta Precursor'
    
    ## Cingulin
    indices.temp <- grep(pattern = '^CGN',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Cingulin'
    
    ## Chorionic gonadotropin
    indices.temp <- grep(pattern = '^CGA|^CGB',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Chorionic gonadotropin'
    
    ## Annexin
    indices.temp <- grep(pattern = '^ANXA',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Annexin'
    
    ## Lipoprotein
    indices.temp <- grep(pattern = '^LPL',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Lipoprotein'
    ## EGF
    indices.temp <- grep(pattern = '^EGF',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'EGF'
    ## IGF
    indices.temp <- grep(pattern = '^IGF',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'IGF'
    ## Vitronectin
    indices.temp <- grep(pattern = '^VTN',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Vitronectin'
    ## Complement
    indices.temp <- grep(pattern = '^C1Q|^C3|^C4',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Complement'
    
    ## HBEGF
    indices.temp <- grep(pattern = '^HBEGF',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'HBEGF'
    ## HGF
    indices.temp <- grep(pattern = '^HGF',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'HGF'
    
    ## Antigen Presentation
    indices.temp <- grep(pattern = '^HLA-|^B2M',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Antigen Presentation'
    
    ## Extracellular Enzyme
    indices.temp <- grep(pattern = '^PLAU|^MMP|^ADAM',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Extracellular Enzyme'
    
    ## Interferon
    indices.temp <- grep(pattern = '^IFN',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Interferon'
    
    ## Calmodulin
    indices.temp <- grep(pattern = '^CALM',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Calmodulin'
    
    ### Collagen / Laminin / Fibronectin ...
    
    ## Decorin
    indices.temp <- grep(pattern = '^BGN',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'BGN'
    
    ## Decorin
    indices.temp <- grep(pattern = '^DCN',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Decorin'
    
    ## Collagens
    indices.temp <- grep(pattern = '^COL',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Collagen'
    ## Laminins
    indices.temp <- grep(pattern = '^LAM',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Laminin'
    ## Fibronectin
    indices.temp <- grep(pattern = '^FN1',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Fibronectin'
    
    ### FGF / BMP / TGFB / VEGFA / PGF / NGF / WNT...
    
    ## WNT
    indices.temp <- grep(pattern = '^WNT|^RSPO',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'WNT'
    ## BMP
    indices.temp <- grep(pattern = 'BMP',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'BMP'
    ## FGF
    indices.temp <- grep(pattern = 'FGF',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'FGF'
    ## TGFB
    indices.temp <- grep(pattern = 'TGFB',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'TGFB'
    
    ## CTGF
    indices.temp <- grep(pattern = '^CTGF',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'CTGF'
    
    ## VEGF
    indices.temp <- grep(pattern = 'VEGF|FIGF', # FIGF = VEGFD!
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'VEGF'
    ## BDNF
    indices.temp <- grep(pattern = 'BDNF',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'BDNF'
    
    ## ANGPT
    indices.temp <- grep(pattern = '^ANGPT',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'ANGPT'
    ## GDF
    indices.temp <- grep(pattern = 'GDF',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'GDF'
    ## NGF
    indices.temp <- grep(pattern = '^NGF|^NTF', #NTF is part of NGF family https://en.wikipedia.org/wiki/Neurotrophin-3',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'NGF'
    ## PDGF
    indices.temp <- grep(pattern = '^PDGF',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'PDGF'
    
    ## Plasmin
    indices.temp <- grep(pattern = '^PLG',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Plasmin'
    
    ## PGF
    indices.temp <- grep(pattern = '^PGF',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'PGF'
    ## BDNF
    indices.temp <- grep(pattern = '^BDNF',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'BDNF'
    
    ## Versican
    indices.temp <- grep(pattern = '^VCAN',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Versican'
    
    ## TNC
    indices.temp <- grep(pattern = '^TNC',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'TNC'
    
    ## CSF
    indices.temp <- grep(pattern = '^CSF',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'CSF'
    
    ### Spatial
    
    ## Semaphorins
    indices.temp <- grep(pattern = '^SEMA',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Semaphorin'
    
     
    
    ## Neurexophilin
    indices.temp <- grep(pattern = '^NXP',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Neurexophilin'
    ## Agouti-signaling protein
    indices.temp <- grep(pattern = '^ASIP',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Agouti-signaling protein'
    
    ## SLIT
    indices.temp <- grep(pattern = 'SLIT',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'SLIT'
    
    ## Neuropeptide Y
    indices.temp <- grep(pattern = 'NPY',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Neuropeptide Y'
    
    ### CC ...
    ## TNF
    indices.temp <- grep(pattern = '^TNF|^LTA|^LTB', #LTA = TNF-B, LTB = TNF-C
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'TNF'
    ## CCL
    indices.temp <- grep(pattern = 'CCL',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'CCL'
    ## CXC
    indices.temp <- grep(pattern = '^CXC|^CX3',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'CXC'
    
    ## GDNF
    indices.temp <- grep(pattern = '^GDNF',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'GDNF'
    
    ## Neuroligin
    indices.temp <- grep(pattern = '^NLGN',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Neuroligin'
    
    ## IL
    indices.temp <- grep(pattern = '^IL',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Interleukin'
    
    ## Urocortin
    indices.temp <- grep(pattern = '^UCN',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Urocortin'
    
    ## Endothelin
    indices.temp <- grep(pattern = '^EDN',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Endothelin'
    
    ## Adrenomedullin
    indices.temp <- grep(pattern = '^ADM',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Adrenomedullin'
    
    ## Agouti-related peptide
    indices.temp <- grep(pattern = '^AGRP',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Agouti-related Peptide'
    
    ## Agrin
    indices.temp <- grep(pattern = '^AGRN',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Agrin'
    
    ## NEGF
    indices.temp <- grep(pattern = '^PTN|^MDK',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'NEGF'
    
    ## CNTF
    indices.temp <- grep(pattern = '^CNTF',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'CNTF'
    
    ## Relaxin
    indices.temp <- grep(pattern = '^RLN',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Relaxin'
    
    ## Hedgehog
    indices.temp <- grep(pattern = '^DHH|^SHH|^IHH',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Hedgehog'
    
    ## Cadherin
    indices.temp <- grep(pattern = '^CDH',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Cadherin'
    
    ## Vimentin
    indices.temp <- grep(pattern = '^VIM',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Vimentin'
    
    ## Tachykinin
    indices.temp <- grep(pattern = '^TAC',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Tachykinin'
    ## Contactin
    indices.temp <- grep(pattern = '^CNTN',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Contactin'
    
    ## Cell Adhesion
    indices.temp <- grep(pattern = '^VCAM|^NCAM|^ICAM',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Cell Adhesion Molecule'
    
    
    ## Apolipoprotein 
    indices.temp <- grep(pattern = '^APO',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Apolipoprotein'
    
    ## Heat-shock Protein 
    indices.temp <- grep(pattern = '^HSP',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Heat-shock Protein '
    
    ## VIP 
    indices.temp <- grep(pattern = '^VIP',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'VIP'
    
    ## Kininogen 
    indices.temp <- grep(pattern = '^KNG',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Kininogen'
    
    
    ## Dickkopf 
    indices.temp <- grep(pattern = '^DKK',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Dickkopf (WNT-inhibitor)'
    
    ## Inhibin 
    indices.temp <- grep(pattern = '^INH',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Inhibin'
    
    ## Oncostatin M 
    indices.temp <- grep(pattern = '^OSM',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Oncostatin M'
    
    ## Cortistatin 
    indices.temp <- grep(pattern = '^CORT',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Cortistatin'
    
    ## Somatostatin 
    indices.temp <- grep(pattern = '^SST',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Somatostatin'
    
    ## Reticulon
    indices.temp <- grep(pattern = '^RTN',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Reticulon'
    
    ## Lactoferrin
    indices.temp <- grep(pattern = '^LTF',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Lactoferrin'
    
    ## Prosaposin
    indices.temp <- grep(pattern = '^PSAP',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Prosaposin'
    
    ## G Protein
    indices.temp <- grep(pattern = '^GNA',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'G Protein Complex'
    
    ## Tryptophan hydroxylase 
    indices.temp <- grep(pattern = '^TPH1',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Tryptophan hydroxylase'
    
    ## Secretin
    indices.temp <- grep(pattern = '^SCT',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Secretin'
    
    ## Netrin 
    indices.temp <- grep(pattern = '^NTN',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Netrin'
    
    ## Serpin 
    indices.temp <- grep(pattern = '^SERP',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Serpin'
    
    ## Thrombospondin 
    indices.temp <- grep(pattern = '^THBS',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Thrombospondin'
    
    ## VWF 
    indices.temp <- grep(pattern = '^VWF',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'VWF'
    
    ## Ephrin 
    indices.temp <- grep(pattern = '^EFN',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Ephrin'
    
    ## CYR61 
    indices.temp <- grep(pattern = '^CYR61',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'CYR61'
    
    ## Renin-Angiotensin
    indices.temp <- grep(pattern = '^ACE|^AGT|^REN',
                         x = feature.metadata$LIGAND,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Renin-Angiotensin'
    
    
    ###### RECEPTOR.FAMILY ######
    ## CXC Family Receptor
    indices.temp <- grep(pattern = '^CXC',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'CXC Family Receptor'
    
    ## CC Family Receptor
    indices.temp <- grep(pattern = '^CCR',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'CC Family Receptor'
    

    ## Neurexin
    indices.temp <- grep(pattern = '^NRX',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Neurexin'
    
    
    ## EGF Family Receptor
    indices.temp <- grep(pattern = '^ERBB|EGFR',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'EGF Family Receptor'
    
    ## Integrins
    indices.temp <- grep(pattern = 'ITG',
                    x = feature.metadata$RECEPTOR,
                    ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Integrin'
    
    ## NOTCH
    indices.temp <- grep(pattern = 'NOTCH',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'NOTCH Family Receptor'
    
    ## TLR
    indices.temp <- grep(pattern = 'TLR',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Toll-Like Receptor'
    
    ## Lipid Handling
    indices.temp <- grep(pattern = 'LDL',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Lipid Handling'
    
    
    ## Reticulon Family Receptor
    indices.temp <- grep(pattern = 'RTN',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Reticulon Family Receptor'
    
    ## Ephrin Family
    indices.temp <- grep(pattern = '^EPH',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Ephrin Family Receptor'
    
    ## WNT Family
    indices.temp <- grep(pattern = '^LGR|FZD|RYK',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'WNT Family Receptor'
    
    ## BMP Family
    indices.temp <- grep(pattern = '^BMPR',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'BMP Family Receptor'
    
    ## Collagen
    indices.temp <- grep(pattern = '^COL',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Collagen'
    
    ## FGF Family
    indices.temp <- grep(pattern = '^FGFR',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'FGF Family Receptor'
    ## NGF Family
    indices.temp <- grep(pattern = '^NGF|^NTRK',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'NGF Family Receptor'
    
    ## GDNF Family Receptor
    indices.temp <- grep(pattern = '^GFR', #GFRA = GDNF receptor https://en.wikipedia.org/wiki/GFRA1
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'GDNF Family Receptor'
    
    ## EGF Family
    indices.temp <- grep(pattern = '^EGF|^EREG|^AREG',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'EGF Family Receptor'
    ## PDGF Family
    indices.temp <- grep(pattern = '^PDGF',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'PDGF Family Receptor'
    ## TGFBR
    indices.temp <- grep(pattern = '^TGFBR',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'TGFB Family Receptor'
    
    ## ANGPT
    indices.temp <- grep(pattern = '^TIE|^TEK',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'ANGPT Family Receptor'
    ## Plexin
    indices.temp <- grep(pattern = '^PLXN',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Plexin'
    
    ## Selectin
    indices.temp <- grep(pattern = '^SEL',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Selectin'
    
    ## Activin
    indices.temp <- grep(pattern = '^ACVR',
                         x = feature.metadata$RECEPTOR,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Activin'
    
    ###### LIGAND.CATEGORY    ##### 
    ## Matrix / Growth Factor / Cytokines / Metabolic / Direct Contact / Spatial Guidance / Bioelectrical
    
    ## Protease Inhibition
    indices.temp <- grep(pattern = 'Serpin',
                         x = feature.metadata$LIGAND.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.CATEGORY <- 'Protease Inhibition'
    
    ## Kinin-Kallikrein System
    indices.temp <- grep(pattern = 'Kininogen',
                         x = feature.metadata$LIGAND.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.CATEGORY <- 'Kinin-Kallikrein System'
    
    ## Cell-Cell Contact
    indices.temp <- grep(pattern = 'Cadherin|Cell Adhesion|Neuroligin|Contactin|Neurexophilin',
                         x = feature.metadata$LIGAND.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.CATEGORY <- 'Cell-Cell Contact'
    
    ## Spatial Guidance
    indices.temp <- grep(pattern = 'Semaphorin|SLIT|Netrin',
                         x = feature.metadata$LIGAND.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.CATEGORY <- 'Spatial Guidance'
    
    ## Metabolic
    indices.temp <- grep(pattern = 'Apolipoprotein|Lipoprotein|Tryptophan hydroxylase',
                         x = feature.metadata$LIGAND.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.CATEGORY <- 'Metabolic'
    
    ## Calcium
    indices.temp <- grep(pattern = 'Calmodulin',
                         x = feature.metadata$LIGAND.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.CATEGORY <- 'Calcium'
    
    ## Complement
    indices.temp <- grep(pattern = 'Complement|VWF',
                         x = feature.metadata$LIGAND.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.CATEGORY <- 'Complement'
    
    ## Matrix Type
    indices.temp <- grep(pattern = 'Collagen|Laminin|Fibronectin|Vitronectin|BGN|Decorin|TNC|Agrin|CTGF|CYR61|Versican',
                         x = feature.metadata$LIGAND.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.CATEGORY <- 'Matrix'
    
    ## Paracrine Type
    indices.temp <- grep(pattern = 'Adrenomedullin|Endothelin|Agouti-related|Renin-Angiotensin|VIP|Urocortin|Relaxin',
                         x = feature.metadata$LIGAND.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.CATEGORY <- 'Paracrine'
    
    ## Antigen Presentation Type
    indices.temp <- grep(pattern = 'Antigen Presentation',
                         x = feature.metadata$LIGAND.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.CATEGORY <- 'Antigen Presentation'
   
    ## Enzyme Type
    indices.temp <- grep(pattern = 'Extracellular Enzyme|Plasmin',
                         x = feature.metadata$LIGAND.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.CATEGORY <- 'Extracellular Enzyme'
    
    ## GF Type
    indices.temp <- grep(pattern = 'WNT|CSF|PDGF|VEGF|NGF|PGF|BMP|TGFB|FGF|IGF|EGF|GDF|ANGPT|BDNF|HGF|GDNF|NEGF|CNTF',
                         x = feature.metadata$LIGAND.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.CATEGORY <- 'Growth Factor'
    
    ## Cytokine Category
    indices.temp <- grep(pattern = 'Interleukin|Interferon|CCL|CXC|TNF|Oncostatin M',
                         x = feature.metadata$LIGAND.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$LIGAND.CATEGORY <- 'Cytokine'
    
    ###### RECEPTOR.CATEGORY    ##### 
    
    ## Plexin
    indices.temp <- grep(pattern = 'Plexin',
                         x = feature.metadata$RECEPTOR.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.CATEGORY <- 'Plexin'
    
    ## Activin
    indices.temp <- grep(pattern = 'Activin',
                         x = feature.metadata$RECEPTOR.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.CATEGORY <- 'Activin'
    
    ## Toll-Like Interactions
    indices.temp <- grep(pattern = 'Toll-Like Receptor',
                         x = feature.metadata$RECEPTOR.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.CATEGORY <- 'Toll-Like Receptor'
    
    ## Ephrin Interactions
    indices.temp <- grep(pattern = 'Ephrin',
                         x = feature.metadata$RECEPTOR.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.CATEGORY <- 'Ephrin Family Receptor'
    
    ## Integrin Interactions
    indices.temp <- grep(pattern = 'Integrin',
                         x = feature.metadata$RECEPTOR.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.CATEGORY <- 'Integrin'
    
    ## Growth Factor Receptor 
    indices.temp <- grep(pattern = 'WNT|CSF|PDGF|VEGF|NGF|PGF|BMP|TGFB|FGF|IGF|EGF|GDF|ANGPT|BDNF|HGF|GDNF|NEGF|CNTF',
                         x = feature.metadata$RECEPTOR.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.CATEGORY <- 'Growth Factor Receptor'
    
    ## Cytokine Receptor 
    indices.temp <- grep(pattern = '^CXC Family|^CC Family',
                         x = feature.metadata$RECEPTOR.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.CATEGORY <- 'Cytokine Receptor'
    ## Metabolic
    indices.temp <- grep(pattern = 'Lipid',
                         x = feature.metadata$RECEPTOR.FAMILY,
                         ignore.case = FALSE)
    feature.metadata[indices.temp,]$RECEPTOR.CATEGORY <- 'Metabolic'
    ## temp
    temp <- feature.metadata[,c('LIGAND.CATEGORY',
                                'RECEPTOR.CATEGORY',
                                "LIGAND.FAMILY",
                                "RECEPTOR.FAMILY")]
    temp <- !is.na(temp)
    
    feature.metadata$DEFINITION.LEVEL <- as.character(rowSums(temp))

    
    
    
    