Receptor Categorization Reference Markdown
Version: 1.0.0
Author: MSBR
Date Last Modified: 2024-12-31

# FANTOM5 Receptors

%%%%%%%%%%

## Cytokine Receptors

### Chemokine Family Receptors
#### CCL Family Receptor
##### '^CCR'
#### CXC Family Receptor
##### '^CXC'
#### C Chemokine (XCL) Family Receptors
##### '^XCR'
#### Colony Stimulating Factor (CSF) Family Receptors
##### '^CSF'

### Interferon Family Receptors
#### IFN Family Receptors
##### '^INF'

### Interleukin Family Receptors
#### CLCF Family
##### '^CRLF'
#### IL Family Receptor Subunis
##### '^IL'
#### Oncostatin M
##### '^OSMR'
#### Leukemia Inhibitory Factor
##### '^LIFR'

### Tumor Necrosis Factor Family Receptors
#### TNF Family
##### '^LTBR$'
#### TNFSF Family Receptor Subunit
##### '^TNFR|^CD27$'

### Lymphoid Interaction Receptors
#### HLA Complex Molecules
##### '^LILRB|^KIR|^CD1A$|^CD1B$|^CD74$'
#### CD40LG
##### '^CD40'
#### BTLA
##### '^VTCN|^CD79A'

### Pattern Recognition Receptors
#### Toll-Like Receptors
##### 'TLR'

### Miscellaneous (Receptors)
#### "Apoptosis inhibitor expressed by macrophages"
##### '^CD5$'
#### Calprotectin Family Receptor
##### NA
#### Cathelicidin Family Receptor
##### '^FPR'
#### Lactotransferrin Family Receptor
##### NA
#### Annexin Family Receptor
##### '^DYSF'

%%%%%%%%%%

## Growth Factor Receptors

### Canonical Growth Factor Receptors
#### Epidermal Growth Factors (EGF) Family Receptor
##### '^ERBB|EGFR'
#### Fibroblast Growth Factors (FGF) Family Receptor
##### '^FGFR'
#### Heparin-binding EGF-like Growth factor (HBEGF Family)
##### '^CD9$'
#### Hepatocyte Growth Factor (HGF Family)
##### '^MET$'
#### Insulin-like Growth Factor (IGF) Family Receptor
##### '^IGF'
#### Neurite Promoting Growth Factor (NEGF Family)
##### '^PTPR'
#### Nerve Growth Factor (NGF) Family Receptor
##### '^NGF|^NTRK'
#### Neuregulins (NRG Family)
##### NA
#### Platelet Derived Growth Factor (PDGF) Family Receptor
##### '^PDGF'
#### Placental Growth Factor (PGF Family)
##### 
#### Transforming Growth Factor Beta (TGFB) Family Receptor
##### '^TGFBR'
#### Vascular Endothelial Growth Factor (VEGF Family)
##### 
#### Angiopoietin (ANGPT) Family Receptor
##### '^TIE|^TEK'
#### Brain-derived Neurotrophic Factor (BDNF Family)
#####
#### Bone Morphogenetic Protein (BMP) Family Receptor
##### '^BMPR'
#### Apelin
##### '^APLNR$'
#### Ciliary Neurotrophic factor
##### '^CNTFR'
#### Growth Differentiation Factors (GDF Family)
##### 
#### Stem Cell Factor (KITLG) Receptor
##### 
#### Glial cell line-derived neurotrophic factors (GDNF) Family Receptor
##### '^RET|^GFR'

### Canonical Morphogen Receptors
#### NOTCH Family Receptors
##### 'NOTCH'
#### WNT Family Receptors
##### '^LGR|FZD|RYK'
#### Dickkopf Family (WNT-inhibitor) Family Receptors
##### '^KREMEN'
#### Hedgehog Family Receptors
##### '^PTCH|^HH|^CDON'
#### Inhibin Family Receptors
##### '^ACVR|^BAMBI'

%%%%%%%%%%

## Endocrine Molecule Receptors
  ## Parathyroid Hormone Receptor
  indices.temp <- grep(pattern = '^PTH',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Parathyroid Hormone Receptor'
  ## Prolactin Family Receptor
  indices.temp <- grep(pattern = '^PRL',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Prolactin Family Receptor'

  ## Glucagon-like peptide receptor
  indices.temp <- grep(pattern = '^GLP',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Glucagon-like peptide receptor'
### Calcium Homeostasis
#### Calcitonin Family Receptor
##### '^CALCR'
#### Calreticulin
#####

### Metabolic Homeostasis
#### Growth Hormone Family Receptor
##### '^GHR'
#### Ghrelin
##### 
#### Pancreatic Polypeptide
##### 
#### "Insulin-like"
##### 
#### Somatostatin
##### 
#### Agouti-signaling protein
##### 

### Lipid Homeostasis
#### Adiponectin Family Receptor
##### '^ADIPOR'

#### Lipoprotein Family Receptor
##### 'LDL'

### Cardiovascular & Osmotic Homeostasis
#### Adrenomedullin Family Receptor
##### '^RAMP'
#### Endothelin Family Receptor
##### '^EDNR'
#### Natriuretic Family Receptor
##### '^NPR'
#### Relaxin Family Receptor
##### '^RXFP'
#### Renin-Angiotensin System Receptor
##### '^AGTR'
#### Secretin Family Receptor
##### '^SCTR$'
#### Vasopressin Family Receptor
##### '^AVPR'
#### VIP Family Receptor
##### '^VIPR'

### Reproductive Homeostasis
#### Oxytocin Receptor
##### '^OXTR'
#### Sex Hormone Binding Protein
##### 
#### Human chorionic gonadotropin (hCG)
##### 

### Stress Homeostasis
#### Urocortin
##### 
#### CRH
##### 

%%%%%%%%%%

## Matrix Molecule Receptors


### Heparan sulfate proteoglycans (HSPGs)
#### Glypicans
##### '^GPC'
#### Syndecans
##### '^SYN'
### Integrins
#### Integrins
##### '^ITG'

### Matricellular Family Receptors
#### Agrin (AGRN)
##### '^MUSK$'
#### Biglycan (BGN)
##### NA
#### Connective Tissue Growth Factor (CTGF)
##### NA
#### CTHRC1
##### NA
#### CYR61
##### NA
#### Decorins
##### NA
#### ECM1
##### NA
#### Fibulins
##### NA
#### Fibrillins
##### NA
#### Lecticans
##### NA
#### Microfibrillar-associated proteins
##### NA
#### Myocilin
##### NA
#### Nephronectins
##### NA
#### Nidogens 
##### NA
#### Osteopontins
##### NA
#### Prosaposins
##### '^GPR37|^CELSR1$'
#### Reelin
##### NA
#### Tenascins
##### 
#### Thrombospondins
##### '^CD36$|^CD47$'
#### Vitronectins
##### NA
#### ZP3
##### NA

### Protease Family Receptors
#### Urokinase (PLAU)
##### '^PLAUR$'
#### Metalloprotease Family
##### NA
#### ADAM Family
##### NA
#### Tissue factor pathway inhibitor (TFPI) 
##### NA
#### Transglutaminase Family
##### '^GPR56$'
#### Presenilin Family
##### '^NCSTN'
#### Granzymes
##### NA
#### Plasmin Family Receptor
##### '^PLGRKT$'

### Protease Inhibition
#### Serpins
##### 
#### Globulins
##### 
#### Tissue inhibitors of metalloproteinases (TIMPs)
##### 
### Circulatory Factors
#### Kinin-Kallikrein System
##### 
#### Complement Molecules
##### 
#### Transferrins
##### 
#### Coagulation Factors
##### 
### Miscellaneous/Uncategorized
#### LYPD3
##### 
#### Secretory Proteins
##### 

  
%%%%%%%%%%

## Spatial-Guidance Molecule Receptors

### Cell-Adhesion Molecules
#### Amyloid Beta Precursor
##### NA
#### Cadherin Family Receptors
##### '^CDH'
#### CAM Family Receptor Subunit
##### 'CAM'
#### Contactins
##### '^CNTN'
#### Neurexins
##### '^NRXN'
#### Lectin Family Receptors
##### '^SEL|^CD93'

### ECM-like Molecules
#### Netrin Family Receptors
##### '^UNC5'
#### SLIT Family Receptors
##### '^ROBO'

### Growth & Migration Inhibitors
#### Repulsive guidance molecules (RGM Family)
##### '^NEO'
#### Reticulon (Nogo) Family Receptors
##### '^RTN'

### Context Signaling Molecules
#### Ephrin Family Receptors
##### ^EPH'
#### Semaphorin Family Receptors
##### ^PLXN'

%%%%%%%%%%

## Neuropeptide Receptors

### Central Melanocortin System
#### Agouti-Related Peptide
##### 
#### NPY
##### 
#### POMC
##### 
### Neuropeptide, Unsorted
  ## Neuropeptide Receptor
  indices.temp <- grep(pattern = '^NPFFR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Neuropeptide Receptor'
#### Pro-melanin stimulating hormone (PMCH)
##### 
#### Orexin (HCRT)
##### 
#### Neurotensin (NTS)
##### 
#### Neuromedins
##### 
#### NPW
##### 
#### NPB 
##### 
#### CNTF
##### 
#### Cortistatin
##### 
#### Tachykinin
##### 

## Pleiotrophic Receptors
### Transmembrane Glycoproteins
#### Neuropilins
##### '^NRP'

%%%%

## Non-Receptors (Cytosolic Molecules)

### Intracellular (Structural)
#### Afadin
##### 
#### Vimentin
##### 
#### Cingulin
##### 
#### Scaffold Proteins
##### 
#### Keratin-associated
##### 
### Intracellular (Signaling)
#### G Protein Complex Molecules
##### 
#### LDL-receptor-related protein-associated protein (LRPAP)
##### 
#### Phosphatidylinositol-glycan biosynthesis class F protein (PIGF)
##### 
#### Heat-shock Proteins
#####
#### GTPase Signaling
##### 
### Intracellular (Enzymes)
#### Tryptophan Hydroxylase
##### 
#### Enzyme, Unsorted
##### 
### Intracellular (Uncategorized)
#### Intracellular (Uncategorized)
##### 

%%%%%%%%%%%%%%%%%%


DefineReceptorFamily <- function(feature.metadata){
  ###### RECEPTOR.FAMILY ######


  ## Asialoglycoprotein receptor
  indices.temp <- grep(pattern = '^ASGR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Asialoglycoprotein receptor'


  ## Prokineticin family receptor
  indices.temp <- grep(pattern = '^PROK',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Prokineticin family receptor'
  ## Guanylate cyclase
  indices.temp <- grep(pattern = '^GUCY',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Guanylate cyclase'

  ## TAM Family Kinase
  indices.temp <- grep(pattern = '^AXL$|^MERTK$|^TYRO3$',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'TAM Family Kinase'
  ## Formyl peptide receptor
  indices.temp <- grep(pattern = '^FPR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Formyl peptide receptor'
  ## RAMP
  indices.temp <- grep(pattern = '^RAMP',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'RAMP'


  ## ABC Transporter
  indices.temp <- grep(pattern = '^ABC',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'ABC Transporter'

  ## Neuropeptide Receptor
  indices.temp <- grep(pattern = '^NPY|^GALR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'NPY'
  
  ## Solute Carrier Family
  indices.temp <- grep(pattern = '^SLC',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Solute Carrier Family'


  ## Opiod Receptor Subunit
  indices.temp <- grep(pattern = '^OPR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Opiod Receptor Subunit'


  ## G-Protein Receptor
  indices.temp <- grep(pattern = '^GPR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'G-Protein Receptor'


  ## LRP Family
  indices.temp <- grep(pattern = '^LRP',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'LRP Family'
  ## Insulin Receptor
  indices.temp <- grep(pattern = '^INSR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Insulin Receptor'



  ## Syndecan
  indices.temp <- grep(pattern = '^SDC',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Syndecan'




  ## Collagen
  indices.temp <- grep(pattern = '^COL',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Collagen'



  ## Activin
  indices.temp <- grep(pattern = '^ACVR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Activin'

