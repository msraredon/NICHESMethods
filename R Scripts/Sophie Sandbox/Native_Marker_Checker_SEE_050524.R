
# Canonical Gene Marker Test - SEE - 05.05.2024

# This function allows you to run through multiple sets of known native (canonical) markers
# For me, SE, this is helpful when I am trying to refine my cluster annotations
# Our engineered systems are complex and so comparing against native data is always helpful.

# Set working directory to be true
setwd("/Users/sophieedelstein/Desktop/SE Single Cell/Organoids_BASC_Mono/Native Lung Atlas")
# Load native object
load("lung.combined.clean.classed.annotated.final.2022-07-24.Robj")

# Now grabbing my engineered sample object (organoids)
setwd("/Users/sophieedelstein/Desktop/SE Single Cell/Organoids_BASC_Mono")
load("sub_basc_3_obj_SEE_05072024_CC.Robj")

# Checking presence of known (canonical) markers by cell type (for loop to run through)
# Define gene groups - you can change these based on your own marker classes of interest
ATI <- c("Pdpn", "Ager", "Col4a4", "Col4a3", "Akap5", "Clic5", "Sema3e", "Pla2g1b", "Aqp5", "Hopx")
ATII <- c("Sftpc", "Napsa", "Lyz2", "Defb4", "Scgb1a1", "S100g", "Abca3", "Sftpd", "Sftpb")
Basal <- c("Tp63", "Krt5", "Krt14")
AlvMac <- c("Mrc1", "Arg1", "Slc39a2", "Prodh2")
InstMac <- c("Cd163", "C1qb", "C1qa", "Trem2", "Pf4", "Clec10a")
BASC <- c("Nrgn", "Plac8", "Gng13")
Tuft <- c("Dclk1", "Rgs13", "Trpm5", "Avil")
EMT <- c("Twist1", "Twist2", "Fn1", "Tagln", "Vim", "Cdh2")
DevEpi <- c("Sox9", "Sox2", "Nkx2-1")
Secretory <- c("Scgb1a1", "Rarres1","Enpp3","Wfdc2","Chia","")
Ciliated <- c("Tppp3", "Sntn", "Riiad1","Wipf3","Dynlrb2")
Cycling <- c("Top2a", "Ube2c", "Mki67","Nusap1")

Native_Mark_Checker = function(gene_group, folder_name, object) {
  # Set the working directory to the specified folder (CHANGE AS NECESSARY)
  setwd("/Users/sophieedelstein/Desktop/SE Single Cell/Organoids_BASC_Mono")
  for (gene in gene_group) {
    # Create violin plot
    vln_plot <- VlnPlot(object, gene)
    # Create feature plot
    feature_plot <- FeaturePlot(object, gene)
    
    # Save plots as PDF files
    vln_filename <- paste(folder_name, "_VlnPlot_", gene, ".pdf", sep = "")
    pdf(file = vln_filename, width = 8, height = 6)  # Adjust width and height as needed
    print(vln_plot)
    dev.off()
    
    feature_filename <- paste(folder_name, "_FeaturePlot_", gene, ".pdf", sep = "")
    pdf(file = feature_filename, width = 8, height = 6)  # Adjust width and height as needed
    print(feature_plot)
    dev.off()
  }
  
}
# Call the function for a specific gene group and save plots to a specific folder
Native_Mark_Checker(ATI, "ATI", sub_basc_3_obj) # change the object based on what your sample is
Native_Mark_Checker(ATII, "ATII", sub_basc_3_obj)
Native_Mark_Checker(Basal, "Basal", sub_basc_3_obj)
Native_Mark_Checker(AlvMac, "AlvMac", sub_basc_3_obj)
Native_Mark_Checker(InstMac, "InstMac", sub_basc_3_obj)
Native_Mark_Checker(BASC, "BASC", sub_basc_3_obj)
Native_Mark_Checker(Tuft, "Tuft", sub_basc_3_obj)
Native_Mark_Checker(EMT, "EMT", sub_basc_3_obj)
Native_Mark_Checker(DevEpi, "DevEpi", sub_basc_3_obj)
Native_Mark_Checker(Cycling, "Cycling", sub_basc_3_obj)
Native_Mark_Checker(Secretory, "Secretory", sub_basc_3_obj)
Native_Mark_Checker(Ciliated, "Ciliated", sub_basc_3_obj)


################################# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####################################
