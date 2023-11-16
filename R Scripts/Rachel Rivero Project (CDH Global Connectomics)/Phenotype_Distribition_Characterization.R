# Miscellany Phenotype Investigations

# Set wd
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Rachel/E17.E19.Merged.Explorations.2023-10-09")

# Load packages
require(Seurat)
require(ggplot2)
require(cowplot)
require(dplyr)
require(ggthemes)
require(ggdark)

# Load previously-saved, already-annotated data
load('frl.combined.annotated.classed.2023-10-09.Robj')


# Look at Cell Type Distribution vs. Time (Number of Cells)
dat <- data.frame(table(frl.combined$prelim.celltypes,frl.combined$Timepoint))
png('Cell Type Distribution vs. Time (Number of Cells).png',width = 7,height = 5,res = 300,units = 'in')
ggplot(dat,
       aes(x=Var1,y=Freq,fill=Var2))+
  geom_bar(stat='identity',
           position = 'dodge')+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ggtitle('Cell Type Distribution vs. Time')+
  xlab('Cell Type')+
  ylab('Number of Cells')+ 
  scale_fill_discrete(name = "Timepoint")
dev.off()

# Look at Cell Type Distribution vs. Time (Percentage of Tissue)
dat <- table(frl.combined$prelim.celltypes,frl.combined$Timepoint)
# Normalize to total number of cells in each timepoint
cells.per.timepoint <- colSums(dat)
dat[,1] <- dat[,1]/cells.per.timepoint[1]
dat[,2] <- dat[,2]/cells.per.timepoint[2]

# Check math
colSums(dat) # each column should equal 1

# Conver to data frame for ggplot
dat <- data.frame(dat)

# Convert fractions to percentages
dat$Freq <- dat$Freq*100

# make plot
png('Cell Type Distribution vs. Time (Percentage of Tissue).png',width = 7,height = 5,res = 300,units = 'in')
ggplot(dat,
       aes(x=Var1,y=Freq,fill=Var2))+
  geom_bar(stat='identity',
           position = 'dodge')+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ggtitle('Cell Type Distribution vs. Time')+
  xlab('Cell Type')+
  ylab('Percentage of Tissue')+ 
  scale_fill_discrete(name = "Timepoint")
dev.off()

# Now let's do where we split by Sample, to confirm that this pattern is present across multiple tissues

# Look at Cell Type Distribution vs. Time (Percentage of Tissue) BY SAMPLE
dat <- table(frl.combined$prelim.celltypes,frl.combined$Sample)
# Normalize to total number of cells in each timepoint
cells.per.sample <- colSums(dat)
for(i in 1:ncol(dat)){
dat[,i] <- dat[,i]/cells.per.sample[i]
}
# Check math
colSums(dat) # each column should equal 1

# Conver to data frame for ggplot
dat <- data.frame(dat)

# Convert fractions to percentages
dat$Freq <- dat$Freq*100

# Organize metadata for plotting
dat$Var2 <- factor(dat$Var2,levels = c('frl1','frl2','cdhfrl1','cdhfrl2','NRL19','CDH19'))
# make plot
png('Cell Type Distribution vs. Time (Percentage of Tissue) BY SAMPLE.png',width = 10,height = 5,res = 300,units = 'in')
ggplot(dat,
       aes(x=Var1,y=Freq,fill=Var2))+
  geom_bar(stat='identity',
           position = 'dodge')+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ggtitle('Cell Type Distribution vs. Time (By Sample)')+
  xlab('Cell Type')+
  ylab('Percentage of Tissue')+ 
  scale_fill_manual(values = c('#63D2FF','#7D80DA','#E6AF2E','#F24236','#6A8D73','#000022'),name = "Sample")
dev.off()

