#set working directory
setwd("C:/Users/beverlie/Google Drive/Proteomics Analysis on Human Monocytes")
#setwd on mac
setwd("~/Google Drive/Proteomics Analysis on Human Monocytes")

#load library
install.packages ("tidyverse")
library(tidyverse)

#install packages for heatmap visualization
install.packages("BiocManager")
library(BiocManager)
install.packages("gplots")
library(gplots)

library(RSvgDevice)

#read data into dataframe
proteomics_alldonors <- read.csv("HumanMonocyteDonor_ProteomicsAnalysis.csv", header = T)

############################generate heatmap of comparisons dd(LPS.LPSvVeh) and d(LPSvVeh) with D188

colnames(df)

#new df with only d and dd values for all donors
prot_alldonors <- proteomics_alldonors %>% select (c(2, 16:23)) %>%   
  distinct(gene, .keep_all = TRUE) %>% 
  column_to_rownames("gene")

#make simple heatmap
heatmap(as.matrix(prot_alldonors), scale="none")

#use baseR to save directly to folder
#Open an svg file
svg("Proteomics_AllDonors_Heatmap.svg")
#Generate image to put in open file
heatmap.2(as.matrix(prot_alldonors), scale="none", col = bluered(100), trace = "none", density.info="none", cexRow = 0.5, cexCol = 1.5, keysize = 0.75)
#Close file
dev.off()

############################generate heatmap of comparisons dd(LPS.LPSvVeh) and d(LPSvVeh) WITHOUT D188

colnames(proteomics_alldonors)

#new df with only d and dd values for 3 donors
prot_threedonors <- proteomics_alldonors %>% select (c(2, 16:19, 22:23)) %>%   
  distinct(gene, .keep_all = TRUE) %>% 
  column_to_rownames("gene")

#make simple heatmap
heatmap(as.matrix(prot_threedonors), scale="none")

#use baseR to save directly to folder
#Open an svg file
svg("Proteomics_ThreeDonors_Heatmap2.svg")
#Generate image to put in open file
heatmap.2(as.matrix(prot_threedonors), scale="none", col = bluered(100), trace = "none", density.info="none", cexRow = 0.5, cexCol = 1.5, keysize = 0.75)
#Close file
dev.off()






