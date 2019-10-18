##Goal1: to generate unique gene lists from the comparison of  1xLPSvVeh & 2xLPSvVeh. 
##Goal2: input unique gene lists into NA for unsupervised networks and generate pathways
##Goal3: input the same unique gene list into SIGORA for enriched pathways dependent on the gene list
##Goal4: compare and identify common pathways and query the biological relevance for each condition

#setwd
setwd("C:/Users/beverlie/Google Drive/Proteomics Analysis on Human Monocytes/Unique gene and pathway lists/Proteomics/PropsedOmics_Proteomics")
#setwd for mac
setwd("~/Google Drive/Proteomics Analysis on Human Monocytes/Unique gene and pathway lists/Proteomics/PropsedOmics_Proteomics")

#load tidyverse library the includes dplyr
install.packages ("tidyverse")
library(tidyverse)
library(dplyr)
install.packages("VennDiagram")
library(VennDiagram)
install.packages("sigora")
library(sigora)

#read proteomic dataframes with paired ttest values included
LPS <- read_csv("Pairedttest_LPSvVeh_Proteomics.csv")
LPSLPS <- read_csv("Pairedttest_LPSLPSvVeh_Proteomics.csv")

#remove na values
filter.LPS <- na.omit(LPS)
filter.LPSLPS <- na.omit(LPSLPS)

#select columns: gene, FC and pvalues
filter.LPS <- filter.LPS %>% 
  select(c(2, 12:13)) 

filter.LPS$p.adj <- p.adjust(filter.LPS$d_pairedttest, method = "BH")


filter.LPSLPS <- filter.LPSLPS %>% select(c(2, 12:13))
filter.LPSLPS$p.adj <- p.adjust(filter.LPSLPS$dd_pairedttest, method = "BH")

#cut off paired ttest value <0.05
filtered.LPS <- filter.LPS %>% 
  filter(p.adj < 0.1)
filtered.LPSLPS <- filter.LPSLPS %>%
  filter(p.adj < 0.1)

#find unique genes for each comparison, no FC cutoff 
unique.LPS <- anti_join(filtered.LPS, filtered.LPSLPS, by = "gene" )
unique.LPSLPS <- anti_join(filtered.LPSLPS, filtered.LPS, by = "gene")

#find gene count of shared genes, no FC cutoff 
shared.genes <- merge(filtered.LPS, filtered.LPSLPS, by= "gene", all = FALSE)

#save csv file and look up pathways on NetworkAnalyst
write.csv(unique.LPS, file = "UniqueProteins_LPSvVehicle_genes.csv")
##on NA 1st order, nodes=3418, edges=6026, seeds=216
##on NA zero order, nodes=52, edges=59, seeds=52
##on minimum, nodes=515, edges=1765, seeds=216
write.csv(unique.LPSLPS, file= "UniqueProteins_LPSLPSvVehicle_genes.csv")
##on NA 1st order, nodes=10279, edges=77126, seeds=2412
##on NA zero order, nodes=2132, edges=17670, seeds=2132
##on minimum, nodes=2132, edges=17670, seeds=2132

#check number of genes
nrow(unique.LPS)
nrow(unique.LPSLPS)

#add unique DE protein lists and category names, no FC cutoff 
vennplot <- draw.pairwise.venn(area1 = 430,area2 = 2723,cross.area = 192,category = c("LPSvVehicle", "LPSLPSvVehicle"))

#save svg and manipulate on inkscape
ggsave(vennplot, file = "proteomics_vennplot.svg", device = "svg")

#########################################script from Travis to generate pathway list from SIGORA
#sigora with no FC cutoff
sigora_unique.LPS <- sigora(reaH, level = 4, queryList = unique.LPS$gene) %>%
  pluck("summary_results") 
#filter(Bonferroni <= 0.001)
##with such low gene counts, removed the significance filter to include all pathways

sigora_unique.LPSLPS <- sigora(reaH, level = 4, queryList = unique.LPSLPS$gene) %>%
  pluck("summary_results") %>%
  filter(Bonferroni <= 0.001)
##found 80 enriched LPSLPSvVehicle pathways from Sigora

#save SIGORA lists as .csv
write.csv(sigora_unique.LPS, file = "SigoraUnique_LPS.csv")
write.csv(sigora_unique.LPSLPS, file = "SigoraUnique_LPSLPS.csv")

############################################find shared pathway lists from NA and Sigora
NA.LPS <- read_csv("Pairedttest_LPSvVeh_Proteomics.csv")





