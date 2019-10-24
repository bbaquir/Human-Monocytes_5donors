#do pathway enrichment with proteomic datasets akin to RNA-Seq
#create pathway lists for all treatment comparisons

setwd("C:/Users/beverlie/Desktop/Human-Monocytes_5donors/HumanMonocyte_Proteomics/")

install.packages("tidyverse")
install.packages("sigora")

##############################################LPS VS VEHICLE##########################################
LPSvVehProteomics <- read_csv("DEProteins_LPSvsVehicle.csv")

#filter pvalues <0.05 & <0.01
LPSvVehProt_0.05 <- LPSvVehProteomics %>% filter (pvalue_ttest1 <0.05) %>% select (c(2, 6:7))
LPSvVehProt_0.01 <- LPSvVehProteomics %>% filter (pvalue_ttest1 <0.01) %>% select (c(2, 6:7))

#sigora pathways from pvalues <0.05 & <0.01  proteins
LPSvVehProt_0.05 <- sigora(reaH, level = 4, queryList =LPSvVehProt_0.05$Gene) %>% 
  pluck("summary_results") %>% 
  filter(Bonferroni <=0.001)
LPSvVehProt_0.01 <- sigora(reaH, level = 4, queryList =LPSvVehProt_0.01$Gene) %>% 
  pluck("summary_results") %>% 
  filter(Bonferroni <=0.001)

#save pathway lists
write.csv(LPSvVehProt_0.05, file= "LPSvsVehicle_0.05_SigoraPathways.csv")
write.csv(LPSvVehProt_0.01, file= "LPSvsVehicle_0.01_SigoraPathways.csv")

#######################################LPSLPS VS VEHICLE#########################

LPSLPSvVehProteomics <- read_csv("DEProteins_LPSLPSvsVehicle.csv")

#filter pvalues <0.05 & <0.01
LPSLPSvVehProt_0.05 <- LPSLPSvVehProteomics %>% filter (pvalue_ttest1 <0.05) %>% select (c(2, 6:7))
LPSLPSvVehProt_0.01 <- LPSLPSvVehProteomics %>% filter (pvalue_ttest1 <0.01) %>% select (c(2, 6:7))

#sigora pathways from pvalues <0.05 & <0.01  proteins
LPSLPSvVehProt_0.05 <- sigora(reaH, level = 4, queryList = LPSLPSvVehProt_0.05$Gene) %>% 
  pluck("summary_results") %>% 
  filter(Bonferroni <=0.001)
LPSLPSvVehProt_0.01 <- sigora(reaH, level = 4, queryList = LPSLPSvVehProt_0.01$Gene) %>% 
  pluck("summary_results") %>% 
  filter(Bonferroni <=0.001)

#save pathway lists
write.csv(LPSLPSvVehProt_0.05, file= "LPSLPSvsVehicle_0.05_SigoraPathways.csv")
write.csv(LPSLPSvVehProt_0.01, file= "LPSLPSvsVehicle_0.01_SigoraPathways.csv")

###################################LPSLPS VS LPS########################################
LPSLPSvLPSProteomics <- read_csv("DEProteins_LPSLPSvsLPS.csv")

#filter pvalues <0.05 & <0.01
LPSLPSvLPSProt_0.05 <- LPSLPSvLPSProteomics %>% filter (pvalue_ttest1 <0.05) %>% select (c(2, 14:15))
LPSLPSvLPSProt_0.01 <- LPSLPSvLPSProteomics %>% filter (pvalue_ttest1 <0.01) %>% select (c(2, 14:15))

#sigora pathways from pvalues <0.05 & <0.01  proteins
LPSLPSvLPSProt_0.05 <- sigora(reaH, level = 4, queryList = LPSLPSvLPSProt_0.05$Gene) %>% 
  pluck("summary_results") %>% 
  filter(Bonferroni <=0.001)
LPSLPSvLPSProt_0.01 <- sigora(reaH, level = 4, queryList = LPSLPSvLPSProt_0.01$Gene) %>% 
  pluck("summary_results") %>% 
  filter(Bonferroni <=0.001)

#save pathway lists
write.csv(LPSLPSvLPSProt_0.05, file= "LPSLPSvsLPS_0.05_SigoraPathways.csv")
write.csv(LPSLPSvLPSProt_0.01, file= "LPSLPSvsLPS_0.01_SigoraPathways.csv")
