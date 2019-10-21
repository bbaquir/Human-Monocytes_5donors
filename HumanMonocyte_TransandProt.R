setwd("C:/Users/beverlie/Desktop/Human-Monocytes_5donors")

#install packages
install.packages("VennDiagram")
library(VennDiagram)
install.packages("tidyverse")
library(tidyverse)
library(dplyr)

#upload LPS.LPS v Vehicle transcriptomic list
#select only Genes, log2FC and padj values
#convert ensembl ids to gene ids by comparing to biomart_table
RNAseq <- read_csv("DEGenes_LPSLPSvsVehicle_10-2_14-06-18.csv")
RNAseq <- RNAseq %>% 
  select(c(1, 3, 7)) 
biomart_table <- rename(biomart_table, Gene=ensembl_gene_id)

RNAseq_genes <- inner_join(RNAseq, biomart_table, by = "Gene") 
RNAseq_genes <- RNAseq_genes %>% dplyr::select (c(3, 7, 10))
RNAseq_genes <- subset(RNAseq_genes, select=c(3, 1, 2))

#upload LPS.LPS v Vehicle proteomic list
Proteomics <- read_csv("DEProteins_LPSLPSvsVehicle.csv")
  
  
  #paired ttest: t.test(x, y, paired = TRUE, alternative = "two.sided")