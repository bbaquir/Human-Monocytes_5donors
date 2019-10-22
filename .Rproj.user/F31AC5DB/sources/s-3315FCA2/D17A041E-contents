setwd("C:/Users/beverlie/Desktop/Human-Monocytes_5donors")

#install packages
install.packages("VennDiagram")
library(VennDiagram)
install.packages("tidyverse")
library(tidyverse)
library(dplyr)

#upload LPS.LPS v Vehicle transcriptomic list with padj<0.05 cutoff included
#select only Genes, log2FC and padj values
#convert ensembl ids to gene ids by comparing to biomart_table
RNAseq <- read_csv("DEGenes_LPSLPSvsVehicle_10-2_14-06-18.csv")
RNAseq <- RNAseq %>% 
  select(c(1, 3, 7)) 
RNAseq <- rename (RNAseq, Gene=X1)
biomart_table <- rename(biomart_table, Gene=ensembl_gene_id)

RNAseq_genes <- inner_join(RNAseq, biomart_table, by = "Gene") 
RNAseq_genes <- RNAseq_genes %>% dplyr::select (c(2, 3, 4))
RNAseq_genes <- subset(RNAseq_genes, select=c(3, 1, 2))
RNAseq_genes <- rename(RNAseq_genes, Gene=hgnc_symbol)
RNAseq_genes <- na.omit(RNAseq_genes) 

#upload LPS.LPS v Vehicle proteomic list
#select only pvalues <0.05
Proteomic <- read_csv("DEProteins_LPSLPSvsVehicle.csv")
Proteomic_genes <- filter(Proteomic, pvalue_ttest1 <0.05) %>%
  dplyr::select(c(2, 6:7))
Proteomic_genes <- na.omit(Proteomic_genes) 

##############determine shared and unique DE genes btwn RNAseq and Proteomics
##############plot in venn diagram
install.packages("sqldf")
library(sqldf)

#combine and generate new dataframe with both RNAseq and Proteomic lists
RNAseqINProteomic <- sqldf('SELECT * FROM RNAseq_genes INNER JOIN Proteomic_genes ON RNAseq_genes.Gene = Proteomic_genes.Gene')

#check if the number of genes are equal to the new dataframe
length(which(Proteomic_genes$Gene %in% RNAseq_genes$Gene))

#find unique genes for each omics method
UniqueRNAseq_genes <- anti_join(RNAseq_genes, Proteomic_genes, by = "Gene" )
UniqueProteomic_genes <- anti_join(Proteomic_genes, RNAseq_genes, by = "Gene")

#find gene count of shared genes btwn omics methods
SharedGenes <- merge(RNAseq_genes, Proteomic_genes, by= "Gene", all = FALSE)

#add LPSLPSvVehicle Unique DE gene lists btwn omics methods
LPSLPSvVehicleVenn <- draw.pairwise.venn(area1 = 2544,area2 = 2723,cross.area = 291,category = c("RNA-Seq", "Proteomics"))

##venn.diagram(
 # x = list(RNAseq_genes, Proteomic_genes),
  #category.names = c("RNA-Seq" , "Proteomics"), 
  #filename = 'LPSLPSvsVehicleOmics_Venn.png',output=TRUE)
##warning message: Length of logical index must be 1 or 2723, not 0 
  #save svg and manipulate on inkscape

ggsave(LPSLPSvVehicleVenn, 
       file = "LPSLPSvsVehicle_omics_Vennplot.svg", device = "svg")



  
  #paired ttest: t.test(x, y, paired = TRUE, alternative = "two.sided")