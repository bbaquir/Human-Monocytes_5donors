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
#only include Log2FC values >1.5 and <-1.5
#find out how many DE genes are upregulated and downregulated
RNAseq <- read_csv("DEGenes_LPSLPSvsVehicle_10-2_14-06-18.csv")
RNAseq <- RNAseq %>% 
  dplyr::select(c(1, 3, 7)) 
#RNAseq <- rename (RNAseq, Gene=X1)
biomart_table <- rename(biomart_table, Gene=ensembl_gene_id)

RNAseq_genes <- inner_join(RNAseq, biomart_table, by = "Gene") 
RNAseq_genes <- RNAseq_genes %>% filter(log2FoldChange>1.5 | log2FoldChange < -1.5, padj<0.05 ) %>% dplyr::select (c(2:4))
RNAseq_genes <- subset(RNAseq_genes, select=c(3, 1, 2))
RNAseq_genes <- rename(RNAseq_genes, Gene=hgnc_symbol)
RNAseq_genes <- na.omit(RNAseq_genes) 

#counts of upregulated and downregulated genes from RNASeq genes
sum(RNAseq_genes[, c("log2FoldChange")]>0)
sum(RNAseq_genes[, c("log2FoldChange")]<0)

###################################################Enriched pathways
#sigora pathways from pvalues <0.05 & <0.01  proteins
RNAseq_pathways <- sigora(reaH, level = 4, queryList = RNAseq_genes$Gene) %>% 
  pluck("summary_results") %>% 
  filter(Bonferroni <=0.001)

#switch to correct folder prior to saving file
write.csv(RNAseq_genes, file= "LPSLPSvsVehicle_0.05_SigoraPathway.csv")
####################################################

#upload LPS.LPS v Vehicle proteomic list
#select only pvalues <0.05
Proteomic <- read_csv("DEProteins_LPSLPSvsVehicle.csv")
Proteomic_genes <- filter(Proteomic, pvalue_ttest1 <0.05) %>%
  dplyr::select(c(2, 6:7))
Proteomic_genes <- na.omit(Proteomic_genes) 

#counts of upregulated and downregulated genes from Proteomic genes
sum(Proteomic_genes[, c("LPSLPS_Log2FC")]>0)
sum(Proteomic_genes[, c("LPSLPS_Log2FC")]<0)

##############determine shared and unique DE genes btwn RNAseq and Proteomics: LPSLPS vs Vehicle
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



##upload LPS v Vehicle transcriptomic list with padj<0.05 cutoff included
##select only Genes, log2FC and padj values
##convert ensembl ids to gene ids by comparing to biomart_table
##find out how many DE genes are upregulated and downregulated
RNAseq <- read_csv("DEGenes_LPSvsVehicle_10-2_14-06-18.csv")
RNAseq <- RNAseq %>% 
  dplyr::select(c(1, 3, 7)) 
RNAseq <- rename (RNAseq, Gene=X1)
biomart_table <- rename(biomart_table, Gene=ensembl_gene_id)

RNAseq_genes <- inner_join(RNAseq, biomart_table, by = "Gene") 
RNAseq_genes <- RNAseq_genes %>% dplyr::select (c(2:4))
RNAseq_genes <- subset(RNAseq_genes, select=c(3, 1, 2))
RNAseq_genes <- rename(RNAseq_genes, Gene=hgnc_symbol)
RNAseq_genes <- na.omit(RNAseq_genes) 

#counts of upregulated and downregulated genes from RNASeq genes
sum(RNAseq_genes[, c("log2FoldChange")]>0)
sum(RNAseq_genes[, c("log2FoldChange")]<0)

#upload LPS.LPS v Vehicle proteomic list
#select only pvalues <0.05
Proteomic <- read_csv("DEProteins_LPSvsVehicle.csv")
Proteomic_genes <- filter(Proteomic, pvalue_ttest1 <0.05) %>%
  dplyr::select(c(2, 6:7))
Proteomic_genes <- na.omit(Proteomic_genes) 

#counts of upregulated and downregulated genes from Proteomic genes
sum(Proteomic_genes[, c("LPS_Log2FC")]>0)
sum(Proteomic_genes[, c("LPS_Log2FC")]<0)

##############determine shared and unique DE genes btwn RNAseq and Proteomics: LPSLPS vs Vehicle
##############plot in venn diagram
#install.packages("sqldf")
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
LPSvVehicleVenn <- draw.pairwise.venn(area1 = 3949,area2 = 430,cross.area = 169,category = c("RNA-Seq", "Proteomics"))

##venn.diagram(
# x = list(RNAseq_genes, Proteomic_genes),
#category.names = c("RNA-Seq" , "Proteomics"), 
#filename = 'LPSLPSvsVehicleOmics_Venn.png',output=TRUE)
##warning message: Length of logical index must be 1 or 2723, not 0 
#save svg and manipulate on inkscape

ggsave(LPSvVehicleVenn, 
       file = "LPSvsVehicle_omics_Vennplot.svg", device = "svg")


  




###upload LPSLPS v LPS transcriptomic list with padj<0.05 cutoff included
###select only Genes, log2FC and padj values
###convert ensembl ids to gene ids by comparing to biomart_table
###find out how many DE genes are upregulated and downregulated
RNAseq <- read_csv("DEGenes_LPSLPSvsLPS_strict_10-2_14-06-18.csv")
RNAseq <- RNAseq %>% 
  dplyr::select(c(1, 3, 7)) 
RNAseq <- rename (RNAseq, Gene=X1)
biomart_table <- rename(biomart_table, Gene=ensembl_gene_id)

RNAseq_genes <- inner_join(RNAseq, biomart_table, by = "Gene") 
RNAseq_genes <- RNAseq_genes %>% dplyr::select (c(2:4))
RNAseq_genes <- subset(RNAseq_genes, select=c(3, 1, 2))
RNAseq_genes <- rename(RNAseq_genes, Gene=hgnc_symbol)
RNAseq_genes <- na.omit(RNAseq_genes) 

#counts of upregulated and downregulated genes from RNASeq genes
sum(RNAseq_genes[, c("log2FoldChange")]>0)
sum(RNAseq_genes[, c("log2FoldChange")]<0)

#upload LPS.LPS v Vehicle proteomic list
#select only pvalues <0.05
Proteomic <- read_csv("DEProteins_LPSLPSvsLPS.csv")
Proteomic_genes <- filter(Proteomic, pvalue_ttest1 <0.05) %>%
  dplyr::select(c(2, 14:15))
Proteomic_genes <- na.omit(Proteomic_genes) 

#counts of upregulated and downregulated genes from Proteomic genes
sum(Proteomic_genes[, c("Log2FC")]>0)
sum(Proteomic_genes[, c("Log2FC")]<0)

##############determine shared and unique DE genes btwn RNAseq and Proteomics: LPSLPS vs Vehicle
##############plot in venn diagram
#install.packages("sqldf")
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

#add LPSLPSvLPS Unique DE gene lists btwn omics methods
LPSLPSvLPSVenn <- draw.pairwise.venn(area1 = 2220,area2 = 3585,cross.area = 385,category = c("RNA-Seq", "Proteomics"))

##venn.diagram(
# x = list(RNAseq_genes, Proteomic_genes),
#category.names = c("RNA-Seq" , "Proteomics"), 
#filename = 'LPSLPSvsVehicleOmics_Venn.png',output=TRUE)
##warning message: Length of logical index must be 1 or 2723, not 0 
#save svg and manipulate on inkscape

ggsave(LPSLPSvLPSVenn, 
       file = "LPSLPSvsLPS_omics_Vennplot.svg", device = "svg")
  #paired ttest: t.test(x, y, paired = TRUE, alternative = "two.sided")





###############################Arjun helped with generating adj pvalues
##############################used bonferroni correction 
#############################FDR and BH are the same correction 
Proteomic_genes$pvalue_ttest1 %>% p.adjust(method = "BH")
table(Proteomic_genes$pvalue_ttest1 %>% p.adjust(method = "BH") <= 0.05)
table(Proteomic$pvalue_ttest1 %>% p.adjust(method = "BH") <= 0.05)
table(Proteomic$pvalue_ttest1 %>% p.adjust(method = "BH") <= 0.10)
table(Proteomic$pvalue_ttest1 %>% p.adjust(method = "BH") <= 0.01)
table(Proteomic$pvalue_ttest1 %>% p.adjust(method = "BH") <= 0.05)
table(Proteomic$pvalue_ttest1 %>% p.adjust(method = "BH") <= 0.010)
table(Proteomic$pvalue_ttest1 %>% p.adjust(method = "BH") <= 0.10)
table(Proteomic$pvalue_ttest1 %>% p.adjust(method = "BH") <= 0.05)
table(Proteomic$pvalue_ttest1 %>% p.adjust(method = "fdr") <= 0.05)
table(Proteomic$pvalue_ttest1 %>% p.adjust(method = "fdr") <= 0.10)
table(Proteomic$pvalue_ttest1 %>% p.adjust(method = "BH") <= 0.10)
table(Proteomic$pvalue_ttest1 %>% p.adjust(method = "bonferroni") <= 0.10)
table(Proteomic$pvalue_ttest1 %>% p.adjust(method = "bonferroni") <= 0.05)
Proteomic$pvalue_ttest1 %>% p.adjust(method = "bonferroni")
Proteomic$pvalue_ttest1 %>% p.adjust(method = "BH")
Proteomic$adjusted <- Proteomic$pvalue_ttest1 %>% p.adjust(method = "BH")




















