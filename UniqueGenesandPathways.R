##Goal1: to generate unique gene lists from the comparison of  1xLPSvVeh & 2xLPSvVeh. 
##Goal2: input unique gene lists into NA for unsupervised networks and generate pathways
##Goal3: input the same unique gene list into SIGORA for enriched pathways dependent on the gene list
##Goal4: compare and identify common pathways and query the biological relevance for each condition

#setwd
setwd("C:/Users/beverlie/Google Drive/Proteomics Analysis on Human Monocytes/Unique gene and pathway lists/Transcriptomics")

#load tidyverse library the includes dplyr
install.packages ("tidyverse")
library(tidyverse)
library(dplyr)


#read transcriptomic monocyte data into dataframe; both already have padj<0.05
LPS <- read.csv("DEGenes_LPSvsVehicle_10-2_14-06-18_copy.csv")
LLPS <- read.csv("DEGenes_LPSLPSvsVehicle_10-2_14-06-18_copy.csv")

#label "gene" column
names (LLPS) <- c("Gene", "log2FoldChange", "padj")

#find unique genes for each comparison, no FC cutoff yet
unique.LPS <- anti_join(LPS, LLPS, by = "Gene" )
unique.LLPS <- anti_join(LLPS, LPS, by = "Gene")

#find gene count of shared genes, no FC cutoff yet
shared.genes <- merge(LPS, LLPS, by= "Gene", all = FALSE)

#change ensembl gene to hgnc_symbol
unique.LPS <- left_join(unique.LPS, biomart_table, by=c("Gene"="ensembl_gene_id"))
unique.LLPS <- left_join(unique.LLPS, biomart_table, by=c("Gene"="ensembl_gene_id"))

#only include log2FC >1.5|<-1.5
#if this is applied to data, there is no shared genes btwn the two comparisons
sig.unique.LPS <- unique.LPS %>% filter(log2FoldChange < -1.5 | log2FoldChange > 1.5) %>% select ("Gene", "log2FoldChange")
sig.unique.LLPS <- unique.LLPS %>% filter(log2FoldChange < -1.5 | log2FoldChange > 1.5) %>% select ("Gene", "log2FoldChange")

#save csv file and look up pathways on NetworkAnalyst
write.csv(unique.LPS, file = "unique_LPSvVehicle_genes.csv")
##on NA 1st order, nodes=8663, edges=35859, seeds=1758
##on NA zero order, nodes=1003, edges=2548, seeds=1004
write.csv(unique.LLPS, file= "unique_LLPSvVehicle_genes.csv")
##on NA 1st order, nodes=4456, edges=9459, seeds=538
##on NA zero order, nodes=129, edges=164, seeds=129

#check number of genes
nrow(unique.LPS)
nrow(unique.LLPS)
                  
#generate a venn diagram 
install.packages("VennDiagram")
library(VennDiagram)

#add unique DE gene lists and category names, no FC cutoff yet
vennplot <- draw.pairwise.venn(area1 = 4020,area2 = 2679,cross.area = 1958,category = c("LPSvVehicle", "LPSLPSvVehicle"))

#save svg and manipulate on inkscape
ggsave(vennplot, file = "transcript_vennplot.svg", device = "svg")
                   
#script from Travis to generate pathway list from SIGORA

#install and load sigora
install.packages("sigora")
library(sigora)

#sigora with no FC cutoff
sigora_unique.LPS <- sigora(reaH, level = 4, queryList = unique.LPS$Gene) %>%
pluck("summary_results") %>%
  filter(Bonferroni <= 0.001)
sigora_unique.LLPS <- sigora(reaH, level = 4, queryList = unique.LLPS$Gene) %>%
  pluck("summary_results") %>%
  filter(Bonferroni <= 0.001)

#save SIGORA lists as .csv
write.csv(sigora_unique.LPS, file = "sigora_unique_LPS.csv")
write.csv(sigora_unique.LLPS, file = "sigora_unique_LLPS.csv")


##


#found that Platelet Degranulation is a shared LLPS pathway btwn NA_zero order and SIGORA
#downloaded the complete gene list of said pathway from innateDB
#compare these two lists for shared genes

#read innatedb dataframe for LLPS 
PD.innatedb <- read.csv("InnateDB_PlateletDegranulation_GeneList.csv")
colnames(PD.innatedb)

#read innatedb data fraom for LPS
Apop.innatedb <- read.csv("InnateDB_Apoptosis_GeneList.csv")
colnames(Apop.innatedb)

#only include ensembl and gene function
#PD.innatedb.ensembl <- PD.innatedb %>% select (c(4, 16))

#combine only shared genes btwn the two lists: ID's the list input resulting pathway enrichment
shared.PDgenes <- merge(unique.LLPS, PD.innatedb, by.x= "Gene", by.y = "ensembl", all = FALSE)
shared.Apopgenes <- merge(unique.LPS, Apop.innatedb, by.x = "Gene", by.y = "ensembl", all = FALSE)

#change ensembl gene to gene ID
shared.PDgenes <- left_join(shared.PDgenes, biomart_table, by=c("Gene"="ensembl_gene_id"))
shared.Apopgenes <- left_join(shared.Apopgenes, biomart_table, by = c("Gene" = "ensembl_gene_id"))

#save shared gene list
write.csv(shared.PDgenes, file = "LLPS_Shared_PlateletDegranulation_GeneList.csv")
write.csv(shared.Apopgenes, file = "LPS_Shared_Apoptosis_GeneList.csv")
