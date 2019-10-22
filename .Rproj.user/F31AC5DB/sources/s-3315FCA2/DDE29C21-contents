##HELLO TEST
##DUMB TEST
#set working directory
setwd("C:/Users/beverlie/Google Drive/Proteomics Analysis on Human Monocytes")
#setwd on mac
setwd("~/Google Drive/Proteomics Analysis on Human Monocytes")

#load library
install.packages ("tidyverse")
library(tidyverse)



##################################PCA PLOT###########################################
#setwd("C:/Users/beverlie/Google Drive/Proteomics Analysis on Human Monocytes")

#load library
library(tidyverse)

#read data into dataframe
mp <- read.csv("MonocyteProteome.csv")


#look at the first 3 rows
head(mp, n=3)

#select gene and convert column to rownames
mp_light <- mp %>%
  select("gene", -starts_with("description"), starts_with("D1"), starts_with("D2")) %>%
  distinct(gene, .keep_all = TRUE) %>% 
  column_to_rownames("gene")

# Get rid of NA values
mp_light[is.na(mp_light)] <- 0

#mpPCA plot
pca_mp <- prcomp(mp_light, scale=TRUE, center = TRUE)

##with help from Travis
#prcomp summary of SD, proportion of variance, cumulative proportion of variance
percentvars <- summary(pca_mp) %>% 
  unclass() %>% 
  pluck("importance") %>%
  as.data.frame() %>% 
  rownames_to_column("Name") %>% 
  select(Name, PC1, PC2) %>% 
  filter(Name == "Proportion of Variance")

pca_donor <- pca_mp$rotation %>%
  as.data.frame() %>% 
  select(PC1, PC2) %>%
  rownames_to_column("samples") %>% 
  mutate(Donor = str_extract(samples, pattern = "D[0-9]{3}"),
         Treatment = str_replace(samples, pattern = "D[0-9]{3}", replacement = ""))

ggplot(pca_donor, aes(PC1, PC2, color = Donor, shape = Treatment)) +
  geom_point(size = 6) +
  labs(x = paste0("PC1: ", percentvars[1, 2]), y = paste0("PC 2: ", percentvars[1, 3]))

pca.plot <- ggplot(pca_donor, aes(PC1, PC2, color = Donor, shape = Treatment)) +
  geom_point(size = 6) +
  labs(x = paste0("PC1: ", percentvars[1, 2]), y = paste0("PC 2: ", percentvars[1, 3]))

#save PCA, but install and load packages first
install.packages("gridExtra")
library(gridExtra)
#library(grid)
library(tidyverse)
library(svglite)
ggsave ("./prot_pca.svg", pca.plot, width=200, height=100, scale=0.075, dpi=300, limitsize=FALSE)
dev.off()





################################HEATMAP##################################
#read PROTEOMICS data into dataframe
prot <- read.csv("HumanMonocyteDonor_ProteomicsAnalysis.csv", header = T)

#look at the first 3 rows and column names
head(prot, n=3)
colnames(prot)

#remove description 
prot_light <- prot %>% select (c(2, 16:23)) %>%   
  distinct(gene, .keep_all = TRUE) %>% 
  column_to_rownames("gene")

#install packages for heatmap visualization
install.packages("BiocManager")
library(BiocManager)
#install.packages("ComplexHeatmap")

#make simple heatmap
heatmap(as.matrix(prot_light), scale="none")

install.packages("gplots")
library(gplots)
heatmap_plot <- heatmap.2(as.matrix(prot_light), scale="none", col = bluered(100), trace = "none", density.info="none", cexRow = 0.1, cexCol = 0.5, margins = c(5, 20))

heatmap_plot2 <- heatmap.2(as.matrix(prot_light), scale="none", col = bluered(100), trace = "none", density.info="none", cexRow = 0.1, cexCol = 0.5, keysize = 0.75)

#save heatmap as svg: easiest way to do this is exporting and saving image as .svg file rather than the following line of code
par("mar")
par(mar=c(1,1,1,1))
ggsave ("./prot_heatmapnew.svg", heatmap_plot, width=200, height=400, scale=0.075, dpi=300, limitsize=FALSE)

#44 gene heatmap save
ggsave ("./44genehm.svg", combo_signature, width = 200, height = 400, scale = 0.075, dpi = 300, limitsize = FALSE)
##################################Q: do protein and transcript lists agree with each other?#########################




hm_40signature <- heatmap.2(as.matrix(combo_df), scale="none", col = bluered(100), trace = "none", density.info="none", cexRow = 0.5, cexCol = 0.5, breaks = col_breaks)

########################find out what CR genes are found in the protein list##############################

#read PROTEOMICS data into dataframe
prot <- read.csv("HumanMonocyteDonor_ProteomicsAnalysis.csv", header = T)

#look at the first 3 rows
head(prot, n=3)
colnames(prot)

#remove description 
##prot_light <- prot %>% select (c(2, 16:23)) %>%   distinct(gene, .keep_all = TRUE) %>% column_to_rownames("gene")

#select gene, dd values and convert column to rownames
#prot_lighter <- prot %>%
  select (c (2, 17, 19, 21, 23)) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  column_to_rownames("gene")

#select gene, dd values
prot_lighter <- prot %>% select (c(2, 17, 19, 21, 23)) %>% distinct(gene, .keep_all = TRUE)

#take average dd values
prot_lighter$dd.avg = rowMeans(prot_lighter[,2:5])

#keep only genes and dd.avg column
prot_lightest <- prot_lighter %>% select(c(1,6))
  
#upload ETT gene signature
signature <- read.csv("pbmc_mono_etsig_de_genes.csv", header=TRUE)
signature$fc_LPS_LPS_Veh %>% log2()

#change ETT sig to log2FC
sig <- read.csv("ETT_monocyte.csv")
sig$fc_LPS_LPS_Veh %>% log2()
sig$fc_LPS_LPS_Veh[!is.finite(sig$fc_LPS_LPS_Veh)] <- 0
sig2 <- sig %>% rename(gene= exter0l_gene_0me, log2FC = fc_LPS_LPS_Veh)

######sig_newer <- na.omit(sig_new)
sig$fc_LPS_LPS_Veh[!is.finite(sig$fc_LPS_LPS_Veh)] <- 0
sig$fc_LPS_LPS_Veh[!is.finite(sig$fc_LPS_LPS_Veh)] <- 0
sig_light <- sig %>% rename(gene= exter0l_gene_0me, log2FC = fc_LPS_LPS_Veh)
sig$fc_LPS_LPS_Veh %>% log2()

#change column external_gene_name to just gene
library(dplyr)
#compare pbmcs versus donor protein lists
signature_light <- signature %>% rename (gene = external_gene_name) %>% select(c(1,2))

#compare donor monocytes with donor proteim
signature_light <- signature %>% rename (gene = external_gene_name) %>% select(c(1,5)) %>% na.omit()
                                                                                                  
signature_lighter <- signature_light[!duplicated(signature_light$gene),,drop=FALSE]

#combine lists and call the signature genes from the protein list 
install.packages("sqldf")
library(sqldf)

#combine and generate new dataframe with both pbmc and protein lists, just curated
signINprot <- sqldf('SELECT * FROM signature_lighter INNER JOIN prot_lighter ON signature_lighter.gene = prot_lighter.gene')

#combine and generate new dataframe with both monocytes and protein lists, just curated
signINprot <- sqldf('SELECT * FROM signature_light INNER JOIN prot_lighter ON signature_light.gene = prot_lighter.gene')

#check if the number of genes are equal to the new dataframe
length(which(prot_lighter$gene %in% signature_lighter$gene))

#combine and generate new dataframe for pbmcs and protein but with one list of genes
combo_signature <- inner_join(signature_lighter, prot_lighter, by = "gene") %>% column_to_rownames("gene")

#combine and generate new dataframe for monocytes and proteins but with one list of genes
combo_signature <- inner_join(signature_light, prot_lighter, by = "gene") %>% column_to_rownames("gene")


shared.CRgenes <- merge(signature_light, prot_lighter, by= "gene", all = FALSE)
write.csv(shared.CRgenes, file = "Shared_44CRgenes.csv")
###############combine and generate new dataframe for monocytes and proteins with one list of genes in log2FC
combo_sig <- inner_join(sig2, prot_lighter, by = "gene") %>% column_to_rownames("gene")
combo_sig2 <- combo_sig %>% select(c(1:5)) %>% 
  distinct(gene, .keep_all = TRUE)

#average protein dd values and graph correlation with pbmcs
combo_signature$dd.avg = rowMeans(combo_signature[,2:5])
qplot(combo_signature$fc_pbmc, combo_signature$dd.avg)

#average protein dd values and graph correlation with monocytes
combo_signature$dd.avg = rowMeans(combo_signature[,2:5])
qplot(combo_signature$fc_LPS_LPS_Veh, combo_signature$dd.avg)

#combo_singature <- inner_join(signature_lighter, prot_lighter, by.x="gene", by.y="gene")

#plot heatmap comparing pbmcs and protein
hm_49signature <- heatmap.2(as.matrix(combo_signature), scale="none", col = bluered(100), trace = "none", density.info="none", cexRow = 0.5, cexCol = 0.5)

#plot heatmap comparing monocytes and protein
hm_44signature <- heatmap.2(as.matrix(combo_signature), scale="none", col = bluered(100), trace = "none", density.info="none", cexRow = 0.5, cexCol = 0.5)

heatmap.2(as.matrix(combo_sig), scale="none", col = bluered (100), trace = "none", density.info="none", cexRow = 0.5, cexCol = 0.5)

#only proteins listed
prot_lighter2 <- prot_lighter %>% select (c(1:5)) %>%  column_to_rownames("gene")
heatmap.2(as.matrix(prot_lighter2), scale="none", col = bluered (100), trace = "none", density.info="none", cexRow = 0.5, cexCol = 0.5)

####################correlate protein values to transcript expression####################################

#average protein dd values and graph correlation with pbmcs
combo_signature$dd.avg = rowMeans(combo_signature[,2:5])
qplot(combo_signature$fc_pbmc, combo_signature$dd.avg)

#average protein dd values and graph correlation with monocytes
combo_signature$dd.avg = rowMeans(combo_signature[,2:5])
qplot(combo_signature$fc_LPS_LPS_Veh, combo_signature$dd.avg)

#################find enriched pathways in proteome list using sigora#####################################

#install and load sigora
install.packages("sigora")
library(sigora)

#script from Travis, change for protein list
##sigora(reaH, level = 4, queryList = degenes$ensemblGene) %>%
  pluck("summary_results") %>%
  filter(Bonferroni <= 0.001)

#sigora pathways from proteomics
sigora_prot <- sigora(reaH, level = 4, queryList = prot_lightest$gene) %>% 
  pluck("summary_results") %>% 
  filter(Bonferroni <=0.001)
#sigora pathways from significant proteins
sigora_signprot <- sigora(reaH, level = 4, queryList =signprot$gene) %>% 
  pluck("summary_results") %>% 
  filter(Bonferroni <=0.001)

#save sigora_prot as .csv
write.csv(sigora_prot, file = "sigora_avgddproteins.csv")
write.csv(sigora_signprot, file = "sigora_significantproteins.csv")


###REDO SCRIPT TO GENERATE SINGLELPSvVEH SIGORA PATHWAYS FROM SIGNIFICANT PROTEINS
#select gene, d values
dprot_lighter <- prot %>% select (c(2, 16, 18, 20, 22)) %>% distinct(gene, .keep_all = TRUE)

#take average d values
dprot_lighter$d.avg = rowMeans(dprot_lighter[,2:5])

#keep only genes and d.avg column
dprot_lightest <- dprot_lighter %>% select(c(1,6))

#sigora pathways from proteomics
dsigora_prot <- sigora(reaH, level = 4, queryList = dprot_lightest$gene) %>% 
  pluck("summary_results") %>% 
  filter(Bonferroni <=0.001)

#sigora pathways from significant proteins
#sigora_signprot <- sigora(reaH, level = 4, queryList =signprot$gene) %>% 
  pluck("summary_results") %>% 
  filter(Bonferroni <=0.001)

#save sigora_prot as .csv
write.csv(dsigora_prot, file = "sigora_avgdproteins.csv")
#write.csv(sigora_signprot, file = "sigora_significantproteins.csv")

########################convert proteome gene list into ensembl for NA application#####################

prot_ensembl <- read.csv("prot_ensembl.csv")
signprot <- read.csv("LPSLPSvVeh_significantproteins.csv")

#combine and generate new protein dataframe with one list of ensembl genes
prot_ensembl_join <- left_join(prot_lightest, biomart_table, by=c("gene"="hgnc_symbol"))
join <- left_join(signprot, biomart_table, by= c("gene"="hgnc_symbol"))

write.csv(prot_ensembl_join, file = "prot_ensemblID.csv")
write.csv(join, file = "significanensembltproteins.csv")
