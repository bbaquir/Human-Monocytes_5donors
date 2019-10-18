#set working directory
setwd("C:/Users/beverlie/Google Drive/Proteomics Analysis on Human Monocytes")
#setwd on mac
setwd("~/Google Drive/Proteomics Analysis on Human Monocytes")

#load library
install.packages ("tidyverse")
library(tidyverse)


#read data into dataframe
complete.protein.data <- read.csv("MonocyteProteome.csv", header = T)

#################################################GOAL: PCA plots that include D188

#select gene and convert column to rownames
all.donors <- complete.protein.data %>%
  select(c(2, 4:15)) %>%
  distinct(gene, .keep_all = TRUE) %>% 
  column_to_rownames("gene")

#Get rid of NA values (but not needed here)
all.donors[is.na(all.donors)] <- 0

#Get PCA data
all.donors.pca <- prcomp(all.donors, scale=TRUE, center = TRUE)

##with help from Travis
#prcomp summary of SD, proportion of variance, cumulative proportion of variance
percentvars <- summary(all.donors.pca) %>% 
  unclass() %>% 
  pluck("importance") %>%
  as.data.frame() %>% 
  rownames_to_column("Name") %>% 
  select(Name, PC1, PC2) %>% 
  filter(Name == "Proportion of Variance")

all.donors.pca <- all.donors.pca$rotation %>%
  as.data.frame() %>% 
  select(PC1, PC2) %>%
  rownames_to_column("samples") %>% 
  mutate(Donor = str_extract(samples, pattern = "D[0-9]{3}"),
         Treatment = str_replace(samples, pattern = "D[0-9]{3}", replacement = ""))

ggplot(all.donors.pca, aes(PC1, PC2, color = Donor, shape = Treatment)) +
  geom_point(size = 6) +
  labs(x = paste0("PC1: ", percentvars[1, 2]), y = paste0("PC 2: ", percentvars[1, 3]))

all.donors.pca.plot <- ggplot(all.donors.pca, aes(PC1, PC2, color = Donor, shape = Treatment)) +
  geom_point(size = 6) +
  labs(x = paste0("PC1: ", percentvars[1, 2]), y = paste0("PC 2: ", percentvars[1, 3]))

#save PCA, but install and load packages first
install.packages("gridExtra")
library(gridExtra)
#library(grid)
library(tidyverse)
library(svglite)

#save All Donors PCA plot
ggsave ("./AllDonorProteins_pca.svg", all.donors.pca.plot, width=200, height=100, scale=0.075, dpi=300, limitsize=FALSE)


#################################################GOAL: PCA plots that EXCLUDE D188

#select gene and convert column to rownames
three.donors <- complete.protein.data %>%
  select(c(2, 4:9, 13:15)) %>%
  distinct(gene, .keep_all = TRUE) %>% 
  column_to_rownames("gene")

#Get rid of NA values (but not needed here)
three.donors[is.na(three.donors)] <- 0

#Get PCA data
three.donors.pca <- prcomp(three.donors, scale=TRUE, center = TRUE)

##with help from Travis
#prcomp summary of SD, proportion of variance, cumulative proportion of variance
percentvars.three <- summary(three.donors.pca) %>% 
  unclass() %>% 
  pluck("importance") %>%
  as.data.frame() %>% 
  rownames_to_column("Name") %>% 
  select(Name, PC1, PC2) %>% 
  filter(Name == "Proportion of Variance")

three.donors.pca <- three.donors.pca$rotation %>%
  as.data.frame() %>% 
  select(PC1, PC2) %>%
  rownames_to_column("samples") %>% 
  mutate(Donor = str_extract(samples, pattern = "D[0-9]{3}"),
         Treatment = str_replace(samples, pattern = "D[0-9]{3}", replacement = ""))

ggplot(three.donors.pca, aes(PC1, PC2, color = Donor, shape = Treatment)) +
  geom_point(size = 6) +
  labs(x = paste0("PC1: ", percentvars[1, 2]), y = paste0("PC 2: ", percentvars[1, 3]))

three.donors.pca.plot <- ggplot(three.donors.pca, aes(PC1, PC2, color = Donor, shape = Treatment)) +
  geom_point(size = 6) +
  labs(x = paste0("PC1: ", percentvars[1, 2]), y = paste0("PC 2: ", percentvars[1, 3]))

#save PCA, but install and load packages first
install.packages("gridExtra")
library(gridExtra)
#library(grid)
library(tidyverse)
library(svglite)

#save All Donors PCA plot
ggsave ("./ThreeDonorProteins_pca.svg", three.donors.pca.plot, width=200, height=100, scale=0.075, dpi=300, limitsize=FALSE)




