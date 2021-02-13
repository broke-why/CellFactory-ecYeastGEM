# This script is used to compare the predicted genes from ecFSEOF and the genes used in the reference strains.

library(ggplot2)
library(tidyverse)
library(hongR)
library(readxl)

# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}

# data preparation
# input the genes from background strains
gene_background <- read.table("../ComplementaryData/genetic_background.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
# only keep products with gene background information
# as a whole only 56 products have gene background information
gene_background0 <- filter(gene_background, !(overexpression=="" & downregulation=="" & deletion==""))

# input the predicted genes targets without background information
# ecFSEOF
ecFSEOF <- read.table("../results/targetsMatrix_ecFSEOF.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
# mech-validated
mech_validated <- read.table("../results/targetsMatrix_mech_validated.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
# compatible
compatible_gene <- read.table("../results/targetsMatrix_compatible.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)


find_common_gene <- function(s1, s2){
  common <- vector()
  # s1 <- background_OE
  # s2 <- predict_OE
  s10 <- unlist(str_split(s1, ","))
  s10 <- str_trim(s10, side = "both")
  s2 <- str_trim(s2, side = "both")
  common <- intersect(s10,s2)
  if(length(common)){
  return(common)} esle{
    return(NA)
  }
}

# comparison
# the comparsion with start from ecFSEOF, mech_validated and compatible
gene_background0$common_OE <- NA
gene_background0$common_KD <- NA
gene_background0$common_KO <- NA

for (i in 1:nrow(gene_background0)){
  print(i)
  product0 <- gene_background0$results_folder[i]
  product00 <- str_replace(product0, "_targets", "")
  product00 <-  product00 %>% str_replace_all(.,"-","")
  product00 <- paste(product00, "_del_",sep = "")
  # as the name of products need to be unified
  product00 <- str_replace(product00, "^[:digit:]","") # remove the number at the start of string
  
  background_OE <- gene_background0$overexpression[i]
  background_KD <- gene_background0$downregulation[i]
  background_KO <- gene_background0$deletion[i]
  
  # find the predicted genes based on the product name
  col_select <- c(c("genes","shortNames","enzymes","subSystems"),product00)
  predict_tagets <- ecFSEOF[,col_select]
  predict_OE <- predict_tagets$genes[predict_tagets$sadenosyllmethionine_del_==3]
  predict_KD <- predict_tagets$genes[predict_tagets$sadenosyllmethionine_del_==2]
  predict_KO <- predict_tagets$genes[predict_tagets$sadenosyllmethionine_del_==1]
  
  common_OE <- find_common_gene(background_OE, predict_OE)
  common_KD <- find_common_gene(background_KD, predict_KD)
  common_KO <- find_common_gene(background_KO, predict_KO)
  
  gene_background0$common_OE[i] <- common_OE
  gene_background0$common_KD[i] <- common_KD
  gene_background0$common_KO[i] <- common_KO
}




















