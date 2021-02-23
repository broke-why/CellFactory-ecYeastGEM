# This script is used to compare the predicted genes from ecFSEOF and the genes used in the reference strains.
library(ggplot2)
library(tidyverse)
library(readxl)
library(pheatmap)
library(viridis)

# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  path <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}

# input the predicted genes targets without background information
# ecFSEOF
# ecFSEOF <- read.table("../results/targetsMatrix_ecFSEOF.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
# mech-validated
# mech_validated <- read.table("../results/targetsMatrix_mech_validated.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
# compatible
compatible_gene <- read.table("../results/targetsMatrix_compatible.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)

# product family information
chemicals_info <- read_excel("../ComplementaryData/chemicals_info.xlsx")
chemicals_info$Name0 <- str_replace_all(chemicals_info$ecModel, ".mat", "")
chemicals_info$Name0 <- str_replace_all(chemicals_info$Name0, "^ec", "")
chemicals_info$Name0 <- str_to_lower(chemicals_info$Name0)


# unique product families
product_families <- unique(chemicals_info$class)

# analyze the result based on product families

# Note two metabolite need check from file chemical_info -- sucretca_del_; salidroside_del_
# maybe the result does not contain these two products


for (x in product_families){
  print(x)
  products_one_family <- chemicals_info$ecModel[chemicals_info$class==x]
  # name standardization
  
  products_one_family <- products_one_family %>%
    str_replace_all(., "\\.mat", "") %>%
    str_replace_all(., "-", "") %>%
    str_replace_all(., "_", "") %>%
    str_replace_all(., ",", "") %>%
    str_replace_all(., "\\(", "") %>%
    str_replace_all(., "\\)", "") %>%
    paste(., "_del_", sep = "") %>% # as the name of products need to be unified
    str_replace_all(., "^ec", "") %>%
    str_replace_all(., "[:digit:]", "")  %>%
    str_to_lower(.)# remove the number at the start of string %>%
  
  # extract gene target
  # check which metaboite does not have result
  all_product <- colnames(compatible_gene)[-(c(1:4))]
  metabolite_need_check <- setdiff(products_one_family, all_product)
  print("metabolite name need check:---")
  print(metabolite_need_check)
  metabolite_with_gene_targes <- intersect(products_one_family, all_product)
  colnames <- c("genes", metabolite_with_gene_targes)
  gene_target <- compatible_gene[, colnames]
  
  
  # analyse the common genes from one product families
  
  gene_target_value <- compatible_gene[, metabolite_with_gene_targes]
  # OE
  gene_target_OE <- gene_target[apply(gene_target_value == 3, 1, all), ]
  # knock-down
  gene_target_KD <- gene_target[apply(gene_target_value == 2, 1, all), ]
  # knock-out
  gene_target_KO <- gene_target[apply(gene_target_value == 1, 1, all), ]
  

}


















