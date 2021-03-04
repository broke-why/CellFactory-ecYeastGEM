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
input_prediction <- read.table("../results/targetsMatrix_ecFSEOF.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
# mech-validated
# input_prediction <- read.table("../results/targetsMatrix_mech_validated.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
# compatible
# input_prediction <- read.table("../results/targetsMatrix_compatible.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)


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

'
for (x in product_families){
  print(x)
  
  # x <- "amino acid"
  # x <- "protein"   
  # x <- "flavonoid" # no over expression
  # x <- "terpene"
  # x <- "fatty acids and lipids" # for test
 
  
  # calculation
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
  all_product <- colnames(input_prediction)[-(c(1:4))]
  metabolite_need_check <- setdiff(products_one_family, all_product)
  print("metabolite name need check:---")
  print(metabolite_need_check)
  metabolite_with_gene_targes <- intersect(products_one_family, all_product)
  colnames <- c("genes", metabolite_with_gene_targes)
  gene_target <- input_prediction[, colnames]
  
  
  # analyse the common genes from one product families
  
  gene_target_value <- input_prediction[, metabolite_with_gene_targes]
  
  # option 1 calculate the common gene targets
  # OE
  gene_target_OE <- gene_target[apply(gene_target_value == 3, 1, all), ]
  # knock-down
  gene_target_KD <- gene_target[apply(gene_target_value == 2, 1, all), ]
  # knock-out
  gene_target_KO <- gene_target[apply(gene_target_value == 1, 1, all), ]
  
  # option 2 calculate the top common genes targets
  gene_top_rank <- data.frame(genes=gene_target[, c("genes")], stringsAsFactors = FALSE)
  gene_top_rank$OE_occurance <- NA
  gene_top_rank$KD_occurance <- NA
  gene_top_rank$KO_occurance <- NA
  gene_top_rank$products_num <- length(metabolite_with_gene_targes)
  for(i in 1:nrow(gene_top_rank)){
    print(gene_top_rank$genes[i])
    ss <- gene_target[i,]
    OE_num <- length(which(ss==3))
    KD_num <- length(which(ss==2))
    KO_num <- length(which(ss==1))
    gene_top_rank$OE_occurance[i] <- OE_num
    gene_top_rank$KD_occurance[i] <- KD_num
    gene_top_rank$KO_occurance[i] <- KO_num
  }
  # for the enrichment analysis
  gene_select <- filter(gene_top_rank, OE_occurance >= length(metabolite_with_gene_targes)*0.5) #only choose the gene targets overexpressed in more than half of products in one families
  gene_select0 <- paste0(gene_select$genes, collapse = ";")
}'

# this function is used to run the above pipeline
# this function could find the common gene targets used for OE, KD, KO over half of products in one family
find_common_gene_target <- function(product_class, input_prediction0 = input_prediction) {
  x <- product_class
  # calculation
  products_one_family <- chemicals_info$ecModel[chemicals_info$class == x]
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
    str_replace_all(., "[:digit:]", "") %>%
    str_to_lower(.) # remove the number at the start of string %>%

  # extract gene target
  # check which metaboite does not have result
  all_product <- colnames(input_prediction0)[-(c(1:4))]
  metabolite_need_check <- setdiff(products_one_family, all_product)
  print("metabolite name need check:---")
  print(metabolite_need_check)
  metabolite_with_gene_targes <- intersect(products_one_family, all_product)
  colnames <- c("genes", metabolite_with_gene_targes)
  gene_target <- input_prediction0[, colnames]


  # analyse the common genes from one product families
  gene_target_value <- input_prediction0[, metabolite_with_gene_targes]

  # option 1 calculate the common gene targets
  # OE
  gene_target_OE <- gene_target[apply(gene_target_value == 3, 1, all), ]
  # knock-down
  gene_target_KD <- gene_target[apply(gene_target_value == 2, 1, all), ]
  # knock-out
  gene_target_KO <- gene_target[apply(gene_target_value == 1, 1, all), ]

  # option 2 calculate the top common genes targets
  gene_top_rank <- data.frame(genes = gene_target[, c("genes")], stringsAsFactors = FALSE)
  gene_top_rank$OE_occurance <- NA
  gene_top_rank$KD_occurance <- NA
  gene_top_rank$KO_occurance <- NA
  gene_top_rank$products_num <- length(metabolite_with_gene_targes)
  for (i in 1:nrow(gene_top_rank)) {
    print(gene_top_rank$genes[i])
    ss <- gene_target[i, ]
    OE_num <- length(which(ss == 3))
    KD_num <- length(which(ss == 2))
    KO_num <- length(which(ss == 1))
    gene_top_rank$OE_occurance[i] <- OE_num
    gene_top_rank$KD_occurance[i] <- KD_num
    gene_top_rank$KO_occurance[i] <- KO_num
  }

  # for the enrichment analysis
  gene_select_OE <- filter(gene_top_rank, OE_occurance >= length(metabolite_with_gene_targes) * 0.5) # only choose the gene targets overexpressed in more than half of products in one families
  gene_select_KD <- filter(gene_top_rank, KD_occurance >= length(metabolite_with_gene_targes) * 0.5) # only choose the gene targets overexpressed in more than half of products in one families
  gene_select_KO <- filter(gene_top_rank, KO_occurance >= length(metabolite_with_gene_targes) * 0.5) # only choose the gene targets overexpressed in more than half of products in one families

  return(list(gene_select_OE, gene_select_KD, gene_select_KO))

}




# test
# ecFSEOF
input_prediction <- read.table("../results/targetsMatrix_ecFSEOF.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
# mech-validated
# input_prediction <- read.table("../results/targetsMatrix_mech_validated.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
# compatible
# input_prediction <- read.table("../results/targetsMatrix_compatible.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
output00 <- find_common_gene_target(product_class = "amino acid", input_prediction0 = input_prediction)
OE_out <- output00[[1]]
KD_out <- output00[[2]]
KO_out <- output00[[3]]


output00 <- find_common_gene_target(product_class = "fatty acids and lipids", input_prediction0 = input_prediction)
OE_out <- output00[[1]]
KD_out <- output00[[2]]
KO_out <- output00[[3]]


output00 <- find_common_gene_target(product_class = "protein", input_prediction0 = input_prediction)
OE_out <- output00[[1]]
KD_out <- output00[[2]]
KO_out <- output00[[3]]


output00 <- find_common_gene_target(product_class = "organic acid", input_prediction0 = input_prediction)
OE_out <- output00[[1]]
KD_out <- output00[[2]]
KO_out <- output00[[3]]




# ecFSEOF
# input_prediction <- read.table("../results/targetsMatrix_ecFSEOF.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
# mech-validated
# input_prediction <- read.table("../results/targetsMatrix_mech_validated.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
# compatible
input_prediction <- read.table("../results/targetsMatrix_compatible.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
output00 <- find_common_gene_target(product_class = "amino acid", input_prediction0 = input_prediction)
OE_out <- output00[[1]]
KD_out <- output00[[2]]
KO_out <- output00[[3]]


output00 <- find_common_gene_target(product_class = "fatty acids and lipids", input_prediction0 = input_prediction)
OE_out <- output00[[1]]
KD_out <- output00[[2]]
KO_out <- output00[[3]]


output00 <- find_common_gene_target(product_class = "protein", input_prediction0 = input_prediction)
OE_out <- output00[[1]]
KD_out <- output00[[2]]
KO_out <- output00[[3]]


output00 <- find_common_gene_target(product_class = "organic acid", input_prediction0 = input_prediction)
OE_out <- output00[[1]]
KD_out <- output00[[2]]
KO_out <- output00[[3]]

