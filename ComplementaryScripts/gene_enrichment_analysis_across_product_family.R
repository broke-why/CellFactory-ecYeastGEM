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
  #print("metabolite name need check:---")
  #print(metabolite_need_check)
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
    #print(gene_top_rank$genes[i])
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

input_prediction <- read.table("../results/targetsMatrix_ecFSEOF.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
newDF <- data.frame()
families <- unique(chemicals_info$class)
for (family in families){
  output00 <- find_common_gene_target(product_class = family, input_prediction0 = input_prediction)
  OE_out <- output00[[1]]
  KD_out <- output00[[2]]
  KO_out <- output00[[3]]
  newDF <- rbind(newDF,cbind(family,nrow(OE_out),nrow(KD_out),nrow(KO_out)))
}
colnames(newDF) <- c("family","OE","KD","KO")
write.table(newDF,'../results/targetsOverlap_FSEOF.txt',sep='\t',quote=FALSE,row.names =FALSE)



input_prediction <- read.table("../results/targetsMatrix_mech_validated.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
newDF <- data.frame()
families <- unique(chemicals_info$class)
for (family in families){
  output00 <- find_common_gene_target(product_class = family, input_prediction0 = input_prediction)
  OE_out <- output00[[1]]
  KD_out <- output00[[2]]
  KO_out <- output00[[3]]
  newDF <- rbind(newDF,cbind(family,nrow(OE_out),nrow(KD_out),nrow(KO_out)))
}
colnames(newDF) <- c("family","OE","KD","KO")
write.table(newDF,'../results/targetsOverlap_mechVal.txt',sep='\t',quote=FALSE,row.names =FALSE)

input_prediction <- read.table("../results/targetsMatrix_compatible.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
newDF <- data.frame()
families <- unique(chemicals_info$class)
for (family in families){
  output00 <- find_common_gene_target(product_class = family, input_prediction0 = input_prediction)
  OE_out <- output00[[1]]
  KD_out <- output00[[2]]
  KO_out <- output00[[3]]
  newDF <- rbind(newDF,cbind(family,nrow(OE_out),nrow(KD_out),nrow(KO_out)))
}
colnames(newDF) <- c("family","OE","KD","KO")
write.table(newDF,'../results/targetsOverlap_compatible.txt',sep='\t',quote=FALSE,row.names =FALSE)


# here maybe we use the targetsMatrix_mech_validated to do the enrichment analysis
# output over-expressed gene targets
input_prediction <- read.table("../results/targetsMatrix_mech_validated.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
#newDF <- data.frame()
families <- unique(chemicals_info$class)

out_gene_list = data.frame(families=families, stringsAsFactors = FALSE)
gene_list_OE= c()
for (family in families){
  output00 <- find_common_gene_target(product_class = family, input_prediction0 = input_prediction)
  OE_out <- output00[[1]]
  gene_list <- paste0(OE_out$genes, collapse = ",")
  print(family)
  print(gene_list)
  gene_list_OE <- c(gene_list_OE, gene_list)
}
out_gene_list$families <- families
out_gene_list$OE_genelist <- gene_list_OE
write.table(out_gene_list, "../results/targetsMatrix_mech_validated_OE_geneList.txt", row.names = FALSE, sep = "\t")

# combine knock out and knock down together and then output the gene list
out_gene_list2 = data.frame(families=families, stringsAsFactors = FALSE)
gene_list_KD_KO= c()
for (family in families){
  output00 <- find_common_gene_target(product_class = family, input_prediction0 = input_prediction)
  
  KD_out <- output00[[2]]
  KO_out <- output00[[3]]
  combine_out <- unique(c(KD_out$genes, KO_out$genes))
  
  gene_list <- paste0(combine_out, collapse = ",")
  print(family)
  print(gene_list)
  gene_list_KD_KO <- c(gene_list_KD_KO, gene_list)
}
out_gene_list2$families <- families
out_gene_list2$KO_KD_genelist <- gene_list_KD_KO
write.table(out_gene_list2, "../results/targetsMatrix_mech_validated_KD_KO_geneList.txt", row.names = FALSE, sep = "\t")


# do the enrichment analysis manually based on online toobox- DAVID https://david.ncifcrf.gov/summary.jsp


# 1 analyze the result from OE
result_dir <- "../results/gene_enrichment_analysis_for_overexpressed_genes"
subfile <- list.files(result_dir)

families1 <-families[families!="other"]


for (file0 in families1){
  print(file0)
  
  # for the test
  #file0 <- "fatty acids and lipids"
  
  file1 <- paste(file0, ".txt", sep = "")
  file_input <- paste(result_dir, file1, sep = "/")
  result_analysis <- readLines(file(file_input))
  result_analysis <-  result_analysis[!str_detect(result_analysis, "Enrichment Score")]
  result_analysis <-  result_analysis[str_detect(result_analysis, "\t")]
  file00 <- paste("update_", file1, sep = "")
  write_lines(result_analysis, paste(result_dir, file00, sep = "/"))
  
  newfile <- read.table(paste(result_dir, file00, sep = "/"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  newfile$Count <- as.numeric(newfile$Count)
  newfile$PValue <- as.numeric(newfile$Benjamini)
  newfile <- newfile[newfile$Category != "Category",]
  
  # plot the bar plot
  # plot
  # initial filteration
  newfile <- filter(newfile, Count >=4 & PValue <= 0.05)
  # remove duplicated term
  newfile <- newfile[!duplicated(newfile$Term),]
  # order
  Factor <-  newfile[order(newfile$Count,decreasing = TRUE),]
  newfile$Term <-factor(newfile$Term, levels=Factor$Term)
  
  # based on keggg
  newfile_kegg <- newfile[str_detect(newfile$Term,"sce"),] # for the kegg pathway maybe remove some very general pathways
  newfile_go <- newfile[str_detect(newfile$Term,"GO:"),]
  # plot
  fileName <- paste('../results/plots/geneSet_KEGG_enrichment_OE_',file0,'.png',sep='')
  png(fileName,width=700, height=600)
  p1 <- ggplot(data=newfile_kegg , aes(x=Term, y=Count)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    theme(legend.position = c(0.85, 0.2)) +
    theme(axis.text=element_text(size=10, family="Arial"),
          axis.title=element_text(size=12,family="Arial"),
          legend.text = element_text(size=10, family="Arial")) +
    ggtitle(file0) #+
  plot(p1)
  dev.off()
  
  fileName <- paste('../results/plots/geneSet_GO_enrichment_OE_',file0,'.png',sep='')
  png(fileName,width=700, height=600)
  p2 <- ggplot(data=newfile_go , aes(x=Term, y=Count)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    theme(legend.position = c(0.85, 0.2)) +
    theme(axis.text=element_text(size=10, family="Arial"),
          axis.title=element_text(size=12,family="Arial"),
          legend.text = element_text(size=10, family="Arial")) +
    ggtitle(file0) #+
  plot(p2)
  dev.off()
}





# 2 analyze the result from KD + KO
result_dir <- "../results/gene_enrichment_analysis_for_KD_KO_genes"
subfile <- list.files(result_dir)
subfile0 <- str_replace_all(subfile, ".txt", "")
families1 <-intersect(families, subfile0)


for (file0 in families1){
  print(file0)
  
  # for the test
  #file0 <- "fatty acids and lipids"
  
  file1 <- paste(file0, ".txt", sep = "")
  file_input <- paste(result_dir, file1, sep = "/")
  result_analysis <- readLines(file(file_input))
  result_analysis <-  result_analysis[!str_detect(result_analysis, "Enrichment Score")]
  result_analysis <-  result_analysis[str_detect(result_analysis, "\t")]
  file00 <- paste("update_", file1, sep = "")
  write_lines(result_analysis, paste(result_dir, file00, sep = "/"))
  
  newfile <- read.table(paste(result_dir, file00, sep = "/"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  newfile$Count <- as.numeric(newfile$Count)
  newfile$PValue <- as.numeric(newfile$Benjamini)
  newfile <- newfile[newfile$Category != "Category",]
  
  # plot the bar plot
  # plot
  # initial filteration
  newfile <- filter(newfile, Count >=4 & PValue <= 0.05)
  # remove duplicated term
  newfile <- newfile[!duplicated(newfile$Term),]
  # order
  Factor <-  newfile[order(newfile$Count,decreasing = TRUE),]
  newfile$Term <-factor(newfile$Term, levels=Factor$Term)
  
  # based on keggg
  newfile_kegg <- newfile[str_detect(newfile$Term,"sce"),] # for the kegg pathway maybe remove some very general pathways
  newfile_go <- newfile[str_detect(newfile$Term,"GO:"),]
  # plot
  fileName <- paste('../results/plots/geneSet_KEGG_enrichment_KD_KO_',file0,'.png',sep='')
  png(fileName,width=700, height=600)
  p1 <- ggplot(data=newfile_kegg , aes(x=Term, y=Count)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    theme(legend.position = c(0.85, 0.2)) +
    theme(axis.text=element_text(size=10, family="Arial"),
          axis.title=element_text(size=12,family="Arial"),
          legend.text = element_text(size=10, family="Arial")) +
    ggtitle(file0) #+
  plot(p1)
  dev.off()
  
  fileName <- paste('../results/plots/geneSet_GO_enrichment_KD_KO_',file0,'.png',sep='')
  png(fileName,width=700, height=600)
  p2 <- ggplot(data=newfile_go , aes(x=Term, y=Count)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    theme(legend.position = c(0.85, 0.2)) +
    theme(axis.text=element_text(size=10, family="Arial"),
          axis.title=element_text(size=12,family="Arial"),
          legend.text = element_text(size=10, family="Arial")) +
    ggtitle(file0) #+
  plot(p2)
  dev.off()
}














