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
colors <- cividis(11)
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
#Presence threshold for defining overlaping genes
threshold <- 0.5

# analyze the result based on product families
# Note two metabolite need check from file chemical_info -- sucretca_del_; salidroside_del_
# maybe the result does not contain these two products

# this function is used to run the above pipeline
# this function could find the common gene targets used for OE, KD, KO over half of products in one family
find_common_gene_target <- function(product_class, input_prediction0 = input_prediction,threshold) {
  x <- product_class
  # calculation
  products_one_family <- chemicals_info$ecModel[chemicals_info$class == x]
  # name standardization
  products_one_family <- products_one_family %>%
    str_replace_all(., "\\.mat", "") %>% str_replace_all(., "-", "") %>%
    str_replace_all(., "_", "") %>% str_replace_all(., ",", "") %>%
    str_replace_all(., "\\(", "") %>% str_replace_all(., "\\)", "") %>%
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
  gene_select_OE <- filter(gene_top_rank, OE_occurance >= length(metabolite_with_gene_targes) * threshold) # only choose the gene targets overexpressed in more than half of products in one families
  gene_select_KD <- filter(gene_top_rank, KD_occurance >= length(metabolite_with_gene_targes) * threshold) # only choose the gene targets overexpressed in more than half of products in one families
  gene_select_KO <- filter(gene_top_rank, KO_occurance >= length(metabolite_with_gene_targes) * threshold) # only choose the gene targets overexpressed in more than half of products in one families
  return(list(gene_select_OE, gene_select_KD, gene_select_KO))
}

targets_levels <- c('ecFSEOF','mech_validated','compatible')
levels_abrv    <- c('FSEOF','mechVal','compatible')

for (i in 1:length(targets_levels)){
  inputName <- paste('../results/targetsMatrix_',targets_levels[i],'.txt',sep='')
  input_prediction <- read.table(inputName, sep="\t", header = TRUE, stringsAsFactors = FALSE)
  newDF <- data.frame()
  families <- unique(chemicals_info$class)
  for (family in families){
    output00 <- find_common_gene_target(product_class = family, input_prediction0 = input_prediction,threshold)
    OE_out <- output00[[1]]
    KD_out <- output00[[2]]
    KO_out <- output00[[3]]
    newDF <- rbind(newDF,cbind(family,nrow(OE_out),nrow(KD_out),nrow(KO_out)))
  }
  colnames(newDF) <- c("family","OE","KD","KO")
  outputName <- paste('../results/targetsOverlap_',levels_abrv[i],'.txt',sep='')
  write.table(newDF,outputName,sep='\t',quote=FALSE,row.names =FALSE)
}
# here maybe we use the targetsMatrix_mech_validated to do the enrichment analysis
# output over-expressed gene targets
input_prediction <- read.table("../results/targetsMatrix_mech_validated.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
#newDF <- data.frame()
families <- unique(chemicals_info$class)

out_gene_list = data.frame(families=families, stringsAsFactors = FALSE)
gene_list_OE= c()
for (family in families){
  output00 <- find_common_gene_target(product_class = family, input_prediction0 = input_prediction,threshold)
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
  output00 <- find_common_gene_target(product_class = family, input_prediction0 = input_prediction,threshold)
  
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
analysisType <- c('OE','KD_KO')
sources <- c('KEGG','GO')
for (directionality in analysisType){
  # 1 analyze the result from OE
  result_dir <- paste('../results/gene_enrichment_',directionality,'_genes',sep='')
  subfile    <- list.files(result_dir)
  families1  <-families[families!="other"]
  if (directionality=='KD_KO'){
    families1  <-families1[families1!="protein"]
  }
  for (file0 in families1){
    print(file0)
    file_input <- paste(result_dir, paste(file0, ".txt", sep = ""), sep = "/")
    result_analysis <- readLines(file(file_input))
    result_analysis <-  result_analysis[!str_detect(result_analysis, "Enrichment Score")]
    result_analysis <-  result_analysis[str_detect(result_analysis, "\t")]
    file00 <- paste("update_",file0, ".txt", sep = "")
    write_lines(result_analysis, paste(result_dir, file00, sep = "/"))
    
    newfile <- read.table(paste(result_dir, file00, sep = "/"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    newfile$Count <- as.numeric(newfile$Count)#/as.numeric(newfile$List.Total)
    newfile$PValue <- as.numeric(newfile$Benjamini)
    newfile <- newfile[newfile$Category != "Category",]
    # initial filteration
    newfile <- filter(newfile, Count >1,PValue <= 0.05)#filter(newfile, Count >=4 & PValue <= 0.05)
    # remove duplicated term
    newfile <- newfile[!duplicated(newfile$Term),]
    # order
    Factor            <- newfile[order(newfile$Count,decreasing = TRUE),]
    newfile_go        <- Factor[grep('GO:',Factor$Term),]
    newfile_go$Term   <- substr(newfile_go$Term,12,nchar(newfile_go$Term))
    newfile_go$Term   <- factor(newfile_go$Term, levels=unique(newfile_go$Term))
    newfile_kegg      <- Factor[grep('sce',Factor$Term),]
    newfile_kegg$Term <- substr(newfile_kegg$Term,10,nchar(newfile_kegg$Term))
    newfile_kegg$Term <- factor(newfile_kegg$Term, levels=unique(newfile_kegg$Term))
    # plot
    for (dataSource in sources){
      if (dataSource=='KEGG'){dataF <- newfile_kegg}
      else{dataF <- newfile_go}
      if (directionality=='OE'){colour <- colors[11]}
      else{colour <- colors[3]}
      dataF$Count <- as.numeric(dataF$Count)/as.numeric(dataF$List.Total)
      
      fileName <- paste('../results/plots/geneSet_',dataSource, '_enrichment_',directionality,'_',file0,'.png',sep='')
      png(fileName,width=700, height=600)
      p1 <- ggplot(data=dataF , aes(x=Term, y=Count)) +
        geom_bar(stat="identity",fill=colour) +
        theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
        theme(panel.background = element_rect(fill = "NA")) +
        theme(panel.grid.major = element_line(colour = "grey90")) +
        theme(axis.line = element_line(size = 1, colour = "black")) +
        theme(legend.position = c(0.85, 0.2)) +
        theme(axis.text.y=element_text(size=18, family="Arial"),
              axis.text.x=element_text(size=12, family="Arial"),
              axis.title=element_blank(),
              legend.text = element_text(size=10, family="Arial"),
              plot.title = element_text(size=18, family="Arial")) +
        ylim(c(0,1)) +
        ggtitle(file0) #+
      plot(p1)
      dev.off()
    }
  }
}