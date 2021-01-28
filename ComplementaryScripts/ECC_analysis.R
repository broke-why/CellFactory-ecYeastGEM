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


# merge all the ECC result together
# datafile
ECC_dir <- "../results/ECC/"
all_strain <- list.files(ECC_dir)
All_gene <- vector()

for (x in all_strain){
  print(x)
  #x <- "serine_ECCs.txt"
  ss0 <- paste(ECC_dir, x, sep = "")
  print(ss0)
  ECC0 <- read.table(ss0, header = TRUE, stringsAsFactors = FALSE)
  # remove the rows with zero ECC in low and high conditions
  ECC1 <- filter(ECC0, CC_lowGlc !=0 | CC_highGlc !=0)
  gene0 <- ECC1$genes
  All_gene <- c(All_gene, gene0)
}

All_gene_unique <- unique(All_gene)



# creat a dataframe
# low ECC
ECC_low_df <- data.frame(gene=All_gene_unique, stringsAsFactors = FALSE)
for (x in all_strain){
  print(x)
  #x <- "serine_ECCs.txt"
  ss0 <- paste(ECC_dir, x, sep = "")
  print(ss0)
  ECC0 <- read.table(ss0, header = TRUE, stringsAsFactors = FALSE)
  # remove the rows with zero ECC in low and high conditions
  ECC1 <- filter(ECC0, CC_lowGlc !=0 | CC_highGlc !=0)
  # product name
  pp <- str_replace(x, "_ECCs.txt", "")
  ECC_low_df[[pp]] <- getSingleReactionFormula(ECC1$CC_lowGlc, ECC1$genes, ECC_low_df$gene)
  ECC_low_df[[pp]] <- as.numeric(ECC_low_df[[pp]])
}



# High ECC
ECC_high_df <- data.frame(gene=All_gene_unique, stringsAsFactors = FALSE)
for (x in all_strain){
  print(x)
  #x <- "serine_ECCs.txt"
  ss0 <- paste(ECC_dir, x, sep = "")
  print(ss0)
  ECC0 <- read.table(ss0, header = TRUE, stringsAsFactors = FALSE)
  # remove the rows with zero ECC in low and high conditions
  ECC1 <- filter(ECC0, CC_lowGlc !=0 | CC_highGlc !=0)
  # product name
  pp <- str_replace(x, "_ECCs.txt", "")
  ECC_high_df[[pp]] <- getSingleReactionFormula(ECC1$CC_lowGlc, ECC1$genes, ECC_high_df$gene)
  ECC_high_df[[pp]] <- as.numeric(ECC_high_df[[pp]])
}


# cluster analysis
ECC_low_df1 <- as.data.frame(t(ECC_low_df), stringsAsFactors = FALSE)
colnames(ECC_low_df1) <- ECC_low_df1[1,]
ECC_low_df1 <-  ECC_low_df1[-c(1),]
ECC_low_df1[] <- sapply(ECC_low_df1, as.numeric)
ECC_low_df1[is.na(ECC_low_df1)] <- 0
# add product annotation information
product_df <- data.frame(product=rownames(ECC_low_df1), stringsAsFactors = FALSE)
chemicals_info <- read_excel("../ComplementaryData/chemicals_info.xlsx")
chemicals_info$Name0 <- str_replace_all(chemicals_info$ecModel, ".mat", "")
chemicals_info$Name0 <- str_replace_all(chemicals_info$Name0, "^ec", "")
chemicals_info$Name0 <- str_to_lower(chemicals_info$Name0)

product_df$class <- getSingleReactionFormula(chemicals_info$class,chemicals_info$Name0,product_df$product)
# it is found some products are not grouped
# also need a unique name of product
ECC_low_df2 <- ECC_low_df1
ECC_low_df2$product <- product_df$product
ECC_low_df2$class <- product_df$class


library(ggfortify)
autoplot(prcomp(ECC_low_df1), data = ECC_low_df2, colour = 'class') +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=20, family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(legend.text=element_text(size=15),
        legend.title =element_text(size=15))

