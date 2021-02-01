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

# PCA analysis
# it seems there is no good classification based on product families?
library(ggfortify)
autoplot(prcomp(ECC_low_df1), data = ECC_low_df2, colour = 'class') +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=20, family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(legend.text=element_text(size=15),
        legend.title =element_text(size=15))



# other statistical analysis
InfluencedProducts <- c()
for ( gene in All_gene_unique ){
  print(gene)
  product_info <- ECC_low_df2[, gene]
  product_info <- product_info[product_info > 0 | product_info < 0]
  product_num <- length(product_info)
  InfluencedProducts <- c(InfluencedProducts, product_num)
}

product_num_gene <- data.frame(gene=All_gene_unique, products=InfluencedProducts, stringsAsFactors = FALSE)





# plot
# choose the top 10 genes which could affect most products in ECC analysis
product_num_gene_top_10 <- product_num_gene[product_num_gene$products > 52,]
Factor <- product_num_gene_top_10
Factor <- Factor[order(Factor$products,decreasing = TRUE),]
product_num_gene_top_10$gene <-factor(product_num_gene_top_10$gene, levels=Factor$gene)

ggplot(data=product_num_gene_top_10, aes(x=gene, y=products)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=10, family="Arial"),
        axis.title=element_text(size=12,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') #+
#theme(panel.background = element_rect(fill = "white", color="black", size = 1)) 

# general analysis
product_num_gene$group <- NA
product_num_gene$group[product_num_gene$products==1] <- "a. 1 product"
product_num_gene$group[product_num_gene$products > 1 &  product_num_gene$products <=3] <- "b. 2-3 products"
product_num_gene$group[product_num_gene$products > 3 &  product_num_gene$products <=5] <- "c. 4-5 products"
product_num_gene$group[product_num_gene$products > 5 &  product_num_gene$products <=11] <- "d. 6-11 products"
product_num_gene$group[product_num_gene$products > 11 &  product_num_gene$products <= 19] <- "e. 12-19 products"
product_num_gene$group[product_num_gene$products > 19] <- "over 20 products"

product_num_gene1 <- product_num_gene[!is.na(product_num_gene$group),]

product_num_df <- as.data.frame(table(product_num_gene1$group), stringsAsFactors = FALSE)
colnames(product_num_df) <- c("Group", "Number")
ggplot(data=product_num_df, aes(x=Group, y=Number)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=10, family="Arial"),
        axis.title=element_text(size=12,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') #+
#theme(panel.background = element_rect(fill = "white", color="black", size = 1)) 

# heatmap
# choose the top 10 genes and prepare a heatmap
top_10_genes <- as.character(product_num_gene_top_10$gene)
ECC_low_df1_top10 <- ECC_low_df1[, top_10_genes]
# randome choose 20 products
ECC_low_df1_top10_20products <- ECC_low_df1_top10[sample(nrow(ECC_low_df1_top10), 20), ]
library(pheatmap)
pheatmap(ECC_low_df1_top10_20products,
         method = c("pearson"),
         clustering_method = "complete",
         treeheight_row = 40,
         treeheight_col = 40,
         cluster_row = FALSE,
         cluster_col = FALSE,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize = 12,
         color = colorRampPalette(c("white", "SandyBrown", "firebrick3"))(100))



