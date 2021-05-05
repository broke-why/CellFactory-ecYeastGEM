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
source('plotVennDiagram.R')

prod_capacity <- read.table("../results/production_capacity_comparison.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
prod_capacity 

colBioY <- grepl('bioY',colnames(prod_capacity))
colproY <- grepl('proY',colnames(prod_capacity))
colproR <- grepl('proR',colnames(prod_capacity))

df1 <- prod_capacity[,colBioY]
df1 <- df1/df1[,1]
rownames(df1) <- prod_capacity$chemical
colnames(df1) <- c('WT','data','pred')
fileName <- '../results/plots/bioYield_comparison.png'
png(fileName,width=1000, height=950)
p <- pheatmap(df1,color = cividis(11),cluster_cols = F,cluster_rows = T, show_rownames = TRUE,scale='none',fontsize = 20)
dev.off()

df2 <- prod_capacity[,colproY]
df2 <- df2/df2[,1]
rownames(df2) <- prod_capacity$chemical
colnames(df2) <- c('WT','data','pred')
fileName <- '../results/plots/productYield_comparison.png'
png(fileName,width=1000, height=950)
p <- pheatmap(df2,color = cividis(11),cluster_cols = F,cluster_rows = T, show_rownames = TRUE,scale='none',fontsize = 20)
dev.off()

df3 <- prod_capacity[,colproR]
df3 <- df3/df3[,1]
rownames(df3) <- prod_capacity$chemical
colnames(df3) <- c('WT','data','pred')
fileName <- '../results/plots/productRate_comparison.png'
png(fileName,width=1000, height=950)
p <- pheatmap(df3,color = cividis(11),cluster_cols = F,cluster_rows = T, show_rownames = TRUE,scale='none',fontsize = 20)
dev.off()
threshold <- 1 - (1E-4)
temp_pred <- rowMeans(data.frame(df1[,3],df2[,3],df3[,3]))
temp_data <- rowMeans(data.frame(df1[,2],df2[,2],df3[,2]))
metric  <- (temp_pred>threshold) 
indexes <-metric>threshold & temp_pred>temp_data
topChemicals <- data.frame(temp_data,temp_pred)
rownames(topChemicals) <- prod_capacity$chemical
topChemicals <- topChemicals[indexes,]

fileName <- '../results/plots/productYield_Rate_comparison.png'
png(fileName,width=1000, height=950)
p <- pheatmap(topChemicals,color = cividis(11),cluster_cols = F,cluster_rows = T, show_rownames = TRUE,breaks=seq(0, 2,0.2),scale='none',fontsize = 20)
dev.off()

indexes <-temp_pred>1
temp_pred <- data.frame(prod_capacity$chemical[indexes])
indexes <-temp_data>1
temp_data <- data.frame(prod_capacity$chemical[indexes])
fileName <- '../results/plots/topChemicals_pred_vs_data.png'
png(fileName,width=500, height=500)
plotVennDiagram(list(temp_pred[,1],temp_data[,1]),c('pred','exp'),c('blue','red'),c(2.5,2.5,2.5),2,FALSE)
dev.off()