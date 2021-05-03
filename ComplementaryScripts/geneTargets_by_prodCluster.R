library(ggplot2)
library(scales)
library(viridis)
library(reshape)
library(RColorBrewer)
library(pheatmap)
library(Rtsne)
library(dplyr)
library(tidyr)
library(cluster)
library(ggfortify)
library(plotly)
library(htmlwidgets)
library(readr)
library(reshape2)

# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}
#Create dir for results
dir.create('../results/cluster_strains')
families <- c('amino acid','alkaloid','organic acid','protein','alcohol','terpene','fatty acids and lipids','flavonoid','aromatic','bioamine')
codes    <- c('_AA','_alk','_oAc','_pro','_alc','_ter','_FA','_fla','_aro','_bioAm')
#Load chemicals cluster info
filename      <- '../ComplementaryData/product_clusters.txt'
prod_clusters <- read_delim(filename,"\t", escape_double = FALSE, na = "NA",trim_ws = TRUE)
clusters      <- unique(prod_clusters$cluster)
#reformat
prod_clusters$ecModel<- gsub('-','_',prod_clusters$ecModel)
prod_clusters$ecModel<- gsub(',','_',prod_clusters$ecModel)
prod_clusters$ecModel<- substr(prod_clusters$ecModel,3,(nchar(prod_clusters$ecModel)-4))
prod_clusters$ecModel<- gsub('_','',prod_clusters$ecModel)
prod_clusters$ecModel<-tolower(prod_clusters$ecModel)
prod_clusters$ecModel<-gsub('[0-9]+', '', prod_clusters$ecModel)
#
allTargetsMat <- read.csv('../results/targetsMatrix_mech_validated.txt',sep='\t',stringsAsFactors = TRUE)
targetsMat    <- allTargetsMat
rownames(targetsMat) <- targetsMat$shortNames
targetsMat <- targetsMat[rowSums(targetsMat[,5:ncol(targetsMat)])>0,]
newDF <- targetsMat[,5:ncol(targetsMat)]
newDF <- as.data.frame(t(newDF))
names <- data.frame(names=rownames(newDF))
names <- names %>% separate(names, c("chemical", "family"), "_fam_")
rownames(newDF) <- names$chemical
actions <- c('KO','KD','OE')
strains_summary <- data.frame()
for (cluster in clusters){
  #cluster <- cluster[1:2]
  clusterData <- prod_clusters[prod_clusters$cluster==cluster,]
  indexes   <- match(clusterData$ecModel,rownames(newDF))
  clust_mat <- newDF[indexes,]
  nChems    <- nrow(clusterData)
  newStrain <- data.frame(stringsAsFactors = FALSE)
  products <- rownames(clust_mat)
  products <- paste(products,collapse = ' // ')
  newRow <- data.frame(cluster,length(indexes),products)
  for (i in 1:3){
    tempMat <- clust_mat
    tempMat[tempMat != i] <- 0
    tempMat[tempMat == i] <- 1
    indexes2  <- which(colSums(tempMat)==nChems)
    genes     <- colnames(clust_mat)[indexes2]
    action    <- rep(actions[i],length(indexes2))
    newStrain <- rbind(newStrain,cbind(genes,action))
    newRow <- cbind(newRow,length(indexes2))
  }
  strains_summary <- rbind(strains_summary,newRow)
  fileName <- paste('../results/cluster_strains/',cluster,'_strain.txt',sep='')
  write.table(newStrain,fileName,quote=FALSE,row.names = FALSE,sep = '\t')
}
colnames(strains_summary) <- c('cluster','n_prods','chemicals','KOs','KDs','OEs')
fileName <- '../results/cluster_strains/cluster_strains_summary.txt'
write.table(strains_summary,fileName,quote=FALSE,row.names = FALSE,sep = '\t')
strains_summary <- strains_summary[which(strains_summary$n_prods>2),]
newDF <- data.frame(OE=strains_summary$OEs,KD=strains_summary$KDs,KO=strains_summary$KOs)
rownames(newDF) <- gsub('cluster_','c',strains_summary$cluster)
legLabels <- c('OE','KD','KO')
newDF <- t(newDF)
maxLim <- 20# max(newDF)
minLim <- 0
newDF <- rbind(rep(maxLim,ncol(newDF)),rep(0,ncol(newDF)),newDF)
newDF <- as.data.frame(newDF,stringsAsFactors = FALSE)

#plot spider plot
colors_border = c(rgb(0.8,0.6,0,0.8), rgb(0.4,0.4,0.40,0.8),rgb(0.1,0,0.8,0.8))
colors_in     = c(rgb(0.8,0.6,0,0.2), rgb(0.4,0.4,0.40,0.2),rgb(0.1,0,0.8,0.2))
plotName <- '../results/plots/panGenes_by_clusters.png'
png(plotName,width = 550, height = 500)
radarchart( newDF , axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(minLim,maxLim*100,(maxLim-minLim)/4), cglwd=1.5,
            #custom labels
            vlcex=2, calcex = 1.5)
#plot(p)
legend(x=1.1, y=1.1, legend = legLabels, bty = "n", pch=20 , col=colors_in , text.col = "black", cex=1.5, pt.cex=3)
dev.off()





