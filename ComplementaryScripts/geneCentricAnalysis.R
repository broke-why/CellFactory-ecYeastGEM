#Gene-centric analysis of metabolic targets
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

# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}

#Theme for plotting
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),    #strip major gridlines
    panel.grid.minor = element_blank(),    #strip minor gridlines
    plot.title=element_text(size=26, face="bold"),
    legend.title = element_text(size=24,face="bold"),
    legend.text = element_text(size=18)
  )

#load matrix indicating the relation between gene targets and chemical products
allTargetsMat <- read.csv('../results/targetsMatrix_compatible.txt',sep='\t',stringsAsFactors = TRUE)
allTargetsMat <- read.csv('../results/targetsMatrix_mech_validated.txt',sep='\t',stringsAsFactors = TRUE)

direction  <- c('OE','KD','KO')
colors     <- cividis(11) 
threshold <- 10

productsNmr <- c()

#threshold <- 10
chemicals <- list(c(),c(),c())
for (action in direction){
  
  if (action=='OE'){
    targetsMat <- allTargetsMat
    targetsMat[targetsMat == 1] <- 0
    targetsMat[targetsMat == 2] <- 0
    targetsMat[targetsMat == 3] <- 1
    barColor <- colors[10]
    j <- 1
  }  
  if (action=='KD'){
    targetsMat <- allTargetsMat
    targetsMat[targetsMat == 1] <- 0
    targetsMat[targetsMat == 2] <- 1
    targetsMat[targetsMat == 3] <- 0
    barColor <- colors[6]
    j <- 2
  }
  if (action=='KO'){
    targetsMat <- allTargetsMat
    targetsMat[targetsMat == 1] <- 1
    targetsMat[targetsMat == 2] <- 0
    targetsMat[targetsMat == 3] <- 0
    barColor <- colors[2]
    j <- 3
  }
  
  matrix           <- targetsMat[,5:ncol(targetsMat)]
  rownames(matrix) <- targetsMat$shortNames
  colnames(matrix) <- colnames(targetsMat)[5:ncol(targetsMat)]
  matrix           <- matrix[rowSums(matrix)>0,]
  #Keep those genes that are targets
  matrix$suma <- rowSums(matrix)
  matrix <- matrix[order(matrix$suma,decreasing=TRUE),]
  top10  <- matrix[1:threshold,]
  sumas <- colSums(top10)
  newDF <- data.frame(genes = rownames(matrix),occurrence = matrix$suma)
  newDF <- newDF[1:threshold,]
  #get chemicals related to the top genes
  topGenes <- top10[1:threshold,]
  for (i in 1:threshold){
    colsumas <- colSums(topGenes)
    chemicals[[j]] <- c(chemicals[[j]],colnames(topGenes)[colsumas==threshold])
  }
  #generate bar plots for top genes
  p <- ggplot(newDF,aes(x= reorder(genes,-occurrence),occurrence))+geom_bar(stat ="identity",fill=barColor)+  
    theme_bw(base_size = 2*12)+xlab('Gene targets') +
    ylab('Occurences')#+ scale_y_continuous(breaks = (seq(0,5,by = 0.5)))
  plotTitle <- paste('../results/plots/top_', action,'_genes.png',sep='')
  png(plotTitle,width = 850, height = 600)
  plot(p)
  dev.off()
  #Generate stacked bar plots for top represented genes, indicating products family
  chemProds <- data.frame(product = colnames(matrix)[1:(ncol(matrix)-1)])
  chemProds <- chemProds %>% separate(product, c("product", "family"), "_fam_")
  geneNames <- c()
  gene_fam  <- c()
  for (k in 1:nrow(top10)){
    indexes <- which(top10[k,1:nrow(chemProds)]>0)
    if (length(indexes)>=1){
    for (l in indexes){
      geneNames <- c(geneNames,rownames(top10)[k])
      gene_fam  <- c(gene_fam,chemProds$family[l])
    }
    }
  }
  topGenesDF <- data.frame(geneNames,Family=gene_fam,stringsAsFactors = TRUE)
  topGenesDF$geneNames <- factor(topGenesDF$geneNames,levels = rownames(top10),ordered = TRUE)
  #plot
  getPalette  <- colorRampPalette(brewer.pal(5, "Set1"))
  colourCount <- length(unique(topGenesDF$Family))
  p <- ggplot(topGenesDF, aes(fill=Family, x=(geneNames))) + 
    geom_bar(position="stack", stat="count") + theme_bw(base_size = 2*12)+
    xlab('') + ylab('Number of chemicals') + scale_fill_manual(values = getPalette(colourCount))
  png(paste('../results/plots/topGenes_chemFam_',action,'.png',sep=''),width = 920, height = 600)
  plot(p)
  dev.off()
  #Get histogram for product-gene specificity
  matrix  <- targetsMat[,5:ncol(targetsMat)]
  sumas   <- rowSums(matrix)[rowSums(matrix)>0]
  geneNmr <- c()
  newDF <- data.frame(chemicals = rownames(matrix)[rowSums(matrix)>0],productNmr =sumas)
  p <- ggplot(newDF, aes(x=productNmr)) + geom_histogram(binwidth=1,fill=barColor)
  p <-  p + theme_bw(base_size = 2*12)+xlab('Number of products') + ylab('Number of gene targets')
  png(paste('../results/plots/target_specificity_',action,'.png',sep=''),width = 650, height = 600)
  plot(p)
  dev.off()
}
