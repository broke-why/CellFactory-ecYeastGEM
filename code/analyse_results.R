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
#library(qgraph)
library(reshape2)


# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}
families <- c('amino acid','alkaloid','organic acid','protein','alcohol','terpene','fatty acids and lipids','flavonoid','aromatic','bioamine')
codes    <- c('_AA','_alk','_oAc','_pro','_alc','_ter','_FA','_fla','_aro','_bioAm')
#Load targets summary
filename        <- paste('../results/production_targets/targets_summary.txt',sep='')
targets_summary <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
#Get pie chart for families of chemicals
classes <- unique(targets_summary$chemClass)
counts  <- c()
for (i in 1:length(classes))
{
  counts <- c(counts,sum(targets_summary$chemClass==classes[i]))
}
#Get data frame for plotting
df <- data.frame(classes,counts)
df <- df[order(-counts),]
df$classes <- factor(df$classes,levels<-df$classes) 
df$perc <- round(df$counts*100/sum(df$counts),digits=1)
#Get pallete of colors
colourCount <- length(unique(df$classes))
getPalette <- colorRampPalette(brewer.pal(colourCount,"Paired"))

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
  
p <- ggplot(df, aes(x='',y=perc,fill=classes)) +
     geom_bar(colour = 'black',width = 1, stat = 'identity')
     pie <- p + coord_polar("y", start=0,direction=1)
     pie <- pie + blank_theme 
     pie <- pie +scale_fill_manual(values = getPalette(colourCount))#scale_fill_viridis(discrete = T, option = "E")
     pie <- pie + labs(fill = 'Chemical families')
     png('../results/plots/chemicalFamilies.png',width = 800, height = 800)
     plot(pie)
     dev.off()


#Analyse targets matrix
allTargetsMat <- read.csv('../results/targetsMatrix_compatible.txt',sep='\t',stringsAsFactors = TRUE)
allTargetsMat <- read.csv('../results/targetsMatrix_mech_validated.txt',sep='\t',stringsAsFactors = TRUE)
targetsMat    <- allTargetsMat
targetsMat <- targetsMat[rowSums(targetsMat[,5:ncol(targetsMat)])>0,]
#substitute values in matrix
targetsMat[targetsMat == 2] <- 20
targetsMat[targetsMat == 1] <- 10
targetsMat[targetsMat == 3] <- 30
targetsMat[targetsMat == 0] <- 1
targetsMat[targetsMat == 20] <- 0.5
targetsMat[targetsMat == 10] <- 0
targetsMat[targetsMat == 30] <- 2
targetsMat <- targetsMat[rowSums(targetsMat[,5:ncol(targetsMat)])>0,]
#hclust heatmap
plotTitle <- paste('../results/plots/targetsMatrix_ALL.png',sep='')
png(plotTitle,width = 5000, height = 13000)
rownames(targetsMat) <- targetsMat$shortNames
pheatmap(targetsMat[,5:ncol(targetsMat)],cluster_cols = T,cluster_rows = T, show_rownames = TRUE,scale='none',fontsize = 28)
dev.off()

#Dimensionality reduction
newDF <- targetsMat[,5:ncol(targetsMat)]
#newDF <- newDF[rowSums(newDF)>0,]
newDF <- as.data.frame(t(newDF))
newDF$extra <- rownames(newDF)
newDF<- newDF %>% separate(extra, c("chemical", "family"), "_fam_")
#newDF$family <- factor(newDF$family)
#Remove duplicate vectors
idxs <- which(!duplicated(newDF[,1:(ncol(newDF)-3)]))
newDF <- newDF[idxs,]
#Add product family info
newDF <- newDF[order(newDF$family),]
famLvls <- as.numeric(unique(factor(newDF$family)))
famLvls <- (unique(factor(newDF$family)))
famLvls <- famLvls[order(famLvls)]
newDF$family <- factor(newDF$family,levels = famLvls)
#PCA 
PCAdata   <- prcomp(newDF[,1:(ncol(newDF)-3)], center = TRUE,scale = FALSE,retx=TRUE)
p         <- autoplot(PCAdata,data = newDF,colour = 'family',size = 4)
plotTitle <- paste('../results/plots/PCA_allTargets.png',sep='')
png(plotTitle,width = 600, height = 600)
plot(p)
dev.off()
#tSNE
set.seed(18) # Set a seed if you want reproducible results
perplxty <- nrow(newDF)/length(famLvls)

tsne_out  <- Rtsne(newDF[,1:(ncol(newDF)-3)],dims=3,perplexity=perplxty) # Run TSNE
tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2],z = tsne_out$Y[,3],newDF$chemical,newDF$family)
#Define color palleter
colourCount2 <- length(unique(tsne_plot$family))
getPalette2  <- colorRampPalette(brewer.pal(colourCount, "Paired"))

colnames(tsne_plot)[ncol(tsne_plot)]<- 'family'
p <- plot_ly(x=tsne_plot$x, y=tsne_plot$y, z=tsne_plot$z,text =tsne_plot$newDF.chemical, type="scatter3d", mode="markers", color=tsne_plot$family,colors = getPalette2(colourCount))
saveWidget(p, "../results/plots/tsne_ALL.html", selfcontained = F, libdir = "lib")

filename  <- paste('../results/production_targets/targets_summary.txt',sep='')
targetsDF <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
for (j in 1:length(families)){
  chemClass <- families[j] 
  famCode   <- codes[j]
  targets_summary <- targetsDF
  if (nchar(chemClass)>1){
    targets_summary <- targets_summary[targets_summary$chemClass==chemClass,]
  }
  nCompounds <- nrow(targets_summary)
  if (nCompounds>=4){
    vector1 <- c(targets_summary$cand_1_OE,targets_summary$cand_2_OE,targets_summary$cand_3_OE)
    vector2 <- c(rep('step1',nCompounds),rep('step2',nCompounds),rep('step3',nCompounds))
    vector3 <- c(rep('OE',length(vector1)))
    df <- data.frame(vector1,vector2,vector3)
    
    vector1 <- c(targets_summary$cand_1_dR,targets_summary$cand_2_dR,targets_summary$cand_3_dR)
    vector2 <- c(rep('step1',nCompounds),rep('step2',nCompounds),rep('step3',nCompounds))
    vector3 <- c(rep('KD',length(vector1)))
    df2 <- data.frame(vector1,vector2,vector3)
    
    vector1 <- c(targets_summary$cand_1_del,targets_summary$cand_2_del,targets_summary$cand_3_del)
    vector2 <- c(rep('step1',nCompounds),rep('step2',nCompounds),rep('step3',nCompounds))
    vector3 <- c(rep('KO',length(vector1)))
    df3 <- data.frame(vector1,vector2,vector3)
    df <- rbind(df,df2,df3)
    colors <- cividis(11)
    #Generate box plots with all targets (KOs, KDs and OEs)
    p <-ggplot(df, aes(x=vector2, y=vector1, fill=vector3)) + geom_boxplot()
    p <- p + theme_bw(base_size = 2*12) + xlab('') +
    ylab('# of targets')+ylim(c(0,100))+labs(fill = 'Modification') +
    scale_fill_manual(values = c(colors[11],colors[6],colors[3]))
    plotTitle <- paste('../results/plots/targetDistributions',famCode,'.png',sep='')
    png(plotTitle,width = 600, height = 600)
    plot(p)
    dev.off()
  } 
}
