library(ggplot2)
library(scales)
library(viridis)
library(reshape)
library(RColorBrewer)
library(pheatmap)
library(Rtsne)
library(tibble)
library(dplyr)
library(tidyr)
library(cluster)
library(ggfortify)
library(plotly)
library(htmlwidgets)
library(reshape2)


# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}
families <- c('amino acid','alkaloid','organic acid','protein','alcohol','terpene','fatty acids and lipids','flavonoid','aromatic','bioamine','stilbenoids')
codes    <- c('_AA','_alk','_oAc','_pro','_alc','_ter','_FAL','_fla','_aro','_bio','_stb')
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
#Analyse targets matrix
allTargetsMat <- read.csv('../results/production_targets/targetsMatrix_L3_discrete.txt',sep='\t',stringsAsFactors = TRUE)
#allTargetsMat <- read.csv('../results/production_targets/targetsMatrix_mech_validated.txt',sep='\t',stringsAsFactors = TRUE)
targetsMat    <- allTargetsMat
targetsMat <- targetsMat[rowSums(targetsMat[,5:ncol(targetsMat)])>0,]
cases <- ncol(targetsMat) - 4
idxs  <- which(rowSums(targetsMat[,5:ncol(targetsMat)])!=cases)
reduced_matrix <- targetsMat[idxs,]
#substitute values in matrix
# targetsMat[targetsMat == 2] <- 20
# targetsMat[targetsMat == 1] <- 10
# targetsMat[targetsMat == 3] <- 30
# targetsMat[targetsMat == 0] <- 1
# targetsMat[targetsMat == 20] <- 0.5
# targetsMat[targetsMat == 10] <- 0
# targetsMat[targetsMat == 30] <- 2
#hclust heatmap
plotTitle <- paste('../results/plots/targetsMatrix_changing.png',sep='')
png(plotTitle,width = 5000, height = 13000)
rownames(reduced_matrix) <- reduced_matrix$shortNames
pheatmap(reduced_matrix[,5:ncol(reduced_matrix)],cluster_cols = T,cluster_rows = T, show_rownames = TRUE,scale='none',fontsize = 28)
dev.off()

#Dimensionality reduction
newDF <- targetsMat[,5:ncol(targetsMat)]
genes <- targetsMat$shortNames
idxs  <- which(rowSums(targetsMat[,5:ncol(targetsMat)])!=cases)
#newDF <- newDF[rowSums(newDF)>0,]
newDF <- newDF[idxs,]
extra <- colnames(newDF)
genes <- genes[idxs]
rownames(newDF) <- genes
newDF <- as.data.frame(t(newDF))
newDF$extra <- extra
newDF<- newDF %>% separate(extra, c("chemical", "family"), "_fam_")
#Remove duplicate vectors
idxs       <- which(!duplicated(newDF[,1:(ncol(newDF)-3)]))
duplicates <- which(duplicated(newDF[,1:(ncol(newDF)-3)]))
duplicates <- newDF[duplicates,]
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
medianTargets <- c()
for (j in 1:length(families)){
  chemClass <- families[j] 
  famCode   <- 'ALL'#codes[j]
  targets_summary <- targetsDF
  #if (nchar(chemClass)>1){
  #  targets_summary <- targets_summary[targets_summary$chemClass==chemClass,]
  #}
  nCompounds <- nrow(targets_summary)
  if (nCompounds>=1){
    
    vector1 <- c(targets_summary$cand_1_OE,targets_summary$cand_2_OE,targets_summary$cand_3_OE)
    vector2 <- c(rep('Lvl1',nCompounds),rep('Lvl2',nCompounds),rep('Lvl3',nCompounds))
    vector3 <- c(rep('OE',length(vector1)))
    df      <- data.frame(vector1,vector2,vector3)
    values  <- cbind(mean(targets_summary$cand_1_OE),mean(targets_summary$cand_2_OE),mean(targets_summary$cand_3_OE))
    values  <- round(values,2)
    str     <- gsub('_','',famCode)
    type    <- 'OE'
    OErow   <- data.frame(str,values,type)
    
    vector1 <- c(targets_summary$cand_1_dR,targets_summary$cand_2_dR,targets_summary$cand_3_dR)
    vector2 <- c(rep('Lvl1',nCompounds),rep('Lvl2',nCompounds),rep('Lvl3',nCompounds))
    vector3 <- c(rep('KD',length(vector1)))
    df2     <- data.frame(vector1,vector2,vector3)
    values  <- cbind(mean(targets_summary$cand_1_dR),mean(targets_summary$cand_2_dR),mean(targets_summary$cand_3_dR))
    values  <- round(values,2)
    str     <- gsub('_','',famCode)
    type <- 'KD'
    KDrow   <- data.frame(str,values,type)
    
    vector1 <- c(targets_summary$cand_1_del,targets_summary$cand_2_del,targets_summary$cand_3_del)
    vector2 <- c(rep('Lvl1',nCompounds),rep('Lvl2',nCompounds),rep('Lvl3',nCompounds))
    vector3 <- c(rep('KO',length(vector1)))
    df3     <- data.frame(vector1,vector2,vector3)
    df      <- rbind(df,df2,df3)
    values  <- cbind(mean(targets_summary$cand_1_del),mean(targets_summary$cand_2_del),mean(targets_summary$cand_3_del))
    values  <- round(values,2)
    str     <- gsub('_','',famCode)
    type <- 'KO'
    KOrow   <- data.frame(str,values,type)
    
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
    medianTargets <- rbind(medianTargets,OErow,KDrow,KOrow)
  } 
}
colnames(medianTargets) <- c('family','Lvl1','Lvl2','LVl3','type')

actions <- c('OE','KD','KO')

for (action in actions){
temp <- medianTargets[which(medianTargets$type==action),1:(ncol(medianTargets)-1)]
temp<-melt(temp,id.vars="family")

p <- ggplot(temp,aes(x=family,y=value,fill=factor(variable)))+
  geom_bar(stat="identity",position="dodge")+theme_bw(base_size = 2*12) +
  scale_fill_manual(name='',values = c(colors[10],colors[6],colors[3]))+
  xlab("Family")+ylab("Mean number of targets per product")
plotTitle <- paste('../results/plots/meanTargets_AllFams',action,'.pdf',sep='')
pdf(plotTitle,width = 8.5, height = 7)
plot(p)
dev.off()
}
write.table(medianTargets,'../results/processed_results/average_targets_number.txt',quote=FALSE,sep='\t',row.names = FALSE)


