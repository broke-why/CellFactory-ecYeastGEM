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
families <- c('amino acid','alkaloid','organic acid','protein','alcohol','terpene','fatty acids and lipids','flavonoid','aromatic','bioamine')
codes    <- c('_AA','_alk','_oAc','_pro','_alc','_ter','_FA','_fla','_aro','_bioAm')
pathWays <- c('Oxidative phosphorylation','Glycolysis','TCA cycle','Pentose phosphate pathway')
#load enzyme info
filename <- paste('../complementaryData/enzymeTable.txt',sep='')
enzTable <- read.csv(filename,sep = '\t',stringsAsFactors = FALSE)
#load KEGG pathways info
filename <- paste('../complementaryData/keggPathways.txt',sep='')
keggDF <- read.csv(filename,sep = '\t',stringsAsFactors = FALSE)
#Load targets summary
filename        <- paste('../results/targets_summary.txt',sep='')
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
getPalette <- colorRampPalette(brewer.pal(5, "Set1"))
colourCount <- length(unique(df$classes))

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
  
p <- ggplot(df, aes(x='',y=perc,fill=classes))+
  geom_bar(width = 1, stat = 'identity')
pie <- p + coord_polar("y", start=0,direction=1)
pie <- pie + blank_theme 
pie <- pie +scale_fill_manual(values = getPalette(colourCount))#scale_fill_viridis(discrete = T, option = "E")
pie <- pie + labs(fill = 'Chemical families')
png('../results/plots/chemicalFamilies.png',width = 800, height = 800)
plot(pie)
dev.off()


#Get heatmap for targets matrix
allTargetsMat <- read.csv('../results/targetsMatrix_compatible.txt',sep='\t',stringsAsFactors = TRUE)
#allTargetsMat <- read.csv('../results/targetsMatrix_mechValidated.txt',sep='\t',stringsAsFactors = TRUE)

targetsMat    <- allTargetsMat
targetsMat[targetsMat == 2] <- 0
targetsMat[targetsMat == 1] <- 0
targetsMat[targetsMat == 3] <- 1
targetsMat <- targetsMat[rowSums(targetsMat[,5:ncol(targetsMat)])>0,]


plotTitle <- paste('../results/plots/targetsMatrix_OE.png',sep='')
png(plotTitle,width = 5000, height = 13000)
rownames(targetsMat) <- targetsMat$shortNames
pheatmap(targetsMat[,5:ncol(targetsMat)],cluster_cols = T,cluster_rows = T, show_rownames = TRUE,scale='none',fontsize = 28)
dev.off()

#Dimensionality reduction
newDF <- targetsMat[,5:ncol(targetsMat)]
newDF <- newDF[rowSums(newDF)>0,]
newDF <- as.data.frame(t(newDF))
newDF$extra <- rownames(newDF)
newDF<- newDF %>% separate(extra, c("chemical", "family"), "_fam_")
#newDF$family <- factor(newDF$family)
idxs <- which(!duplicated(newDF[,1:(ncol(newDF)-3)]))
newDF <- newDF[idxs,]
newDF <- newDF[order(newDF$family),]
famLvls <- as.numeric(unique(factor(newDF$family)))
famLvls <- (unique(factor(newDF$family)))

famLvls <- famLvls[order(famLvls)]
newDF$family <- factor(newDF$family,levels = famLvls)

PCAdata  <- prcomp(newDF[,1:(ncol(newDF)-3)], center = TRUE,scale = FALSE,retx=TRUE)
p        <- autoplot(PCAdata,data = newDF,colour = 'family',size = 4)
p

set.seed(18) # Set a seed if you want reproducible results
tsne_out  <- Rtsne(newDF[,1:(ncol(newDF)-3)],dims=3,perplexity=5)#nrow(newDF)/length(famLvls)) # Run TSNE
tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2],z = tsne_out$Y[,3],newDF$chemical,newDF$family)
colnames(tsne_plot)[ncol(tsne_plot)]<- 'family'
plot_ly(x=tsne_plot$x, y=tsne_plot$y, z=tsne_plot$z,text =tsne_plot$newDF.chemical, type="scatter3d", mode="markers", color=tsne_plot$family)
# fig <- fig %>% add_annotations(x = tsne_plot$x,
#                                y = tsne_plot$y,
#                                text = tsne_plot$newDF.chemical,
#                                xref = "x",
#                                yref = "y")
# fig
#ggplot(tsne_plot) + geom_point(aes(x=x, y=y,color=family))
## Show the objects in the 2D tsne representation

filename  <- paste('../results/targets_summary.txt',sep='')
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
    plotTitle <- paste('../results/plots/allTargets',famCode,'.png',sep='')
    png(plotTitle,width = 600, height = 600)
    plot(p)
    dev.off()
    
    #Get all deletion targets for chemical family
    filename    <- paste('../results/all_deletions.txt',sep='')
    all_deletions <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
    if (nchar(chemClass)>1){
      all_deletions <- all_deletions[all_deletions$chem_class_del==chemClass,]
    }
    uniqueDels <- unique(all_deletions$del_targets)
    counts <- c()
    for (i in 1:length(uniqueDels)){counts <- c(counts,sum(all_deletions$del_targets==uniqueDels[i]))}
    df <- data.frame(uniqueDels,counts,stringsAsFactors = FALSE)
    df$counts <- df$counts/nCompounds
    dt <- df[order(-counts),]
    dt$uniqueDels <- factor(dt$uniqueDel, levels<- dt$uniqueDel)
    dt <- dt[1:10,]
    #Plot top targets for deletions as barplot
    p <- ggplot(data=dt, aes(x=uniqueDels,y=(counts))) + geom_bar(stat='identity',fill=colors[2])+ 
      theme_bw(base_size = 2*12)+xlab('Gene targets') +
      ylab('Relative frequency')+ylim(c(0,1))#+ scale_y_continuous(breaks = (seq(0,5,by = 0.5)))
    plotTitle <- paste('../results/plots/topKOs',famCode,'.png',sep='')
    png(plotTitle,width = 950, height = 600)
    plot(p)
    dev.off()
    
    #Get all deletion targets for chemical family
    filename      <- paste('../results/all_downRegs.txt',sep='')
    all_deletions <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
    if (nchar(chemClass)>1){
      all_deletions <- all_deletions[all_deletions$chem_class_dR==chemClass,]
    }
    if (nrow(all_deletions)>0){
      uniqueDels <- unique(all_deletions$dR_targets)
      counts <- c()
      for (i in 1:length(uniqueDels)){counts <- c(counts,sum(all_deletions$dR_targets==uniqueDels[i]))}
      df <- data.frame(uniqueDels,counts,stringsAsFactors = FALSE)
      df$counts <- df$counts/nCompounds
      dt <- df[order(-counts),]
      dt$uniqueDels <- factor(dt$uniqueDel, levels<- dt$uniqueDel)
      dt <- dt[1:10,]
      #Plot top targets for deletions as barplot
      p <- ggplot(data=dt, aes(x=uniqueDels,y=(counts))) + geom_bar(stat='identity',fill=colors[6])+ 
        theme_bw(base_size = 2*12)+xlab('Gene targets') +
        ylab('Relative frequency')+ylim(c(0,1))#+ scale_y_continuous(breaks = (seq(0,5,by = 0.5)))
      plotTitle <- paste('../results/plots/topKDs',famCode,'.png',sep='')
      png(plotTitle,width = 950, height = 600)
      plot(p)
      dev.off()
    }
    
    #repeat for OE targets
    filename    <- paste('../results/all_OEs.txt',sep='')
    all_OEs <- read.delim(filename,sep='\t',stringsAsFactors = FALSE)
    if (nchar(chemClass)>1){
      all_OEs <- all_OEs[all_OEs$chem_class_OE==chemClass,]
    }
    uniqueOEs <- unique(all_OEs$OE_targets)
    counts <- c()
    for (i in 1:length(uniqueOEs)){counts <- c(counts,sum(all_OEs$OE_targets==uniqueOEs[i]))}
    df <- data.frame(uniqueOEs,counts,stringsAsFactors = FALSE)
    df$counts <- df$counts/nCompounds
    dt <- df[order(-counts),]
    dt$uniqueOEs <- factor(dt$uniqueOEs, levels<- dt$uniqueOEs)
    dt <- dt[1:10,]
    
    p <- ggplot(data=dt, aes(x=uniqueOEs,y=(counts))) + geom_bar(stat='identity',fill=colors[10])+ 
      theme_bw(base_size = 2*12)+xlab('Gene targets') +
      ylab('Relative frequency')+ylim(c(0,1))#+ scale_y_continuous(breaks = (seq(0,5,by = 0.5)))
    plotTitle <-paste('../results/plots/topOEs',famCode,'.png',sep='')
    png(plotTitle,width = 950, height = 600)
    plot(p)
    dev.off()
    #Get Pathways DF
    KOs <- c()
    OEs <- c()
    for (pathway in pathWays){
      numero <- length(grep(pathway,all_deletions$subSystems_del))
      KOs <- c(KOs,numero)
      numero <- length(grep(pathway,all_OEs$subSystems_OE))
      OEs <- c(OEs,numero)
    }
    uniqueChem <- unique(all_OEs$chem_comp_OE)
    # counts  <- c()
    # actions <- c()
    # pathVec <- c()
    # for (k in 1:length(pathWays)){
    #   for (l in 1:length(uniqueChem)){
    #     #isolate data for compound
    #     idxs    <- grep(uniqueChem[l],all_OEs$chem_comp_OE)
    #     tempOEs <- all_OEs[idxs,]
    #     idxs    <- grep(uniqueChem[l],all_deletions$chem_comp_del)
    #     tempDel <- all_deletions[idxs,]
    #     #get specific targets
    #     targetsOE  <- tempOEs$OE_targets
    #     targetsDel <- tempDel$del_targets
    #     numero <- length(grep(pathWays[l],tempOEs$subSystems_OE))/length(grep(pathWays[l],enzTable$subSystems))
    #     counts  <- c(counts,numero) 
    #     actions <- c(actions,'OEs')
    #     
    #     numero <- length(grep(pathWays[l],tempDel$subSystems_del))/length(grep(pathWays[l],enzTable$subSystems))
    #     counts  <- c(counts,numero) 
    #     actions <- c(actions,'KOs')
    #   }
    #   pathway <- rep(pathWays[k],length(uniqueChem))
    #   pathVec <- c(pathVec,pathway)
    # }
    # pathDF  <- data.frame(counts,actions,pathVec)
    # pathCodes <- c('OxPhos','Glyc','TCA','PPP')
    # for (i in 1:length(pathWays)){
    #   pathDF$pathVec <- gsub(pathWays[i],pathCodes[i],pathDF$pathVec)
    # }
   # p<-ggplot(pathDF, aes(x=pathVec, y=counts, fill=actions)) +
    #   geom_boxplot()
    # p <- p+theme_bw(base_size = 2*12) + xlab('') +
    #   labs(fill = 'Modification')+
    #   ylab('Genes fraction (targets)')+ylim(c(0,0.15))
    # plotTitle <- paste('../results/plots/PathwayEnrichment_',famCode,'.png',sep='')
    # png(plotTitle,width = 600, height = 600)
    # plot(p)
    # dev.off()
  } 
}
