#stoich_distance analysis
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
library(readr)

if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}

families <- c('amino acid','alkaloid','organic acid','protein','alcohol','terpene','fatty acids and lipids','flavonoid','aromatic','bioamine')
codes    <- c('_AA','_alk','_oAc','_pro','_alc','_ter','_FA','_fla','_aro','_bioAm')
#Load chemicals info
filename        <- '../data/chemicals_info.txt'
chemicals_info <- read_delim(filename,"\t", escape_double = FALSE, na = "NA",trim_ws = TRUE)
#reformat
chemicals_info$ecModel<- gsub('-','_',chemicals_info$ecModel)
#chemicals_info$ecModel<- gsub('(S)_reticuline','S_reticuline',chemicals_info$ecModel)
chemicals_info$ecModel<- gsub(',','_',chemicals_info$ecModel)
#chemicals_info$ecModel<- gsub('ecGlutamate.mat','ecGlutamate',chemicals_info$ecModel)
#chemicals_info$ecModel<- gsub('.mat','',chemicals_info$ecModel)
data <- read.csv('../results/met_precursors_turnovers_allChemicals.txt',sep='\t',stringsAsFactors = TRUE)
chemicals <- data$Var1
data <- data[,2:ncol(data)]
vector <- data[1,]
data <- scale(data,center=F,scale =t(vector))
chemicals <- c('G6P','R5P','G3P','PEP','CoA','SUC','F6P','E4P','3PG','PYR','OXO','OAC')
rownames(data) <- chemicals
#Isolate numeric values
t.data <- t(data)
#get a heatmap
colores <- cividis(100)
plotName <- '../results/met_turnover/met_turnover_heatmap.pdf'
pdf(plotName,width = 6, height = 20)

pheatmap(t.data,
                  method = c("pearson"),
                  clustering_method = "complete",
                  cluster_row = T,
                  cluster_col = T,
                  show_rownames = T,
                  show_colnames = T,
                  legend = T,
                  fontsize = 12,
                  color = colores)

dev.off()
df <- t(data)
#t.data <- t.data[,2:ncol(t.data)]
PCAdata   <- prcomp(t.data, center = FALSE,scale = FALSE,retx=TRUE)
t.data <- as.data.frame(t.data)


p         <- autoplot(PCAdata,data = t.data)
plotTitle <- paste('../results/met_turnover/PCA_allChems.png',sep='')
png(plotTitle,width = 600, height = 600)
plot(p)
dev.off()


