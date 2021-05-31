#production_capabilities
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

filename        <- paste('../results/production_capabilities/prodCapabilities_allChemicals.txt',sep='')
df <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
df$yield2 <- df$prodRate_ec/df$cFlux
df[, "maxYield"] <- apply(df[, c(8,12)], 1, max)
df$maxYield    <- 1/df$maxYield
df$protYield <- df$Pburden/df$prodRate_ec

p <- ggplot(df, aes(x=maxYield,y=protYield,color=family,shape=type)) +
  geom_point() + scale_x_log10() + scale_y_log10() 
p <- p + theme_bw() 
#png('../results/plots/chemicalFamilies.png',width = 800, height = 800)
plot(p)
#dev.off()