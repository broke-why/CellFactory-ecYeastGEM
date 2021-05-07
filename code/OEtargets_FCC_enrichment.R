library(ggplot2)
library(tidyverse)

# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}
source('plotVennDiagram.R')

# merge all the ECC result together
# datafile
ECC_dir <- "../results/ECC/"
res_dir <- "../results/"
all_product <- list.files(ECC_dir)
All_gene <- vector()
dir.create('../results/OE_FCC_enrichment')
df <- data.frame(stringsAsFactors = FALSE)
chem_info <- read.table('../data/chemicals_info.txt', header = TRUE, sep='\t',stringsAsFactors = FALSE,quote="")
chem_info$Name <- tolower(chem_info$Name)
for (x in all_product){
  print(x)
  #x <- "serine_ECCs.txt"
  chemical <- gsub('_ECCs.txt','_targets/compatible_genes_results.txt',x)
  fileName <-paste(res_dir,chemical,sep='')
  candidates <- read.table(fileName, header = TRUE, sep='\t',stringsAsFactors = FALSE,quote="")
  ss0 <- paste(ECC_dir, x, sep = "")
  print(ss0)
  ECCs <- read.table(ss0, header = TRUE, stringsAsFactors = FALSE)
  # remove the rows with zero ECC in low and high conditions
  ECCs <- filter(ECCs,CC_highGlc !=0)
  
  targets_ECC <- ECCs$names
  targets_OEs <- candidates$shortNames[candidates$actions=='OE']
  chemical <- gsub('_targets/compatible_genes_results.txt','',chemical)
  fileName <- paste('../results/OE_FCC_enrichment/',chemical,'.png',sep='')
  png(fileName,width=500, height=500)
  intersection <- plotVennDiagram(list(targets_ECC,targets_OEs),c('ECCs','OEs'),c('blue','red'),c(2.5,2.5,2.5),2,FALSE)
  dev.off()
  ratio_OEs <- length(intersection)/length(targets_OEs)
  ratio_ECC <- length(intersection)/length(targets_ECC)
  index <- grep(chemical,chem_info$Name,ignore.case=TRUE)
  if (length(index)>0){family <- chem_info$class[index]}
  else{family<- NA}
    
  
  df <- rbind(df,cbind(chemical,family,ratio_OEs,ratio_ECC))
}
write.table(df,'../results/OE_enrichment.txt',sep='\t',quote=FALSE)
matVals <- df[,3:4]
rownames(matVals)<- df$chemical
colnames(matVals)<- c('OEs','FCCs')
fileName <- '../results/plots/OE_enrichment.png'
png(fileName,width=800, height=(nrow(temp)/22)*750)
p <- pheatmap(matVals,color = cividis(11),cluster_cols = F,cluster_rows = T, show_rownames = TRUE,scale='none',fontsize = 20)
dev.off()

mean(as.double(df[,3]))
median(as.double(df[,3]))


