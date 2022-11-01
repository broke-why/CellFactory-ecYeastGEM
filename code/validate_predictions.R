library(matrixStats)
# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}
#families <- c('amino acid','alkaloid','organic acid','protein','alcohol','terpene','fatty acids and lipids','flavonoid','aromatic','bioamine','stilbenoids')
#codes    <- c('_AA','_alk','_oAc','_pro','_alc','_ter','_FAL','_fla','_aro','_bio','_stb')
#Load targets summary
filename        <- paste('../data/genetic_background.txt',sep='')
genetic_background <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
#classes <- unique(targets_summary$chemClass)
#Analyse targets matrix
allTargetsMat <- read.csv('../results/production_targets/targetsMatrix_L3.txt',sep='\t',stringsAsFactors = TRUE)
targetsMat <- allTargetsMat
cases      <- ncol(targetsMat) - 4
targetsMat <- targetsMat[rowSds(as.matrix(targetsMat[,5:ncol(targetsMat)]))!=0,]
