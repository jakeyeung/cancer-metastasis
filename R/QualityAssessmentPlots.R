##########################################################
# SETTINGS #
###########################################################
WorkDir <- "/Data/Raunak/STAT540_CRPC"
RawDataDir <- "/Data/Raunak/rotproj_2012_prostate_marray/RawData/GSE3988_CRPC"
AnalysisDir <- file.path(WorkDir, "Analysis")
###########################################################

setwd() <- RawDataDir
library("limma")

#Read the target (sample file list) file to process
targets <- readTargets("target_GPL6480.txt", path=AnalysisDir, row.names="GSMID", sep="\t")

#Read channel intensity values from the data files specified in the target file
x <- read.maimages(targets, source="agilent.median", green.only=TRUE)

#Boxplot of log2 intensities for all the arrays to observe any suspicious background intensities.
pdf('BoxplotRawIntensities.pdf')
boxplot(data.frame(log2(x$E)), 
        xlab = "Samples", 
        ylab = "log2 Intensities", 
        main = "Boxplot of Background Intensities")
dev.off()


#MA plots to determine if normalization is necessary
pdf('MAPlotsAll.pdf')
for(i in 1:ncol(x)){
  plotMA(x, array = i)
}
dev.off()

  

