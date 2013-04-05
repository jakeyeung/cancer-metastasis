## Stat 540 Project## 
## Quantile Normalization ## 
setwd("~/Desktop/Grad/Stat540/Project")
###########################################################
# SETTINGS #
###########################################################
WorkDir <- "/Users/melissalee89/Desktop/Grad/Stat540/Project" #path of your working directory
RawDataDir <- "/Users/melissalee89/Desktop/Grad/Stat540/Project/Data/samples" #path of the decompressed RawFiles
AnalysisDir <- file.path(WorkDir, "Data") #Path of the Analysis directory
###########################################################

library(preprocessCore) # normalize.quantiles()
library("limma") #load limma

#Read the target (sample file list) file to process
targets <- readTargets("target_GPL6480.txt", path=AnalysisDir, row.names="GSMID", sep="\t")

#Read channel intensity values (median intensity value) from the data files specified in the target file
#if you want to extract the mean intensity rather than median the use source="agilent.mean"
x <- read.maimages(targets, source="agilent.median", path =RawDataDir, green.only=TRUE)

nrow(x)
#45015

logx <- log2(x$E) # log transform the data 
str(logx)

# quantile normalize the data
normDat <- normalize.quantiles(as.matrix(logx))
normDat <- as.data.frame(normDat) # make it a data frame
names(normDat) <- colnames(x) 

# MA plots pre-normalization 
#pdf("MAplots.pdf")
for(i in 1:ncol(logx)){
  plotMA(logx, array = i)
}
#dev.off()

# MA plots after normalization
pdf("Plots/qnormalizedMAplots.pdf")
for(i in 1:ncol(normDat)){
  plotMA(normDat, array = i)
 }
dev.off()

# boxplots after normalization
pdf("boxplots_afterNorm.pdf")
boxplot(normDat, 
        xlab = "Samples", 
        ylab = "log2 Intensities", 
        main = "Boxplot of Background Intensities After Quantile Normalization")
dev.off()

# lowess normalization
norm2Dat <- normalizeCyclicLoess(normDat)
head(norm2Dat)

pdf("Plots/boxplots_afterNorm2.pdf")
boxplot(norm2Dat, 
        xlab = "Samples", 
        ylab = "log2 Intensities", 
        main = "Boxplot of Background Intensities After \n Quantile and Lowess Normalization")
dev.off()

pdf("Plots/qnormalized2MAplots.pdf")
for(i in 1:ncol(normDat)){
  plotMA(norm2Dat, array = i)
}
dev.off()

write.table(norm2Dat, "NormalizedData.txt", sep = "\t", 
            row.names = T, col.names = NA)