#Stat 540 project

#Quantile normalization and plots for second list of targets
#Loess normalization for both targets
#Taking mean expression for duplicate probes for both target lists
#Merging data from both platforms to give an object (with NA's where probe was not present in a platform)

library("limma")
library(preprocessCore)
library(plyr)
setwd("Project")

WorkDir <- "Working" #path of your working directory
RawDataDir <- "Working/GSE35988_RAW" #path of the decompressed RawFiles
AnalysisDir <- file.path(WorkDir, "Analysis") #Path of the Analysis directory

targets6480 <- readTargets("target_GPL6480.txt", path=AnalysisDir, row.names="GSMID", sep="\t")
targets6848 <- readTargets("target_GPL6848.txt", path=AnalysisDir, row.names="GSMID", sep="\t") 
x <- read.maimages(targets6480, path = RawDataDir, source="agilent.median", green.only=TRUE)
y <- read.maimages(targets6848, path = RawDataDir, source = "agilent.median", green.only = TRUE)

#Taking log2, quantile normalization
logx <- log2(x$E)
normDat <- as.data.frame(normalize.quantiles(as.matrix(logx)))
names(normDat) <- colnames(x)

logy <- log2(y$E)
normDat6848 <- as.data.frame(normalize.quantiles(as.matrix(logy)))
names(normDat6848) <- colnames(y)

#Loess normalization on quantile normalized expression, and adding probe names as a new column to avoid duplicate names issue with rownames later
norm2Datx <- normalizeCyclicLoess(normDat)
norm2Datx <- as.data.frame(norm2Datx)
norm2Datx <- cbind(probe = x$genes$ProbeName, norm2Datx)

norm2Daty <- normalizeCyclicLoess(normDat6848)
norm2Daty <- as.data.frame(norm2Daty)
norm2Daty <- cbind(probe = y$genes$ProbeName, norm2Daty)


#Preparing pre- and post-normalization plots for second targets file
  #Boxplot raw data
pdf('BoxplotRawIntensities_targets2.pdf')
boxplot(data.frame(log2(x$E)), 
        xlab = "Samples", 
        ylab = "log2 Intensities", 
        main = "Boxplot of Background Intensities (6848)")
dev.off()


plotnorm2Daty <- norm2Daty[,-1] # an object with rownames removed for plotting purposes
  #boxplot quantile normalized
pdf("boxplots_afterNorm_targets2.pdf")
boxplot(normDat6848, 
        xlab = "Samples", 
        ylab = "log2 Intensities", 
        main = "Boxplot of Background Intensities After Quantile Normalization (6848)")
dev.off()

  #boxplot loess normalized
pdf("boxplots_afterNorm2_targets2.pdf")
boxplot(plotnorm2Daty, 
        xlab = "Samples", 
        ylab = "log2 Intensities", 
        main = "Boxplot of Background Intensities After \n Quantile and Lowess Normalization (6848)")
dev.off()

  #MA plot before normalized
pdf("MAplot_beforeNorm_targets2.pdf")
for(i in 1:ncol(logy)){
  plotMA(logy, array = i)
}
dev.off()
  #MA plot after normalized
pdf("qnormalizedMAplots_targets2.pdf")
for(i in 1:ncol(normDat6848)){
  plotMA(normDat6848, array = i)
}
dev.off()


#Taking mean of all duplicate probes (using loess normalized data)
#This is pretty computationally intensive! Unless someone knows a beter way to do this, it might take a while if anyone runs this script

probeMeanx <- ddply(norm2Datx,"probe",numcolwise(mean))
probeMeany <- ddply(norm2Daty,"probe",numcolwise(mean))
rownames(probeMeanx) <- (probeMeanx$probe)
rownames(probeMeany) <- (probeMeany$probe)

#Confirming removed duplicates
table(duplicated(probeDatx$probe))
table(duplicated(probeMeanx$probe))

table(duplicated(probeDaty$probe))
table(duplicated(probeMeany$probe))


#merging data from both platforms
mergeTargets <- merge(probeMeanx, probeMeany, by=0, all=TRUE)
table(mergeTargets$Row.names == mergeTargets$probe.x) #Confirming that first column correctly contains all probe names for first platform
table(mergeTargets$Row.names == mergeTargets$probe.y) # Confirming that first column correctly contains all probe names for second platform

#removing redundant probe name columns, setting rownames properly
merged <- mergeTargets[,-2]
which(colnames(merged) == "probe.y") #(finding column number for list of probe names from second platform)
merged <- merged[,-90]

rownames(merged) <- merged[,1]
merged <- merged[,-1]
#### "merged" is combined quantile normalized and loess normalized log2 probe-level expression from each platform, with NA's where a probe is not present for a platform.

# Removing reduntant probe columns in unmerged data, writing merged and unmerged (second platform) data to files

write.table(merged, file = "merged_normalized_probelevel.txt", sep = "\t", row.names = TRUE, col.names = NA)

probeMeany <- probeMeany[,-1]
write.table(probeMeany, file = "normalized_data_targets2.txt", sep = "\t", row.names = TRUE, col.names = NA)
