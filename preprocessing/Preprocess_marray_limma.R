#Preprocess the Agilent GeneExpr Raw Data using limma

###########################################################
# SETTINGS #
###########################################################
WorkDir <- "/Dropbox/STAT540_CRPC" #path of your working directory
RawDataDir <- "/Data/Raunak/rotproj_2012_prostate_marray/RawData/GSE3988_CRPC" #path of the decompressed RawFiles
AnalysisDir <- file.path(WorkDir, "Analysis") #Path of the Analysis directory
###########################################################

setwd() <- RawDataDir #This step can be omitted
library("limma") #load limma

#Read the target (sample file list) file to process
targets <- readTargets("target_GPL6480.txt", path=AnalysisDir, row.names="GSMID", sep="\t")

#Read channel intensity values (median intensity value) from the data files specified in the target file
#if you want to extract the mean intensity rather than median the use source="agilent.mean"
x <- read.maimages(targets, source="agilent.median", green.only=TRUE)

#To view structure of the data
str(x)

#To view the summary of the intensity distribution of samples
summary(x$E)

#The value of the median intensity are stored in
x$E

#To get the Gene level annotation for each probe
Genes <- x$genes$GeneName


#Now we will have to look at the data quality and use the appropriate background correction & normalization steps.