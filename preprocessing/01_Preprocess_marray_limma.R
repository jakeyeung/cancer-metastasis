##########################################################
# SETTINGS #
###########################################################
WorkDir <- "/Data/Raunak/STAT540_CRPC"
RawDataDir <- "/Data/Raunak/rotproj_2012_prostate_marray/RawData/GSE3988_CRPC"
AnalysisDir <- file.path(WorkDir, "Analysis")
###########################################################

setwd(RawDataDir)
library("limma")
library("preprocessCore")
library("plyr")

#Read the target (sample file list) file to process
targets6480 <- readTargets("target_GPL6480.txt", path=AnalysisDir, row.names="GSMID", sep="\t")
targets6848 <- readTargets("target_GPL6848.txt", path=AnalysisDir, row.names="GSMID", sep="\t") 

#Read channel intensity values from the data files specified in the target file
x <- read.maimages(targets6480, source="agilent.median", green.only=TRUE)
y <- read.maimages(targets6848, source = "agilent.median", green.only = TRUE)

###########
#To view structure of the data
str(x)

#To view the summary of the intensity distribution of samples
summary(x$E)

#To get the Gene level annotation for each probe
Genes <- x$genes$GeneName
################

#Taking log2, quantile normalization
dQnormDat.x <- as.data.frame(normalize.quantiles(as.matrix(x$E)))
names(dQnormDat.x) <- colnames(x)

dQnormDat.y <- as.data.frame(normalize.quantiles(as.matrix(y$E)))
names(dQnormDat.y) <- colnames(y)

#Lowess Normalization and Putting the GeneNames as the first column
LdQnormDat.x <- normalizeCyclicLoess(dQnormDat.x)
LdQnormDat.x <- as.data.frame(LdQnormDat.x)
LdQnormDat.x <- cbind(Gene = x$genes$GeneName, LdQnormDat.x)

LdQnormDat.y <- normalizeCyclicLoess(dQnormDat.y)
LdQnormDat.y <- as.data.frame(LdQnormDat.y)
LdQnormDat.y <- cbind(Gene = y$genes$GeneName, LdQnormDat.y)

#Remove the unwanted rows 
removeProbes <- c("NegativeControl","BrightCorner","GE_BrightCorner","DarkCorner")
LdQnormDat.x_rp <- subset(LdQnormDat.x, !(LdQnormDat.x$Gene %in% removeProbes))
LdQnormDat.y_rp <- subset(LdQnormDat.y, !(LdQnormDat.y$Gene %in% removeProbes))

#Getting the common list of genes between the two platforms
CommonGenes <- sort(intersect(LdQnormDat.x_rp$Gene, LdQnormDat.y_rp$Gene), decreasing=FALSE)

#Subset of the common genes
x_CG <- subset(LdQnormDat.x_rp, LdQnormDat.x_rp$Gene %in% CommonGenes)
y_CG <- subset(LdQnormDat.y_rp, LdQnormDat.y_rp$Gene %in% CommonGenes)

x_CG <- x_CG[order(x_CG$Gene, decreasing=FALSE),]
y_CG <- y_CG[order(y_CG$Gene, decreasing=FALSE),]

rownames(x_CG) <- c(1:nrow(x_CG))
rownames(y_CG) <- c(1:nrow(y_CG))

#Get the Unique genes for each platform
UniqueGenes.x <- sort(unique(as.character(x_CG$Gene)), decreasing=FALSE)
UniqueGenes.y <- sort(unique(as.character(y_CG$Gene)), decreasing=FALSE)

#Create an empty matrix for GeneExpression
GeneExpr.x <- matrix(nrow=length(UniqueGenes.x), ncol=ncol(x_CG)-1) 
GeneExpr.y <- matrix(nrow=length(UniqueGenes.y), ncol=ncol(y_CG)-1) 

# filling entries in gene expr matrix
for(gene in UniqueGenes.x) 
{
	geneIndex <- which(UniqueGenes.x == gene)
	GeneExpr_subset <- subset(x_CG, as.character(x_CG$Gene) == gene)[,-1]
	GeneExpr.x[geneIndex,] <- apply(GeneExpr_subset, 2, mean, na.rm=T)
	cat("X. PROCESS COMPLETE :", gene, "\n", sep="\t")
}
rownames(GeneExpr.x) <- UniqueGenes.x
colnames(GeneExpr.x) <- colnames(x_CG[2:ncol(x_CG)])

for(gene in UniqueGenes.y) 
{
	geneIndex <- which(UniqueGenes.y == gene)
	GeneExpr_subset <- subset(y_CG, as.character(y_CG$Gene) == gene)[,-1]
	GeneExpr.y[geneIndex,] <- apply(GeneExpr_subset, 2, mean, na.rm=T)
	cat("Y. PROCESS COMPLETE :", gene, "\n", sep="\t")
}
rownames(GeneExpr.y) <- UniqueGenes.y
colnames(GeneExpr.y) <- colnames(y_CG[2:ncol(y_CG)])

#Merging the two platforms 
GeneExpr <- cbind(GeneExpr.x,GeneExpr.y)
GeneExprlog <- log2(GeneExpr)
GeneExprlog[GeneExprlog == "NaN"] <- 0

#Load Design Table
DesignTbl <- read.table(file.path(AnalysisDir,"DesignTable.txt"), sep="\t", row.names=NULL, header=TRUE)

Gname <- rep(NA, length(colnames(GeneExprlog)))

for(ctr in 1:length(colnames(GeneExprlog)))
{
	GnameIndex <- which(as.character(DesignTbl$GSM_SNAME) == as.character(colnames(GeneExprlog))[ctr])
	Gname[ctr] <- as.character(DesignTbl$GSMID)[GnameIndex]
}

colnames(GeneExpr) <- Gname
colnames(GeneExprlog) <- Gname

write.table(GeneExprlog, file=file.path(AnalysisDir,"Grasso_GeneExpr_AllSamples_log2.txt"), sep="\t", row.names=TRUE, col.names=NA)

#Define Groups
NormalIndex <- which(DesignTbl$SampleType == "Normal", arr.in=TRUE)
NormalSamples <- as.character(DesignTbl$GSMID[NormalIndex])

PrimaryIndex <- which(DesignTbl$SampleType == "PrimaryTumor", arr.in=TRUE)
PrimarySamples <- as.character(DesignTbl$GSMID[PrimaryIndex])

MetIndex <- which(DesignTbl$SampleType == "MetastaticTumor", arr.in=TRUE)
MetSamples <- as.character(DesignTbl$GSMID[MetIndex])

#Generate Sepatate Matrix for the Defined Groups
GeneExpr_Normal <- GeneExprlog[,NormalSamples]
GeneExpr_Primary <- GeneExprlog[,PrimarySamples]
GeneExpr_Met <- GeneExprlog[,MetSamples]

write.table(GeneExpr_Normal, file=file.path(AnalysisDir,"Grasso_GeneExpr_Normal.txt"), sep="\t", row.names=TRUE, col.names=NA)
write.table(GeneExpr_Primary, file=file.path(AnalysisDir,"Grasso_GeneExpr_Primary.txt"), sep="\t", row.names=TRUE, col.names=NA)
write.table(GeneExpr_Met, file=file.path(AnalysisDir,"Grasso_GeneExpr_Metastasis.txt"), sep="\t", row.names=TRUE, col.names=NA)
