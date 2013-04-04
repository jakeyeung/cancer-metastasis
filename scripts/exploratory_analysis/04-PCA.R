### PCA ###
source('~/Desktop/Grad/Stat540/Project/Code/ComBat.R')
library(directlabels)
library(lattice)

# load scaled dataset 
geneZ <- read.table("~/Dropbox/STAT 540 Project/Scripts/Z-score Gene data/geneZ.txt", header=T, quote="\"")

# PCA analysis
pca <- prcomp(geneZ)
pdf("PCA_plots.pdf")
plot(pca)
biplot(pca)
summary(pca) # PC1 accounts for 84% of variance 

rotations <- cbind(pca$rotation,  sample= colnames(geneZ))
rotations <-  as.data.frame(rotations)
rotations$PC1 <- as.numeric(as.character(rotations$PC1))
rotations$PC2 <- as.numeric(as.character(rotations$PC2))

p=xyplot(PC2~PC1,data=rotations, group=sample, pch = 20, cex = 1.5)
p
dev.off()

## COMBAT to correct batch effect  
# formatting design table 
DesignTable <- read.delim("~/Dropbox/STAT 540 Project/MetaData/DesignTable.txt")
DesignTable$RunDate <- factor(DesignTable$RunDate, levels(DesignTable$RunDate)[c(2, 4, 3, 1)])

DesignTable$Batch <- (DesignTable$RunDate)
levels(DesignTable$Batch) <- c(1,2,3,4)
colnames(DesignTable) <- c("Array name", "GSM_SNAME",    "Sample name" ,  "GPLID",  
                           "Covariate 1",   "Covariate 2",
                          "Covariate 3", "RunDate" ,     "Batch"  )
combatDesign <- DesignTable[,c(1,3, 9, 5, 6, 7)] 
write.table(combatDesign, "CombatDesignTable.txt", sep='\t',  row.names=F, col.names =T)
write.table(gExpDat, "gExpDat.txt", sep='\t',  row.names=F, col.names =T)

ComBat("/Users/melissalee89/Desktop/Grad/Stat540/Project/Data/gExpDat.txt",
       "/Users/melissalee89/Desktop/Grad/Stat540/Project/Data/CombatDesignTable.txt")

