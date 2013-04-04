######################
#  Stat 540 Project  #
######################

setwd("~/Dropbox/STAT 540 Project/Analysis/UpdatedPlots")
library(RColorBrewer)

# import data after correcting for batch effect
gExpDat <- read.delim("~/Dropbox/STAT 540 Project/ProcessedData/SeriesMatrixFile/Grasso_NormPriMet.txt", row.names = 1)
str(gExpDat, list.len = 6)
gDes <- read.delim("~/Dropbox/STAT 540 Project/MetaData/DesignTable.txt")
str(gDes)
(n <- nrow(gDes)) # 122
sum(is.na(gExpDat)) # 0 

table(gDes$RunDate)
# 07-Feb-2007 07-Sep-2005 15-Aug-2006 6-June-2006 
# 88          12          12          10 

# converting to date class
gDes$Date <- as.Date(gDes$RunDate, "%d-%b-%Y")
table(gDes$Date)
# 2005-09-07 2006-06-06 2006-08-15 2007-02-07 
# 12         10         12         88 

# converting to factor class
gDes$RunDate <- factor(gDes$RunDate, levels(gDes$RunDate)[c(2, 4, 3, 1)])
        
# order by date
(mOrd <- with(gDes, order(RunDate, SampleType, TissueType, DiseaseState)))
gDes <- gDes[mOrd,]
gExpDat <- gExpDat[, mOrd]

gDes$dayCode <- sapply(unclass(gDes$RunDate), function(i) {
      foo <- rep("-", nlevels(gDes$RunDate) * 2)
      foo[(i * 2 - 1):(i * 2)] <- i
      paste(foo, collapse = "")
  })
gDes$dayCode <- paste(gDes$dayCode, colnames(gExpDat))

mCols <- colorRampPalette(rev(brewer.pal(n = 9, "Greys")))

# Heatmap showing batch effect 
pdf("HeatmapAfterCorrectingBatchEffect.pdf")
heatmap(cor(gExpDat), Rowv = NA, Colv = NA, symm = TRUE, revC = TRUE, labRow = gDes$dayCode, 
          col = mCols(256), margins = c(10, 10), 
        RowSideColor = brewer.pal(11, "RdGy")[c(4, 7)][unclass(gDes$SampleType)],
        ColSideColor = brewer.pal(11, "RdGy")[c(4, 7)][unclass(gDes$SampleType)])
dev.off()

# PCA Analysis
pdf("PCAplotsAfterCorrectingBatchEffect.pdf")
pcs <- prcomp(gExpDat, center = F, scale = F)
plot(pcs)
gDes <- cbind(gDes, pcs$rotation[gDes$GSMID, 1:10])
colnames(gDes)
plot(gDes[, c(1,5,6,11,12,13)], pch = 19, cex = 0.8)

# plot data on first two PCs, colored by sample type

par(xpd=TRUE) # turn off clipping to plot legend outside
plot(gDes[, c("PC1", "PC2")], bg = gDes$SampleType, pch = 21, cex = 1.5, col = "black")
legend(list(x = -0.06, y = 0.17), as.character(levels(gDes$SampleType)), pch = 21, 
       pt.bg = c(1, 2, 3), cex = 0.8)
dev.off()

# compute pairwise distances
pr.dis = dist(t(gExpDat), method = "euclidean")
pr.nms = with(gDes, paste(gExpDat, SampleType, sep = "_"))

# compute hierarchical clustering using different linkage types
pr.hc.s <- hclust(pr.dis, method = "single")
pr.hc.c <- hclust(pr.dis, method = "complete")
pr.hc.a <- hclust(pr.dis, method = "average")
pr.hc.w <- hclust(pr.dis, method = "ward")

# plot them
pdf("dendogramsAfterBatchCorrection.pdf")
par(mar = c(0, 4, 4, 2))
par(mfrow = c(2, 2))

plot(pr.hc.s, labels = FALSE, main = "Single", xlab = "")
plot(pr.hc.c, labels = FALSE, main = "Complete", xlab = "")
plot(pr.hc.a, labels = FALSE, main = "Average", xlab = "")
plot(pr.hc.w, labels = FALSE, main = "Ward", xlab = "")
dev.off()
# Objects in columns

set.seed(31)
samplecenters <- 5
pr.km <- kmeans(t(gExpDat), centers = samplecenters, nstart = 50)

# We can look at the within sum of squares of each cluster
pr.km$withinss
library(xtable)
pr.kmTable = data.frame(SampleType = gDes$SampleType, cluster = pr.km$cluster)
prTable = xtable(with(pr.kmTable, table(SampleType, cluster)), caption = "Number of samples from each sample type within each cluster")

align(prTable) <- "lccccc"
print(prTable, caption.placement = "top")
