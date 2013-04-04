######################
#  Stat 540 Project  #
######################

# checking for batch effect 

setwd("~/Dropbox/STAT 540 Project")

# import quantile, lowess normalized, logged data
gExpDat <- read.delim("~/Dropbox/STAT 540 Project/ProcessedData/Grasso_GeneExpr_AllSamples_log2.txt", row.names = 1)
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

library(RColorBrewer)
mCols <- colorRampPalette(rev(brewer.pal(n = 9, "Greys")))

# Heatmap showing batch effect 
pdf("HeatmapExaminingBatchEffect.pdf")
heatmap(cor(gExpDat), Rowv = NA, Colv = NA, symm = TRUE, revC = TRUE, labRow = gDes$dayCode, 
          col = mCols(256), margins = c(10, 10), 
        RowSideColor = brewer.pal(11, "RdGy")[c(4, 7)][unclass(gDes$SampleType)],
        ColSideColor = brewer.pal(11, "RdGy")[c(4, 7)][unclass(gDes$SampleType)])
dev.off()
