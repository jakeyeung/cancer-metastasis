# 26 March 2013
# Jake Yeung
# Limma for top 40 genes found in OptDis

# Load libraries ----------------------------------------------------------

rm(list=ls())
library(limma)
library(lattice)

# Load data and set directory ---------------------------------------------


# Limma on subset of data (one run date) 
dataDirectory <- "D:/Dropbox/Dropbox/Courses/STAT 540/STAT 540 Project"
writeDirectory <- "D:/Dropbox/Dropbox/Courses/STAT 540 personal"
setwd(dataDirectory)

# import Series Matrix File
gExpDat <- read.delim("./ProcessedData/SeriesMatrixFile/Grasso_NormPriMet.txt", row.names = 1)
str(gExpDat, list.len = 6)
gDes <- read.delim("./MetaData/DesignTable.txt")
str(gDes)
(n <- nrow(gDes)) # 122
sum(is.na(gExpDat)) # 0 
gExpDat_sort <- gExpDat[, as.vector(gDes$GSMID)]

# Import top hits from 40 subnetworks
topOptDis <- read.delim("./Analysis/OptDis/PriAllMet_Grasso_SNGenesTop40.txt")
topOptDisList <- as.vector(unlist(topOptDis))

# Import HPRD genes 
PPIgenes <- read.delim("./ProcessedData/PPI_BackgroundGenes_SN.txt")
PPIList <- as.vector(unlist(PPIgenes))


# Subset data for subnetwork ----------------------------------------------
gExpDatSubset <- gExpDat_sort[PPIList, ]
gExpDatSubset <- na.omit(gExpDatSubset)
gExpDatSubset_OptDis <- gExpDat_sort[topOptDisList, ]


# Primary vs. Metastatic ------------------------------------------------------

gDes$SampleType <- factor(gDes$SampleType, 
                          levels = c("PrimaryTumor", "MetastaticTumor", "Normal"))
DesMat <- model.matrix(~SampleType, gDes)
Fit <- lmFit(gExpDatSubset, DesMat)
ebFit <- eBayes(Fit)
colnames(coef(ebFit))
hits <- topTable(ebFit, coef = grep("Metastatic", colnames(coef(ebFit))), 
                 adjust.method = "BH",
                 number = Inf)

# Plot top and bottom hits -----------------------------------------------------------

tophits <- hits$ID[c(1, 2, 3, 4, nrow(hits), nrow(hits)-1, nrow(hits)-2, nrow(hits)-3)]
tophit_data <- gExpDat_sort[tophits, ]
tophit_list <- as.list(t(gExpDat_sort[tophits, ]))
tophit_data.long <- data.frame(gDes, 
                               ID=factor(rep(rownames(tophit_data), 
                                             each=ncol(tophit_data))),
                               geneExp=unlist(tophit_list))
tophit_data.long$SampleType <- factor(gDes$SampleType, 
                                      levels = c("Normal", "PrimaryTumor", "MetastaticTumor"))

# Plot as pdf file
setwd(writeDirectory)
pdf("TopandBotHitsStripplotPrimaryvsAll.pdf")
stripplot(geneExp ~ SampleType | ID, tophit_data.long,
          auto.key = TRUE, grid = TRUE, type = c("p", "a"))
dev.off()

# Write table of hits as .txt file
write.table(hits, "topTable_Pri_vs_Met_SN.txt", sep="\t", col.names=TRUE,
            row.names=TRUE, quote=FALSE)


# Plot log-log of P-value vs. Rank ----------------------------------------------------------

hits$rank <- seq(1:nrow(hits))
hits$log.q.val <- -log10(hits$adj.P.Val)
OptDisIndex <- which(hits$ID %in% topOptDisList)
hits_optdis <- hits[OptDisIndex, ]
colors_index <- c("red", "blue")
colors <- rep(colors_index[1], nrow(hits))
colors[OptDisIndex] <- colors_index[2]
symbols <- rep(1, nrow(hits))
symbols[OptDisIndex] <- 15
# symbols <- factor(symbols)
xlab="Rank"
ylab="-log(q-value)"
group1="PPI Network Genes in Grasso dataset from t-test"
group2="Genes Discovered by OptDis"
main="Differentially expressed genes (Primary vs. Metastatic) \n t-tests vs. subnetwork approach"

# Plot as pdf file
pdf("log_log_pval_rank_bargraph_legend.pdf")
df.bar <- barplot(hits$log.q.val, border="blue",
          ylim=c(0, 40),
          ylab="-log(BH p-value)",
          xlab="Gene Rank",
          space=0.001, main=main)
rownames(df.bar) <- hits$rank
df.bar <- df.bar[hits_optdis$rank, ]
points(df.bar, hits_optdis$log.q.val, col="red", cex=1, pch=8)
legend("topright", c(group1, group2), 
       col=c("blue", "red"),
       pch=c(15, 8))
dev.off()


pdf("Log_log_pval_rank.pdf")
plot(hits$rank, hits$log.q.val, log="x",
     type="p", col=colors, pch=symbols,
     xlab=xlab, ylab=ylab,
     cex=2,
     main="Significant genes between primary and metastatic samples")
legend("topright", c(group1, group2), 
       col=colors_index,
       pch=as.numeric(levels(factor(symbols))))
dev.off()



