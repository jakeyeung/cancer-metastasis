setwd("/Users/dorjipelzom/Desktop/courses/Winter2013/STAT540/FinalProject/Data/")
library(lattice)
library("limma")
gene<-read.table("Grasso_GeneExpr_AllSamples_log2.txt", header=T)
design<-read.table("DesignTable.txt", header=T)

gene_t<-t(gene[,-1])
heatmap(gene_t[1:50,], Rowv = NA, Colv = NA, scale = "none", margins = c(5, 8))

View(gene_t)
View(design)
nrow(design)
nrow(gene_t)
names(design)

gene_dat<-as.data.frame(gene_t)
gene_dat$stage<-design$SampleType
head(gene_dat)
summary(gene_dat$stage)

normal<-subset(gene_dat, gene_dat$stage=="Normal")
Meta<-subset(gene_dat, gene_dat$stage=="MetastaticTumor")
PrimaryTumor<-subset(gene_dat, gene_dat$stage=="PrimaryTumor")

dim(normal)
normal_m<-as.matrix(normal[,-26380])
nor_cr<-cor(t(normal_m))
meta_m<-as.matrix(Meta[,-26380])
meta_cr<-cor(t(meta_m))
pri_m<-as.matrix(PrimaryTumor[,-26380])
pri_cr<-cor(t(pri_m))
pdf("heatmap.pdf")
heatmap(nor_cr, Rowv = NA, Colv = NA, scale = "none", margins = c(5, 8), main="Normal Samples")
heatmap(meta_cr, Rowv = NA, Colv = NA, scale = "none", margins = c(5, 8), main="Metastatic Tumor")
heatmap(pri_cr, Rowv = NA, Colv = NA, scale = "none", margins = c(5, 8), main="Primary Tumor")
dev.off()

#random correlation of the samples 
set.seed(624)
dat <- sample(1:ncol(gene), size = 5)
pairDat <- subset(gene, select = dat)
str(pairDat)
pairs(pairDat)


#Z-score conversion
str(gene)
geneMat <- as.matrix(gene[,-1])
mean<-apply(geneMat, 2, mean)
sd<-apply(geneMat, 2, sd)
geneZ<-(geneMat-mean)/sd

rn<-gene[,1]

rownames(geneZ)<-rn
write.table(geneZ, "geneZ.txt")
