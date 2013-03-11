# Jake Yeung
# batch_effect_heatmap.R
# Looks at batch effects using heatmap

rm(list=ls())
T# Load library and data
library(lattice)
library(ggplot2)
dataDirectory <- "D:/Dropbox/Dropbox/Courses/STAT 540/STAT 540 Project/ProcessedData"
designDirectory <- "D:/Dropbox/Dropbox/Courses/STAT 540/STAT 540 Project/MetaData"
setwd(dataDirectory)
exprs_normal <- read.table("Grasso_GeneExpr_Normal.txt", 
                           header=TRUE, sep="\t")
exprs_metastasis <- read.table("Grasso_GeneExpr_Metastasis.txt", 
                               header=TRUE, sep="\t") 
exprs_primary <- read.table("Grasso_GeneExpr_Primary.txt",
                            header=TRUE, sep="\t")
design <- read.table(paste(designDirectory, "/DesignTable.txt", sep=""),
                     header=TRUE, sep="\t")


# Function definitions
MakeRownames <- function(data){
# Takes a dataframe and makes rownames the first column of the dataframe,
# then delete first column of dataframe.
    rownames(data) <- data[, 1]
    data <- data[, -1]
}


# Make rownames the first column in each table. 
exprs_normal <- MakeRownames(exprs_normal)
exprs_metastasis <- MakeRownames(exprs_metastasis)
exprs_primary <- MakeRownames(exprs_primary)


# Do correlation matrix to see if any batch effects. 
# First, make correlation matrix for normal, mets and primary
# Second, name rownames and column names to include dates.
# Third, plot heatmap to see if any batch effects. 
# 1.
normal_cor <- cor(exprs_normal)
normal_metastasis <- cor(exprs_metastasis)
normal_primary <- cor(exprs_primary)
# 2. 
heatmap(normal_cor)


