##########################################################
# SETTINGS #
###########################################################
workDir <- "/Data/Raunak/STAT540_CRPC"
Analysis0Dir <- file.path(workDir, "Analysis")
AnalysisDir <- file.path(workDir, "Analysis/Pathway")

NetworkDir <- "/Data/Raunak/DataInt_DriverNet/InteractionNetwork"
networkfile <- "HPRDNetworkWithComplexes_GeneID_nonDuplicateInt.txt"

BatchDir <- file.path(workDir, "Outputs/Results/MarkerDiscovery/PriAllMet_Grasso-Chandran")
ModuleDir <- file.path(BatchDir, "Progression/INFORGAIN/modules")

ExprSummaryfile <- "PriAllMet_Grasso_ExprSummary.txt"
pathwayfile <- "PriAllMet_IPAPath_Grasso.txt"
SNPath_file <- "PriAllMet_SNinPathway_Grasso.txt"
OutGraphfile <- "PriAllMet_SNPathNet_Grasso.pdf"

title <- "CRPC Grasso: Primary -> Metastasis"

#CommonGenes <- c("AR","CNN1","EZH2","MEIS1","MEIS2","MRVI1","PBX1","PLN","PTN","SMAD4","SORBS1","VCL") #Common Genes
###########################################################
library('Matrix')
library("igraph")
library("stringr")


######## PPI NETWORK #########
#Load Network
network <- read.table(file.path(NetworkDir, networkfile ), sep="\t", header=FALSE, stringsAsFactors=FALSE, na.strings = "null")
colnames(network) <- c("A", "B")

#Create Graph
ggPPI <- graph.data.frame(network, directed=FALSE)

#delete edges from the graph which has self loop
gPPI <- delete.edges(ggPPI, which(is.loop(ggPPI)))

####################################

#Load Expression Summary
ExprSummary <- read.table(file.path(Analysis0Dir, ExprSummaryfile), sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings = "null")
ExprSummary$Color <- ifelse(ExprSummary$FC > 0, "red", "green")
colnames(ExprSummary)[1] <- "Gene"

######### PATHWAY #########
#Get IPA Pathway Enrichment file
Path <- read.table(file.path(AnalysisDir,pathwayfile), sep="\t", header=TRUE)

cutoff <- 1.00
Path_passvalIndex <- which(Path$X.log.B.H.p.value. >= cutoff)
Path_pass <- Path[Path_passvalIndex,]

Path_pass_order <- Path_pass[order(Path_pass$X.log.B.H.p.value., decreasing=TRUE),]

PathGene_split <- strsplit(as.character(Path_pass_order$Molecules), ",")
TotalGenes <- unlist(strsplit(as.character(Path_pass_order$Molecules), ","))
UniqueGenes <- unique(TotalGenes)

for(ctr in 1:length(PathGene_split))
{
  Path_pass_order$NumOfGenes[ctr] <- length(PathGene_split[[ctr]])
}

Path_pass_order$RatioGenes <- Path_pass_order$NumOfGenes/min(Path_pass_order$NumOfGenes)
Path_pass_order$HeatColor <- heat.colors(nrow(Path_pass_order))

#Compute Pathway Relation Pn -> Pm
PathGeneRelation <- matrix(0, ncol=2, nrow=length(TotalGenes))
rowcount <- 1
for(pathIndex in 1:nrow(Path_pass))
{
	pathGenelist <- PathGene_split[[pathIndex]]
	
	for(geneIndex in 1:length(pathGenelist))
	{
			PathGeneRelation[rowcount,1] <- as.character(Path_pass_order$IngenuityCanonicalPathways)[pathIndex]
			PathGeneRelation[rowcount,2] <- pathGenelist[geneIndex]
			rowcount <- rowcount + 1
	}
	cat("Complete for Pathway :", pathIndex, "of", length(Path_pass_order$IngenuityCanonicalPathways), as.character(Path_pass_order$IngenuityCanonicalPathways)[pathIndex],"\n", sep="\t")
}
colnames(PathGeneRelation) <- c("Pathway","Gene")
PathGeneRelation.df <- as.data.frame(PathGeneRelation)
####################################

####### SUBNETWORKS ################
#Get list of files containing Sub-Network Genes
list.of.files <- list.files(ModuleDir, full.names=TRUE)
filenames <- list.of.files[grep("modules/[0-9]{1,2}", list.of.files)]

#Get FileNames
f_stringlist <- str_split(str_extract(filenames, "modules/[0-9]*"), "/", n=2)
modules.df <- as.data.frame(do.call(rbind, f_stringlist))
modules <- sort(as.numeric(modules.df$V2), decreasing=FALSE)

#Get Genes in each SubNetworks, Direction and aggregate them
SNgenes_all <- data.frame(row.names=NULL)
GenesinSN <- data.frame(row.names=NULL)
genecount <- 1

for(ctr in 1:40)
{
	file <- modules[ctr]
	
	#Load each module file
	SNgenes_file <- read.table(file.path(ModuleDir,file), header=FALSE, sep="\t")
	colnames(SNgenes_file) <- c("Gene","Value")
	SNgenes <- as.character(SNgenes_file$Gene)
	
	SNnumber <- paste("SN", file, sep="-")
	
	#Assign Colors to each gene according to Direction of Expression
	GeneColor <- rep(NA,length(SNgenes))
	
	for(count in 1:length(SNgenes))
	{
		SNgenes_all[genecount,1] <- SNgenes[count]
		SNgenes_all[genecount,2] <- ExprSummary$Color[which(ExprSummary$Gene == SNgenes[count])]
		GeneColor[count] <- ExprSummary$Color[which(ExprSummary$Gene == SNgenes[count])]
		SNgenes_all[genecount,3] <- SNnumber
		genecount <- genecount + 1
	}
	
	SN.dat <- data.frame(cbind(SNgenes_file, GeneColor))
	GenesinSN[ctr,1] <- SNnumber
	GenesinSN[ctr,2] <- paste(as.character(SN.dat$Gene),collapse=",")
}
colnames(SNgenes_all) <- c("Gene","Color","SNnumber")
colnames(GenesinSN) <- c("SNnumber","Genes")

#Get SubNetworks mapped to each of the gene
Genes <- unique(SNgenes_all$Gene)
vec.sn <- rep(NA, length(Genes))
for(i in 1:length(Genes))
{
	gIndex <- which(SNgenes_all$Gene == Genes[i])
	a <- SNgenes_all[gIndex,]
	vec.sn[i] <- paste(a$SNnumber, collapse=",")
}

SNinGenes <- data.frame(Genes=Genes, SN=vec.sn)

######### SUB-NETWORK + PATHWAY ############

### Only a fraction of genes in subnetworks define the significant list of pathways
gPathway <- as.character(unique(PathGeneRelation.df$Gene))
PathGeneRelation.df$SN <- rep(NA, nrow(PathGeneRelation.df))

#Arrange Subnetworks mapped to each genes from PathGeneRelation,df$Gene
for(i in 1:length(gPathway))
{
	v.sn <- as.character(SNinGenes$SN[which(SNinGenes$Genes == gPathway[i])])
	PathGeneRelation.df$SN[which(as.character(PathGeneRelation.df$Gene) == gPathway[i])] <- v.sn
	cat(i, gPathway[i], "\n", sep="\t")
}	

#Get SubNetworks related to each of the pathway
SNinPathway <- data.frame(Pathways=Path_pass_order$IngenuityCanonicalPathways, Genes=Path_pass_order$Molecules)
SNinPathway$SubNetwork <- rep(NA, nrow(SNinPathway))

for(i in 1:nrow(SNinPathway))
{
	v.sn <- PathGeneRelation.df$SN[which(as.character(PathGeneRelation.df$Pathway) == as.character(SNinPathway$Pathways)[i])]
	SNinPathway$SubNetwork[i] <-  paste(v.sn, collapse=",")
}	
###################################################################

############## SUB-NETWORK + PATHWAY + GENES ######################
#### Get All the genes in the SubNetworks for each Pathway
SNinPathway$SubNetworkGenes <- rep(NA, nrow(SNinPathway))

SN_split <- strsplit(as.character(SNinPathway$SubNetwork), ",")

for(p in 1:nrow(SNinPathway))
{
	vec.SNgenes <- rep(NA, length(SN_split[[p]]))
	for(s in 1:length(SN_split[[p]]))
	{
		 vec.SNgenes[s] <- GenesinSN$Genes[which(as.character(GenesinSN$SNnumber) == SN_split[[p]][s])]
	}
	SNinPathway$SubNetworkGenes[p] <- paste(unique(strsplit(paste(vec.SNgenes, collapse=","), ",")[[1]]), collapse=",")
	#cat(p,"\n")
}

write.table(SNinPathway, file.path(AnalysisDir, SNPath_file), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
###################################################################

###############################################
##### Plot Network Graph for each Pathway ##### 
###############################################

SNgene_split <- strsplit(as.character(SNinPathway$SubNetworkGenes), ",")

pdf(file.path(AnalysisDir, OutGraphfile))

for(p in 1:nrow(SNinPathway))
{
	vec.genes <- SNgene_split[[p]]
	
	#create a subgraph of genes
	sub.gPPI <- subgraph(gPPI, vec.genes)
	
	sub.gPPI.color <- rep(NA, length(vec.genes))
	y <- 1
	for(gene in as.character(V(sub.gPPI)$name))
	{
		sub.gPPI.color[y] <- SNgenes_all$Color[which(SNgenes_all$Gene == gene)]
		y <- y + 1
	}
	
	V(sub.gPPI)$label <- V(sub.gPPI)$name
	V(sub.gPPI)$label.cex <- .8
	V(sub.gPPI)$color <- sub.gPPI.color
	V(sub.gPPI)$size <- 8
	
	#### For Assigning FrameColor to CommonGenes ####
	#CGIndex <- which(as.character(V(sub.gPPI)$name) %in% CommonGenes)
	#fcol <- rep(NA,length(V(sub.gPPI)))
	#fcol[CGIndex] <- "black"
	#V(sub.gPPI)$frame.color <- fcol
	#################################################

	V(sub.gPPI)$frame.color <- NA
	
	plot(sub.gPPI, main = as.character(SNinPathway$Pathways)[p], layout=layout.fruchterman.reingold)
}
dev.off()

