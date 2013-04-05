##########################################################
# SETTINGS #
###########################################################
WorkDir <- "/Data/Raunak/STAT540_CRPC/"
AnalysisDir <- file.path(WorkDir, "Analysis")
BatchDir <- file.path(WorkDir, "Outputs/Results/MarkerDiscovery/PriAllMet_Grasso-Chandran")
InfogainDir <- file.path(BatchDir, "Progression/INFORGAIN")
ModuleDir <- file.path(BatchDir, "Progression/INFORGAIN/modules")
NetworkDir <- "/Data/Raunak/DataInt_DriverNet/InteractionNetwork"
batch <- "PriAllMet_SNGraph_Grasso"
Bigtitle <- "CRPC Grasso: Primary -> Metastasis, SubNetworks 1 to 40"
#CommonGenes <- c() #Common Genes
###########################################################
library("stringr")
library("igraph")

#Load Expression Summary
ExprSummary <- read.table(file.path(AnalysisDir,"PriAllMet_Grasso_ExprSummary.txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings = "null")
ExprSummary$Color <- ifelse(ExprSummary$FC > 0, "red", "green")
colnames(ExprSummary)[1] <- "Gene"

#Load Network
network <-read.table(file.path(NetworkDir,"HPRDNetworkWithComplexes_GeneID_nonDuplicateInt.txt"), sep="\t", header=FALSE, stringsAsFactors=FALSE, na.strings = "null")
colnames(network) <- c("A", "B")

#note: Invoke the following in command line
## grep -o 'Weighted Density:'.*'[0-9].[0-9]*' modules.txt > modules_density.txt

#Get Module Score file
ModuleDensity <- read.table(file.path(InfogainDir , "modules_density.txt"), sep="\t", header=FALSE, stringsAsFactors=FALSE, na.strings = "null")

#Create Graph
gg <- graph.data.frame(network, directed=FALSE)

#delete edges from the graph which has self loop
g <- delete.edges(gg, which(is.loop(gg)))

#Get list of files containing Sub-Network Genes
list.of.files <- list.files(ModuleDir, full.names=TRUE)
filenames <- list.of.files[grep("modules/[0-9]{1,2}", list.of.files)]

#Get FileNames
f_stringlist <- str_split(str_extract(filenames, "modules/[0-9]*"), "/", n=2)
modules.df <- as.data.frame(do.call(rbind, f_stringlist))
modules <- sort(as.numeric(modules.df$V2), decreasing=FALSE)

#Get SubNetwork Genes for each module and plot graph
pdf(file.path(InfogainDir, paste(batch,"pdf", sep=".")))

SNgenes_all <- data.frame(row.names=NULL)
genecount <- 1

for(ctr in 1:40)
{
  file <- modules[ctr]
	
	#Load each module file
	SNgenes_file <- read.table(file.path(ModuleDir,file), header=FALSE, sep="\t")
	colnames(SNgenes_file) <- c("Gene","Value")
	SNgenes <- as.character(SNgenes_file$Gene)
	
	#Assign Colors to each gene according to Direction of Expression
	GeneColor <- rep(NA,length(SNgenes))
	
	for(count in 1:length(SNgenes))
	{
		SNgenes_all[genecount,1] <- SNgenes[count]
		SNgenes_all[genecount,2] <- ExprSummary$Color[which(ExprSummary$Gene == SNgenes[count])]
		GeneColor[count] <- ExprSummary$Color[which(ExprSummary$Gene == SNgenes[count])]
		
		genecount <- genecount + 1
	}
	
	SN.dat <- data.frame(cbind(SNgenes_file, GeneColor))
	sub.g <- subgraph(g, SNgenes)
	
	GeneColor.c <- rep(NA,length(SNgenes))
	for(i in 1:nrow(SN.dat))
	{
		gIndex <- which(as.character(SN.dat$Gene) == as.character(V(sub.g)$name)[i])
		GeneColor.c[i] <- as.character(SN.dat$GeneColor)[gIndex]
	}

	V(sub.g)$label <- V(sub.g)$name
	V(sub.g)$color <- GeneColor.c
	V(sub.g)$frame.color <- NA
	
	title <- paste("SubNetwork-", file, sep="")
	subtitle <- ModuleDensity$V1[ctr]
	plot(sub.g, main=title , sub=subtitle, layout=layout.fruchterman.reingold)
		
	cat("Plot Generated :", title,"\n", sep="\t")	
}
colnames(SNgenes_all) <- c("Gene","Color")
SNgenes_all_unique <- SNgenes_all[!duplicated(SNgenes_all),]

#Plot for all the genes in the SubNetworks
sub.g_all <- subgraph(g, SNgenes_all_unique$Gene)

color_value <- rep(NA, length(V(sub.g_all)$name))
for(v.index in 1:length(V(sub.g_all)$name))
{
	color_value[v.index] <- SNgenes_all_unique$Color[which(SNgenes_all_unique$Gene == V(sub.g_all)$name[v.index])]
}	
	
V(sub.g_all)$label <- V(sub.g_all)$name
V(sub.g_all)$label.cex <- .3
V(sub.g_all)$color <- color_value
V(sub.g_all)$size <- 6

#### For Assigning FrameColor to CommonGenes ####
#CGIndex <- which(as.character(V(sub.g_all)$name) %in% CommonGenes)
#fcol <- rep(NA,length(V(sub.g_all)))
#fcol[CGIndex] <- "black"
#V(sub.g_all)$frame.color <- fcol
##############################################

#### If CommonGenes are not to be given different frame color
V(sub.g_all)$frame.color <- NA

E(sub.g_all)$width <- .6
E(sub.g_all)$color <- "grey"

plot(sub.g_all, main=Bigtitle , layout=layout.fruchterman.reingold)
cat("Plot Generated :", file.path(ModuleDir , paste(batch,"pdf", sep=".")), sep=" ")
dev.off()


write.table(SNgenes_all_unique$Gene, file.path(InfogainDir, "PriAllMet_Grasso_SNGenesTop40.txt"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(SNgenes_all_unique, file.path(AnalysisDir, "Top50SNgenes_PriAllMet_Grasso.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
