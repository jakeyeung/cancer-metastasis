##########################################################
# SETTINGS #
###########################################################
workDir <- "/Data/Raunak/STAT540_CRPC"
processedDataDir <- file.path(workDir, "ProcessedData")
outputDir <- file.path(workDir, "Outputs")

scriptDir <- "/Data/kwang/bdvtools/TestFramework"
snDir <- "/Data/kwang/bdvtools/Subnetwork"


#. Set parameters and load libraries .#
options(width = 100, digits=8)
source(file.path("/Data/kwang/bdvtools/ArrayProcessing", "visualizeBatch.R"))
source(file.path(scriptDir, "runAnalysis.R"))
source(file.path(scriptDir, "utilities.R"))

###########################################################
#. LOAD PROCESSED DATA
###########################################################

Chandran_expr <- read.table(file.path(processedDataDir, "PriAllMets_Chandran.txt"), sep="\t", row.names=1, header=TRUE, stringsAsFactors=F)
Grasso_expr <- read.table(file.path(processedDataDir, "PriAllMets_Grasso.txt"), sep="\t", row.names=1, header=TRUE, stringsAsFactors=F)

Train_expr <- Grasso_expr
Test_expr <-Chandran_expr

#External Validation
source(file.path(scriptDir, "runAnalysis.R"))
runExternalValidation_SN(trainExpr=Train_expr, testExpr=Test_expr, analysis="MarkerDiscovery", batch="PriAllMet_Grasso-Chandran", endpoint="Progression", 
  network="HPRD", optSN="OptDis", optActivity="Single", scale="Brief", platform="U133Plus2", sameplatform=T, ftrSelection=F, nSNs=50, supervisedLearning="Classification")
plotPerformance_CustomSN(analysis="MarkerDiscovery", batch="PriAllMet_Grasso-Chandran", endpoint="Progression", optSN="OptDis", optActivity="Single", nSNs=50)

#External Cross-Validation (Rep=5; Fold=3)
MarkerModuleDir <- file.path(workDir, "Outputs/Results/MarkerDiscovery/PriAllMet_Grasso-Chandran/Progression/INFORGAIN/modules")

runInternalValidation_SN_Index(exprMatrix=Test_expr, nReps=5, nFolds=3, analysis="CV_R5_N3", batch="PriAllMet_Grasso-Chandran", endpoint="Progression", 
	network="HPRD", optSN="OptDis", optActivity="Single", scale="Brief", platform="U133Plus2", createFolds=T, nSNs=50, snModuleDir=MarkerModuleDir)


###########################################################
#. DATA PREPARATION FOR IPA
###########################################################
AnalysisDir <- file.path(workDir, "Analysis")
bglist_file <- "/Data/Raunak/rotproj_2012_prostate_marray/Analysis/IPA_bg_gene/IPA_BackgroundGenes_SN.txt"
ModuleDir <- file.path(workDir, "Outputs/Results/MarkerDiscovery/PriAllMet_Grasso-Chandran/Progression/INFORGAIN/modules")

# Get gene list for IPA Analysis
expr <- Train_expr
classLabel <- as.character(expr[1,])
expr.matrix <- data.matrix(expr[-1,])

# Find median, stdev in each group
expr_subset <- expr.matrix
median_A <- apply(expr_subset[,classLabel==unique(classLabel)[1]], 1, median)
median_B <- apply(expr_subset[,classLabel==unique(classLabel)[2]], 1, median)
sd_A <- apply(expr_subset[,classLabel==unique(classLabel)[1]], 1, sd)
sd_B <- apply(expr_subset[,classLabel==unique(classLabel)[2]], 1, sd)
fc <- ifelse(median_B > median_A, 2^median_B/2^median_A, -2^median_A/2^median_B)

#conduct t-test
ttestFun <- function(x, labels) {
	result <- t.test(x~as.factor(labels), na.rm=T, paired=F, var.equal=F)
	return(result$p.value)
}
pVal <- apply(expr_subset, 1, ttestFun, classLabel)
pVal_BH <- p.adjust(pVal, method="BH")
expr_summary <- data.frame(Median_Norm=median_A, SD_Norm=sd_A, Median_Pri=median_B, SD_Pri=sd_B, FC=fc, pVal=pVal, pVal_BH=pVal_BH)

write.table(expr_summary, file=file.path(AnalysisDir,"PriAllMet_Grasso_ExprSummary.txt"), sep="\t", col.names=NA)


#Top 50 Marker List

retrieveMultipleModuleGeneListForIPA(bgGeneListFile=bglist_file, moduleDir=ModuleDir, nModules=50)

ipa_list <- read.table(file.path(ModuleDir,"ipa_40.txt"), sep="\t", stringsAsFactors=F, header=T, row.names=NULL)
ipa_list_color <- ipa_list

geneList <- rownames(expr_subset)

for (index in match(geneList, ipa_list_color[,1])) {
	if (!is.na(index)) {
		ipa_list_color[index,]$DummyExpr <- expr_summary[ipa_list_color[index,1],]$FC
		if (is.infinite(ipa_list_color[index,]$DummyExpr) | is.nan(ipa_list_color[index,]$DummyExpr))
			ipa_list_color[index,]$DummyExpr <- 0
	}
}

write.table(ipa_list_color, file=file.path(AnalysisDir,"PriAllMet_Grasso_IPA40.txt"), sep="\t", col.names=T, row.names=F, quote=F)






