##########################################################
# SETTINGS #
###########################################################
workDir <- "/Data/Raunak/STAT540_CRPC"
AnalysisDir <- file.path(workDir, "Analysis/Pathway")
###########################################################

pathwayfile <- "PriAllMet_IPAPath_Grasso.txt"
OutBarGraphfile <- "PriAllMet_IPApathway_Grasso.pdf"
title <- "CRPC Grasso: Pathway Enrichment, Primary -> Metastasis"

Path <- read.table(file.path(AnalysisDir,pathwayfile), sep="\t", header=TRUE, as.is=c("IngenuityCanonicalPathways", "X.log.p.value.", "X.log.B.H.p.value.", "Ratio", "Molecules"))

cutoff <- 1.00
Path_passvalIndex <- which(Path$X.log.B.H.p.value. >= cutoff)
Path_pass <- Path[Path_passvalIndex,]
Path_pass_order <- Path_pass[order(Path_pass$X.log.B.H.p.value., decreasing=FALSE),]

####### barplot
pdf(file.path(AnalysisDir, OutBarGraphfile))
par(mar=c(5,15,4,2))
barplot(Path_pass_order$X.log.B.H.p.value., names = Path_pass_order$IngenuityCanonicalPathways, main=title,
  xlab="-log(BHpvalue)", horiz = TRUE, space=1, cex.main = 0.8, cex.names=0.3, las=2, col="brown")
dev.off()
