##########################################################
# SETTINGS #
###########################################################
WorkDir <- "/Data/Raunak/STAT540_CRPC"
RawDataDir <- "/Data/Raunak/rotproj_2012_prostate_marray/RawData/GSE3988_CRPC"
AnalysisDir <- file.path(WorkDir, "Analysis")
###########################################################

setwd(RawDataDir)
library("limma")

#Read the target (sample file list) file to process
targets6480 <- readTargets("target_GPL6480.txt", path=AnalysisDir, row.names="GSMID", sep="\t")
targets6848 <- readTargets("target_GPL6848.txt", path=AnalysisDir, row.names="GSMID", sep="\t")

#Merge both the targets
targets <- rbind(targets6480,targets6848) 

#Create an empty dataframe to store GSMID and dates
dat <- data.frame(row.names=NULL)

#Get GSMID and dates from each of the rawfile
for(count in 1:nrow(targets))
{
	RawFile <- paste(RawDataDir,targets$FileName[count], sep="/")
	
	#Get GSMID
	id <- targets$GSMID[count]
	dat[count,1] <- strsplit(id, " ")[[1]][1]
	
	#Read the RawFile (skip 1st line, take 2nd line as header and read only 1 row/line after header)
	metadat <- read.table(RawFile, sep="\t", header=TRUE, skip=1, nrow=1)
	
	#Get Date_Time and split data & time
	date_time <- as.character(metadat$Protocol_date)
	date <- strsplit(date_time, " ")[[1]][1]
	dat[count,2] <- date
	
	cat("PROCESS COMPLETE FOR : ",  targets$GSMID[count],"\n", sep="\t")
}
colnames(dat) <- c("GSMID","RunDate")

#Load DesignTable
DesignTbl <- read.table(file.path(AnalysisDir,"DesignTable.txt"), sep="\t", header=TRUE, row.names=1)

#Arrange the RunDates in the same order of samples as in the DesignTable
RunDate <- rep(NA, nrow(dat))

for(ctr in 1:nrow(dat))
{
	rIndex <- which(as.character(DesignTbl$GSMID) == dat$GSMID[ctr])
	RunDate[ctr] <- dat$RunDate[rIndex]
}

#Merge the RunDates to the DesignTable
des <- cbind(DesignTbl, RunDate)
write.table(des, file=file.path(AnalysisDir,"DesignTable.txt"), sep="\t", row.names=TRUE, col.names=NA)

