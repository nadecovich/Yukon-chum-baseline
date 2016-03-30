#Tile: WASC Baseline
#Date:
#Name:Nick DeCovich
#   rm(list=ls(all=TRUE))
#This function ??????????


#Run this to load the workspace
load("V:/Analysis/3_AYK/Chum/Yukon/2016 Baseline/2016baselineData")
##save.image("V:/Analysis/3_AYK/Chum/Yukon/2016 Baseline/2016baselineData")
#Set the working directory
setwd("V:/Analysis/3_AYK/Chum/Yukon/2016 Baseline")
#Source all .GCL functions
source("C:/Users/nadecovich/Documents/R/Functions.GCL.r")
[1] "AllPossiblePhenotypes.GCL"       "AttributesToIDs.GCL"             "CheckDupWithinSilly.GCL"         "CombineConflictsWithPlateID.GCL"
[5] "CombineLoci.GCL"                 "ConfusionMatrices.GCL"           "CreateBaseline.GCL"              "CreateControlFile.GCL"          
[9] "CreateLocusControl.GCL"          "CreateMixture.GCL"               "CreateSPAMcontrol.GCL"           "CustomCombineBAYESOutput.GCL"   
[13] "CustomCombineHWLEROutput.GCL"    "DupCheckBetweenSillys.GCL"       "FailureRate.GCL"                 "FindAlternateSpecies.GCL"       
[17] "FisherCompute.GCL"               "FishersTest.GCL"                 "FreqPop.GCL"                     "gcl2Genepop.GCL"                
[21] "gcl2Nexus.GCL"                   "GenepopTOgcl.GCL"                "HoFisFstTable.GCL"               "HWLERControlFile.GCL"           
[25] "LeaveOneOutDist.GCL"             "LinkageCorrelationJAGS.GCL"      "LOKI2R.GCL"                      "LOKI2R_GAPS.GCL"                
[29] "MultiChainInits.GCL"             "PairwiseFstTree.GCL"             "PoolCollections.GCL"             "PostPriorDiffAlleleFreq.GCL"    
[33] "Prior.GCL"                       "ProofTest.GCL"                   "RandomInits.GCL"                 "ReadBiomarkQC.GCL"              
[37] "ReadProjectLOKI2R.GCL"           "RemoveAlternateSpecies.GCL"      "RemoveDups.GCL"                  "RemoveIDs.GCL"                  
[41] "RemoveIndMissLoci.GCL"           "RepeatedProofTest.GCL"           "SampSizeByLocus.GCL"             "StratifiedEstimator.GCL"        
[45] "treeColor.GCL"             




####ReadLOKI.GCL(sillyvec=YukonChinookSillys,markersuite="ChinookYukon2014_41SNPs")

#
species <- "chum"
markersuite <- "WAK_Chum_192SNPs"

project <- "CM045"
projectID <- 2210

username <- "nadecovich"
password <- "Terje0623*"
CreateLocusControl.GCL(markersuite = markersuite, username = username, password = password)

# Read in Project data
ReadProjectLOKI2R.GCL(projectID = projectID, username = username, password = password)

# Verify SILLYs imported
objects(pattern = "\\.gcl")

# Create a list of project sillys
ProjectSillys <- unlist(strsplit(objects(pattern = "\\.gcl"),  split = "\\.gcl"))

loci <- LocusControl$locusnames

##### Get rid of markers not scored by lab
### [1] "Oke_F036"     "Oke_hepcid5"  "Oke_mitofer1" "Oke_mitofer2" "Oke_RAD11690" "Oke_RAD15073" "Oke_RAD3948"  "Oke_RAD6703" 
##### These are part of the WASSIP set and are part of a linked pair and should be dropped
### Oke_pgap-92   Oke_gdh1-62
#### These are mtDNA and will probably be dealt with seperately 
###   "Oke_Cr30"        "Oke_Cr386" "Oke_ND3-69" "Oke_12776_MT" "Oke_13594_MT"
write.table(loci,file= "Output/locuscontrol.xls",sep="\t", col.names=NA, row.names=T)
loci182 <- loci[-c(22,30,35,42,43,51,61,72,104,116)]

### Narrow to only list of project sillies
ProjectSillys <- scan("clipboard",what='')

### Group vec for project sillies

GroupVec42  <- scan("clipboard",what='')
GroupVec42  <- as.numeric(GroupVec42)
GroupNames7 <- scan("clipboard",what='')
Group7Colors <- scan("clipboard",what='')


#Create folders
MixtureFolders <- c("Output","BAYES","SPAM")
for (x in MixtureFolders){
  dir <- paste(getwd(),"/",x,sep='')
  dir.create(dir)}
#Subfolders for Bayes and SPAM
files <- c("Mixture","Baseline","Output","Control")
for (x in files){
  dir<- paste(getwd(),"/BAYES/",x,sep='')
  dir.create(dir)
  dir<- paste(getwd(),"/SPAM/",x,sep='')
  dir.create(dir)}  

#Create Mixture Output folders
for (x in Mixtures){
  dir<- paste(getwd(),"/BAYES/Output/",x,sep='')
  dir.create(dir)}

#Dump .gcl objects for safe keeping
dump(paste(ProjectSillys ,".gcl",sep=''),"Output/WASC.gcl_objects.txt")


#Check sample sizes before removing individuals with >20% missing locus scores.
ColSize <- sapply(paste(ProjectSillys,".gcl",sep=''), function(x) get(x)$n)
write.table(ColSize,file= "Output/OriginalColSize.xls",sep="\t", col.names=NA, row.names=T)
missloci<-RemoveIndMissLoci.GCL(sillyvec=ProjectSillys ,proportion=0.8)
ColSize2 <- sapply(paste(ProjectSillys ,".gcl",sep=''), function(x) get(x)$n)
ColSizePostMiss <- ColSize-ColSize2 #checking to see how many individuals were removed.
write.table(ColSize2,file= "Output/ColSizePostMissLoci.xls",sep="\t", col.names=NA, row.names=T)



#Check for duplicate fish
dupcheck<-CheckDupWithinSilly.GCL(sillyvec=ProjectSillys,loci=loci,quantile=NULL,minproportion=0.95)


##Paste duplicate info here
Check dup check object for dups. Nothing removed as of 2/29/16. Will revisit

#Remove duplicate individuals
removedDups<-RemoveDups.GCL(dupcheck)
unlist(removedDups[!removedDups=="Nothing Removed"])#Gets the individuals removed from workspace
##  Paste removed info here





ColSizePostDup <- sapply(paste(ProjectSillys,".gcl",sep=''), function(x) get(x)$n)
write.table(ColSizePostDup,file= "Output/ColSizePostDup.xls",sep="\t", col.names=NA, row.names=T)


## generate GENEPOP format data file for HWE and Linkage check
gcl2Genepop.GCL(sillyvec=ProjectSillys,loci=DipLoci,path="V:/Analysis/3_AYK/Chum/WASC/Genepop/WASCchum",VialNums = FALSE)

HWE=ReadGenepopHWE.GCL(file="V:/Analysis/3_AYK/Chum/WASC/Genepop/WASCchum.P")
write.table((HWE$SummaryPValues),"Output/HWEprepooling.xls",sep='\t',col.names=NA, row.names=T)

## Dump allele 2 frequencies
AlleleCounts <- FreqPop.GCL(sillyvec=ProjectSillys,loci=loci182)
str(AlleleCounts)
Freqs <- AlleleCounts[,,"Allele 2"]/(AlleleCounts[,,"Allele 2"]
                                     + AlleleCounts[,,"Allele 1"])
str(Freqs)
write.table(Freqs,"Output/frequencies.xls",sep='\t',col.names=NA, row.names=T)

## Use Fishers test to check homogeneity of temporal samples before pooling
## Check temporal colletions for pooling
TemporalTest1=list(ProjectSillys[11:12],
                   ProjectSillys[14:16],
                   ProjectSillys[18:19],
                   ProjectSillys[21:22],
                   ProjectSillys[29:30],
                   ProjectSillys[31:32],
                   ProjectSillys[35:37],
                   ProjectSillys[38:39])


FisherTemporalResults <- FishersTest.GCL(freq=AlleleCounts, loci=loci182,tests=TemporalTest1)
sink(file="Output/FisherTemporalResults.xls")
print(FisherTemporalResults)
sink()



## 

lapply(1:length(TemporalTest1),function(x){PoolCollections.GCL(TemporalTest1[[x]], loci=loci182, IDs = NULL,
                                                               newname = paste(TemporalTest1[[x]],collapse = "."))})


ls(pattern=".gcl")

SillysPooled <- scan("clipboard",what='')
GroupVec32  <- scan("clipboard",what='')
GroupVec32  <- as.numeric(GroupVec32)
GroupNames7 <- scan("clipboard",what='')
Group7Colors <- scan("clipboard",what='')

PooledColSize <- sapply(paste(SillysPooled,".gcl",sep=''), function(x) get(x)$n)

## Create freq plots using pooled pops
source("V:/Analysis/R files/Scripts/DEV/FreqFisPlot4SNPs.GCL.r")

FreqFisPlot4SNPs.GCL(sillyvec=SillysPooled, loci=loci182, groupvec=GroupVec32, alpha=0.05,groupcol=NULL,file="V:/Analysis/3_AYK/Chum/WASC/Output/WASCpopfreq.pdf")

