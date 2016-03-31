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
markersuite <- "ChumGolden2011_96SNPs"

##project <- "CM045"
##projectID <- 2210

username <- "nadecovich"
password <- "Terje0623*"
CreateLocusControl.GCL(markersuite = markersuite, username = username, password = password)

# Read in Project data
ReadProjectLOKI2R.GCL(projectID = projectID, username = username, password = password)

# Verify SILLYs imported
objects(pattern = "\\.gcl")

# Create a list of project sillys
ProjectSillys <- unlist(strsplit(objects(pattern = "\\.gcl"),  split = "\\.gcl"))



ProjectSillys  <- scan("clipboard",what='')

LOKI2R.GCL(sillyvec = ProjectSillys, username = username, password = password)




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
dump(paste(ProjectSillys ,".gcl",sep=''),"Output/Yukon2016.gcl_objects.txt")


#Check sample sizes before removing individuals with >20% missing locus scores.
ColSize <- sapply(paste(ProjectSillys,".gcl",sep=''), function(x) get(x)$n)
write.table(ColSize,file= "Output/OriginalColSize.xls",sep="\t", col.names=NA, row.names=T)
missloci<-RemoveIndMissLoci.GCL(sillyvec=ProjectSillys ,proportion=0.8)
ColSize2 <- sapply(paste(ProjectSillys ,".gcl",sep=''), function(x) get(x)$n)
ColSizePostMiss <- ColSize-ColSize2 #checking to see how many individuals were removed.
write.table(ColSize2,file= "Output/ColSizePostMissLoci.xls",sep="\t", col.names=NA, row.names=T)



#Check for duplicate fish
dupcheck<-CheckDupWithinSilly.GCL(sillyvec=ProjectSillys,loci=loci,quantile=NULL,minproportion=0.95)



#Remove duplicate individuals
removedDups<-RemoveDups.GCL(dupcheck)
unlist(removedDups[!removedDups=="Nothing Removed"])#Gets the individuals removed from workspace
##  Paste removed info here





ColSizePostDup <- sapply(paste(ProjectSillys,".gcl",sep=''), function(x) get(x)$n)
write.table(ColSizePostDup,file= "Output/ColSizePostDup.xls",sep="\t", col.names=NA, row.names=T)

## Create new locus list based on linkage work done in WASSIP and 2013 Yukon analysis

### Keep Oke_pgap-111 (drop Oke_pgap-92)
### Keep Oke_gdh1-191 (drop Oke_gdh1-62)
### Drop mtDNA
loci.91 <- LocusControl$locusnames[-c(14,15,24,38,41)]

## generate GENEPOP format data file for HWE and Linkage check
gcl2Genepop.GCL(sillyvec=ProjectSillys,loci=DipLoci,path="V:/Analysis/3_AYK/Chum/WASC/Genepop/WASCchum",VialNums = FALSE)

HWE=ReadGenepopHWE.GCL(file="V:/Analysis/3_AYK/Chum/WASC/Genepop/WASCchum.P")
write.table((HWE$SummaryPValues),"Output/HWEprepooling.xls",sep='\t',col.names=NA, row.names=T)

## Dump allele 2 frequencies
AlleleCounts <- FreqPop.GCL(sillyvec=ProjectSillys,loci=loci.91)
str(AlleleCounts)
Freqs <- AlleleCounts[,,"Allele 2"]/(AlleleCounts[,,"Allele 2"]
                                     + AlleleCounts[,,"Allele 1"])
str(Freqs)
write.table(Freqs,"Output/frequencies.xls",sep='\t',col.names=NA, row.names=T)

## Use Fishers test to check homogeneity of temporal samples before pooling
## Check temporal colletions for pooling
TemporalTest=list(ProjectSillys[3:4],
                  ProjectSillys[9:13],
                  ProjectSillys[16:17],
                  ProjectSillys[26:28],
                  ProjectSillys[29:32],
                  ProjectSillys[33:35],
                  ProjectSillys[38:39])
                  


TemporalTest2=list(ProjectSillys[42:44],
                  ProjectSillys[45:47],
                  ProjectSillys[52:53],
                  ProjectSillys[55:56],
                  ProjectSillys[58:63],
                  ProjectSillys[65:67],
                  ProjectSillys[70:71],
                  ProjectSillys[72:73],
                  ProjectSillys[74:75]
)



FisherTemporalResults <- FishersTest.GCL(freq=AlleleCounts, loci=loci.91,tests=TemporalTest2)
sink(file="Output/FisherTemporalResults2.xls")
print(FisherTemporalResults)
sink()



## 

lapply(1:length(TemporalTest2),function(x){PoolCollections.GCL(TemporalTest2[[x]], loci=loci, IDs = NULL,
                                                               newname = paste(TemporalTest2[[x]],collapse = "."))})

##Forgot to do Tatchun
PoolCollections.GCL(collections = ProjectSillys[c(72:73)], loci=loci, IDs = NULL,
                    newname = paste(ProjectSillys[c(72:73)],collapse = "."))

ls(pattern=".gcl")

SillysPooled <- scan("clipboard",what='')
GroupVec3 <- scan("clipboard",what='')
GroupVec3 <- as.numeric(GroupVec3)
GroupNames3 <- scan("clipboard",what='')
Group3Colors <- scan("clipboard",what='')

PooledColSize <- sapply(paste(SillysPooled,".gcl",sep=''), function(x) get(x)$n)

## Create freq plots using pooled pops
source("V:/Analysis/R files/Scripts/DEV/FreqFisPlot4SNPs.GCL.r")

FreqFisPlot4SNPs.GCL(sillyvec=SillysPooled, loci=loci182, groupvec=GroupVec32, alpha=0.05,groupcol=NULL,file="V:/Analysis/3_AYK/Chum/WASC/Output/WASCpopfreq.pdf")


## Convert genind to a genpop list for dist.genpop function and create distance matrix
gcl2Genepop.GCL(sillyvec=SillysPooled,loci=loci.91,path="Genepop/YukonChum.gen",VialNums = FALSE)
library("ape")
genind = read.genepop(file="Genepop/YukonChum.gen")
genpop <- genind2genpop(x=genind)
AdegenetNei <- dist.genpop(genpop, method = 2, diag = T, upper = T)## methods 1=Nei's (1972), 2=CSE, 3=Reynold's, 4=Rogers', 5=Provesti's

#####################################################################################################################

NJofAdegenetNei318tree <- nj(as.matrix(AdegenetNei318))
NJofAdegenetNei318tree$tip.label <- names131
NJofAdegenetNei318tree$edge.length=pmax(0,NJofAdegenetNei318tree$edge.length) #Get rid of negative branches
plot(NJofAdegenetNei318tree,font=2,cex=.3)
axis(1)  ## Adds scale to bottom of plot



#######################################################################################################################
SillysPooled <- scan("clipboard",what='')
GroupVec3 <- scan("clipboard",what='')
GroupVec3 <- as.numeric(GroupVec3)
GroupNames3 <- scan("clipboard",what='')
Group3Colors <- scan("clipboard",what='')

WAK_Likelihood_Profile <- LeaveOneOutDist.GCL(sillyvec=WAKPop60, loci=loci.319[-c(305)], groupvec=GroupVec60)

WAK_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist=WAK_Likelihood_Profile, groupnames=GroupNames60, groupvec=GroupVec60, sillyvec=WAKPop60)

sapply(rownames(WAK_Confusion[[1]]), function(row) {barplot(height = WAK_Confusion[[1]][row, ], col = WAKColors, main = row, ylim = c(0, 1))})


############## Multi Dimensional Scaling (MDS) Plots:
####These are from AYK juvenile work. 


####### Create MDS Plots ##############



x=as.vector(cmdscale(as.matrix(AdegenetNei)[1:45,1:45],k=3)[,1])   
y=as.vector(cmdscale(as.matrix(AdegenetNei)[1:45,1:45],k=3)[,2])
z=as.vector(cmdscale(as.matrix(AdegenetNei)[1:45,1:45],k=3)[,3])

library('rgl')

plot3d(x,y,z+abs(range(z)[1]),xlab='',ylab='',zlab='',aspect=F,col=Group3Colors,size=.75,type='s',box=T,axes=F,top=T,cex=4)
plot3d(x,y,z+abs(range(z)[1]),aspect=F,col="Black",size=100,type='h',box=T,axes=T,top=F,add=T) #adds pins to spheres
axes3d(edges="bbox",'x-')
texts3d(x,y,z+abs(range(z)[1]),adj=c(1.2,-.5),text=labels,font=2,cex=.6,add=T,top=T,axes=F)#adds numbers to points(adj moves the numbers around the points)
rgl.snapshot(filename="MDS/NeiYukonchum.png",fmt="png",top=TRUE)  #this saves a snapshot view of your MDS plot


labels= scan('clipboard',what='')



