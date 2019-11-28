# TODO Note: all region part are OUT-TO-DATED

# I. Current data (loaded from 34.R) ----
# 1. General data
# Source: /home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/Data/FS101013_Datasets_V3_26R.RData
# Created in: 02_COPY_UploadDataSets.R. Updated in: 10_UpdateFS101013File.R. Updated in: 26_UpdateDATAINFO_24062019.R
# Content: dataAllelesList, dataAnnotatedList, dataDistancesList2, dataInfo(WEEK UPDATED), dataNameLociList,
# numLociList, dirPathDatasets, fileAllelesList, fileDistancesList, fileInfo, numIsolates,
# dataCoreGenome

# 2. Core genome data
# lengthCore, distancesCore, distances1643, dataAllelesCore
# Remove if more than 200 NA's in the core genome
#plot(table(colSums(is.na(dataAllelesList[[1]]))))
#plot(table(colSums(is.na(dataAllelesList[[1]][dataCoreGenome$columnsFrom1643,]))))

# 3. Duplicated information (more than one isolate per patient)
# Source: /home/laura/Dropbox/Laura/PhD_Year2/06_MixedModels/Data/FS101013_DataInfoDup.RData
# Created in: 06_MixedModels/12.R 16.02
# Content: dataInfoDup, pairRelatedIsolates, distanceMatrixDuplicates

# 4. ?? Distance of all neighbours per isolate
# Source: /home/laura/Dropbox/Laura/PhD_Year2/06_MixedModels/RCode_201901/09_13022019_FS101013Neighbours.RData
# Created in: 06_MixedModels/09.R/6.
# Content: neighboursTotal

# 5. MCMC input (for time/genome)
# Source: /home/laura/Dropbox/Laura/PhD_Year2/06_MixedModels/RCode_201903/15_Input_MCMCCorrected04032019.RData
# Created in: 14b.R/1. on the 04.03.2019
# Content: numSequences, numWeeks, y, distanceMatrixGroups, groupSequencesCore
rm(numWeeks, y, numSequences)
suppressWarnings(groupSequencesCore[, numWeek := NULL])
suppressWarnings(groupSequencesCore[, weekId := NULL])

# 6 and 7 are OUT-TO-DATE 
# 6. Location info
# a. Load UK sectors with location and RU classification
# Source: /home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201906/21a_04062019_SectorsLocationAndRU.RData
# Created in: 21a.R/1.
# Content: ukSectors, RUClass

# b. Upload dataInfo copy with postcodes
# Source: /home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201907/21a_10062019_DataInfoWithPostCodes.RData
# Created in: 21a.R/2.
# Content: dataInfoPost, cat(review)

# c. Shapefiles and areas information
# Source: /home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201906/21a_30052019_FilesForMapping.RData
# Created in: 21a.R/3.
# Content: area.pointsAll, area.postcodes, data.postcodesDT (DT to link postcodes with polygons in area.postcodes),
# neighbouringRegions (orden given by data.postcodesDT$rnInt)

# 7. Cases data for all models (TO BE UPDATED!!!!)
# TODO
# Source: /home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201907/22_17062019_BasicDataForModels.RData
# Created in: 22.R/H.
# Content: weeksForModels, regionsForModels, casesForModels, matrixForGMRF

# 8.
# Source: /home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201909/33a_03092019_DataOXTW_Perfect_Depr.RData
# Created in: 33a.R
# Content: dataLSOA, dataMSOA, dataInfoPostSOA, chosenAreasOX, chosenAreasTW

# Extra: Add duplicates information
setkey(groupSequencesCore, columnInDataAlleles)
setkey(dataInfoPostSOA, columnInDataAlleles)
groupSequencesCore[dataInfoPostSOA, isDuplicate := isDuplicate]

# II. New auxiliar data ----
# 8. Partition for k's and distance matrix for chosen clusters
if(1){
  # Source: /home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201905/20_2a_14052019_ClustersKInfo.RData
  # Created in: 20.R on the 14.05.2019 as in the following block of code
  # Content: hClustOut, clusterTable, heighCutLow, heighCutHigh, numClust, distanceMatrixClusters, numSequences (corrected),
  #          numBlockClust, groupSequencesCore (new column: groupDupId)
  # Main outputs: clusterTable, groupSequencesCore (modified)
  load("/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38_II_27102019_ClustersKInfo.RData")
  
  clusterTable[, clusterId := clusterLowId]
  
  setkey(groupSequencesCore, groupDupId)
  setkey(clusterTable, groupDupId)
  groupSequencesCore[clusterTable, c("clusterHighId", "clusterLowId") := .(clusterHighId, clusterLowId)]
  
  clusterInfo <- clusterTable[, .N, .(clusterLowId, clusterHighId)]
  
}else{
  #load("/home/laura/Dropbox/Laura/PhD_Year2/06_MixedModels/RCode_201903/15_Input_MCMCCorrected04032019.RData")
  # Created on Year2/06/RCode_201903/14a.R/2.a. on the 04.03.2019
  # numSequences, numWeeks, y, distanceMatrixGroups, groupSequencesCore
  
  # Modified to remove duplicates:
  groupSequencesCore[isDuplicate == FALSE, groupDupId := frank(groupId, ties.method = "dense")]
  tempcols <- groupSequencesCore[isDuplicate == FALSE, sort(unique(groupId))]
  
  numSequences <- length(tempcols)
  tempDistanceMatrix <- distanceMatrixGroups[tempcols, tempcols]
  
  # Create dendogram
  distObj <- as.dist(tempDistanceMatrix)
  hClustOut <- fastcluster::hclust(distObj, method = "average")
  plot(hClustOut)
  
  # Explore effect of cut height
  tempS <- sapply(seq(5,1000,5), function(cut) length(unique(cutree(hClustOut, h = cut))))
  plot(seq(5,1000,5), tempS)
  
  # Cut tree
  # Tradeoff between few points and loss of info or many points but hard to fit
  # CHOOSE heighCutLow for background
  heighCutLow <- 50
  hClustCutLow <- cutree(hClustOut, h = heighCutLow)
  numClustLow <- length(unique(hClustCutLow))
  cat("Num. of clusters: ", numClustLow)
  plot(hClustOut)
  rect.hclust(hClustOut , h = heighCutLow, border = 2:6)
  sum(sort(unique(hClustCutLow)) != 1:numClustLow) == 0 # check indexes are coherent / ordered
  
  # CHOOSE heighCutLow for outbreak indicators
  heighCutHigh <- 300
  hClustCutHigh <- cutree(hClustOut, h = heighCutHigh)
  numClustHigh <- length(unique(hClustCutHigh))
  cat("Num. of clusters: ", numClustHigh)
  plot(hClustOut)
  rect.hclust(hClustOut , h = heighCutHigh, border = 2:6)
  sum(sort(unique(hClustCutHigh)) != 1:numClustHigh) == 0 # check indexes are coherent / ordered
  
  # Store cluster indices in a data table
  # Note that hClustCut has the same ordering as distanceMatrixGroups
  clusterTable <- data.table(groupDupId = 1:numSequences, clusterLowId = hClustCutLow, clusterHighId = hClustCutHigh)
  
  # Create distance matrix of clusters (ordered according to indices in hClustCutLow)
  tempMeltedDist <- data.table(melt(tempDistanceMatrix))[X1 != X2]
  setnames(tempMeltedDist, c("X1", "X2", "value"), c("groupFrom", "groupTo", "groupDist"))
  setkey(tempMeltedDist, groupFrom)
  setkey(clusterTable, groupDupId)
  tempMeltedDist[clusterTable, clusterFrom := clusterLowId]
  setkey(tempMeltedDist, groupTo)
  setkey(clusterTable, groupDupId)
  tempMeltedDist[clusterTable, clusterTo := clusterLowId]
  
  clusDistTable <- tempMeltedDist[, mean(groupDist), by = .(clusterFrom, clusterTo)]
  setnames(clusDistTable, c("V1"), c("clustDist"))
  distanceMatrixClusters <- matrix(0, nrow = numClustLow, ncol = numClustLow)
  distanceMatrixClusters[as.matrix(clusDistTable[, .(clusterFrom, clusterTo)])] <- clusDistTable$clustDist
  distanceMatrixClusters[1:10,1:10] # ordered as clusterId in clusterTable
  
  diag(distanceMatrixClusters) <- 0
  isSymmetric(distanceMatrixClusters)
  
  numClust <- numClustLow
  numBlockClust <- numClustHigh
  
  # Save
  #OLDsave(hClustOut, clusterTable, heighCut, numClust, distanceMatrixClusters, file = "07_MixedModelsP2/RCode_201905/20_2a_14052019_ClustersKInfo.RData")
  #save(hClustOut, clusterTable, heighCutLow, heighCutHigh, numClust, distanceMatrixClusters, numSequences, groupSequencesCore, numBlockClust,
  # # numSequences corrected
  #     file = "07_MixedModelsP2/RCode_201911/38_II_27102019_ClustersKInfo.RData")
}

# 17.07.2019 :D CORRECTION on -
#nrow(clusDistTable) == sum(clusterTable[, .N, by = clusterId]$N != 1) + 232806
#sum(diag(distanceMatrixClusters) == 0)
#diag(distanceMatrixClusters) <- 0

# No. (14.05.2019) Create matrices of indices for G. Create matrix of time distanceMatrixTime
if(0){
  load("07_MixedModelsP2/20_2a_14052019_ClustersKInfo.RData") # hClustOut, clusterTable, heighCut, numClust, distanceMatrixClusters
  
  indexK <- matrix(hClustCut, nrow = numWeeks, ncol = numSequences, byrow = TRUE)
  
  intervalAgg <- 4*6
  tempJ <- floor(((1:numWeeks) - 1)/intervalAgg) + 1
  indexJ <- matrix(tempJ, nrow = numWeeks, ncol = numSequences, byrow = FALSE)
  
  sizePartitionK <- numClust
  sizePartitionJ <- max(tempJ)
  
  # matrix/vector to go from a small matrix (sizePartitionJ, sizePartitionK) to a large one (numWeeks, numSequences)
  indexJK <- cbind(c(indexJ), c(indexK))
  indexJKrow <- apply(indexJK, 1, function(x) sizePartitionJ*(x[[2]] - 1) + x[[1]])
  sizePartitionJ*sizePartitionK == max(indexJKrow)
  
  distanceMatrixTime <- abs(matrix(1:sizePartitionJ, nrow = sizePartitionJ, ncol = sizePartitionJ, byrow = TRUE) -
                              matrix(1:sizePartitionJ, nrow = sizePartitionJ, ncol = sizePartitionJ, byrow = FALSE))
  
  # Save
  #save(sizePartitionK, sizePartitionJ, indexK, indexJ, indexJK, indexJKrow, distanceMatrixTime,
  #     file = "07_MixedModelsP2/RCode_201905/20_2a_14052019_AuxForG.RData")
}

# III. Input for MCMC (for time/space/genome) ----
# Content: weeksForModels, regionsForModels, casesForModels, matrixForGMRF
# Content: hClustOut, clusterTable, heighCut, numClust, distanceMatrixClusters

# Idea: Create matrix y and auxiliar tables explaining the content of y
# Output:   weeksForModels, regionsForModels, casesForModels
#           distanceMatrixMCMC matrixForGMRF y pop
#           iToGroups jToGroups kToGroups
#           numRegions numWeeks numSequences numWeeksGroups numRegionsGroups numSequenceGroups dimBeta

# III.1. Input values ----
lengthBlockPeriod <- 3
minWeek <- min(dataInfoPostSOA$numWeek_corrected)
maxWeek <- max(dataInfoPostSOA$numWeek_corrected)
numPeriods <- maxWeek - minWeek + 1
numBlockPeriods <- floor((maxWeek - minWeek + 1)/lengthBlockPeriod)
numWeeks <- numPeriods # to keep notation of MCMC in other files
numBlocksSequences <- numBlockClust

# III.2. Create weeksForModels (as weeksInfo in 36a.R/2.) ----
weeksForModels <- data.table(weekId = 1:numPeriods, key = "weekId")
# ..Add outbreak interval id
weeksForModels[order(weekId), row := 1:numWeeks]
weeksForModels[order(weekId), groupId := ceiling(row/lengthBlockPeriod)] # TODO check
# ..better recalculating dates than merging datasets
tempMinDate <- dataInfoPostSOA[numWeek_corrected == minWeek, min(received_date_nextFriday)]
weeksForModels[, numWeek := minWeek + weekId - 1]
weeksForModels[, received_date_nextFriday := tempMinDate + 7*(weekId - 1)]
weeksForModels[, received_date_nextFriday_group := tempMinDate + 7*lengthBlockPeriod*(groupId - 1)]
dataInfoPostSOA[numWeek_corrected == minWeek + numPeriods - 1, max(received_date_nextFriday)] ==
  weeksForModels[, max(received_date_nextFriday)] # useless?

# III.3. Create regionsForModels ----
regionsForModels <- data.table(regionId = 1, col = 1, group = 1, estimatedPop = 1)
matrixForGMRF <- matrix(1, nrow = 1, ncol = 1)

# III.4. Create genomeForModels ----
# ..we have to guarantee that no clusters are empty because of duplicates

genomeForModels <- clusterInfo[, .(clusterLowId, clusterHighId)]
genomeForModels[, col := 1:numClust]

# III.5. Create casesForModels (as casesInfo in 36b.R/5.) ----
# TODO remove last period if too small: not here because maybe genetic clusters are left empty
casesForModels <- dataInfoPostSOA[isDuplicate == FALSE,# & numWeek_corrected <= minWeek + lengthBlockPeriod*numBlockPeriods - 1,
                             .(id, columnInDataAlleles, numWeek_corrected, idInDataLSOA, received_date_nextFriday, LSOA11CD, MSOA11CD)]
setnames(casesForModels, "id", "idInDataInfo")
casesForModels[, idInWeek := numWeek_corrected] # messy
cat("Number of cases:", nrow(casesForModels))
# Create and add row/column numbers
setkey(casesForModels, idInWeek)
setkey(weeksForModels, weekId)
casesForModels[weeksForModels, dim1Cases := row]
#setkey(casesForModels, idInDataLSOA)
#setkey(regionsForModels, idInDataLSOA)
#casesInfo[regionsForModels, colCases := col]
casesForModels[, dim2Cases := 1]
setkey(casesForModels, columnInDataAlleles)
setkey(groupSequencesCore, columnInDataAlleles)
casesForModels[groupSequencesCore, c("clusterHighId", "clusterLowId") := .(clusterHighId, clusterLowId)]
setkey(casesForModels, clusterLowId)
setkey(genomeForModels, clusterLowId)
casesForModels[genomeForModels, dim3Cases := col]
# .Check
sum(is.na(casesForModels$dim1Cases)) == 0
sum(is.na(casesForModels$dim3Cases)) == 0
range(casesForModels$dim1Cases)
range(casesForModels$dim3Cases)
# TODO store

# III.6. MCMC input ----
numSequences <- numClust
numWeeks <- nrow(weeksForModels)
numRegions <- nrow(regionsForModels)

iToGroups <- weeksForModels[order(row), groupId] # length of numWeeks, max value is numWeeksGroups
jToGroups <- regionsForModels[order(col), group]
kToGroups <- genomeForModels[order(col), clusterHighId]
numWeeksGroups <- max(iToGroups)
numRegionsGroups <- max(jToGroups)
numSequenceGroups <- max(kToGroups)

dimBeta <- 3 # dimension where B is gonna change (2 or 3)
numBeta <- (dimBeta == 2)*numRegionsGroups + (dimBeta == 3)*numSequenceGroups

#allToGroups <- c(sapply(jToGroups, function(y) sapply(iToGroups, function(x) (y-1)*numWeeksGroups + x)))
#numWeeksGroups*numRegionsGroups == max(allToGroups)
## large array containing small array indices and transformed in vector byrow. Only used twice for building matrixX (28_00.R/4 and 28_08.R).
## careful: allToGroups follows the R convention of indices filled by row, e.g. matrix(1:6, ncol = 2)[3]
## e.g. sapply(c(1,1,1,2), function(y) sapply(c(1,1,1,2,2), function(x) (y-1)*2 + x))
## e.g. matrix(matrix(21:24, ncol = 2)[c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2,3,3,3,4,4)], ncol = 4)
## e.g. matrix(matrix(21:24, ncol = 2)[matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2,3,3,3,4,4), ncol = 4)], ncol = 4)

# e.g.  5,4,3 groups 2,2,3
#       ttt <- c(sapply(c(1,2,2), function(z) sapply(c(1,1,2,2), function(y) sapply(c(1,1,2,3,3), function(x) (z-1)*2*3 + (y-1)*3 + x))))
#       array(ttt, dim = c(5,4,3))
allToGroups <- c(sapply(kToGroups,
                        function(z) sapply(jToGroups,
                                           function(y) sapply(iToGroups,
                                                              function(x) (z-1)*numRegionsGroups*numWeeksGroups + (y-1)*numWeeksGroups + x))))
numWeeksGroups*numRegionsGroups*numSequenceGroups == max(allToGroups)
numWeeks*numRegions*numSequences == length(allToGroups)

temp <- casesForModels[, .N, by = .(dim1Cases, dim2Cases, dim3Cases)]
y <- array(0, dim = c(numWeeks, numRegions, numSequences)) # large array
y[as.matrix(temp[,.(dim1Cases, dim2Cases, dim3Cases)])] <- temp$N
sum(y) == nrow(casesForModels)
pop <- regionsForModels[order(col), estimatedPop] # lengh: numRegions
distanceMatrixMCMC <- distanceMatrixClusters
sum(y)

# BORRAR
#temp <- casesForModels[, .N, by = .(rowCases, colCases)]
#y <- array(0, dim = c(numWeeks, numSequencesLocations)) # large array
#y[as.matrix(temp[,.(rowCases, colCases)])] <- temp$N
#sum(y) == nrow(casesForModels)
#pop <- regionsForModels[order(col), estimatedPop] # lengh: numRegions
#distanceMatrixMCMC <- distanceMatrixClusters

if(simulatedData){
  load("07_MixedModelsP2/RCode_201907/28_Files/28_13_26062019_SimulatedData.RData") # yy, simulatedParam
  y <- yy
}


