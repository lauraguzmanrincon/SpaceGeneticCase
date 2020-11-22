
# V2: Substantial changes that might be reverted
# Comments TODO are real
#installed.packages("ggdendro")

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
# Recall: column order in distanceMatrixGroups is given by groupSequencesCore$groupId
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
# Source: /home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201907/22_17062019_BasicDataForModels.RData
# Created in: 22.R/H.
# Content: weeksForModels, regionsForModels, casesForModels, matrixForGMRF

# 9.
# Source: /home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201909/33a_03092019_DataOXTW_Perfect_Depr.RData
# Created in: 33a.R
# Content: dataLSOA, dataMSOA, dataInfoPostSOA, chosenAreasOX, chosenAreasTW

# Extra: Add duplicates information
setkey(groupSequencesCore, columnInDataAlleles)
setkey(dataInfoPostSOA, columnInDataAlleles)
groupSequencesCore[dataInfoPostSOA, isDuplicate := isDuplicate]

# Extra: create unique chosenArea # 14.04.2020
# Note we only need chosenLSOA, if S is not in the model.
if(typeSpatial == "OX"){
  chosenAreas <- chosenAreasOX
}else if(typeSpatial == "TW"){
  chosenAreas <- chosenAreasTW
}else if(typeSpatial == "ALL"){
  # Note that "ALL" should be used only if S is not in the model, since we rely on epiclustR-Input files
  chosenAreas <- list(LSOACDs = c(chosenAreasOX$LSOACDs, chosenAreasTW$LSOACDs))
}

# II. New other data ----
# 10. Function to load spatial data from 36b files
#' Load files stored previously. since it is inside a fn, other files are not stored.
#' (The time parameter does not matter)
#' typ: {string} OX or TW
#' size: {string} 20/40/60/MSOA are the only ones covered by this function
getSpatialFn <- function(typ, size){
  # Load
  load(paste0("07_MixedModelsP2/36b_epiclustRInput/36b_25092019_DataForEpiclustR_",typ,"_sp",size,"_tm","1",".RData")) # dataST, spatialInfo, ...
  # Reconstruct nb (looks inefficient but it's actually better)
  temp <- data.table(melt(dataST$nb[, 2:ncol(dataST$nb)]))[value != 0]
  matrixNeigh <- matrix(0, nrow = nrow(dataST$nb), ncol = nrow(dataST$nb))
  matrixNeigh[as.matrix(temp[, .(X1, value)])] <- 1
  if(isSymmetric(matrixNeigh) == FALSE) warning("Error in construction the GMRF matrix")
  if(sum(matrixNeigh) != sum(temp$value != 0)) warning("Error in construction the GMRF matrix")
  return(list(spatialInfo = spatialInfo, nb = dataST$nb, matrix = matrixNeigh))
}

# 10. Modify groupSequencesCore to indicate/remove duplicates, create colsDistanceMatrixGroupsNoDups
# Content: groupSequencesCore$groupDupId colsDistanceMatrixGroupsNoDups
# *Recall: column order in distanceMatrixGroups is given by groupSequencesCore$groupId
#OXADD groupSequencesCore[isDuplicate == FALSE, groupDupId := frank(groupId, ties.method = "dense")]
#OXADD colsDistanceMatrixGroupsNoDups <- groupSequencesCore[isDuplicate == FALSE, sort(unique(groupId))]

# Adjustment OXADD: 18.12.2019 :)
# Only include data in the spatial region
# NOTE: groupDupId is like a groupChosenId, but it's too complicated to modify it
colsDistanceMatrixGroupsNoDups <- groupSequencesCore[isDuplicate == FALSE &
                                                       columnInDataAlleles %in% dataInfoPostSOA[isDuplicate == FALSE & LSOA11CD %in% chosenAreas$LSOACDs, columnInDataAlleles],
                                                     sort(unique(groupId))]
groupSequencesCore[isDuplicate == FALSE & columnInDataAlleles %in% dataInfoPostSOA[isDuplicate == FALSE & LSOA11CD %in% chosenAreas$LSOACDs, columnInDataAlleles],
                   groupDupId := frank(groupId, ties.method = "dense")] # 

# 8. Partition for k's and distance matrix for chosen clusters
# UPLOAD Unless hasn't been constructed. Then CHANGE 1 to 0.
if(1){
  # Source: /home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38_II_27102019_ClustersKInfo.RData
  # Created in: 38V2.R on the 18.11.2019 as in the following block of code
  # Content: hClustOut ddata cat(readme) exploreCutFn clusterInfoFn
  # NOTE that hClustCut has the same ordering as colsDistanceMatrixGroupsNoDups and groupToClusterTable$groupDupId
  #load("/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38V2_II_18112019_ClustersKInfo.RData")
  if(typeSpatial == "OX"){
    load("/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38V2_II_18122019_ClustersKInfo_OX.RData") # Adjustment OXADD
    # hClustOut, ddata, readme, exploreCutFn, clusterInfoFn, colsDistanceMatrixGroupsNoDups*
    # *added on the 05.12.2019 (forgotten to store link of genomes by mimstake (?))
  }else if(typeSpatial == "TW"){
    load("/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_202004/38V2_II_14042020_ClustersKInfo_TW.RData")
    # hClustOut, ddata, readme, exploreCutFn, clusterInfoFn, colsDistanceMatrixGroupsNoDups
  }else if(typeSpatial == "ALL"){
    load("/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_202004/38V2_II_14042020_ClustersKInfo_ALL.RData")
    # hClustOut, ddata, readme, exploreCutFn, clusterInfoFn, colsDistanceMatrixGroupsNoDups
  }
}else{
  # Already uploaded!
  #load("/home/laura/Dropbox/Laura/PhD_Year2/06_MixedModels/RCode_201903/15_Input_MCMCCorrected04032019.RData")
  # Input: distanceMatrixGroups, groupSequencesCore, rm(numWeeks, y, numSequences)
  # Output: hClustOut exploreCutFn clusterInfoFn
  
  # Create dendogram
  # TODO why average?
  tempDistanceMatrix <- distanceMatrixGroups[colsDistanceMatrixGroupsNoDups, colsDistanceMatrixGroupsNoDups]
  distObj <- as.dist(tempDistanceMatrix)
  hClustOut <- fastcluster::hclust(distObj, method = "average")
  
  dhc <- as.dendrogram(hClustOut)
  ddata <- data.table(ggdendro::segment(ggdendro::dendro_data(dhc, type = "rectangle")))
  readme <- c("hClustObj built using hclust with average method\n", "ddata created using package ggdendro",
              "*Recall: column order in distanceMatrixGroups is given by groupSequencesCore$groupId")
  
  # Create functions
  #' Explores how cutting the tree looks like for different pairs of heighCutLow, heighCutHigh
  #' Example: exploreCutFn(heighCutLow = 50, heighCutHigh = 300, hClustOut)
  exploreCutFn <- function(heighCutLow, heighCutHigh, hClustOut){
    # Tradeoff between few points and loss of info or many points but hard to fit
    # Low cut
    hClustCutLow <- cutree(hClustOut, h = heighCutLow)
    numClustLow <- length(unique(hClustCutLow))
    cat("Num. of clusters: ", numClustLow, "\n")
    if(!sum(sort(unique(hClustCutLow)) != 1:numClustLow) == 0) stop("Unidentified error") # check indexes are coherent / ordered
    
    # High cut
    hClustCutHigh <- cutree(hClustOut, h = heighCutHigh)
    numClustHigh <- length(unique(hClustCutHigh))
    cat("Num. of clusters: ", numClustHigh, "\n")
    if(!sum(sort(unique(hClustCutHigh)) != 1:numClustHigh) == 0) stop("Unidentified error") # check indexes are coherent / ordered
    
    # Plots
    par(mfrow = c(2,1))
    plot(hClustOut)
    rect.hclust(hClustOut , h = heighCutLow, border = 2:6)
    plot(hClustOut)
    rect.hclust(hClustOut , h = heighCutHigh, border = 2:6)
  }
  
  #' Cut the hClustOut tree according to heighCutLow and heighCutHigh and return all info corresponding to the cut
  #' Requires: groupSequencesCore$groupId exists
  #' Output: clustInfoList(clusterInfo, groupToClusterTable, distanceMatrixClusters, numClustLow, numClustHigh)
  #' Example: hola <- clusterInfoFn(heighCutLow = 50, heighCutHigh = 300, hClustOut, groupSequencesCore, distanceMatrixGroups, colsDistanceMatrixGroupsNoDups)
  clusterInfoFn <- function(heighCutLow, heighCutHigh, hClustOut, groupSequencesCore, distanceMatrixGroups, colsDistanceMatrixGroupsNoDups){
    numGroupsNoDups <- length(colsDistanceMatrixGroupsNoDups)
    
    # Cut tree
    hClustCutLow <- cutree(hClustOut, h = heighCutLow)
    numClustLow <- length(unique(hClustCutLow))
    hClustCutHigh <- cutree(hClustOut, h = heighCutHigh)
    numClustHigh <- length(unique(hClustCutHigh))
    
    # Store cluster indices in a data table
    # Note that hClustCut has the same ordering as colsDistanceMatrixGroupsNoDups
    groupToClusterTable <- data.table(groupDupId = 1:numGroupsNoDups, clusterLowId = hClustCutLow, clusterHighId = hClustCutHigh)
    
    # Create distance matrix of clusters (ordered according to indices in hClustCutLow)
    tempDistanceMatrix <- distanceMatrixGroups[colsDistanceMatrixGroupsNoDups, colsDistanceMatrixGroupsNoDups]
    tempMeltedDist <- data.table(melt(tempDistanceMatrix))[X1 != X2]
    setnames(tempMeltedDist, c("X1", "X2", "value"), c("groupFrom", "groupTo", "groupDist"))
    setkey(tempMeltedDist, groupFrom)
    setkey(groupToClusterTable, groupDupId)
    tempMeltedDist[groupToClusterTable, clusterFrom := clusterLowId]
    setkey(tempMeltedDist, groupTo)
    setkey(groupToClusterTable, groupDupId)
    tempMeltedDist[groupToClusterTable, clusterTo := clusterLowId]
    
    clusDistTable <- tempMeltedDist[, mean(groupDist), by = .(clusterFrom, clusterTo)]
    setnames(clusDistTable, c("V1"), c("clustDist"))
    distanceMatrixClusters <- matrix(0, nrow = numClustLow, ncol = numClustLow)
    distanceMatrixClusters[as.matrix(clusDistTable[, .(clusterFrom, clusterTo)])] <- clusDistTable$clustDist
    
    diag(distanceMatrixClusters) <- 0
    isSymmetric(distanceMatrixClusters)
    
    # Create main output
    groupToClusterTable[, clusterId := clusterLowId]
    clusterInfo <- groupToClusterTable[, .N, .(clusterId, clusterLowId, clusterHighId)]
    
    return(list(clusterInfo = clusterInfo, groupToClusterTable = groupToClusterTable,
                distanceMatrixClusters = distanceMatrixClusters, numClustLow = numClustLow, numClustHigh = numClustHigh))
  }
  #save(hClustOut, ddata, readme, exploreCutFn, clusterInfoFn, colsDistanceMatrixGroupsNoDups, file = "07_MixedModelsP2/RCode_201911/38V2_II_18112019_ClustersKInfo.RData")
  #save(hClustOut, ddata, readme, exploreCutFn, clusterInfoFn, colsDistanceMatrixGroupsNoDups, file = "07_MixedModelsP2/RCode_201911/38V2_II_18122019_ClustersKInfo_OX.RData")
  #save(hClustOut, ddata, readme, exploreCutFn, clusterInfoFn, colsDistanceMatrixGroupsNoDups, file = "07_MixedModelsP2/RCode_202004/38V2_II_14042020_ClustersKInfo_TW.RData")
  #save(hClustOut, ddata, readme, exploreCutFn, clusterInfoFn, colsDistanceMatrixGroupsNoDups, file = "07_MixedModelsP2/RCode_202004/38V2_II_14042020_ClustersKInfo_ALL.RData")
}

# 11. Choose heighCutLow, heighCutHigh (exploration)
if(0){
  # Plot dendogram
  # TODO plot colouring of cut
  ggplot(ddata) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    coord_flip() + scale_y_reverse(expand = c(0.2, 0)) + theme_laura()
  
  # Explore effect of cut height
  tempS <- sapply(seq(5,1000,5), function(cut) length(unique(cutree(hClustOut, h = cut))))
  plot(seq(5,1000,5), tempS)
  exploreCutFn(heighCutLow = 10, heighCutHigh = 50, hClustOut) # 50 300
  # The chosen heighCutLow and heighCutHigh must be modified in 38_00.R/1.
}

# III. Input for MCMC (for time/space/genome) ----
# Idea: Create matrix y and auxiliar tables explaining the content of y
# Output:   weeksForModels, regionsForModels, casesForModels
#           distanceMatrixMCMC matrixForGMRF y pop
#           iToGroups jToGroups kToGroups
#           numRegions numPeriods(and numWeeks) numSequences numBlockPeriods numBlockRegions numBlocksSequences dimBeta

# __.1. Constants ----
minWeek <- min(dataInfoPostSOA$numWeek_corrected)
maxWeek <- max(dataInfoPostSOA$numWeek_corrected)
numBlockPeriods <- floor((maxWeek - minWeek + 1)/lengthBlockPeriod)
numPeriods <- numBlockPeriods*lengthBlockPeriod

# __.2. Create weeksForModels (as weeksInfo in 36a.R/2.) ----
weeksForModels <- data.table(weekId = 1:numPeriods, key = "weekId")
weeksForModels[order(weekId), row := 1:numPeriods]
weeksForModels[order(weekId), groupId := ceiling(row/lengthBlockPeriod)]
# (better recalculating dates than merging datasets)
tempMinDate <- dataInfoPostSOA[numWeek_corrected == minWeek, min(received_date_nextFriday)]
weeksForModels[, numWeek := minWeek + weekId - 1]
weeksForModels[, received_date_nextFriday := tempMinDate + 7*(weekId - 1)]
weeksForModels[, received_date_nextFriday_group := tempMinDate + 7*lengthBlockPeriod*(groupId - 1)]
dataInfoPostSOA[numWeek_corrected == minWeek + numPeriods - 1, max(received_date_nextFriday)] ==
  weeksForModels[, max(received_date_nextFriday)] # useless?

# __.3. Create regionsForModels ----
# NOTE e.g. in 36 files we call it spatialInfo. Do not get confused
if(!2 %in% dimToInclude & typeSpatial == "ALL"){
  # TODO check is ok!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  temp <- getSpatialFn("OX", numClustSpatial) # TRICK to avoid problems with ALL
}else{
  temp <- getSpatialFn(typeSpatial, numClustSpatial) # taken from 36b
}
regionsForModels <- temp$spatialInfo
matrixForGMRF_nb <- temp$nb
matrixForGMRF <- temp$matrix
numRegions <- nrow(regionsForModels)
numBlockRegions <- length(unique(regionsForModels$regionId))

# __.4. Create genomeForModels ----
# "we have to guarantee that no clusters are empty because of duplicates" No, empty clusters are also imformative
# Pre-adjustments: create clusters
clusterInfoList <- clusterInfoFn(heighCutLow, heighCutHigh, hClustOut, groupSequencesCore, distanceMatrixGroups, colsDistanceMatrixGroupsNoDups) # ~ 4sec
setkey(groupSequencesCore, groupDupId)
setkey(clusterInfoList$groupToClusterTable, groupDupId)
groupSequencesCore[clusterInfoList$groupToClusterTable, c("clusterHighId", "clusterLowId") := .(clusterHighId, clusterLowId)]

numSequences <- clusterInfoList$numClustLow
numBlocksSequences <- clusterInfoList$numClustHigh
genomeForModels <- clusterInfoList$clusterInfo[, .(clusterLowId, clusterHighId)]
genomeForModels[, col := 1:numSequences]
distanceMatrixClusters <- clusterInfoList$distanceMatrixClusters

# __.5. Create casesForModels (as casesInfo in 36b.R/5.) ----
casesForModels <- dataInfoPostSOA[isDuplicate == FALSE &
                                    columnInDataAlleles %in% dataInfoPostSOA[isDuplicate == FALSE & LSOA11CD %in% chosenAreas$LSOACDs, columnInDataAlleles], # Adjustment OXADD
                                  # & numWeek_corrected <= minWeek + lengthBlockPeriod*numBlockPeriods - 1,
                             .(id, columnInDataAlleles, numWeek_corrected, idInDataLSOA, received_date_nextFriday, LSOA11CD, MSOA11CD)]
setnames(casesForModels, "id", "idInDataInfo")
cat("Number of cases:", nrow(casesForModels))
# Create and add row/column numbers
# (First be sure the week number is a nice key)
#casesForModels[, idInWeek := numWeek_corrected] # better key but messy
sum(weeksForModels[, .N, numWeek]$N != 1) == 0
setkey(casesForModels, numWeek_corrected)
setkey(weeksForModels, numWeek)
casesForModels[weeksForModels, c("idInWeek", "dim1Cases") := .(weekId, row)]
setkey(casesForModels, idInDataLSOA)
setkey(regionsForModels, idInDataLSOA)
casesForModels[regionsForModels, dim2Cases := col]
setkey(casesForModels, columnInDataAlleles)
setkey(groupSequencesCore, columnInDataAlleles)
casesForModels[groupSequencesCore, c("groupDupId", "clusterHighId", "clusterLowId") := .(groupDupId, clusterHighId, clusterLowId)]
setkey(casesForModels, clusterLowId)
setkey(genomeForModels, clusterLowId)
casesForModels[genomeForModels, dim3Cases := col]
# Population
pop <- regionsForModels[order(col), estimPop2015]
# .Check
sum(is.na(casesForModels$dim1Cases)) == 0 # can be false
sum(is.na(casesForModels$dim2Cases)) == 0 # can be false
sum(is.na(casesForModels$dim3Cases)) == 0 # can be false??
range(casesForModels$dim1Cases, na.rm = T)
range(casesForModels$dim2Cases, na.rm = T)
range(casesForModels$dim3Cases, na.rm = T)

# __.6. Mapping to groups ----
iToGroups <- weeksForModels[order(row), groupId] # length of numPeriods, max value is numBlockPeriods
jToGroups <- regionsForModels[order(col), region]
kToGroups <- genomeForModels[order(col), clusterHighId]

# __.7. Collapse non-included dimensions ----
# Not a clean way to reuse same code
if (length(dimToInclude) == 3) stop("Operation cannot be supported. Choose maximum two dimensions.")
if(!1 %in% dimToInclude){
  casesForModels[, dim1Cases := 1]
  weeksForModels[, c("row", "groupId") := .(1, 1)]
  numPeriods <- 1
  numBlockPeriods <- 1
  iToGroups <- 1
}else{
  casesForModels <- casesForModels[!is.na(dim1Cases)] # TODO dangerous. Better one? test?
}
if(!2 %in% dimToInclude){
  casesForModels[, dim2Cases := 1]
  regionsForModels[, c("col", "regionId", "region") := .(1, 1, 1)]
  matrixForGMRF <- matrix(1, nrow = 1, ncol = 1)
  matrixForGMRF_nb <- matrix(0, nrow = 1, ncol = 1)
  numRegions <- 1
  numBlockRegions <- 1
  pop <- c(1)
  jToGroups <- 1
}else{
  casesForModels <- casesForModels[!is.na(dim2Cases)] # TODO dangerous. Better one? test?
}
if(!3 %in% dimToInclude){
  casesForModels[, dim3Cases := 1]
  groupSequencesCore[, c("clusterLowId", "clusterHighId") := .(1, 1)]
  genomeForModels <- data.table(col = 1, clusterLowId = 1, clusterHighId = 1)
  numSequences <- 1
  numBlocksSequences <- 1
  kToGroups <- 1
}else{
  # Nothing to do supposedly
  # Adjustment OXADD
  casesForModels <- casesForModels[!is.na(dim3Cases)]
}

# III.7. MCMC input ----
numBlockPeriods == max(iToGroups)
numBlockRegions == max(jToGroups)
numBlocksSequences == max(kToGroups)

numBeta <- (dimBeta == 2)*numBlockRegions + (dimBeta == 3)*numBlocksSequences

#allToGroups <- c(sapply(jToGroups, function(y) sapply(iToGroups, function(x) (y-1)*numBlockPeriods + x)))
#numBlockPeriods*numBlockRegions == max(allToGroups)
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
                                                              function(x) (z-1)*numBlockRegions*numBlockPeriods + (y-1)*numBlockPeriods + x))))
numBlockPeriods*numBlockRegions*numBlocksSequences == max(allToGroups)
numPeriods*numRegions*numSequences == length(allToGroups)

temp <- casesForModels[, .N, by = .(dim1Cases, dim2Cases, dim3Cases)]
y <- array(0, dim = c(numPeriods, numRegions, numSequences)) # large array
y[as.matrix(temp[,.(dim1Cases, dim2Cases, dim3Cases)])] <- temp$N
sum(y) == nrow(casesForModels)
sum(y)

if(simulatedData){
  load(file = "07_MixedModelsP2/RCode_201911/38_Files/38_01_21112019_MCMCInputTinis.RData")
  # TODO update
  #load("UPDATE 07_MixedModelsP2/RCode_201907/28_Files/28_13_26062019_SimulatedData.RData") # yy, simulatedParam
  load("07_MixedModelsP2/RCode_201911/38_Files/38_13_21112019_SimulatedData.RData") # yy, simulatedParam
  y <- yy
}

# Adjustments to keep notation of MCMC in the other files
numWeeks <- numPeriods
distanceMatrixMCMC <- distanceMatrixClusters
numWeeksGroups <- numBlockPeriods
numRegionsGroups <- numBlockRegions
numSequenceGroups <- numBlocksSequences
# New one
numDims <- c(numPeriods, numRegions, numSequences)
numBlockDims <- c(numBlockPeriods, numBlockRegions, numBlocksSequences)
listDims <- list(iToGroups, jToGroups, kToGroups)

# ----
inputForData <- list(simulatedData = simulatedData, dimToInclude = dimToInclude, lengthBlockPeriod = lengthBlockPeriod, dimBeta = dimBeta,
                     heighCutLow = heighCutLow, heighCutHigh = heighCutHigh, typeSpatial = typeSpatial, numClustSpatial = numClustSpatial)
groupSequencesCore_chosen <- groupSequencesCore # 18.12.2019 :)
groupSequencesCore_chosen[, groupChosenId := groupDupId]
save(inputForData,
     simulatedData, dimToInclude, lengthBlockPeriod, dimBeta, heighCutLow, heighCutHigh, typeSpatial, numClustSpatial,
     weeksForModels, regionsForModels, genomeForModels, casesForModels,
     groupSequencesCore_chosen,
     numSequences, numPeriods, numRegions, numBlockPeriods, numBlockRegions, numBlocksSequences,
     iToGroups, jToGroups, kToGroups,
     numBeta, pop,
     numWeeks, numWeeksGroups, numRegionsGroups, numSequenceGroups,
     numDims, numBlockDims, listDims, allToGroups,
     y, distanceMatrixMCMC, matrixForGMRF_nb, matrixForGMRF,
     #file = "007_MixedModelsP2/RCode_201911/38_Files/38_01_21112019_MCMCInputTinis.RData")
     #file = "007_MixedModelsP2/RCode_201911/38_Files/38_01_21112019_SimulInputTinis.RData")
     #file = "007_MixedModelsP2/RCode_201912/38_01_18122019_MCMCInput_GOX.RData")
     #file = "007_MixedModelsP2/RCode_201912/38_01_06012020_MCMCInput_GOX2.RData")
     #file = "007_MixedModelsP2/RCode_201912/38_01_06012020_MCMCInput_GOX3.RData")
     #file = "007_MixedModelsP2/RCode_202002/38_01_05032020_MCMCInput_TGOX.RData")
     #file = "007_MixedModelsP2/RCode_202002/38_01_29032020_MCMCInput_TGOX.RData")
     file = nameToStoreInput)
     #file = "007_MixedModelsP2/RCode_202004/38_01_14042020_MCMCInput_SGTW.RData")
     #file = "007_MixedModelsP2/RCode_202004/38_01_14042020_MCMCInput_TGALL.RData")
     
     



