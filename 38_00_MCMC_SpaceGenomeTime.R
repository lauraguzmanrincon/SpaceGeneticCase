# ------------------------------ #
# MCMC
# Data: FS101013 V2
# Created on 26.10.2019
# Adapted from the 28 files, created on 24.06.2019
# Changes:  --.--.2019:
#           26.02.2020: Add block update for G
# ------------------------------ #

runTinis <- 1 # !!!!!!!!!!!!!!!!!!!!!
simulatedData <- 0
typeModel <- "TG2"

if(runTinis == 0){
  if(1){
    # 0. Load Workspace ----
    source("/home/laura/Dropbox/Laura/PhD_Year2/04_AboutApproaches/RCode/2018_04/00_Functions.R")
    setwd("/home/laura/Dropbox/Laura/PhD_Year3/")
    source("07_MixedModelsP2/RCode_201908/34_HeaderFor201908.R") # load workspace for August 2019
    load("07_MixedModelsP2/RCode_201909/33a_03092019_DataOXTW_Perfect_Depr.RData") # dataLSOA, dataMSOA, dataInfoPostSOA, chosenAreasOX, chosenAreasTW
    dirFiles <- "07_MixedModelsP2/RCode_201911/38_Files/"
    dirOutFiles <- "07_MixedModelsP2/RCode_202002/"
  }
  if(1){
    # 1. Create/load data ----
    # TODO organise
    dimToInclude <- sort(c(1,3)) # T S G # Never do the three of them!
    lengthBlockPeriod <- 4 # should be 1 if 1 is not in dimToInclude!!!
    dimBeta <- 3 # dimension where B is gonna change (2 or 3 or 123) # if 123 numBeta must be ignored TODO check
    heighCutLow <- 10 #   50  10   NO1
    heighCutHigh <- 50 #  300 50  NO25
    typeSpatial <- "OX"
    numClustSpatial <- "MSOA"
    stGeneticGroups <- c(21, 353, 828, 257, 206, 464, 48, 45) # TG2 only # groups not in list are merged in one group # assumed that there are left groups
    numSTGroups <- 1 # Do NOT modify. For SG and TG.
    if(simulatedData == 0){
      # Construct data
      #source(paste0(dirFiles, "38_01_MCMC_CreateLoadData_V2.R")) # ~ 8 sec
      #source(paste0(dirFiles, "38_01b_MCMC_CreateLoadData_TG2.R")) # ~ 8 sec # TG2 only
      # Or load (not formal)
      #NOTE:length(unique(casesForModels$clusterLowId)) == numSequences
      if(typeModel == "SG"){
        # CURRENT ONE for SG:
        load("07_MixedModelsP2/RCode_201912/38_01_18122019_MCMCInput_GOX.RData") # CONFIG: inputForData OR notsim,dims23,interval3?,beta3,cuts10/50,regionOXMSOA
        # used for 38_00_18122019, 38_00_02012020_newprior, 38_00_04012020_newprior, 38_00_06012020_newprior[*]
        #load("07_MixedModelsP2/RCode_201912/38_01_06012020_MCMCInput_GOX2.RData") # CONFIG: inputForData OR notsim,dims23,interval3,beta3,cuts1/25,regionOXMSOA
        # used for 38_00_06012020_newradii
        #load("07_MixedModelsP2/RCode_201912/38_01_06012020_MCMCInput_GOX3.RData") # CONFIG: inputForData OR notsim,dims23,interval3,beta3,cuts50/300,regionOXMSOA
        # used for 38_00_06012020_newradii2
      }else if(typeModel == "TG"){
        # CURRENT ONE for TG:
        #load("07_MixedModelsP2/RCode_202002/38_01_05032020_MCMCInput_TGOX.RData") # CONFIG: inputForData OR notsim,dims13,interval1,beta3,cuts10/50,regionOXMSOA
        load("07_MixedModelsP2/RCode_202002/38_01_29032020_MCMCInput_TGOX.RData") # CONFIG: inputForData OR notsim,dims13,interval4,beta3,cuts10/50,regionOXMSOA
      }else if(typeModel == "TG2"){
        # CURRENT ONE for TG2:
        #load("07_MixedModelsP2/RCode_202003/38_01_30032020_MCMCInput_TG2OX.RData") # CONFIG: inputForData OR notsim,dims13,interval1,beta3,cuts10/50,regionOXMSOA
        load("07_MixedModelsP2/RCode_202003/38_01_30032020_MCMCInput_TG2OX_I4.RData") # CONFIG: inputForData OR notsim,dims13,interval4,beta3,cuts10/50,regionOXMSOA
      }
      length(unique(casesForModels$clusterLowId)) == numSequences
      load("/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38V2_II_18122019_ClustersKInfo_OX.RData") # add clustering for G block update! 26.02.2020
      # hClustOut, ddata, readme, exploreCutFn, clusterInfoFn, colsDistanceMatrixGroupsNoDups*
    }else{
      # note simulatedData is stored as 0 in the first file!
      load("07_MixedModelsP2/RCode_201911/38_Files/38_01_21112019_MCMCInputTinis.RData") # CONFIG: notsim, dims23, interval3, beta3, cuts50/300, regionOXMSOA
      #load("07_MixedModelsP2/RCode_201911/38_Files/38_13_21112019_SimulatedData.RData") # yy, simulatedParam 24112019
      #load("07_MixedModelsP2/RCode_201911/38_Files/38_13_fastSim27112019_SimulatedData.RData") # yy, simulatedParam 24112019 (Skype meeting Simon 25.11.2019)
      load("07_MixedModelsP2/RCode_201911/38_Files/38_13_03122019_SimulatedData.RData") # yy, simulatedParam
      y <- yy
      simulatedData <- 1
    }
  }
}else{
  dirFiles <- "FilesFromLocal/38_Files/"
  dirOutFiles <- dirFiles
  source(paste0(dirFiles, "38_14_LoadForTinis.R"))
}

#nameOutput <- "38_00_21112019_MCMCOutputTinis.RData"
#nameOutput <- "38_00_21112019_SimulOutputTinis_a.RData"
#nameOutput <- "38_00_21112019_SimulOutputTinis_r.RData"
#nameOutput <- "38_00_21112019_SimulOutputTinis_g.RData"
#nameOutput <- "38_00_21112019_SimulOutputTinis_b.RData"
#nameOutput <- "38_00_21112019_SimulOutputTinis_x.RData"
#nameOutput <- "38_00_21112019_SimulOutputTinis_p.RData"
#nameOutput <- "38_00_21112019_SimulOutputTinis_tr.RData"
#nameOutput <- "38_00_21112019_SimulOutputTinis_ltg.RData"
#nameOutput <- "38_00_03122019_MCMCOutputDellToSim.RData"
#nameOutput <- "38_00_0312201902_MCMCOutputTinis.RData"
#nameOutput <- "38_00_1812201901_MCMCOutput_GOX.RData"
# Trying new priors and parameters
nameOutput <- "38_00_0201202001_MCMCOutput_newprior.RData" # tauG: 2,1 rho: 10,0.5 r: 10/50*
nameOutput <- "38_00_0401202001_MCMCOutput_newprior.RData" # tauG: 1,1 rho: 10,0.5 r: 10/50
nameOutput <- "38_00_0601202001_MCMCOutput_newprior.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50
nameOutput <- "38_00_0601202001_MCMCOutput_newprior_SE.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50* SEkernel
nameOutput <- "38_00_0601202001_MCMCOutput_newprior_largeL.RData" # tauG: 1,1 rho: 10,0.2 r: 10/50
nameOutput <- "38_00_0601202001_MCMCOutput_newradii.RData" # tauG: 1,0.01 rho: 10,0.5 r: 1/25 BORRAR!! (does not behave well)
nameOutput <- "38_00_0601202001_MCMCOutput_newradii2.RData" # tauG: 1,0.01 rho: 10,0.5 r: 50/300
# Chosen ones for Tinis (see Notes 06.01.2020, p.2)
#nameOutput <- "38_00_0601202001_MCMCOutputTinis_MAT.RData" # tauG: 1,0.01 rho: 10,0.2 r: 10/50* # 798311
#nameOutput <- "38_00_0601202001_MCMCOutputTinis_SE.RData" # tauG: 1,0.01 rho: 10,0.2 r: 10/50* # 798312 # WRONG:SE
#nameOutput <- "38_00_0601202001_MCMCOutputTinis_MAT32.RData" # tauG: 1,0.01 rho: 10,0.2 r: 10/50* # 798341
# Try new beta on Dell
nameOutput <- "38_00_0701202001_MCMCOutputTinis_blockBeta.RData" # tauG: 1,0.01 rho: 10,0.2 r: 10/50 kernel: 3/2 # too slow beta(dimBeta = 123)
#dimBeta <- 123
nameOutput <- "38_00_temp_MAT32.RData" # tauG: 1,0.01 rho: 10,0.2 r: 10/50 kernel: 3/2
nameOutput <- "temp_SE_Tinis.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 0 WRONG:SE
nameOutput <- "38_00_temp_SE_Dell.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 0 # with SE correction
nameOutput <- "38_00_temp_SE_Dell_noconstBORRAR.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 0 # plus constraint correction
nameOutput <- "38_00_2601202002_Tinis_SE_Constraint.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 0
#nameOutput <- "38_00_2601202001_Tinis_MAT12_Constraint.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 1/2
#nameOutput <- "38_00_2601202001_Tinis_MAT32_Constraint.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 3/2
# Running with G block update on Tinis (change for different chains...). Also, several revisions to improve mixing...
nameOutput <- "38_00_2602202001_Tinis_SE_GBlock.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 0
nameOutput <- "38_00_2602202001_Tinis_MAT12_GBlock_rho.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 1/2
nameOutput <- "38_00_0203202002_Tinis_MAT12_GBlock_rho_two.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 1/2
nameOutput <- "38_00_0303202001_Tinis_MAT12_GBlock_rho_two.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 1/2
nameOutput <- "38_00_0303202001_Tinis_MAT32_GBlock_rho_two.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 3/2
nameOutput <- "38_00_0303202001_Tinis_SE_GBlock_rho_two.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 0
# Running TG 05.03.2020
nameOutput <- "38_00_0503202001_TinisTG_MAT12.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 0.5
nameOutput <- "38_00_0503202001_TinisTG_MAT12_500it.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 0.5
# Running SG for all cutting points
nameOutput <- "38_00_1403202001_TinisSG_MAT12_1000it_cuts.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 0.5
# Running TG with S block update to improve mixing 17.03.2020 + other intervals for weeks 29.03.2020
nameOutput <- "38_00_1703202001_TinisTG_MAT12_1000it_Sblock.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 0.5 interval:1
nameOutput <- "38_00_2903202001_TinisTG_MAT12_5000it_Sblock.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 0.5 interval:1
nameOutput <- "38_00_2903202001_TinisTG_MAT12_1000it_Interval4.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 0.5 interval:4
# TG2
nameOutput <- "38_00_3003202001_TinisTG2_MAT12_1000it.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 0.5 interval:1 # 860901
nameOutput <- "38_00_3003202001_TinisTG2_MAT12_400it.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 0.5 interval:1 # 860917
nameOutput <- "38_00_3003202001_TinisTG2_MAT12_1000it_I4.RData" # tauG: 1,0.01 rho: 10,0.5 r: 10/50 kernel: 0.5 interval:4 # 861062
#nameOutput <- "BORRAR"

# 2. Set-up environment ----
source(paste0(dirFiles, "38_02_AuxiliarFn.R"))
# Functions: all plot...
source(paste0(dirFiles, "38_03_ConstantsFn.R"))
# Functions: initConstants, initConfig, initParam, createStorage, adaptive

# 3. Auxiliar functions ----
source(paste0(dirFiles, "38_11_Prior_Likelihoods.R"))

# 4. Pre-iterations ----
cat("Sim", simulatedData, "/n")
constants <- initConstants(aB = 2, bB = 1, aG = 1, bG = 0.01, aL = 10, bL = 0.5) # B large, L larger, tauG specific
#constants <- initConstants(aP = 1, bP = 50, aB = 2, bB = 1) # for fastSim27112019
#constants <- initConstants(aP = 1, bP = 5, aB = 2, bB = 1) # for 03122019
#config <- initConfig(numIterations = 1000, burnIn = 50) # estandar
#config <- initConfig(numIterations = 50, burnIn = 10) # quick exploration
#config <- initConfig(numIterations = 1000, burnIn = 50, ifAUpdate = 1, ifRUpdate = 0, ifSUpdate = 0, ifGUpdate = 0, ifBUpdate = 0, ifXUpdate = 0,
#                     ifPUpdate = 0, ifTauRUpdate = 0, ifTauSUpdate = 0, ifLTauGUpdate = 0)
config <- initConfig(numIterations = 1000, burnIn = 50, maternParameter = 0.5,
                     ifAUpdate = 1, ifRUpdate = 1, ifGUpdate = 1, ifBUpdate = 1, ifXUpdate = 1,
                     ifPUpdate = 1, ifTauRUpdate = 1, ifLTauGUpdate = 1)
# TODO create optional constructParam() function
out <- initParam(configUpdates = config) # Run after initConfig
storage <- createStorage(config)
adaptive <- adaptiveConfig(FALSE)

#simulatedParam$p <- mean(simulatedParam$X) # !!!!!!!!!!!!!!!!!!!!!!!! BORRAR

parameters <- out$parameters
sigmaJumps <- out$sigmaJumps
cat("A", parameters$a, "/n")

#parameters$a <- 0 # 19.01.2020 constrain correction
#config$ifAUpdate <- 0 # 19.01.2020 constrain correction

# Requires: distanceMatrixMCMC pop matrixForGMRF y
# iToGroups jToGroups kToGroups
# numRegions numWeeks numSequences numWeeksGroups numRegionsGroups numSequenceGroups dimBeta

# TODO The following will eventually be in a function

# Auxiliar matrices
sqrDistanceMatrix <- distanceMatrixMCMC^2
if(1 %in% dimToInclude & typeModel == "TG"){
  seasonalCoefficientMatrixSqr <- seasonalCoefficientFn(numWeeks)
}else if(1 %in% dimToInclude & typeModel == "TG2"){
  seasonalCoefficientMatrixSqr <- matrix(0, nrow = numPeriods, ncol = numPeriods)
  for(scm in 1:numSTGroups){
    seasonalCoefficientMatrixSqr[(scm - 1)*numPeriodsPerST + (1:numPeriodsPerST), (scm - 1)*numPeriodsPerST + (1:numPeriodsPerST)] <- seasonalCoefficientFn(numPeriodsPerST)
  }
  numSTGroups*sum(abs(seasonalCoefficientFn(numPeriodsPerST))) == sum(abs(seasonalCoefficientMatrixSqr))
}
matrixPop <- aperm(array(pop, dim = c(numRegions, numWeeks, numSequences)), perm = c(2,1,3))
#neighLength <- rowSums(matrixForGMRF)
neighLength <- matrixForGMRF_nb[,1]
AMatrix <- diag(1, numSequences, numSequences) - matrix(1/numSequences, numSequences, numSequences) # 22.01.2020

# Initial covariance matrix and other initial matrices
# TODO c++??
#deltaMatrix <- deltaMatrixFn(parameters$l, maternParameter) # + 0.5*diag(numSequences)
#invDeltaMatrix <- solve(deltaMatrix)
# SPAM deltaMatrix <- as.spam(deltaMatrixFn(parameters$l, distanceMatrixMCMC, sqrDistanceMatrix, config$maternParameter) + 0.01*diag(numSequences))
deltaMatrix <- deltaMatrixFn(parameters$l, distanceMatrixMCMC, sqrDistanceMatrix, config$maternParameter) + 0.01*diag(numSequences) # ??? 03.03.2020 17.45
#invDeltaMatrix <- solve(deltaMatrix) # TODO works always as a matrix?
invDeltaMatrix <- AMatrix%*%solve(deltaMatrix)%*%AMatrix # 22.01.2020 new precision matrix, no spam (see 43.R)
matrixS <- array(parameters$S, dim = c(numWeeks, numRegions, numSequences))
matrixR <- aperm(array(parameters$R, dim = c(numRegions, numWeeks, numSequences)), perm = c(2,1,3))
if(dimBeta == 2){
  matrixB <- aperm(array(parameters$B[jToGroups], dim = c(numRegions, numWeeks, numSequences)), perm = c(2,1,3))
}else if(dimBeta == 3){
  matrixB <- aperm(array(parameters$B[kToGroups], dim = c(numSequences, numWeeks, numRegions)), perm = c(2,3,1))
}else if(dimBeta == 123){
  matrixB <- array(parameters$B[allToGroups], dim = c(numWeeks, numRegions, numSequences))
}
matrixG <- aperm(array(parameters$G, dim = c(numSequences, numWeeks, numRegions)), perm = c(2,3,1))
#matrixX <- replicate(numSequences, matrix(parameters$X[allToGroups], nrow = numWeeks, ncol = numRegions, byrow = FALSE))
matrixX <- array(parameters$X[allToGroups], dim = c(numWeeks, numRegions, numSequences))
#e.g. replicate(3, matrix(1:10, 5)). See 28_01.R. Careful with indices conventions in R
matrixRR <- matrixForGMRF*parameters$R
# TODO test all matrixSomething
# "If you try to add a vector to an array, then the usual vector recycling rules apply, but the dimension of the results is taken from the array".
matrixSexp <- exp(matrixS)
matrixRexp <- exp(matrixR)
matrixGexp <- exp(matrixG)
#matrixBexp <- exp(matrixB)
matrixXBexp <- exp(matrixX*matrixB)

# 5. Iterations ----
counter <- Sys.time()
numIterations <- config$numIterations
cat("HOLA", sum(y, na.rm = T), "\n")
for (it in 1:numIterations) {# 1:numIterations  (numIterations + 1):(2*numIterations)
  cat(it, "\n")
  x <- sqrt(2)
  #x <- 2/sqrt(it)
  
  # PARAMETERS
  # Update a
  source(paste0(dirFiles, "38_04_UpdateA.R"))
  
  # Update R
  source(paste0(dirFiles, "38_05_UpdateR.R"))
  
  # Update S
  source(paste0(dirFiles, "38_06_UpdateS.R")) # !
  
  # Update G
  source(paste0(dirFiles, "38_07_UpdateG.R"))
  
  # Update X
  source(paste0(dirFiles, "38_08_UpdateX.R"))
  
  # Update B
  source(paste0(dirFiles, "38_09_UpdateB.R"))
  
  # HYPERPARAMETERS
  # Update hyperparameters
  source(paste0(dirFiles, "38_10_UpdateHyperparameters.R"))
  
  # Add constrains
  #parameters$a <- parameters$a + mean(parameters$G)
  parameters$a <- parameters$a + mean(parameters$S)
  parameters$a <- parameters$a + mean(parameters$R)
  #parameters$G <- parameters$G - mean(parameters$G)
  parameters$S <- parameters$S - mean(parameters$S)
  parameters$R <- parameters$R - mean(parameters$R)
  
  # Store
  storage$parameters$a[it] <- parameters$a
  storage$parameters$G[,it] <- parameters$G
  storage$parameters$S[,it] <- parameters$S
  if(dimBeta == 123){
    storage$parameters$B <- storage$parameters$B + parameters$B # only the mean unfortunately
  }else{
    storage$parameters$B[,it] <- parameters$B
  }
  storage$parameters$R[,it] <- parameters$R
  storage$parameters$X[[it]] <- which(parameters$X == 1, arr.ind = TRUE)
  storage$parameters$tau.R[it] <- parameters$tau.R
  storage$parameters$tau.S[it] <- parameters$tau.S
  storage$parameters$tau.G[it] <- parameters$tau.G
  storage$parameters$p[it] <- parameters$p
  storage$parameters$l[it] <- parameters$l
  storage$parameters$cutHeight[it] <- parameters$cutHeight # 27.02.2020
  storage$parameters$sBlockSize[,it] <- parameters$sBlockSize # 17.03.2020
  
  # Store expected number of cases
  # TODO test
  # TODO it has not been implemented really because the storage capacity is insane if we wanna do it for every cell... better manually afterwards
  #if(1 %in% dimToInclude){ ???
  #  # Recall: matrixXBexp already exists
  #  matrixG <- aperm(array(parameters$G, dim = c(numSequences, numWeeks, numRegions)), perm = c(2,3,1))
  #  matrixS <- array(parameters$S, dim = c(numWeeks, numRegions, numSequences))
  #  matrixR <- aperm(array(parameters$R, dim = c(numRegions, numWeeks, numSequences)), perm = c(2,1,3))
  #  storage$parameters$scases[,it] <- apply(matrixPop*exp(parameters$a + matrixS + matrixR + matrixG), 1, sum)
  #  storage$parameters$ecases[,it] <- apply(matrixPop*exp(parameters$a + matrixS + matrixR + matrixG)*matrixXBexp, 1, sum)
  #}
  
  cat(sprintf("%.*f", 2, lll()), "\n")
  cat("\n")
}
runningTime <- counter - Sys.time()
print(runningTime)
#config$numIterations <- 220
save(inputForData, constants, config, out, storage, adaptive, it, runningTime,
     file = paste0(dirOutFiles, nameOutput)) # keep similar name with Input!
cat("Output stored successfully... in", paste0(dirOutFiles, nameOutput))

# TODO create function to extend/cut size of variables if numIterations changes (keep running after a certain point for instance)

#temp <- lapply(1:it, function(x) cbind(x, storage$parameters$X[[x]]))
#temp2 <- as.data.table(do.call(rbind, temp)) # Excpect a warning for each iteration without X's
#setnames(temp2, c("x","row","col"), c("it","week","region"))
#storage$parameters$X <- temp2[, .N, .(week, region)]
#save(inputForData, constants, config, out, storage, adaptive, it, runningTime,
#     file = paste0(dirFiles, "38_00_21112019_MCMCOutputTinis_short.RData")) # keep similar name with Input!
#cat("Output (shorter verson) stored successfully...")

