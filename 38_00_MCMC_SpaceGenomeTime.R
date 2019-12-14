# ------------------------------ #
# MCMC
# Data: FS101013 V2
# Created on 26.10.2019
# Adapted from the 28 files, created on 24.06.2019
# Changes: --.--.2019:
# ------------------------------ #

runTinis <- 1 # !!!!!!!!!!!!!!!!!!!!!
simulatedData <- 0

if(runTinis == 0){
  if(1){
    # 0. Load Workspace ----
    source("/home/laura/Dropbox/Laura/PhD_Year2/04_AboutApproaches/RCode/2018_04/00_Functions.R")
    setwd("/home/laura/Dropbox/Laura/PhD_Year3/")
    source("07_MixedModelsP2/RCode_201908/34_HeaderFor201908.R") # load workspace for August 2019
    load("07_MixedModelsP2/RCode_201909/33a_03092019_DataOXTW_Perfect_Depr.RData") # dataLSOA, dataMSOA, dataInfoPostSOA, chosenAreasOX, chosenAreasTW
    dirFiles <- "07_MixedModelsP2/RCode_201911/38_Files/"
  }
  if(1){
    # 1. Create/load data ----
    dimToInclude <- sort(c(2,3)) # T S G # Never do the three of them!
    lengthBlockPeriod <- 3
    dimBeta <- 3 # dimension where B is gonna change (2 or 3)
    heighCutLow <- 50
    heighCutHigh <- 300
    typeSpatial <- "OX"
    numClustSpatial <- "MSOA"
    if(simulatedData == 0){
      # Construct data
      #source(paste0(dirFiles, "38_01_MCMC_CreateLoadData_V2.R")) # ~ 8 sec
      # Or load (not formal)
      load("07_MixedModelsP2/RCode_201911/38_Files/38_01_21112019_MCMCInputTinis.RData") # CONFIG: notsim, dims23, interval3, beta3, cuts50/300, regionOXMSOA
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
nameOutput <- "38_00_0312201902_MCMCOutputTinis.RData"

# 2. Set-up environment ----
source(paste0(dirFiles, "38_02_AuxiliarFn.R"))
# Functions: all plot...
source(paste0(dirFiles, "38_03_ConstantsFn.R"))
# Functions: initConstants, initConfig, initParam, createStorage, adaptive

# 3. Auxiliar functions ----
source(paste0(dirFiles, "38_11_Prior_Likelihoods.R"))

# 4. Pre-iterations ----
cat("Sim", simulatedData, "/n")
constants <- initConstants(aP = 1, bP = numBlockDims[dimBeta] - 1, aB = 2, bB = 1)
#constants <- initConstants(aP = 1, bP = 50, aB = 2, bB = 1) # for fastSim27112019
#constants <- initConstants(aP = 1, bP = 5, aB = 2, bB = 1) # for 03122019
#config <- initConfig(numIterations = 1000, burnIn = 50) # estandar
#config <- initConfig(numIterations = 50, burnIn = 10) # quick exploration
#config <- initConfig(numIterations = 1000, burnIn = 50, ifAUpdate = 1, ifRUpdate = 0, ifSUpdate = 0, ifGUpdate = 0, ifBUpdate = 0, ifXUpdate = 0,
#                     ifPUpdate = 0, ifTauRUpdate = 0, ifTauSUpdate = 0, ifLTauGUpdate = 0)
config <- initConfig(numIterations = 5000, burnIn = 50,
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

# Requires: distanceMatrixMCMC pop matrixForGMRF y
# iToGroups jToGroups kToGroups
# numRegions numWeeks numSequences numWeeksGroups numRegionsGroups numSequenceGroups dimBeta

# TODO The following will eventually be in a function

# Auxiliar matrices
sqrDistanceMatrix <- distanceMatrixMCMC^2
if(1 %in% dimToInclude){
  seasonalCoefficientMatrix <- cbind(rep(0,numWeeks - 2), rep(0,numWeeks - 2), diag(numWeeks - 2)) +
    cbind(rep(0,numWeeks - 2), -2*diag(numWeeks - 2), rep(0,numWeeks - 2)) +
    cbind(diag(numWeeks - 2), rep(0,numWeeks - 2), rep(0,numWeeks - 2))
  seasonalCoefficientMatrixSqr <- t(seasonalCoefficientMatrix)%*%seasonalCoefficientMatrix
}
matrixPop <- aperm(array(pop, dim = c(numRegions, numWeeks, numSequences)), perm = c(2,1,3))
#neighLength <- rowSums(matrixForGMRF)
neighLength <- matrixForGMRF_nb[,1]

# Initial covariance matrix and other initial matrices
# TODO c++??
#deltaMatrix <- deltaMatrixFn(parameters$l, maternParameter) # + 0.5*diag(numSequences)
#invDeltaMatrix <- solve(deltaMatrix)
deltaMatrix <- as.spam(deltaMatrixFn(parameters$l, distanceMatrixMCMC, config$maternParameter) + 0.01*diag(numSequences))
invDeltaMatrix <- solve(deltaMatrix) # TODO works always as a matrix?
matrixS <- array(parameters$S, dim = c(numWeeks, numRegions, numSequences))
matrixR <- aperm(array(parameters$R, dim = c(numRegions, numWeeks, numSequences)), perm = c(2,1,3))
if(dimBeta == 2){
  matrixB <- aperm(array(parameters$B[jToGroups], dim = c(numRegions, numWeeks, numSequences)), perm = c(2,1,3))
}else if(dimBeta == 3){
  matrixB <- aperm(array(parameters$B[kToGroups], dim = c(numSequences, numWeeks, numRegions)), perm = c(2,3,1))
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
cat("HOLA", sum(y), "\n")
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
  parameters$G <- parameters$G - mean(parameters$G)
  parameters$S <- parameters$S - mean(parameters$S)
  parameters$R <- parameters$R - mean(parameters$R)
  
  # Store
  storage$parameters$a[it] <- parameters$a
  storage$parameters$G[,it] <- parameters$G
  storage$parameters$S[,it] <- parameters$S
  storage$parameters$B[,it] <- parameters$B
  storage$parameters$R[,it] <- parameters$R
  storage$parameters$X[[it]] <- which(parameters$X == 1, arr.ind = TRUE)
  storage$parameters$tau.R[it] <- parameters$tau.R
  storage$parameters$tau.S[it] <- parameters$tau.S
  storage$parameters$tau.G[it] <- parameters$tau.G
  storage$parameters$p[it] <- parameters$p
  storage$parameters$l[it] <- parameters$l
  
  # Store expected number of cases
  # TODO test
  if(1 %in% dimToInclude){
    matrixG <- aperm(array(parameters$G, dim = c(numSequences, numWeeks, numRegions)), perm = c(2,3,1))
    matrixS <- array(parameters$S, dim = c(numWeeks, numRegions, numSequences))
    matrixR <- aperm(array(parameters$R, dim = c(numRegions, numWeeks, numSequences)), perm = c(2,1,3))
    storage$parameters$scases[,it] <- apply(matrixPop*exp(parameters$a + matrixS + matrixR + matrixG), 1, sum)
    storage$parameters$ecases[,it] <- apply(matrixPop*exp(parameters$a + matrixS + matrixR + matrixG)*matrixXBexp, 1, sum)
  }
  
  cat("\n")
}
runningTime <- counter - Sys.time()
print(runningTime)
save(inputForData, constants, config, out, storage, adaptive, it, runningTime,
     file = paste0(dirFiles, nameOutput)) # keep similar name with Input!
cat("Output stored successfully...")

#temp <- lapply(1:it, function(x) cbind(x, storage$parameters$X[[x]]))
#temp2 <- as.data.table(do.call(rbind, temp)) # Excpect a warning for each iteration without X's
#setnames(temp2, c("x","row","col"), c("it","week","region"))
#storage$parameters$X <- temp2[, .N, .(week, region)]
#save(inputForData, constants, config, out, storage, adaptive, it, runningTime,
#     file = paste0(dirFiles, "38_00_21112019_MCMCOutputTinis_short.RData")) # keep similar name with Input!
#cat("Output (shorter verson) stored successfully...")

