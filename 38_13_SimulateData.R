# 0. Load Workspace ----
source("/home/laura/Dropbox/Laura/PhD_Year2/04_AboutApproaches/RCode/2018_04/00_Functions.R")
setwd("/home/laura/Dropbox/Laura/PhD_Year3/")
#source("07_MixedModelsP2/RCode_201911/25_HeaderFor201907.R") # load workspace for July 2019

# Create simulation after running MCMC and storing output
if(0){
  #old
  simulatedData <- 0
  source("07_MixedModelsP2/RCode_201907/28_Files/28_01_MCMC_CreateLoadData.R")
  load("07_MixedModelsP2/RCode_201907/28_Files/28_00_26062019_ForSimulation.RData") # constants, config, out, storage, adaptive, it
  #RUN 28_01.R/III
}else{
  #load("07_MixedModelsP2/RCode_201911/38_Files/38_01_21112019_MCMCInputTinis.RData") # CONFIG: notsim, dims23, interval3, beta3, cuts50/300, regionOXMSOA
  load("07_MixedModelsP2/RCode_201911/38_Files/38_01_21112019_MCMCInputTinis.RData") # CONFIG: notsim, dims23, interval3, beta3, cuts50/300, regionOXMSOA
  #load("/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_201911/38_00_21112019_MCMCOutputTinis.RData")
  # constants, config, out, storage, adaptive, it, runningTime
  #load("/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_201911/38_00_27112019fastSim_MCMCOutputTinis.RData")
  load("/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_201911/38_00_03122019_MCMCOutputDellToSim.RData")
}

finalIterations <- 1001:5000
finalIterations <- (config$burnIn + 1):config$numIterations

outputA <- storage$parameters$a[finalIterations]
outputG <- storage$parameters$G[,finalIterations]
outputS <- storage$parameters$S[,finalIterations]
outputS <- matrix(storage$parameters$S[,finalIterations], nrow = 1) # TODO
outputR <- storage$parameters$R[,finalIterations]
outputB <- storage$parameters$B[,finalIterations]
outputTauG <- storage$parameters$tau.G[finalIterations]
outputTauS <- storage$parameters$tau.S[finalIterations]
outputTauR <- storage$parameters$tau.R[finalIterations]
outputP <- storage$parameters$p[finalIterations]
outputL <- storage$parameters$l[finalIterations]
temp <- lapply(finalIterations, function(x) cbind(x, storage$parameters$X[[x]]))
outputX <- as.data.table(do.call(rbind, temp)) # Excpect a warning for each iteration without X's
#nrow(outputX)/(numWeeksGroups*numRegionsGroups)
setnames(outputX, c("x","row","col"), c("it","week","region"))

matX <- array(0, dim = c(1, numBlockDims[dimToInclude[1]], numBlockDims[dimToInclude[2]]))
tempX <- outputX[, .N/length(finalIterations), by = .(week, region)]
cat("Choose a threshold accordingly:")
# plot(sort(tempX$V1))
# OLD matX[as.matrix(tempX[, .(1, week, region)])] <- tempX$V1 > 0.6 # 0.7(21112019) NO0.6(03122019)
matX[as.matrix(tempX[, .(1, week, region)])] <- runif(nrow(tempX)) < tempX$V1
cat("Proba:", 100*sum(matX)/(numBlockDims[dimToInclude[1]]*numBlockDims[dimToInclude[2]]), "%")

#y
matrixS <- array(rowMeans(outputS), dim = c(numWeeks, numRegions, numSequences))
matrixR <- aperm(array(rowMeans(outputR), dim = c(numRegions, numWeeks, numSequences)), perm = c(2,1,3))
if(dimBeta == 2){
  matrixB <- aperm(array(rowMeans(outputB)[jToGroups], dim = c(numRegions, numWeeks, numSequences)), perm = c(2,1,3))
}else if(dimBeta == 3){
  matrixB <- aperm(array(rowMeans(outputB)[kToGroups], dim = c(numSequences, numWeeks, numRegions)), perm = c(2,3,1))
}
matrixG <- aperm(array(rowMeans(outputG), dim = c(numSequences, numWeeks, numRegions)), perm = c(2,3,1))
#OLDmatrixX <- replicate(numSequences, matrix(matX[allToGroups], nrow = numWeeks, ncol = numRegions, byrow = FALSE))
matrixX <- array(matX[allToGroups], dim = c(numWeeks, numRegions, numSequences))
matrixSexp <- exp(matrixS)
matrixRexp <- exp(matrixR)
matrixGexp <- exp(matrixG)
matrixXBexp <- exp(matrixX*matrixB)
matrixPop <- aperm(array(pop, dim = c(numRegions, numWeeks, numSequences)), perm = c(2,1,3))

matrixLambda <- matrixPop*exp(mean(outputA))*matrixSexp*matrixRexp*matrixGexp*matrixXBexp
#plot(sort(sample(matrixLambda, 1000)))

yy <- apply(matrixLambda, 1:3, function(x) rpois(1, x))
sum(yy) # 100*1242179/prod(dim(yy)) 1242179
sum(matrixLambda)

simulatedParam <- list(a = mean(outputA), R = rowMeans(outputR), S = rowMeans(outputS), G = rowMeans(outputG),
                       X = matX, B = rowMeans(outputB), tau.R = mean(outputTauR), tau.S = mean(outputTauS),
                       tau.G = mean(outputTauG), p = mean(outputP), l = mean(outputL))
#save(yy, simulatedParam, file = "07_MixedModelsP2/RCode_201907/28_Files/28_13_26062019_SimulatedData.RData")
#save(yy, simulatedParam, file = "007_MixedModelsP2/RCode_201911/38_Files/38_13_21112019_SimulatedData.RData") # 24112019
#save(yy, simulatedParam, file = "007_MixedModelsP2/RCode_201911/38_Files/38_13_fastSim27112019_SimulatedData.RData") # Simons meeting
#save(yy, simulatedParam, file = "07_MixedModelsP2/RCode_201911/38_Files/38_13_03122019_SimulatedData.RData")

#pop[5]*exp(mean(outputA) +mean(outputS[4,]) + mean(outputR[5,]) + mean(outputG[3,]) + mean(outputB[jToGroups[5],])*mean(matX[iToGroups[4],jToGroups[5]]))
#matrixLambda[4,5,3]
#pop[6]*exp(mean(outputA) +mean(outputS[4,]) + mean(outputR[6,]) + mean(outputG[3,]) + mean(outputB[jToGroups[6],])*mean(matX[iToGroups[4],jToGroups[6]]))
#matrixLambda[4,6,3]



# Now: RUN simulation in 28_00.R, change value of simulatedData. Fix whatever you wanna fix in initConfig. Check manually or do plots
burninEnd <- 1 # 1 51 21
finalIterations <- burninEnd:(it - 1)
outputA <- storage$parameters$a[finalIterations]
outputG <- storage$parameters$G[,finalIterations]
outputS <- storage$parameters$S[,finalIterations]
outputR <- storage$parameters$R[,finalIterations]
outputB <- storage$parameters$B[,finalIterations]
outputTauG <- storage$parameters$tau.G[finalIterations]
outputTauS <- storage$parameters$tau.S[finalIterations]
outputTauR <- storage$parameters$tau.R[finalIterations]
outputP <- storage$parameters$p[finalIterations]
outputL <- storage$parameters$l[finalIterations]
temp <- lapply(finalIterations, function(x) cbind(x, storage$parameters$X[[x]]))
outputX <- as.data.table(do.call(rbind, temp)) # Excpect a warning for each iteration without X's
setnames(outputX, c("x","row","col"), c("it","week","region"))
matX <- matrix(0, nrow = numWeeksGroups, ncol = numRegionsGroups)
tempX <- outputX[, .N/length(finalIterations), by = .(week, region)]
matX[as.matrix(tempX[, .(week, region)])] <- tempX$V1

plot(mean(outputA), simulatedParam$a)
plot(rowMeans(outputR), simulatedParam$R)
plot(rowMeans(outputS), simulatedParam$S)
plot(rowMeans(outputG), simulatedParam$G)
plot(rowMeans(outputB), simulatedParam$B) # fixed! problem with the index j in lpriorB and in the update of matrixB and matrixXBexp :D
plot(c(matX), c(simulatedParam$X)) # :D

plot(mean(outputTauR), simulatedParam$tau.R) # ??
plot(mean(outputTauS), simulatedParam$tau.S) # ??
plot(mean(outputTauG), simulatedParam$tau.G)
plot(mean(outputP), simulatedParam$p) # fixed? pop was missing in the X update
plot(mean(outputL), simulatedParam$l) # rejecting everyone # BUG with the recent update of deltaMatrixFn!!!

plot(outputL)
plot(outputTauG)

ttt <- replicate(numSequences, matrix(simulatedParam$X[allToGroups], nrow = numWeeks, ncol = numRegions, byrow = FALSE))
sum(ttt)/prod(dim(ttt))

storage$accept$a/(storage$accept$a + storage$reject$a)
summary(storage$accept$G/(storage$accept$G + storage$reject$G))
summary(storage$accept$ScondUp/(storage$accept$ScondUp + storage$reject$ScondUp))
summary(storage$accept$S/(storage$accept$S + storage$reject$S))
summary(storage$accept$RcondUp/(storage$accept$RcondUp + storage$reject$RcondUp))
summary(storage$accept$R/(storage$accept$R + storage$reject$R))
summary(storage$accept$B/(storage$accept$B + storage$reject$B))
storage$accept$l/(storage$accept$l + storage$reject$l)


