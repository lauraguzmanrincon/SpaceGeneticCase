# All sections require 0. 1. 2.
# Never re-run 0. (may overwrite some variables)
# Run to 2. and 3. if MCMC just finished and wanna check output

if(0){
  # 0. Load workspace ----
  source("/home/laura/Documents/PhD_Year2_NoDropbox/PhD_Year2_Dell/04_AboutApproaches/RCode/2018_04/00_Functions.R")
  source("/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201908/34_HeaderFor201908.R") # header included?
  setwd("/home/laura/Dropbox/Laura/PhD_Year3/")
  allPlotsPath <- "07_MixedModelsP2"
  
  uploadPackages("latex2exp") # for TeX
  source(paste0("07_MixedModelsP2/RCode_201911/38_Files/", "38_02_AuxiliarFn.R"))
  load("07_MixedModelsP2/RCode_201909/33a_03092019_DataOXTW_Perfect_Depr.RData") # dataLSOA, dataMSOA, dataInfoPostSOA, chosenAreasOX, chosenAreasTW
  
  # For map and colours
  #uploadPackages(c("rgdal", "scales"))
  tempDir <- "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/Areafiles_201907/"
  tempName <- paste0(tempDir, "15_Local_Authority_Districts_December_2018_Boundaries_GB_BFC")
  tempMap <- rgdal::readOGR(dsn = tempName, layer = "Local_Authority_Districts_December_2018_Boundaries_GB_BFC")
  data.shapeLADS <- sp::spTransform(tempMap, CRS("+proj=longlat +datum=WGS84"))
  
  # REMOVE following after including header before?
  # For groups and clusters
  #load("/home/laura/Dropbox/Laura/PhD_Year2/06_MixedModels/RCode_201903/15_Input_MCMCCorrected04032019.RData")
  # Content: numSequences, numWeeks, y, distanceMatrixGroups, groupSequencesCore
  # Recall: column order in distanceMatrixGroups is given by groupSequencesCore$groupId
  #rm(numWeeks, y, numSequences)
  #suppressWarnings(groupSequencesCore[, numWeek := NULL])
  #suppressWarnings(groupSequencesCore[, weekId := NULL])
  
  # General data
  #load("/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/Data/FS101013_Datasets_V3_26R.RData") # ...
}

# 1. Upload MCMC output file + others ----
#dirFiles <- "FilesFromLocal/38_Files/"
#dirFiles <- "07_MixedModelsP2/RCode_201911/38_Files/"
#dirFiles <- "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_201911/"
# --
#nameFiles <- "21112019"
#nameFilesOut <- "0312201901"
#simulatedDataRes <- 0
#typ <- "OX"
#dirInputFiles <- "07_MixedModelsP2/RCode_201911/38_Files/"
#dirOutputFiles <- "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_201911/"
# --
#nameFiles <- "18122019_MCMCInput_GOX"
#nameFilesOut <- "1812201901_MCMCOutput_GOX"
#dirInputFiles <- "07_MixedModelsP2/RCode_201912/"
#dirOutputFiles <- "07_MixedModelsP2/RCode_201912/"
# -- Short simulation on Dell
#nameFiles <- "18122019_MCMCInput_GOX"
#nameFilesOut <- "0601202001_MCMCOutput_newprior_SE" ######### 0201202001_MCMCOutput_newprior 0601202001_MCMCOutput_newprior_SE
#dirInputFiles <- "07_MixedModelsP2/RCode_201912/"
#dirOutputFiles <- "07_MixedModelsP2/RCode_201912/"
# -- Large simulations on Tinis, three kernels
#nameFiles <- "18122019_MCMCInput_GOX"
#nameFilesOut <- "0601202001_MCMCOutputTinis_MAT" # tauG: 1,0.01 rho: 10,0.2 r: 10/50* # 798311
#nameFilesOut <- "0601202001_MCMCOutputTinis_SE" # tauG: 1,0.01 rho: 10,0.2 r: 10/50* # 798312
#nameFilesOut <- "0601202001_MCMCOutputTinis_MAT32" # tauG: 1,0.01 rho: 10,0.2 r: 10/50* # 798341
#dirInputFiles <- "07_MixedModelsP2/RCode_201912/"
#dirOutputFiles <- "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202001/"
# -- Short simulation with Bijk USELESS
#nameFilesOut <- "0701202001_MCMCOutputTinis_blockBeta" # tauG: 1,0.01 rho: 10,0.2 r: 10/50 kernel: 3/2 # too slow beta # 7983496
#dirInputFiles <- "07_MixedModelsP2/RCode_201912/"
# -- Short simulation with B bug corrected 07.01
#nameFiles <- "18122019_MCMCInput_GOX"
#nameFilesOut <- "temp_MAT32"
#nameFilesOut <- "temp_SE_Tinis"
#dirInputFiles <- "07_MixedModelsP2/RCode_201912/"
#dirOutputFiles <- "07_MixedModelsP2/RCode_202001/"
# -- Large simulations on Tinis, three kernels, constraint correction
nameFiles <- "18122019_MCMCInput_GOX"
nameFilesOut <- "2601202002_Tinis_MAT12_Constraint" # tauG: 1,0.01 rho: 10,0.2 r: 10/50* # 798311
nameFilesOut <- "2601202001_Tinis_SE_Constraint" # tauG: 1,0.01 rho: 10,0.2 r: 10/50* # 798312
nameFilesOut <- "2601202001_Tinis_MAT32_Constraint" # tauG: 1,0.01 rho: 10,0.2 r: 10/50* # 798341
dirInputFiles <- "07_MixedModelsP2/RCode_201912/"
dirOutputFiles <- "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202001/"
# Other parameters
simulatedDataRes <- 0
typ <- "OX"
dirInputFile <- paste0(dirInputFiles, "38_01_", nameFiles, ".RData")
load(dirInputFile) # ...
if(simulatedDataRes == 0){
  dirOutputFile <- paste0(dirOutputFiles, "38_00_", nameFilesOut, ".RData")
  load(dirOutputFile) # constants, config, out, storage, adaptive, it, runningTime
}else{
  #nameOutput <- "38_00_21112019_SimulOutputTinis_a.RData"
  nameOutput <- "38_00_21112019_SimulOutputTinis_r.RData"
  #nameOutput <- "38_00_21112019_SimulOutputTinis_g.RData"
  #nameOutput <- "38_00_21112019_SimulOutputTinis_b.RData"
  #nameOutput <- "38_00_21112019_SimulOutputTinis_x.RData"
  #nameOutput <- "38_00_21112019_SimulOutputTinis_p.RData"
  #nameOutput <- "38_00_21112019_SimulOutputTinis_tr.RData"
  #nameOutput <- "38_00_21112019_SimulOutputTinis_ltg.RData"
  dirOutputFile <- paste0(dirOutputFiles, nameOutput)
  
  #dirOutputFile <- paste0(dirOutputFiles, "38_00_", nameFiles, "_SimulOutputTinis.RData")
  load(dirOutputFile) # constants, config, out, storage, adaptive, it, runningTime
  dirSimOutputFile <- paste0(dirInputFiles, "38_13_", nameFiles, "_SimulatedData.RData")
  load(dirSimOutputFile) # yy, simulatedParam 24112019
}

# Other updates
if(typ == "OX"){
  chosenAreas <- chosenAreasOX
  mapLims <- list(x = c(-1.8,-0.7), y = c(51.44, 52.25), xE = c(-1.3,-1.165), yE = c(51.71, 51.8))
}else{
  chosenAreas <- chosenAreasTW
  mapLims <- list(x = c(-2.66,-1.36), y = c(54.8, 55.8), xE =  c(-1.41,-1.78), yE = c(54.95, 55.1))
}

# Overview
if(0){
  str(inputForData)
  str(constants)
  str(config)
  str(out)
  str(storage, max.level = 2)
  str(adaptive)
  runningTime
}

# 2. Build data from MCMC output for visualisations ----
# a. Traces
#finalIterations <- 1001:5000
finalIterations <- (config$burnIn + 1):config$numIterations
#finalIterations <- 1:75
#finalIterations <- 1:config$numIterations
#finalIterations <- 1:(it - 1)
sizePlot <- 15
simulatedDataRes <- 0

# Check output
#gIndices <- union((freqGUpdate != 1)*which((finalIterations)%%freqGUpdate == 1), (freqGUpdate == 1)*(finalIterations))
outputA <- storage$parameters$a[finalIterations]
outputG <- storage$parameters$G[,finalIterations]
outputS <- array(storage$parameters$S[,finalIterations], dim = c(numPeriods, length(finalIterations))) # just in case is of size 1
outputR <- array(storage$parameters$R[,finalIterations], dim = c(numRegions, length(finalIterations)))
if(dimBeta != 123) outputB <- storage$parameters$B[,finalIterations]
if(dimBeta == 123) outputB <- storage$parameters$B/config$numIterations
outputTauG <- storage$parameters$tau.G[finalIterations]
outputTauS <- storage$parameters$tau.S[finalIterations]
outputTauR <- storage$parameters$tau.R[finalIterations]
outputP <- storage$parameters$p[finalIterations]
outputL <- storage$parameters$l[finalIterations]
temp <- lapply(finalIterations, function(x) cbind(x, storage$parameters$X[[x]]))
outputX <- as.data.table(do.call(rbind, temp)) # Excpect a warning for each iteration without X's
setnames(outputX, c("x"), c("it"))
#setnames(outputX, c("x","row","col"), c("it","week","region"))
nrow(outputX)/(prod(numBlockDims)*length(finalIterations))

# b. Outbreak info
THRESHOLD <- 0.1
outProbs <- outputX[, .N/length(finalIterations), .(row, col)]
outProbs[, id := .I]
setnames(outProbs, "V1", "proba")
outProbs[, sizeB := rowMeans(outputB)[col]]

setkey(casesForModels, dim2Cases)
setkey(regionsForModels, col)
casesForModels[regionsForModels, region := region]
setkeyv(casesForModels, c("region", "clusterHighId"))
setkeyv(outProbs, c("row", "col"))
casesForModels[outProbs, c("idOutbreak", "probaOutbreak") := .(id, proba)]
casesForModels[, sizeOutbreak := .N, idOutbreak]

# 3. MCMC output: quick chain visualisation to overview MCMC ----
# Requires 0. 1. 2.
# Check acceptance rate ----
# Recall: it's 1 for X, tauS, tauG, p
data.table(names = names(storage$accept),
           do.call(rbind, lapply(1:length(storage$accept), function(x) summary(storage$accept[[x]]/(storage$accept[[x]] + storage$reject[[x]])))))

# Plot traces ----
#BORRAR
#sampleSPlot <- sample(1:numWeeks, 1)
#ggplot(data.table(it = finalIterations, value = outputS[sampleSPlot,]), aes(x = it, y = value)) +
#  geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $S_j$")) +
#  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$S[sampleSPlot], NaN))
#sampleRPlot <- sample(1:numRegions, 1)
#ggplot(data.table(it = finalIterations, value = outputR[sampleRPlot,]), aes(x = it, y = value)) +
#  geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $R_j$")) +
#  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$R[sampleRPlot], NaN))
#sampleGPlot <- sample(1:numSequences, 1)
#ggplot(data.table(it = finalIterations, value = outputG[sampleGPlot,]), aes(x = it, y = value)) +
#  geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $G_k$")) +
#  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$G[sampleGPlot], NaN))
#sigmaJumps$G[sampleGPlot]
#storage$accept$G[sampleGPlot]/(storage$accept$G[sampleGPlot] + storage$reject$G[sampleGPlot])
#sum(distanceMatrixMCMC[,sampleGPlot] < 100)
#sum(y[,,sampleGPlot])

sampleSPlot <- sample(1:numWeeks, 1)
sampleRPlot <- sample(1:numRegions, 1)
#sampleRPlot <- 257
#ggplot(data.table(it = finalIterations, value = outputR[sampleRPlot,]), aes(x = it, y = value)) +
#  geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $S_j$")) +
#  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$R[sampleRPlot], NaN))
sampleBPlot <- sample(1:numRegionsGroups, 1)
sampleGPlot <- sample(1:numSequences, 1)
#sampleGPlot <- 28
#ggplot(data.table(it = finalIterations, value = outputG[sampleGPlot,]), aes(x = it, y = value)) +
#  geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $G_k$")) +
#  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$G[sampleGPlot], NaN))

rm(g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)
g0 <- ggplot(data.table(it = finalIterations, value = outputA),
             aes(x = it, y = value)) + geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $\\alpha$")) +
  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$a, NaN))
g1 <- ggplot(data.table(it = finalIterations, value = outputTauG),
             aes(x = it, y = value)) + geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $\\tau_G$"))
g2 <- ggplot(data.table(it = finalIterations, value = outputTauS),
             aes(x = it, y = value)) + geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $\\tau_R$"))
g9 <- ggplot(data.table(it = finalIterations, value = outputTauR),
             aes(x = it, y = value)) + geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $\\tau_U$"))
g3 <- ggplot(data.table(it = finalIterations, value = outputP),
             aes(x = it, y = value)) + geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $p$")) +
  #geom_hline(yintercept = ifelse(simulatedDataRes, sum(simulatedParam$X)/prod(dim(simulatedParam$X)), NaN))
  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$p, NaN))
g4 <- ggplot(data.table(it = finalIterations, value = outputS[sampleSPlot,]), aes(x = it, y = value)) +
  geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $R_j$")) +
  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$S[sampleSPlot], NaN))
g10 <- ggplot(data.table(it = finalIterations, value = outputR[sampleRPlot,]), aes(x = it, y = value)) +
  geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $U_j$")) +
  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$R[sampleRPlot], NaN))
g5 <- ggplot(data.table(it = finalIterations, value = outputB[sampleBPlot,]), aes(x = it, y = value)) +
  geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $B_k$")) +
  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$B[sampleBPlot], NaN))
g6 <- ggplot(data.table(it = finalIterations, value = outputL),
             aes(x = it, y = value)) + geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $\\rho$"))
g7 <- ggplot(data.table(it = finalIterations, value = outputG[sampleGPlot,]), aes(x = it, y = value)) +
  geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $G_k$")) +
  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$G[sampleGPlot], NaN))
g8 <- ggplot(outputX[, .N/prod(numBlockDims), by = .(it)], aes(x = it, y = V1)) +
  geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $X_{jk}$")) +
  geom_hline(yintercept = ifelse(simulatedDataRes, mean(simulatedParam$X), NaN))
if(2 %in% dimToInclude & 3 %in% dimToInclude){
  #multiplot(g4,g10,g7,g5,g0, g2,g9,g1,g6,g3, cols = 2) # Sj Rj Gk Bk a // tauS tauR tauG ro p
  multiplot(g10,g7,g5,g0, g9,g1,g6,g3, cols = 2) # Sj Rj Gk Bk a // tauS tauR tauG ro p
  #g8
  #g3
  if(dimBeta == 123) multiplot(g10,g7,ggplot(),g0, g9,g1,g6,g3, cols = 2) 
}
#savePDF(multiplot(g10,g7,g5,g0, g9,g1,g6,g3, cols = 2), "Plot25112019_01_MCMCoutput_Real", 8, 10)#12)
#savePDF(multiplot(g4,g7,g5,g0,g2,g1,g6,g3, cols = 2), "Plot25112019_01_MCMCoutput_Real", 8, 10)#12)
#savePDF(multiplot(g7,g1,g6, cols = 1), "Plot27112019_01_MCMCoutput_RealG", 5, 7)
#savePDF(multiplot(g10,g7,g5,g0, g9,g1,g6,g3, cols = 2), "Plot27112019_01_MCMCoutput_Real2", 8, 10)#12)

#plot(outputX[, .N/prod(numBlockDims), by = .(it)]$V1, outputP)

# Posterior G R S B ----
#outputToPlot <- outputB
dataToPlotGSB <- rbind(data.table(parameter = 1:nrow(outputG), mean = rowMeans(outputG), sd = apply(outputG, 1, sd),
                                  q025 = apply(outputG, 1, quantile, 0.025),
                                  q975 = apply(outputG, 1, quantile, 0.975),
                                  y = apply(y, 3, sum),
                                  meanExp = apply(outputG, 1, function(x) mean(exp(x))),
                                  label = "G", minQ = NA, maxQ = NA, maxI = NA, meanQ = NA,
                                  real = ifelse(simulatedDataRes, simulatedParam$G, NaN),
                                  jump = out$sigmaJumps$G),
                       data.table(parameter = 1:nrow(outputB), mean = rowMeans(outputB), sd = apply(outputB, 1, sd),
                                  q025 = apply(outputB, 1, quantile, 0.025),
                                  q975 = apply(outputB, 1, quantile, 0.975),
                                  y = 0, #apply(y, 2, sum),
                                  meanExp = apply(outputB, 1, function(x) mean(exp(x))),
                                  label = "B", minQ = qgamma(0.025, shape = constants$aB, rate = constants$bB),
                                  maxQ = qgamma(0.975, shape = constants$aB, rate = constants$bB),
                                  maxI = numSequences, meanQ = qgamma(0.5, shape = constants$aB, rate = constants$bB),
                                  real = ifelse(simulatedDataRes, simulatedParam$B, NaN),
                                  jump = out$sigmaJumps$B))
if(1 %in% dimToInclude){
  dataToPlotGSB <- rbind(dataToPlotGSB, data.table(parameter = 1:nrow(outputS), mean = rowMeans(outputS), sd = apply(outputS, 1, sd),
                                                   q025 = apply(outputS, 1, quantile, 0.025),
                                                   q975 = apply(outputS, 1, quantile, 0.975),
                                                   y = apply(y, 1, sum),
                                                   meanExp = apply(outputS, 1, function(x) mean(exp(x))),
                                                   label = "S", minQ = NA, maxQ = NA, maxI = NA, meanQ = NA,
                                                   real = ifelse(simulatedDataRes, simulatedParam$S, NaN),
                                                   jump = out$sigmaJumps$S))
}
if(2 %in% dimToInclude){
  dataToPlotGSB <- rbind(dataToPlotGSB, data.table(parameter = 1:nrow(outputR), mean = rowMeans(outputR), sd = apply(outputR, 1, sd),
                                                   q025 = apply(outputR, 1, quantile, 0.025),
                                                   q975 = apply(outputR, 1, quantile, 0.975),
                                                   y = apply(y, 2, sum),
                                                   meanExp = apply(outputR, 1, function(x) mean(exp(x))),
                                                   label = "R", minQ = NA, maxQ = NA, maxI = NA, meanQ = NA,
                                                   real = ifelse(simulatedDataRes, simulatedParam$R, NaN),
                                                   jump = out$sigmaJumps$R))
}

#dataToPlot[, parameterName := parameter]
#dataToPlotGSB[label == "G", orderedParam := factor(parameter, levels = unique(dataToPlotGSB[label == "G"][order(mean), parameter]))] # not the best...
# ... ordering x axis
dataToPlotGSB[label != "S", orderedParam := rank(mean), by = .(label)]
dataToPlotGSB[label == "S", orderedParam := as.double(parameter)]
gGSB <- ggplot(dataToPlotGSB, aes(x = orderedParam, y = mean, ymin = mean - sd, ymax = mean + sd)) +
  geom_errorbar(color = "#919994") +
  geom_point(size = 1) + facet_wrap(~ label, ncol = 1, scales = "free", strip.position = "right") +
  #geom_rect(data = dataToPlotGSB[parameter == 1], aes(xmin = 0, xmax = maxI, ymin = minQ, ymax = maxQ), alpha = 0.3, fill = "red") +
  #geom_hline(aes(yintercept = meanQ), color = "red") +
  theme_laura(size = sizePlot) + theme(axis.text.x = element_blank()) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "index (not to be compared)", y = TeX("mean and $(p_{0.025}, p_{0.975})$ interval")) # TeX("mean $\\pm$ ")
#gGSB #savePDF(gGSB, "Plot16042019_MCMCpgsb_Real", 6, 5)
if(simulatedDataRes){
  gGSB + geom_point(aes(y = real), color = "red", size = 1)
}else{
  gGSB
}

if(0){
  # Plots for meeting 27.11.2019
  gGSB <- ggplot(dataToPlotGSB[label == "G"], aes(x = orderedParam, y = mean, ymin = mean - sd, ymax = mean + sd)) +
    geom_errorbar(color = "#919994") +
    geom_point(size = 1) + facet_wrap(~ label, ncol = 1, scales = "free", strip.position = "right") +
    #geom_rect(data = dataToPlotGSB[parameter == 1], aes(xmin = 0, xmax = maxI, ymin = minQ, ymax = maxQ), alpha = 0.3, fill = "red") +
    #geom_hline(aes(yintercept = meanQ), color = "red") +
    theme_laura(size = sizePlot) + theme(axis.text.x = element_blank()) +
    scale_x_continuous(expand = c(0,0)) +
    labs(x = "index (not to be compared)", y = TeX("mean and $(p_{0.025}, p_{0.975})$ interval")) + # TeX("mean $\\pm$ ")
    geom_point(aes(y = real), color = "red", size = 1)
  #savePDF(gGSB, "Plot27112019_MCMCpgsb_Real", 8, 5)
}

# Prior-posterior hyperparameters ----
# a, TauG, TauS, rho, Bk, p
aValToPlot <- plotNorm(mean = constants$ma, sd = constants$sa, xlab = "", returnPlot = 0)
gValToPlot <- plotGamma(shape = constants$aG, rate = constants$bG, xlab = "", returnPlot = 0)
lValToPlot <- plotGamma(shape = constants$aL, rate = constants$bL, xlab = "", returnPlot = 0)
pValToPlot <- plotBeta(shape1 = constants$aP, shape2 = constants$bP, xlab = "", returnPlot = 0)
dataToPlotPos <- rbind(
  data.table(x = gValToPlot$x, y = gValToPlot$y, label = "prior", var = "tauG"),
  data.table(x = outputTauG, y = 0, label = "posterior", var = "tauG"),
  data.table(x = lValToPlot$x, y = lValToPlot$y, label = "prior", var = "rho"),
  data.table(x = outputL, y = 0, label = "posterior", var = "rho"),
  data.table(x = pValToPlot$x, y = pValToPlot$y, label = "prior", var = "p"),
  data.table(x = outputP, y = 0, label = "posterior", var = "p")
  )
if(1 %in% dimToInclude){
  sValToPlot <- plotGamma(shape = constants$aS, rate = constants$bS, xlab = "", returnPlot = 0)
  dataToPlotPos <- rbind(dataToPlotPos,
                         data.table(x = sValToPlot$x, y = sValToPlot$y, label = "prior", var = "tauS"),
                         data.table(x = outputTauS, y = 0, label = "posterior", var="tauS"))
}
if(2 %in% dimToInclude){
  rValToPlot <- plotGamma(shape = constants$aR, rate = constants$bR, xlab = "", returnPlot = 0)
  dataToPlotPos <- rbind(dataToPlotPos,
                         data.table(x = rValToPlot$x, y = rValToPlot$y, label = "prior", var = "tauR"),
                         data.table(x = outputTauR, y = 0, label = "posterior", var="tauR"))
}
dataToPlotPos[, title := paste0(var, " - ", label)]
gPos <- ggplot(dataToPlotPos, aes(x = x, fill = label, color = label)) + facet_wrap(~ title, scales = "free", ncol = 2) +
  geom_density(data = dataToPlotPos[label == "posterior"], alpha = 0.3) +
  geom_line(aes(y = y)) + theme_laura() + labs(x = "", y = "density", color = "", fill = "")
gPos #savePDF(gPos, "Plot25112019_02_MCMCposterior_Real", 8, 7)

# Other exploratory plots
multiplot(plotGamma(constants$aB, constants$bB), ggplot(data.table(rowMeans(outputB)), aes(V1)) + geom_histogram(binwidth = 0.1), cols = 2)
ggplot(data.table(outputTauG, outputL), aes(x = outputTauG, y = outputL)) + geom_path() + geom_point() + theme_laura()
if(config$maternParameter == 1/2){
  plot(1:50, exp(-seq(1,50,1)/mean(outputL))/mean(outputTauG))
}else if(config$maternParameter == 3/2){
  plot(1:50, (1 + (sqrt(3)*seq(1,50,1)/mean(outputL)))*exp(-sqrt(3)*seq(1,50,1)/mean(outputL))/mean(outputTauG))
}else if(config$maternParameter == 0){
  plot(1:50, exp(-seq(1,50,1)^2/(2*mean(outputL)^2))/mean(outputTauG))
}
  
# 4. MCMC output: formal chain visualisation (table, plots) ----
# Requires: 0. 1. 2.
# Plot traces* ----
# copied and adapted from 36d_V2.R/9.Traces
#sampleSPlot <- sample(1:numWeeks, 1)
sampleRPlot <- sample(1:numRegions, 1)
sampleBPlot <- sample(1:numRegionsGroups, 1)
sampleGPlot <- sample(1:numSequences, 1)
sampleGPlot
plotTraces <- rbind(data.table(id = 1, variable = "tau[U]", iteration = 1:length(finalIterations), chain = 1, value = outputTauR),
                    data.table(id = 2, variable = "tau[G]", iteration = 1:length(finalIterations), chain = 1, value = outputTauG),
                    data.table(id = 3, variable = "p", iteration = 1:length(finalIterations), chain = 1, value = outputP),
                    data.table(id = 5, variable = "U[i]", iteration = 1:length(finalIterations), chain = 1, value = outputR[sampleRPlot,]),
                    data.table(id = 4, variable = "G[k]", iteration = 1:length(finalIterations), chain = 1, value = outputG[sampleGPlot,]),
                    #data.table(id = 4, variable = "R[t]", iteration = 1:length(finalIterations), chain = 1, value = outputR[sampleRPlot,]),
                    #data.table(id = 1, variable = "Xjk", melt(epiclustR:::ssapply(modTr, epiclustR:::extract_variable, "X")[17,40,,])),
                    data.table(id = 6, variable = "B[sigma]", iteration = 1:length(finalIterations), chain = 1, value = outputB[sampleBPlot,]))
plotTraces[, variableFc := factor(variable, levels = unique(plotTraces[order(id), variable]))]
# .Plot traces
gtr <- ggplot(plotTraces, aes(x = iteration, y = value, colour = factor(chain))) +
  facet_grid(variableFc ~ ., scales = "free_y", labeller = label_parsed) + geom_path() + theme_laura(16) +
  scale_color_manual(values = c("#003B36", "#E98A15", "#59114D")) + labs(x = "iteration", y = "", colour = "chain") +
  theme(strip.background = element_rect(fill = NA, colour = "gray90"), strip.text.y = element_text(angle = 0)) +
  scale_x_continuous(breaks = c(config$burnIn, seq(1000, config$numIterations, 1000)))
gtr#savePDF(gtr, fileName = "Plot17102019_02_36d_Traces", 8, 7)
#savePDF(gtr, fileName = "Plot05122019_02_38p12c4_0312201901_Tinis_OX_Traces", 8, 7)

# 5. MCMC output: results visualisation ----
# Requires: 0. 1. 2.
# Probability of outbreaks* (... as in 36d.R/10.c. for ChST) ----
# Note that we are plotting only points with prob > 0...!!!
gThr <- ggplot(data.table(x = 1:nrow(outProbs), p = sort(outProbs$proba)), aes(x = x, y = p)) + geom_point() + theme_laura(size = 16) +
  labs(x = TeX("$X_{\\sigma(i),\\xi(k)}$"), y = "probability of outbreak") + scale_x_continuous(breaks = c(1,5000,10000,nrow(outProbs)), labels = rep("", 4)) +
  geom_hline(yintercept = THRESHOLD, linetype = 2, colour = "gray50")
gThr#savePDF(gThr, fileName = "Plot05122019_01_33p12c5_0312201901_Tinis_OX_outP", 6, 4)

#ggplot(outputX[, .N, by = .(row, col)][N>50], aes(x = col, y = row, color = N)) + geom_point() + theme_laura()

# Time proximity within outbreaks (SG Only)* ----
if(0){
  outProbs[proba > 0.2] # region clusterHighId
  regionsForModels[region == 26]
  genomeForModels[clusterHighId == 238]
  casesForModels[dim2Cases %in% regionsForModels[region == 84, col] & dim3Cases %in% genomeForModels[clusterHighId == 395, clusterLowId]]
  casesForModels[dim2Cases %in% regionsForModels[region == 26, col] & dim3Cases %in% genomeForModels[clusterHighId == 238, clusterLowId]]
}

#matToCollapse <- drop(y)
#ttt <- collapseMatrixFn()

dataToPlot <- casesForModels[!is.na(idOutbreak), .(.N, max(numWeek_corrected) - min(numWeek_corrected)), .(idOutbreak, probaOutbreak)]
setnames(dataToPlot, c("N", "V2"), c("sizeOutbreak", "interval"))
#ggplot(dataToPlot, aes(x = probaOutbreak, y = interval, colour = sizeOutbreak == 1)) + geom_point() + theme_laura()
gTOut <- ggplot(dataToPlot[sizeOutbreak > 1], aes(x = probaOutbreak, y = interval)) + geom_point(size = 4) + theme_laura(size = 15) +
  labs(x = TeX("probability of outbreak $\\rho_{\\sigma\\xi}$"), y = TeX("max. temporal distance within block $\\sigma\\xi$ (weeks)"))
gTOut#savePDF(gThr, fileName = "Plot05122019_02_33c8_0312201901_Tinis_OX_outP", 6, 4)
#savePDF(gTOut, fileName = "Plot05122019_05_38p12c4_0312201901_Tinis_OX_Traces", 8, 7)
#savePDF(gTOut, fileName = "Plot19122019_05_38p12c4_1812201901_Local_OX_Traces", 8, 7)
#savePDF(gTOut, fileName = "Plot19022020_04_38p12c4_2601202002_M12_TimeComp", 8, 5.5)

# Plot splitstree TODO ----
distObj <- as.dist(distanceMatrixMCMC)
nnetObj <- neighborNet(distObj)
#save(nnetObj, file = "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_201912/39_04122019_REF21112019_NexusMatrixOX.RData")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager") # like hundred years
if(0) BiocManager::install("treeio")
if(0) BiocManager::install("ggtree")
uploadPackages(c("phangorn", "tidytree"))
uploadPackages(c("ggtree", "treeio")) # treeio useless?

data(yeast)
dm <- dist.ml(yeast)
nnet <- neighborNet(dm) # never do ggplot(nnet)
plot(nnet, "2D")

#nexus2Tree <- treeio::read.mega(nnet)
#nexus2Tree <- ape::read.nexus(nnet)
nexus2df <- tidytree::as_tibble(nnet)
df2tree <- tidytree::as.treedata(nexus2df)
treenet <- as.treedata(nnet)
ggplot(df2tree, aes(x, y)) + geom_tree()

ttt <- ape::read.nexus(file = "/home/laura/Desktop/hhh.nex")
ttt <- treeio::read.tree(file = "/home/laura/Desktop/hhh.nex")

# Plot G using hclust ----
# CHOOSE the clustering used for data in dirInputFile
#load("/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38V2_II_18112019_ClustersKInfo.RData")
load("/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38V2_II_18122019_ClustersKInfo_OX.RData")
# hClustOut, ddata, readme, exploreCutFn, clusterInfoFn, colsDistanceMatrixGroupsNoDups* (forgotten)
cat(readme)

# (very untidy) (groupDupId info is still missing)
# TODO fix/organise
groupSequencesCore
temp$groupToClusterTable[]
genomeForModels
# NOTE that hClustCut has the same ordering as colsDistanceMatrixGroupsNoDups and groupToClusterTable$groupDupId

# Paste risk
genomeInfo <- genomeForModels
range(genomeForModels$col) == c(1, nrow(genomeInfo))
genomeInfo[order(col), riskMean := apply(outputG, 1, mean)]
genomeInfo[order(col), riskSd := apply(outputG, 1, sd)]
genomeInfo[order(col), risk025 := apply(outputG, 1, quantile, 0.025)]
genomeInfo[order(col), risk975 := apply(outputG, 1, quantile, 0.975)]
genomeInfo[order(col), relRisk := apply(exp(outputG), 1, mean)]

groupSequencesCore[, groupDupIdTODO := sapply(groupId, function(x) which(colsDistanceMatrixGroupsNoDups == x)[1])]

dataInfoPostSOA[!is.na(isDuplicate) & LSOA11CD %in% chosenAreas$LSOACDs]
casesForModels

#temp <- clusterInfoFn(heighCutLow = heighCutLow, heighCutHigh = heighCutHigh, hClustOut, groupSequencesCore, distanceMatrixGroups, colsDistanceMatrixGroupsNoDups)
setkey(casesForModels, columnInDataAlleles)
setkey(groupSequencesCore, columnInDataAlleles)
casesForModels[groupSequencesCore, groupDupIdTODO := groupDupIdTODO]
setkey(casesForModels, idInDataInfo)
setkey(dataAnnotatedList[[2]], id)
casesForModels[dataAnnotatedList[[2]], st := clonal_complex..MLST.]

temp <- dcast(casesForModels, formula = clusterLowId ~ st) # expect a warning
meltedST <- data.table(melt(temp, id.vars = "clusterLowId"))[!is.na(clusterLowId)]
sum(meltedST$value) == nrow(casesForModels)

genomeInfo[, orderedParam := rank(relRisk)]
setkey(meltedST, clusterLowId)
setkey(genomeInfo, clusterLowId)
meltedST[genomeInfo, orderedParam := orderedParam]

gG1 <- ggplot(genomeInfo, aes(x = orderedParam, y = riskMean, ymin = risk025, ymax = risk975)) +
  #geom_errorbar(color = "#919994") +
  geom_ribbon(fill = "#CACFD6") +
  geom_point(size = 1) +
  theme_laura(size = sizePlot) + theme(axis.text.x = element_blank()) +
  scale_x_continuous(limits = c(minParam, numDims[3])) +
  labs(x = "", y = TeX("mean and $(p_{0.025}, p_{0.975})$ interval"), subtitle = TeX("Posterior distribution of $G_k$")) # TeX("mean $\\pm$ ")

meltedST[value != 0, valueNA := value]
minParam <- min(meltedST$orderedParam)
orderSTs <- casesForModels[st != "", .N, st][order(N), st]
#ggplot(meltedST[variable != "V1"], aes(x = orderedParam, y = variable, fill = valueNA)) + geom_tile() + theme_laura()
#ggplot(meltedST[variable != "V1"], aes(x = orderedParam, y = variable, size = valueNA)) + geom_point() + theme_laura()
#ggplot(meltedST[variable != "V1" & orderedParam > 1190], aes(x = orderedParam, y = variable, size = valueNA)) + geom_point() + theme_laura()
gG2 <- ggplot(meltedST[variable != "V1"], aes(x = orderedParam, y = factor(variable, levels = orderSTs), size = valueNA)) + geom_point() + theme_laura() +
  labs(x = "", y = "", size = "Num. cases", subtitle = "") + scale_x_continuous(limits = c(minParam, numDims[3])) +
  theme(axis.text.x = element_blank())

#NOdevtools::install_github("wilkelab/cowplot")
#NOlibrary(cowplot) # To align plot areas
devtools::install_github("thomasp85/patchwork")
library(patchwork)
gGT <- gG1 + gG2 + plot_layout(ncol = 1)
#savePDF(gGT, fileName = "Plot05122019_06_38p12c5_0312201901_Tinis_OX_G", 12, 10)
#savePDF(gGT, fileName = "Plot19122019_06_38p12c5_1812201901_Local_OX_G", 12, 10)

plot(hClustOut)
# plot adapted from 39_01.R/II/11.
ggplot(ddata) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + scale_y_reverse(expand = c(0.2, 0)) + theme_laura()

# Time series ----

# Maps* ----
# copied and adapted from 36d.R/10.b.
spatialInfo <- regionsForModels
range(regionsForModels$col) == c(1, nrow(spatialInfo))
spatialInfo[order(col), riskMean := apply(outputR, 1, mean)]
spatialInfo[order(col), riskSd := apply(outputR, 1, sd)]
spatialInfo[order(col), risk025 := apply(outputR, 1, quantile, 0.025)]
spatialInfo[order(col), risk975 := apply(outputR, 1, quantile, 0.975)]
spatialInfo[order(col), relRisk := apply(exp(outputR), 1, mean)]

tempQuants <- quantile(abs(log(spatialInfo$relRisk)), c(1/9,3/9,5/9,7/9,9/9))
# (04.01.2020 correction: better max label and avoid possible NAs by round)
tempMaxCorrection <- ceiling(100*max(spatialInfo$relRisk))/100 # trick not done in Maps to avoid problems with round, and also avoid a large maxima label
spatialInfo[, labels := cut(relRisk,
                            breaks = c(0, sort(round(exp(c(rev(-tempQuants), tempQuants)), digits = 2))[2:(2*length(tempQuants) - 1)], tempMaxCorrection),
                            #breaks = c(0, sort(round(exp(c(rev(-tempQuants), tempQuants)), digits = 2))[2:(2*length(tempQuants))]),
                            include.lowest = F)]

dataToMap <- data.table(fortify(chosenAreas$shpLSOA), key = "id") # or broom::tidy
setkey(spatialInfo, mapRn)
dataToMap[spatialInfo, c("labels", "relRisk") := .(labels, relRisk)]
tempLads <- data.shapeLADS[data.shapeLADS$lad18cd %in% dataLSOA[LSOA11CD %in% chosenAreas$LSOACDs, LAD17CD],]
g4 <- ggplot(dataToMap, aes(x = long, y = lat)) +
  geom_polygon(aes(group = group, fill = labels), alpha = 1, colour = "gray90", size = 0.1) +
  geom_polygon(data = tempLads, aes(group = group), colour = "black", fill = NA, size = 1) +
  annotate(geom = "rect", xmin = mapLims$xE[1], xmax = mapLims$xE[2], ymin = mapLims$yE[1], ymax = mapLims$yE[2],
           fill = "transparent", color = "black", size = 1.5) +
  # CHOOSE:
  annotate(geom = "text", x = c(-1.6, -1.5), y = c(52.15, 51.5), label = c("NORTHAMPTONSHIRE", "OXFORDSHIRE")) + # OX
  #annotate(geom = "text", x = c(-2.45, -1.58), y = c(55.75, 54.9), label = c("NORTHUMBERLAND", "NEWCASTLE,\nNORTH TYNESIDE")) + # TW
  scale_fill_manual(guide = "legend", values = rev(brewer_pal(type = "div", palette = "RdBu")(9))) +
  coord_quickmap(xlim = mapLims$x, ylim = mapLims$y) + theme_void(base_size = 16) +
  labs(x = "", y = "", fill = "Relative risk of\nsporadic cases") +
  guides(fill = guide_legend(reverse = TRUE))
g4#savePDF(g4, fileName = "Plot21082019_04_33c8_2008201901_Tinis_OX_map", 8, 8)
#savePDF(g4, fileName = "Plot05122019_03_38p12c5_0312201901_Tinis_OX_map", 8, 8)
# Enlargement
tempLads2 <- tempLads[!tempLads$lad18nm %in% ifelse(typ == "OX", "Oxford", "Northumberland")]
g5 <- ggplot(dataToMap, aes(x = long, y = lat)) +
  geom_polygon(aes(group = group, fill = labels), alpha = 1, colour = "gray90", size = 0.1) +
  geom_polygon(data = tempLads2, aes(group = group), colour = "black", fill = NA, size = 1) +
  #geom_text(data = chosenSites, aes(x = xCoord, y = yCoord, label = numCases), size = 2) +
  scale_fill_manual(guide = "legend", values = rev(brewer_pal(type = "div", palette = "RdBu")(9))) + # PuBuGn
  coord_quickmap(xlim = mapLims$xE, ylim = mapLims$yE) + theme_void(base_size = 16) +
  guides(fill = FALSE) + theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(subtitle = ifelse(typ == "OX", "Oxford", "Newcastle upon Tyne - North Tyneside"))
g5#savePDF(g5, fileName = "Plot21082019_05_33c8_2008201901_Tinis_OX_map2", 4, 4)
#savePDF(g5, fileName = "Plot05122019_04_38p12c5_0312201901_Tinis_OX_map2", 4, 4)

# .Spatial overview ----
# (05.02.2020)
spatialInfo[order(-relRisk)]
spatialInfo[, .N, .(RUC11CD,RUC11)]
spatialInfo[, .(min(relRisk), mean(relRisk), max(relRisk), mean(riskMean)), .(RUC11CD,RUC11)]
mean(spatialInfo$relRisk)
mean(exp(outputR))
plot(exp(rowMeans(outputR)), rowMeans(exp(outputR)))
sum(exp(rowMeans(outputR)) == rowMeans(exp(outputR)))
mean(exp(rowMeans(outputR)))
mean(rowMeans(exp(outputR)))
plot(apply(outputR, 1, e1071::skewness))
plot(apply(outputR, 1, e1071::skewness), rowMeans(exp(outputR)) - exp(rowMeans(outputR)))
plot(apply(outputR, 1, e1071::skewness), apply(outputR, 1, sd))
plot(spatialInfo$estimPop2015, spatialInfo$relRisk)

# Plot tree - sandbox ----
# Conclusion: Tried mst but unclear how to make a nice plot out of it :/

coo <- matrix(runif(10), ncol = 2)
ttt <- ape::mst(dist(coo))
spdep::plot.mst(x = ttt, coords = coo)

library(igraph)
g <- erdos.renyi.game(10, 3/10)
mst <- minimum.spanning.tree(g)
par(mfrow=c(1,2), mar=c(0,1,0.75,0)) # sub-plots and margins
plot(g , main="Graph")
plot(mst, main = "MST")
igraph::plot.igraph(mst)
igraph::plot.igraph(ttt)

library(fossil)
data(fdata.mat)
fdata.dist<-dino.dist(fdata.mat)
dino.mst<-dino.mst(fdata.dist)

# 17.12
ig <- graph_from_adjacency_matrix(dino.mst, mode = "undirected")
mst <- igraph::minimum.spanning.tree(ig)
igraph::plot.igraph(mst)
ape::as.igraph.phylo(mst, FALSE)
library(ggtree)
ggplot(mst)

distM <- as.dist(distanceMatrixGroups)
dino.mst <- dino.mst(distM)
save(dino.mst, file = "hola.RData")
load("hola.RData")
ig <- graph_from_adjacency_matrix(dino.mst, mode = "undirected")
mst <- igraph::minimum.spanning.tree(ig)
igraph::plot.igraph(mst)
tkplot(mst)

#18.12 :)
# Este si fue!
# https://yulab-smu.github.io/treedata-book/chapter9.html
dino.nj <- ape::nj(distM)
#ggtree(dino.nj) + theme_tree2() + geom_tiplab(align = TRUE, linesize=.5)
p <- ggtree(dino.nj) + theme_tree()
d1 <- data.frame(id = dino.nj$tip.label, val = rnorm(3695, sd = 3))
facet_plot(p, panel = "dot", data = d1, geom = geom_point, aes(x = val), color = 'red3')

# Plot tree ----
distM <- as.dist(distanceMatrixMCMC)
njTree <- ape::nj(distM)
p <- ggtree(njTree) + theme_tree2()
d1 <- data.frame(id = njTree$tip.label, val = rowMeans(outputG))
gt <- facet_plot(p, panel = "G_k", data = d1, geom = geom_point, aes(x = val), color = 'red3', size = 4) +
  theme(strip.background = element_rect(fill = NA, colour = "gray90"), strip.text = element_text(size = 22),
        axis.text = element_text(size = 20))
gt#savePDF(gt, fileName = "Plot19122019_07_38p12c5_1812201901_Local_OX_G_Tree", 12, 10)

# im guessing the order of matrix is same as outputG:
#outputG
#genomeForModels
if(0){
  difOutput <- outer(rowMeans(outputG), rowMeans(outputG), FUN = function(x,y) abs(x - y))
  distanceMatrixMCMC
  dataToPlot <- data.table(cbind(melt(difOutput), melt(distanceMatrixMCMC)))
  names(dataToPlot) <- c("x", "y", "difG", "xx", "yy", "dist")
  dataToPlot[, sum(x != xx) + sum(y != yy)] == 0
  ggplot(dataToPlot[x < y], aes(x = dist, y = difG)) + geom_point() + theme_laura()
}

# Covariance function**41.R/II ----
dataDistKernel <- data.table(dis = 0:50)
if(config$maternParameter == 1/2){
  dataDistKernel[, mean := sapply(dis, function(d) mean(exp(-d/outputL)/outputTauG))] # forgot how to do it in once
  dataDistKernel[, pDn := sapply(dis, function(d) quantile(exp(-d/outputL)/outputTauG, 0.025))]
  dataDistKernel[, pUp := sapply(dis, function(d) quantile(exp(-d/outputL)/outputTauG, 0.975))]
}else if(config$maternParameter == 0){
  dataDistKernel[, mean := sapply(dis, function(d) mean(exp(-d^2/2*outputL^2)/outputTauG))]
  dataDistKernel[, pDn := sapply(dis, function(d) quantile(exp(-d^2/2*outputL^2)/outputTauG, 0.025))]
  dataDistKernel[, pUp := sapply(dis, function(d) quantile(exp(-d^2/2*outputL^2)/outputTauG, 0.975))]
}else if(config$maternParameter == 3/2){
  dataDistKernel[, mean := sapply(dis, function(d) mean((1 + (sqrt(3)*d/outputL))*exp(-sqrt(3)*d/outputL)/outputTauG))]
  dataDistKernel[, pDn := sapply(dis, function(d) quantile((1 + (sqrt(3)*d/outputL))*exp(-sqrt(3)*d/outputL)/outputTauG, 0.025))]
  dataDistKernel[, pUp := sapply(dis, function(d) quantile((1 + (sqrt(3)*d/outputL))*exp(-sqrt(3)*d/outputL)/outputTauG, 0.975))]
}
ggplot(dataDistKernel, aes(x = dis, y = mean)) + geom_ribbon(aes(ymin = pDn, ymax = pUp), fill = "grey70") + geom_line() + theme_laura() +
  labs(x = "distance", y = "covariance")

#install.packages("ggExtra")
gg <- ggplot(data.table(outputTauG, outputL), aes(x = outputTauG, y = outputL)) + geom_point() + theme_laura() + # geom_path()
  labs(x = TeX("$\\tau_G$"), y = TeX("$\\rho$")) # trace and density plot
ggMarginal(gg, type = 'density', colour = '#41959E', fill = '#41959E') # 779997

# Expected number of cases per genome ----
# Long... that's why I wanted it inside the iterations... TODO
# TODO test
# TODO do for dim = 1
subsetIterations <- finalIterations #sample(finalIterations, 200)
subsetBlockG <- c(238, 395, 522) # dim3Cases from outbreak table... More than that would be crazy # temp_MAT32
subsetG <- sort(unique(genomeForModels[clusterHighId %in% subsetBlockG, clusterLowId]))
matrixPop <- aperm(array(pop, dim = c(numRegions, numWeeks, numSequences)), perm = c(2,1,3))
sCasesU <- matrix(0, nrow = numDims[2], ncol = length(subsetIterations))
eCasesU <- matrix(0, nrow = numDims[2], ncol = length(subsetIterations))
sCasesG <- matrix(0, nrow = numDims[3], ncol = length(subsetIterations))
eCasesG <- matrix(0, nrow = numDims[3], ncol = length(subsetIterations))
sCasesAll <- array(0, dim = c(numDims[2], length(subsetG), length(subsetIterations)))
eCasesAll <- array(0, dim = c(numDims[2], length(subsetG), length(subsetIterations)))
for(it in 1:length(subsetIterations)){ # ~ 30s
  matrixG <- aperm(array(outputG[,it], dim = c(numSequences, numWeeks, numRegions)), perm = c(2,3,1))
  matrixS <- array(outputS[,it], dim = c(numWeeks, numRegions, numSequences))
  matrixR <- aperm(array(outputR[,it], dim = c(numRegions, numWeeks, numSequences)), perm = c(2,1,3))
  if(!1 %in% dimToInclude){
    ttt <- cbind(1, storage$parameters$X[it][[1]])
    ttt2 <- array(0, dim = numBlockDims)
    ttt2[ttt] <- 1
    matrixX <- array(ttt2[allToGroups], dim = c(numWeeks, numRegions, numSequences))
  }else{
    # TODO
  }
  if(dimBeta == 2){
    matrixB <- aperm(array(outputB[jToGroups, it], dim = c(numRegions, numWeeks, numSequences)), perm = c(2,1,3))
  }else if(dimBeta == 3){
    matrixB <- aperm(array(outputB[kToGroups, it], dim = c(numSequences, numWeeks, numRegions)), perm = c(2,3,1))
  }
  matrixXBexp <- exp(matrixX*matrixB)
  sCasesU[,it] <- apply(matrixPop*exp(outputA[it] + matrixS + matrixR + matrixG), 2, sum) # 3!!!!!!!!!!!
  eCasesU[,it] <- apply(matrixPop*exp(outputA[it] + matrixS + matrixR + matrixG)*matrixXBexp, 2, sum)
  sCasesG[,it] <- apply(matrixPop*exp(outputA[it] + matrixS + matrixR + matrixG), 3, sum)
  eCasesG[,it] <- apply(matrixPop*exp(outputA[it] + matrixS + matrixR + matrixG)*matrixXBexp, 3, sum)
  sCasesAll[,,it] <- (matrixPop*exp(outputA[it] + matrixS + matrixR + matrixG))[,,subsetG]
  eCasesAll[,,it] <- (matrixPop*exp(outputA[it] + matrixS + matrixR + matrixG)*matrixXBexp)[,,subsetG]
} # Expect warnings
dataToPlot <- rbind(data.table(k = 1:numDims[dimForPlot], cases = sapply(1:numDims[dimForPlot], function(kk) median(sCasesG[kk,])), type = "sporadic"),
                    data.table(k = 1:numDims[dimForPlot], cases = sapply(1:numDims[dimForPlot], function(kk) median(eCasesG[kk,])), type = "total"),
                    data.table(k = 1:numDims[dimForPlot], cases = apply(y, dimForPlot, sum), type = "observed"))
tempOrderPlot <- dataToPlot[type == "sporadic", frank(cases, ties.method = "dense")]
dataToPlot[, orderK := sapply(k, function(x) tempOrderPlot[x])]
ggplot(dataToPlot[order(orderK)], aes(x = orderK, y = cases, colour = factor(type, levels = c("total", "sporadic", "observed")))) + geom_path() + theme_laura() +
  labs(colour = "type") # k 418

if(0){
  plot(sCasesG[418,])
  sum(y[,,418])
  which.max(eCasesG[418,]) # it 290
  
  # for it 290
  (matrixPop*exp(outputA[it] + matrixS + matrixR + matrixG)*matrixXBexp)[,,418]
  which.max((matrixPop*exp(outputA[it] + matrixS + matrixR + matrixG)*matrixXBexp)[,,418]) # region 295
  (matrixPop*exp(outputA[it] + matrixS + matrixR + matrixG)*matrixXBexp)[,295,418]
  matrixXBexp[,295,418]
  matrixX[,295,418]
  matrixB[,295,418]
  # one B is extremely large!
  
  # THEREFORE: change mean to meadian in loop above
}

# Outbreak table* ----
# REQUIRES 5./Expected number of cases per genome with dimForPlot = 3 (just before)

THRESHOLD <- 0.15
outProbs[proba > THRESHOLD]
casesForModels[probaOutbreak > THRESHOLD]

setkey(casesForModels, idInDataInfo)
setkey(dataInfoPostSOA, id)
casesForModels[dataInfoPostSOA, received_date_format := received_date_format]
setkey(casesForModels, LSOA11CD)
setkey(dataLSOA, LSOA11CD)
casesForModels[dataLSOA, LSOA11NM := LSOA11NM]
setkey(casesForModels, idInDataInfo)
setkey(dataAnnotatedList[[2]], id)
casesForModels[dataAnnotatedList[[2]], ccStr := clonal_complex..MLST.]

# -- WRITE down dim3Cases of chosen outbreaks and MODIFY subsetBlockG (opt.)
# sum of the mean is the mean of the sum
eCasesAll_col <- apply(eCasesAll, 1:2, mean)
sCasesAll_col <- apply(sCasesAll, 1:2, mean)
#e.g. eCasesAll_col[regionsForModels[region == 26, col], which(subsetG %in% genomeForModels[clusterHighId == 238, col])]
casesForModels[probaOutbreak > THRESHOLD,
               ecasesBlock := mapply(function(x, y) sum(eCasesAll_col[regionsForModels[region == x, col], which(subsetG %in% genomeForModels[clusterHighId == y, col])]),
                                region, clusterHighId)]
casesForModels[probaOutbreak > THRESHOLD,
               scasesBlock := mapply(function(x, y) sum(sCasesAll_col[regionsForModels[region == x, col], which(subsetG %in% genomeForModels[clusterHighId == y, col])]),
                                region, clusterHighId)]
#casesForModels[, ecasesBlock := sum(ecases), idOutbreak]
#casesForModels[, scasesBlock := sum(scases), idOutbreak]
# --

casesForModels[probaOutbreak > THRESHOLD][order(received_date_format)]
# print
casesForModels[probaOutbreak > THRESHOLD][order(received_date_format)][,.(sizeOutbreak, probaOutbreak, received_date_format, ccStr, LSOA11NM)]
casesForModels[probaOutbreak > THRESHOLD][order(received_date_format)][,.(sizeOutbreak, probaOutbreak, received_date_format, ccStr, LSOA11NM, scasesBlock, ecasesBlock)]

# Distance between clusters histogram* ----
# Only requires the info in dirInputFile
dataToPlot <- data.table(melt(distanceMatrixMCMC))[X1 < X2]
plot(dataToPlot[order(value), value])
ghist <- ggplot(dataToPlot[value < 1000]) + geom_histogram(aes(x = value, y = ..density..), breaks = seq(10,1287,10)) +
  scale_x_continuous(limits = c(10,1000), breaks = c(10,seq(250,1000,250))) + theme_laura(15) +
  theme(panel.grid.major = element_line(colour = "gray", size = 0.3), panel.grid.minor = element_line(colour = "gray", size = 0.1)) + # as in theme_laura2
  labs(x = "distance")
ghist#savePDF(ghist, fileName = "Plot30012020_01_38p12c5_18122019Input_HistClust", 8, 4) # warning

# BORRAR!! ----

plot(colSums(drop(y)), rowMeans(outputG))

plot(table(distanceMatrixMCMC[distanceMatrixMCMC < 50 & distanceMatrixMCMC > 0]))

plot(table(colSums(drop(y))))

plot(1:50, exp(-seq(1,50,1)/mean(outputL))/mean(outputTauG))
plot(1:50, exp(-seq(1,50,1)/mean(outputL))/100)

xInput <- seq(-50,50,1)
xCov <- exp(-abs(xInput)/50)/10 # 10:1 10:10 50:1  50:10
#xCov <- exp(-abs(xInput)/mean(outputL))/mean(outputTauG)
plot(xInput, rnorm(length(xInput), mean = 0, sd = xCov))


# # 3/Prior-posterior hyperparameters/a
gPos2 <- ggplot(dataToPlotPos, aes(x = x, fill = label, color = label)) +
  geom_density(data = dataToPlotPos[label == "posterior"], alpha = 0.3) +
  geom_line(aes(y = y)) + theme_laura() + labs(x = "", y = "density", color = "", fill = "")
gPos2







# TODO expected number of cases
#matrixX <- replicate(numSequences, matrix(matX[allToGroups], nrow = numWeeks, ncol = numRegions, byrow = FALSE))

ggplot(dataToPlotGSB[label == "G"], aes(x = jump, y = mean, ymin = mean - sd, ymax = mean + sd)) + geom_errorbar(color = "#919994") +
  geom_point(size = 1)
ggplot(dataToPlotGSB[label == "G"], aes(x = rank(real), y = mean, ymin = mean - sd, ymax = mean + sd)) + geom_errorbar(color = "#919994") +
  geom_point(size = 1) + geom_point(aes(y = real), color = "red", size = 1)













# Output S and y. Inla output.
if(0){
  dataToPlot2 <- melt(dataToPlot, id.vars = "parameter", measure.vars = c("mean","y"))
  ggplot(dataToPlot2, aes(x = parameter, y = value)) + facet_grid(variable ~ ., scales="free") + geom_point() + theme_laura()
  
  load("../PhD_Year2/06_MixedModels/RCode_201903/15_21032019_OutputInlaOnlyS.RData") # outputInla, dataToPlotInla # Created in 15.R/5.
  ggplot(dataToPlotInla, aes(x = parameterName, y = mean, ymin = mean - sd, ymax = mean + sd)) + geom_errorbar() +
    geom_point() + theme_laura() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))
}

# Output X
if(0){
  ggplot(outputX[, .N, by = .(row, col)], aes(x = col, y = row, fill = N > 500)) + geom_tile() + theme_laura()
  ggplot(outputX[, .N, by = .(row, col)][N>50], aes(x = col, y = row, color = N)) + geom_point() + theme_laura()
}


# Compare to annotated STs
outputG
groupSequencesCore
dataInfo
dataAnnotated <- dataAnnotatedList[[2]]
# add column in outputG to dataAnnotated
setkey(dataAnnotated, id)
setkey(dataInfo, id)
dataAnnotated[dataInfo, columnInDataAlleles := columnInDataAlleles]
setkey(dataAnnotated, columnInDataAlleles)
setkey(groupSequencesCore, columnInDataAlleles)
dataAnnotated[groupSequencesCore, columnInOutput := groupId]
dataAnnotated[, meanG := rowMeans(outputG)[columnInOutput]]
#(why did we group them??)
dataToPlotST <- dataAnnotated[, list(unique(clonal_complex..MLST.), min(meanG)), by = .(columnInOutput)] # TODO improoove
setnames(dataToPlotST, c("V1","V2"), c("ST","mean"))
temp <- dataToPlotST[, mean(mean), by = .(ST)][order(V1), ST]
gST <- ggplot(dataToPlotST, aes(x = factor(ST, levels = temp), y = mean)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) + theme_laura(size = sizePlot) +
  labs(x = "ST type", y = TeX("$G_k$")) + theme(panel.grid.major = element_line(colour = "gray", size = 0.3))
gST#savePDF(gST, "Plot16042019_MCMCST_Real", 8, 6)

# X vs. simulation
xAllValues <- CJ(i = 1:numWeeksGroups, j = 1:numRegionsGroups)
temp <- outputX[, .N, by = .(week, region)]
setkeyv(xAllValues, c("i", "j"))
setkeyv(temp, c("week", "region"))
xAllValues[temp, c("N", "numItOutbr") := .(N, 100*N/(it - 1))]
xAllValues[is.na(N), N := 0]
xAllValues[is.na(numItOutbr), numItOutbr := 0]

simulatedX <- data.table(melt(simulatedParam$X))
setnames(simulatedX, c("X1","X2","value"), c("week","region","isOutbr"))
setkeyv(xAllValues, c("i", "j"))
setkeyv(simulatedX, c("week", "region"))
xAllValues[simulatedX, isOutbr := isOutbr]
sum(simulatedParam$X) == sum(xAllValues$isOutbr)

plot(xAllValues$N, xAllValues$isOutbr)
xAllValues[isOutbr == 1]

