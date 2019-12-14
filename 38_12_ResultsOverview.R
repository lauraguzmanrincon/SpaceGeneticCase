# All sections require 0. 1. 2.
# Never re-run 0. (may overwrite some variables)

if(0){
  # 0. Load workspace ----
  source("/home/laura/Dropbox/Laura/PhD_Year2/04_AboutApproaches/RCode/2018_04/00_Functions.R")
  #source("/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201908/34_HeaderFor201908.R")
  setwd("/home/laura/Dropbox/Laura/PhD_Year3/")
  allPlotsPath <- "07_MixedModelsP2"
  
  uploadPackages("latex2exp") # for TeX
  source(paste0("07_MixedModelsP2/RCode_201911/38_Files/", "38_02_AuxiliarFn.R"))
  load("07_MixedModelsP2/RCode_201909/33a_03092019_DataOXTW_Perfect_Depr.RData") # dataLSOA, dataMSOA, dataInfoPostSOA, chosenAreasOX, chosenAreasTW
  
  # For map and colours
  #uploadPackages(c("rgdal", "scales"))
  tempDir <- "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/Areafiles_201907/"
  tempName <- paste0(tempDir, "15_Local_Authority_Districts_December_2018_Boundaries_GB_BFC")
  tempMap <- readOGR(dsn = tempName, layer = "Local_Authority_Districts_December_2018_Boundaries_GB_BFC")
  data.shapeLADS <- spTransform(tempMap, CRS("+proj=longlat +datum=WGS84"))
  
  # REMOVE following after including header before?
  # For groups and clusters
  load("/home/laura/Dropbox/Laura/PhD_Year2/06_MixedModels/RCode_201903/15_Input_MCMCCorrected04032019.RData")
  # Content: numSequences, numWeeks, y, distanceMatrixGroups, groupSequencesCore
  # Recall: column order in distanceMatrixGroups is given by groupSequencesCore$groupId
  rm(numWeeks, y, numSequences)
  suppressWarnings(groupSequencesCore[, numWeek := NULL])
  suppressWarnings(groupSequencesCore[, weekId := NULL])
  
  # General data
  load("/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/Data/FS101013_Datasets_V3_26R.RData") # ...
}

# 1. Upload MCMC output file + others ----
#dirFiles <- "FilesFromLocal/38_Files/"
#dirFiles <- "07_MixedModelsP2/RCode_201911/38_Files/"
#dirFiles <- "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_201911/"
nameFiles <- "21112019"
nameFilesOut <- "0312201901"
simulatedDataRes <- 0
typ <- "OX"
dirInputFiles <- "07_MixedModelsP2/RCode_201911/38_Files/"
dirOutputFiles <- "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_201911/"
dirInputFile <- paste0(dirInputFiles, "38_01_", nameFiles, "_MCMCInputTinis.RData")
load(dirInputFile) # ...
if(simulatedDataRes == 0){
  dirOutputFile <- paste0(dirOutputFiles, "38_00_", nameFilesOut, "_MCMCOutputTinis.RData")
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

# Other updates
if(typ == "OX"){
  chosenAreas <- chosenAreasOX
  mapLims <- list(x = c(-1.8,-0.7), y = c(51.44, 52.25), xE = c(-1.3,-1.165), yE = c(51.71, 51.8))
}else{
  chosenAreas <- chosenAreasTW
  mapLims <- list(x = c(-2.66,-1.36), y = c(54.8, 55.8), xE =  c(-1.41,-1.78), yE = c(54.95, 55.1))
}

# 2. Build data from MCMC output for visualisations ----
#finalIterations <- 1001:5000
finalIterations <- (config$burnIn + 1):config$numIterations
#finalIterations <- 1:500
#finalIterations <- 1:config$numIterations
#finalIterations <- 1:(it - 1)
sizePlot <- 15
#simulatedDataRes <- 0 # BORRAR

# Check output
#gIndices <- union((freqGUpdate != 1)*which((finalIterations)%%freqGUpdate == 1), (freqGUpdate == 1)*(finalIterations))
outputA <- storage$parameters$a[finalIterations]
outputG <- storage$parameters$G[,finalIterations]
outputS <- array(storage$parameters$S[,finalIterations], dim = c(numPeriods, length(finalIterations))) # just in case is of size 1
outputR <- array(storage$parameters$R[,finalIterations], dim = c(numRegions, length(finalIterations)))
outputB <- storage$parameters$B[,finalIterations]
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

# 3. MCMC output: quick chain visualisation ----
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
  geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $U_j$")) +
  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$S[sampleSPlot], NaN))
g10 <- ggplot(data.table(it = finalIterations, value = outputR[sampleRPlot,]), aes(x = it, y = value)) +
  geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $S_j$")) +
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
  g8
  #g3
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

# 4. MCMC output: formal chain visualisation (table, plots) ----
# Plot traces* ----
# copied and adapted from 36d_V2.R/9.Traces
#sampleSPlot <- sample(1:numWeeks, 1)
sampleRPlot <- sample(1:numRegions, 1)
sampleBPlot <- sample(1:numRegionsGroups, 1)
sampleGPlot <- sample(1:numSequences, 1)
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

# Probability of outbreaks* (... as in 36d.R/10.c. for ChST) ----
# Note that we are plotting only points with prob > 0...!!!
THRESHOLD <- 0.05
outProbs <- outputX[, .N/length(finalIterations), .(row, col)]
outProbs[, id := .I]
#plot(sort(temp$V1))
gThr <- ggplot(data.table(x = 1:nrow(outProbs), p = sort(outProbs$V1)), aes(x = x, y = p)) + geom_point() + theme_laura(size = 16) +
  labs(x = TeX("$X_{\\sigma(i),\\xi(k)}$"), y = "probability of outbreak") + scale_x_continuous(breaks = c(1,5000,10000,nrow(outProbs)), labels = rep("", 4)) +
  geom_hline(yintercept = THRESHOLD, linetype = 2, colour = "gray50")
gThr#savePDF(gThr, fileName = "Plot05122019_01_33p12c5_0312201901_Tinis_OX_outP", 6, 4)

#ggplot(outputX[, .N, by = .(row, col)][N>50], aes(x = col, y = row, color = N)) + geom_point() + theme_laura()

# Time proximity within outbreaks (SG Only) ----
outProbs[V1 > 0.2] # region clusterHighId
regionsForModels[region == 47]
genomeForModels[clusterHighId == 271]
casesForModels[dim2Cases %in% regionsForModels[region == 47, col] & dim3Cases %in% genomeForModels[clusterHighId == 271, clusterLowId]]
casesForModels[dim2Cases %in% regionsForModels[region == 26, col] & dim3Cases %in% genomeForModels[clusterHighId == 234, clusterLowId]]

#matToCollapse <- drop(y)
#ttt <- collapseMatrixFn()

setkey(casesForModels, dim2Cases)
setkey(regionsForModels, col)
casesForModels[regionsForModels, region := region]

setkeyv(casesForModels, c("region", "clusterHighId"))
setkeyv(outProbs, c("row", "col"))
casesForModels[outProbs, c("idOutbreak", "probaOutbreak") := .(id, V1)]

dataToPlot <- casesForModels[!is.na(idOutbreak), .(.N, max(numWeek_corrected) - min(numWeek_corrected)), .(idOutbreak, probaOutbreak)]
setnames(dataToPlot, c("N", "V2"), c("sizeOutbreak", "interval"))
#ggplot(dataToPlot, aes(x = probaOutbreak, y = interval, colour = sizeOutbreak == 1)) + geom_point() + theme_laura()
gTOut <- ggplot(dataToPlot[sizeOutbreak > 1], aes(x = probaOutbreak, y = interval)) + geom_point(size = 4) + theme_laura(size = 16) +
  labs(x = TeX("probability of outbreak $\\rho_{\\sigma,\\xi}$"), y = TeX("max. temporal distance within block $\\sigma\\xi$ (weeks)"))
gTOut#savePDF(gThr, fileName = "Plot05122019_02_33c8_0312201901_Tinis_OX_outP", 6, 4)
#savePDF(gTOut, fileName = "Plot05122019_05_38p12c4_0312201901_Tinis_OX_Traces", 8, 7)

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

# Plot G using hclust* ----
# CHOOSE the clustering used for data in dirInputFile
load("/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38V2_II_18112019_ClustersKInfo.RData")
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
spatialInfo[, labels := cut(relRisk,
                            breaks = c(0, sort(round(exp(c(rev(-tempQuants), tempQuants)), digits = 2))[2:(2*length(tempQuants))]),
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


# BORRAR!! ----

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

