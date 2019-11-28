
if(0){
  # Load workspace
  source("/home/laura/Dropbox/Laura/PhD_Year2/04_AboutApproaches/RCode/2018_04/00_Functions.R")
  #source("/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201908/34_HeaderFor201908.R")
  setwd("/home/laura/Dropbox/Laura/PhD_Year3/")
  allPlotsPath <- "07_MixedModelsP2"
}

# Upload file
#dirFiles <- "FilesFromLocal/38_Files/"
#dirFiles <- "07_MixedModelsP2/RCode_201911/38_Files/"
#dirFiles <- "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_201911/"
nameFiles <- "21112019"
simulatedDataRes <- 1
dirInputFiles <- "07_MixedModelsP2/RCode_201911/38_Files/"
dirOutputFiles <- "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_201911/"
dirInputFile <- paste0(dirInputFiles, "38_01_", nameFiles, "_MCMCInputTinis.RData")
load(dirInputFile) # ...
if(simulatedDataRes == 0){
  dirOutputFile <- paste0(dirOutputFiles, "38_00_", nameFiles, "_MCMCOutputTinis.RData")
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
}

# Take values----
finalIterations <- 1001:5000
finalIterations <- (config$burnIn + 1):config$numIterations
finalIterations <- 1:300
sizePlot <- 15

# Check output
#gIndices <- union((freqGUpdate != 1)*which((finalIterations)%%freqGUpdate == 1), (freqGUpdate == 1)*(finalIterations))
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
setnames(outputX, c("x"), c("it"))
#setnames(outputX, c("x","row","col"), c("it","week","region"))
nrow(outputX)/(prod(numBlockDims)*length(finalIterations))

# Check acceptance rate
# Recall: it's 1 for X, tauS, tauG, p
data.table(names = names(storage$accept),
           do.call(rbind, lapply(1:length(storage$accept), function(x) summary(storage$accept[[x]]/(storage$accept[[x]] + storage$reject[[x]])))))

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
sampleBPlot <- sample(1:numRegionsGroups, 1)
sampleGPlot <- sample(1:numSequences, 1)
rm(g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)
g0 <- ggplot(data.table(it = finalIterations, value = outputA),
             aes(x = it, y = value)) + geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $\\alpha$")) +
  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$a, NaN))
g1 <- ggplot(data.table(it = finalIterations, value = outputTauG),
             aes(x = it, y = value)) + geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $\\tau_G$"))
g2 <- ggplot(data.table(it = finalIterations, value = outputTauS),
             aes(x = it, y = value)) + geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $\\tau_S$"))
g9 <- ggplot(data.table(it = finalIterations, value = outputTauR),
             aes(x = it, y = value)) + geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $\\tau_R$"))
g3 <- ggplot(data.table(it = finalIterations, value = outputP),
             aes(x = it, y = value)) + geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $p$")) +
  #geom_hline(yintercept = ifelse(simulatedDataRes, sum(simulatedParam$X)/prod(dim(simulatedParam$X)), NaN))
  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$p, NaN))
g4 <- ggplot(data.table(it = finalIterations, value = outputS[sampleSPlot,]), aes(x = it, y = value)) +
  geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $S_j$")) +
  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$S[sampleSPlot], NaN))
g10 <- ggplot(data.table(it = finalIterations, value = outputR[sampleRPlot,]), aes(x = it, y = value)) +
  geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $R_j$")) +
  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$R[sampleRPlot], NaN))
g5 <- ggplot(data.table(it = finalIterations, value = outputB[sampleBPlot,]), aes(x = it, y = value)) +
  geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $B_k$")) +
  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$B[sampleBPlot], NaN))
#geom_hline(yintercept = inSimulatedB[1,sampleBPlot])
g6 <- ggplot(data.table(it = finalIterations, value = outputL),
             aes(x = it, y = value)) + geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $\\rho$"))
g7 <- ggplot(data.table(it = finalIterations, value = outputG[sampleGPlot,]), aes(x = it, y = value)) +
  geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $G_k$")) +
  geom_hline(yintercept = ifelse(simulatedDataRes, simulatedParam$G[sampleGPlot], NaN))
g8 <- ggplot(outputX[, .N/prod(numBlockDims), by = .(it)], aes(x = it, y = V1)) +
  geom_point() + theme_laura(size = sizePlot) + labs(x = "iteration", y = TeX("trace of $X_{jk}$"))
#geom_hline(yintercept = sum(inSimulatedX)/(numWeeks*numSequences))
if(2 %in% dimToInclude & 3 %in% dimToInclude){
  #multiplot(g4,g10,g7,g5,g0, g2,g9,g1,g6,g3, cols = 2) # Sj Rj Gk Bk a // tauS tauR tauG ro p
  multiplot(g10,g7,g5,g0, g9,g1,g6,g3, cols = 2) # Sj Rj Gk Bk a // tauS tauR tauG ro p
}
#savePDF(multiplot(g10,g7,g5,g0, g9,g1,g6,g3, cols = 2), "Plot25112019_01_MCMCoutput_Real", 8, 10)#12)
#savePDF(multiplot(g4,g7,g5,g0,g2,g1,g6,g3, cols = 2), "Plot25112019_01_MCMCoutput_Real", 8, 10)#12)

# Output G S B
#outputToPlot <- outputB
dataToPlotGSB <- rbind(data.table(parameter = 1:nrow(outputG), mean = rowMeans(outputG), sd = apply(outputG, 1, sd),
                                  q025 = apply(outputG, 1, quantile, 0.025),
                                  q975 = apply(outputG, 1, quantile, 0.975),
                                  y = apply(y, 3, sum),
                                  meanExp = apply(outputG, 1, function(x) mean(exp(x))),
                                  label = "G", minQ = NA, maxQ = NA, maxI = NA, meanQ = NA,
                                  real = simulatedParam$G,
                                  jump = out$sigmaJumps$G),
                       data.table(parameter = 1:nrow(outputB), mean = rowMeans(outputB), sd = apply(outputB, 1, sd),
                                  q025 = apply(outputB, 1, quantile, 0.025),
                                  q975 = apply(outputB, 1, quantile, 0.975),
                                  y = 0, #apply(y, 2, sum),
                                  meanExp = apply(outputB, 1, function(x) mean(exp(x))),
                                  label = "B", minQ = qgamma(0.025, shape = constants$aB, rate = constants$bB),
                                  maxQ = qgamma(0.975, shape = constants$aB, rate = constants$bB),
                                  maxI = numSequences, meanQ = qgamma(0.5, shape = constants$aB, rate = constants$bB),
                                  real = simulatedParam$B,
                                  jump = out$sigmaJumps$B))
if(1 %in% dimToInclude){
  dataToPlotGSB <- rbind(dataToPlotGSB, data.table(parameter = 1:nrow(outputS), mean = rowMeans(outputS), sd = apply(outputS, 1, sd),
                                                   q025 = apply(outputS, 1, quantile, 0.025),
                                                   q975 = apply(outputS, 1, quantile, 0.975),
                                                   y = apply(y, 1, sum),
                                                   meanExp = apply(outputS, 1, function(x) mean(exp(x))),
                                                   label = "S", minQ = NA, maxQ = NA, maxI = NA, meanQ = NA,
                                                   real = simulatedParam$S,
                                                   jump = out$sigmaJumps$S))
}
if(2 %in% dimToInclude){
  dataToPlotGSB <- rbind(dataToPlotGSB, data.table(parameter = 1:nrow(outputR), mean = rowMeans(outputR), sd = apply(outputR, 1, sd),
                                                   q025 = apply(outputR, 1, quantile, 0.025),
                                                   q975 = apply(outputR, 1, quantile, 0.975),
                                                   y = apply(y, 2, sum),
                                                   meanExp = apply(outputR, 1, function(x) mean(exp(x))),
                                                   label = "R", minQ = NA, maxQ = NA, maxI = NA, meanQ = NA,
                                                   real = simulatedParam$R,
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
}

# Output prior-posterior
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
  ggplot(outputX[, .N, by = .(week, sequence)], aes(x = sequence, y = week, fill = N > 500)) + geom_tile() + theme_laura()
  ggplot(outputX[, .N, by = .(week, sequence)][N>50], aes(x = sequence, y = week, color = N)) + geom_point() + theme_laura()
}

# 15.04.2019
# Output prior-posterior
# a, TauG, TauS, rho, Bk, p
aValToPlot <- seq(qnorm(0.05, mean = constants$ma, sd = constants$sa), qnorm(0.95, mean = constants$ma, sd = constants$sa), length.out = 51)
rValToPlot <- seq(0, qgamma(0.95, shape = constants$aR, rate = constants$bR), length.out = 51)
dataToPlotPos <- rbind(#data.table(x = seq(1,400,5), y = dgamma(seq(1,400,5), shape = a_tau, rate = b_tau), label = "prior", var = "tauG"),
                       #data.table(x = sapply(parameters.stored[burninEnd:storedIterations], function(x) x$tau.G), y = 0, label = "posterior", var="tauG"),
                       #data.table(x = seq(1,400,5), y = dgamma(seq(1,400,5), shape = a_tau, rate = b_tau), label = "prior", var = "tauS"),
                       #data.table(x = sapply(parameters.stored[burninEnd:storedIterations], function(x) x$tau.S), y = 0, label = "posterior", var="tauS"),
                       #data.table(x = seq(1,10,0.5), y = dgamma(seq(1,10,0.5), shape = a_L, rate = b_L), label = "prior", var = "rho"),
                       #data.table(x = sapply(parameters.stored[burninEnd:storedIterations], function(x) x$l), y = 0, label = "posterior", var = "rho"),
                       #data.table(x = seq(0.01, 0.25, 0.01), y = dbeta(seq(0.01, 0.25, 0.01), shape1 = a_p, shape2 = b_p), label = "prior", var = "p"),
                       #data.table(x = sapply(parameters.stored[burninEnd:storedIterations], function(x) x$p), y = 0, label = "posterior", var = "p"),
                       data.table(x = aValToPlot, y = dnorm(aValToPlot, mean = constants$ma, sd = constants$sa), label = "prior", var = "alpha"),
                       data.table(x = outputA, y = 0, label = "posterior", var = "alpha"),
                       data.table(x = rValToPlot, y = dgamma(aValToPlot, shape = constants$ma, rate = constants$sa), label = "prior",var=TeX('$\\tau_R')),
                       data.table(x = outputTauR, y = 0, label = "posterior", var = TeX('$\\tau_R')))
dataToPlotPos[, title := paste0(var, " - ", label)]
gPos <- ggplot(dataToPlotPos, aes(x = x, fill = label, color = label)) + facet_wrap(~ title, scales = "free", ncol = 2) +
  geom_density(data = dataToPlotPos[label == "posterior"], alpha = 0.3) +
  geom_line(aes(y = y)) + theme_laura() + labs(x = "", y = "density", color = "", fill = "")
gPos #savePDF(gPos, "Plot16042019_MCMCposterior_Real", 8, 7)
gPos2 <- ggplot(dataToPlotPos, aes(x = x, fill = label, color = label)) +
  geom_density(data = dataToPlotPos[label == "posterior"], alpha = 0.3) +
  geom_line(aes(y = y)) + theme_laura() + labs(x = "", y = "density", color = "", fill = "")
gPos2

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

