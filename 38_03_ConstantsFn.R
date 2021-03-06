initConstants <- function(ma = 0,
                       sa = 10,
                       aR = 1,
                       bR = 0.1,
                       aS = 1,
                       bS = 0.01,
                       aG = 1,
                       bG = 0.01,
                       aL = 5,
                       bL = 1,
                       M = NULL, # unused
                       aB = 1,
                       bB = 0.5,
                       aP = 1,
                       bP = ifelse(dimBeta != 123, numBlockDims[dimBeta] - 1, prod(numBlockDims)),
                       aP10 = 2,
                       bP10 = 2,
                       aP01 = 1,
                       bP01 = 51,
                       ifPlot = FALSE){
  constants <- list(ma = ma,
                 sa = sa,
                 va = sa^2,
                 aR = aR,
                 bR = bR,
                 aS = aS,
                 bS = bS,
                 aG = aG,
                 bG = bG,
                 aL = aL,
                 bL = bL,
                 M = M,
                 aB = aB,
                 bB = bB,
                 aP = aP,
                 bP = bP,
                 aP10 = aP10,
                 bP10 = bP10,
                 aP01 = aP01,
                 bP01 = bP01)
  
  constants$plots <- list()
  if(ifPlot){
    pA <- plotNorm(mean = ma, sd = sa, xlab = TeX('$\\alpha'))
    pR <- plotGamma(shape = aR, rate = bR, xlab = TeX('$\\tau_R'))
    pS <- plotGamma(shape = aS, rate = bS, xlab = TeX('$\\tau_S'))
    pG <- plotGamma(shape = aG, rate = bG, xlab = TeX('$\\tau_G'))
    pL <- plotGammaTruncated(shape = aG, rate = bG, upp = M, xlab = TeX('$\\tau_\\rho'))
    pB <- plotGamma(shape = aB, rate = bB, xlab = "B")
    pP <- plotBeta(shape1 = aP, shape2 = bP, xlab = "X")
    pP10 <- plotBeta(shape1 = aP10, shape2 = bP10, xlab = "X_10")
    pP01 <- plotBeta(shape1 = aP01, shape2 = bP01, xlab = "X_01")
    constants$plots <- list(pA, pR, pS, pG, pL, pB, pP, p10, p01)
  }
  
  return(constants)
}
#multiplot(plotlist = ttt$plots, cols = 2)

initConfig <- function(numIterations = 100,
                       thinning = 1,
                       burnIn = 50,
                       maternParameter = 1/2,
                       adaptiveMCMC = 0,
                       ifAUpdate = 1,
                       ifRUpdate = 1,
                       ifSUpdate = 1,
                       ifGUpdate = 1,
                       ifBUpdate = 1,
                       ifXUpdate = 1,
                       ifPUpdate = 1,
                       ifTauRUpdate = 1,
                       ifTauSUpdate = 1,
                       ifTauGUpdate = 0, # not available. use 0
                       ifLUpdate = 0, # not available. use 0
                       ifLTauGUpdate = 1){
  config <- list(
    # Algorithm specifications
    numIterations = numIterations,
    thinning = thinning,
    burnIn = burnIn,
    maternParameter = maternParameter,
    adaptiveMCMC = adaptiveMCMC,
    #freqGUpdate = freqGUpdate,
    
    # Turn on/off updates
    ifAUpdate = ifAUpdate,
    ifRUpdate = ifRUpdate,
    ifSUpdate = ifSUpdate,
    ifGUpdate = ifGUpdate,
    ifBUpdate = ifBUpdate,
    ifXUpdate = ifXUpdate,
    ifPUpdate = ifPUpdate,
    ifTauRUpdate = ifTauRUpdate,
    ifTauSUpdate = ifTauSUpdate,
    ifTauGUpdate = ifTauGUpdate,
    ifLUpdate = ifLUpdate,
    ifLTauGUpdate = ifLTauGUpdate && !ifTauGUpdate && !ifLUpdate
  )
  if(ifLTauGUpdate && (ifTauGUpdate || ifLUpdate)){
    warning("ifLTauGUpdate ignored")
  }
  
  # Modify according to the dimensions included
  if(!1 %in% dimToInclude){
    config$ifSUpdate <- 0
    config$ifTauSUpdate <- 0
  }
  if(!2 %in% dimToInclude){
    config$ifRUpdate <- 0
    config$ifTauRUpdate <- 0
  }
  if(!3 %in% dimToInclude){
    config$ifGUpdate <- 0
    config$ifTauGUpdate <- 0
    config$ifLUpdate <- 0
    config$ifLTauGUpdate <- 0
  }
  
  return(config)
}

#' parametersVar1: mean
#' parametersVar2: sd
#' Note: configUpdates should be the output of initConfig() instead of using the default 13.12.2019
initParam <- function(parametersVar1 = list(a = 0, R = 0, S = 0, G = 0, B = 0, tR = 3.5, tS = 3.5, tG = 0.5, p = 1, l = 2.5, p10 = 2, p01 = 1),
                      parametersVar2 = list(a = 1, R = 1, S = 1, G = 1, B = 1, tR = 1, tS = 1, tG = 0.15, p = 100, l = 0.8, p10 = 2, p01 = 51),
                      a.sigma = 1,
                      R.sigma = rep(1, numRegions),
                      S.sigma = rep(1, numWeeks),
                      G.sigma = rep(1, numSequences),
                      B.sigma = NULL,
                      #B.sigma = matrix(1, nrow = numRegionsGroups, ncol = numSequenceGroups),
                      #B.sigma = rep(1, numBeta),
                      l.sigma = 1,
                      p10.sigma = 1,
                      p01.sigma = 1,
                      ifPlot = FALSE,
                      configUpdates = list(ifAUpdate = 1, ifRUpdate = 1, ifSUpdate = 1, ifGUpdate = 1, ifBUpdate = 1, ifXUpdate = 1,
                                           ifPUpdate = 1, ifTauRUpdate = 1, ifTauSUpdate = 1, ifLTauGUpdate = 1)){
  
  if(is.null(B.sigma) & dimBeta == 123){
    B.sigma <- array(1, dim = c(numWeeksGroups, numRegionsGroups, numSequenceGroups))
  }else if(is.null(B.sigma) & !dimBeta == 123){
    B.sigma = rep(1, numBeta)
  }
  
  # Parameters and variables initialisation
  if(simulatedData == 0){
    parameters <- list(a = rnorm(1, mean = parametersVar1$a, sd = parametersVar2$a),
                       R = rnorm(numRegions, mean = parametersVar1$R, sd = parametersVar2$R),
                       S = rnorm(numWeeks, mean = parametersVar1$S, sd = parametersVar2$S),
                       G = rnorm(numSequences, mean = parametersVar1$G, sd = parametersVar2$G),
                       X = array(0, dim = c(numWeeksGroups, numRegionsGroups, numSequenceGroups)),
                       B = exp(rnorm(numBeta, mean = parametersVar1$B, sd = parametersVar2$B)),
                       tau.R = exp(rnorm(1, parametersVar1$tR, parametersVar2$tR)),
                       tau.S = exp(rnorm(1, parametersVar1$tS, parametersVar2$tS)),
                       tau.G = exp(rnorm(1, parametersVar1$tG, parametersVar2$tG)),
                       p = rbeta(1, shape1 = parametersVar1$p, shape2 = parametersVar2$p), #plotBeta(1,100)
                       #l = 3,
                       l = exp(rnorm(1, parametersVar1$l, parametersVar2$l)), # plotLogNorm
                       cutHeight = 0, # 27.02.2020
                       sBlockSize = c(0,0)) # 17.03.2020
    #print(parameters$l)
    parameters$G <- parameters$G - mean(parameters$G) # update G inititalisation 26.02.2020
    if(dimBeta == 123){
      parameters$B <- array(exp(rnorm(prod(numBlockDims), mean = parametersVar1$B, sd = parametersVar2$B)), dim = c(numWeeksGroups, numRegionsGroups, numSequenceGroups))
    }
    if(autocorrPrior == 1){ # autocorrelation X update 13.04.2020
      parameters$p <- NULL
      parameters$p10 <- rbeta(1, shape1 = parametersVar1$p10, shape2 = parametersVar2$p10)
      parameters$p01 <- rbeta(1, shape1 = parametersVar1$p01, shape2 = parametersVar2$p01)
    }
    
    # Modify according to the dimensions included
    if(!1 %in% dimToInclude) parameters$S <- rep(0, numWeeks)
    if(!2 %in% dimToInclude) parameters$R <- rep(0, numRegions)
    if(!3 %in% dimToInclude) parameters$G <- rep(0, numSequences)
    # Note hyperparameters are not modified
  }else{
    # Parameters if using simulated data OLD (stored in inSumilated#) ... or even NEW
    #parameters <- list(a = ifAUpdate*rnorm(1, mean = 0, sd = 1) + (!ifAUpdate)*inSimulatedA[1,1],
    #                   G = ifGUpdate*rnorm(numSequences, mean = -1.1, sd = 0.5) + (!ifGUpdate)*inSimulatedG[1,],
    #                   S = ifSUpdate*rnorm(numWeeks, mean = -3, sd = 0.5) + (!ifSUpdate)*inSimulatedS[,1],
    #                   X = ifXUpdate*matrix(0, nrow = numWeeks, ncol = numSequences) + (!ifXUpdate)*inSimulatedX,
    #                   B = ifBUpdate*rnorm(numBeta, mean = 3, sd = 0.5) + (!ifBUpdate)*inSimulatedB[1,],
    #                   tau.S = 220,
    #                   tau.G = 1.1,
    #                   p = ifPUpdate*0.01 + (!ifPUpdate)*sum(inSimulatedX)/(numWeeks*numSequences),
    #                   l = 3)#exp(rnorm(1, 1, 0.5)))

    # OR
    #parameters <- simulatedParam
    # OR
    parameters <- list(a = configUpdates$ifAUpdate*rnorm(1, mean = parametersVar1$a, sd = parametersVar2$a) + (!configUpdates$ifAUpdate)*simulatedParam$a,
                       R = configUpdates$ifRUpdate*rnorm(numRegions, mean = parametersVar1$R, sd = parametersVar2$R) + (!configUpdates$ifRUpdate)*simulatedParam$R,
                       S = configUpdates$ifSUpdate*rnorm(numWeeks, mean = parametersVar1$S, sd = parametersVar2$S) + (!configUpdates$ifSUpdate)*simulatedParam$S,
                       G = configUpdates$ifGUpdate*rnorm(numSequences, mean = parametersVar1$G, sd = parametersVar2$G) + (!configUpdates$ifGUpdate)*simulatedParam$G,
                       X = configUpdates$ifXUpdate*array(0, dim = c(numWeeksGroups, numRegionsGroups, numSequenceGroups)) + (!configUpdates$ifXUpdate)*simulatedParam$X,
                       B = configUpdates$ifBUpdate*exp(rnorm(numBeta, mean = parametersVar1$B, sd = parametersVar2$B)) + (!configUpdates$ifBUpdate)*simulatedParam$B,
                       tau.R = configUpdates$ifTauRUpdate*exp(rnorm(1, parametersVar1$tR, parametersVar2$tR)) + (!configUpdates$ifTauRUpdate)*simulatedParam$tau.R,
                       tau.S = configUpdates$ifTauSUpdate*exp(rnorm(1, parametersVar1$tS, parametersVar2$tS)) + (!configUpdates$ifTauSUpdate)*simulatedParam$tau.S,
                       tau.G = configUpdates$ifLTauGUpdate*exp(rnorm(1, parametersVar1$tG, parametersVar2$tG)) + (!configUpdates$ifLTauGUpdate)*simulatedParam$tau.G,
                       p = configUpdates$ifPUpdate*rbeta(1, shape1 = parametersVar1$p, shape2 = parametersVar2$p) + (!configUpdates$ifPUpdate)*simulatedParam$p, #plotBeta(1,100)
                       l = configUpdates$ifLTauGUpdate*exp(rnorm(1, parametersVar1$l, parametersVar2$l)) + (!configUpdates$ifLTauGUpdate)*simulatedParam$l) # plotLogNorm
    if(configUpdates$ifGupdate) parameters$G <- parameters$G - mean(parameters$G) # update G inititalisation 26.02.2020
    if(dimBeta == 123){
      parameters$B <- configUpdates$ifBUpdate*array(exp(rnorm(prod(numBlockDims), mean = parametersVar1$B, sd = parametersVar2$B)),
                                                    dim = c(numWeeksGroups, numRegionsGroups, numSequenceGroups)) + (!configUpdates$ifBUpdate)*simulatedParam$B
    }
    if(autocorrPrior == 1){ # autocorrelation X update 13.04.2020
      parameters$p <- NULL
      parameters$p10 <- configUpdates$ifPUpdate*rbeta(1, shape1 = parametersVar1$p10, shape2 = parametersVar2$p10) + (!configUpdates$ifPUpdate)*simulatedParam$p10
      parameters$p01 <- configUpdates$ifPUpdate*rbeta(1, shape1 = parametersVar1$p01, shape2 = parametersVar2$p01) + (!configUpdates$ifPUpdate)*simulatedParam$p01
    }
    
    # ???
    #if(0){
    #  #str(parameters)
    #  #str(simulatedParam)
    #  parameters <- simulatedParam
    #  if(config$ifXUpdate){
    #    parameters$X <- matrix(0, nrow = numWeeksGroups, ncol = numRegionsGroups)
    #  }
    #  if(config$ifLTauGUpdate){
    #    parameters$l <- 1
    #    parameters$tau.G <- 1
    #  }
    #}
  }
  
  sigmaJumps <- list(a = a.sigma,
                     R = R.sigma,
                     S = S.sigma,
                     G = G.sigma,
                     B = B.sigma,
                     l = l.sigma,
                     p10 = p10.sigma,
                     p01 = p01.sigma)
  
  param <- list(parameters = parameters,
                # Jump variance
                sigmaJumps = sigmaJumps)
  
  param$plots <- list()
  if(0*ifPlot){
    # TODO do plots
    pA <- plotNorm(mean = ma, sd = sa, xlab = TeX('$\\alpha'))
    pR <- plotGamma(shape = aR, rate = bR, xlab = TeX('$\\tau_R'))
    pS <- plotGamma(shape = aS, rate = bS, xlab = TeX('$\\tau_S'))
    pG <- plotGamma(shape = aG, rate = bG, xlab = TeX('$\\tau_G'))
    pL <- plotGammaTruncated(shape = aG, rate = bG, upp = M, xlab = TeX('$\\tau_\\rho'))
    pB <- plotGamma(shape = aB, rate = bB, xlab = "B")
    pP <- plotBeta(shape1 = aP, shape2 = bP, xlab = "X")
    constants$plots <- list(pA, pR, pS, pG, pL, pB, pP)
  }
  
  return(param)
}

createStorage <- function(config){
  numIterations <- config$numIterations
  accept <- list(a = 0,
                 R = rep(0, numRegions),
                 RcondUp = rep(0, numRegions),
                 S = rep(0, numWeeks),
                 ScondUp = rep(0, numWeeks),
                 G = rep(0, numSequences),
                 GcondUp = rep(0, numSequences),
                 X = array(0, dim = c(numWeeksGroups, numRegionsGroups, numSequenceGroups)),
                 ##B = matrix(0, nrow = numRegionsGroups, ncol = numSequenceGroups),
                 B = rep(0, numBeta),
                 tau.R = 0,
                 tau.S = 0,
                 tau.G = 0,
                 p = 0,
                 l = 0)
  reject <- list(a = 0,
                 R = rep(0, numRegions),
                 RcondUp = rep(0, numRegions),
                 S = rep(0, numWeeks),
                 ScondUp = rep(0, numWeeks),
                 G = rep(0, numSequences),
                 GcondUp = rep(0, numSequences),
                 X = array(0, dim = c(numWeeksGroups, numRegionsGroups, numSequenceGroups)),
                 ##B = matrix(0, nrow = numRegionsGroups, ncol = numSequenceGroups),
                 B = rep(0, numBeta),
                 tau.R = 0,
                 tau.S = 0,
                 tau.G = 0,
                 p = 0,
                 l = 0)
  #parameters.stored <- vector("list", length = numIterations)
  parameters.stored <- list(a = rep(0, numIterations),
                            R = matrix(0, nrow = numRegions, ncol = numIterations),
                            S = matrix(0, nrow = numWeeks, ncol = numIterations),
                            G = matrix(0, nrow = numSequences, ncol = numIterations),
                            X = vector("list", length = numIterations),
                            #B = array(0, dim = c(numRegionsGroups, numSequenceGroups, numIterations)),
                            B = matrix(0, nrow = numBeta, ncol = numIterations),
                            tau.R = rep(0, numIterations),
                            tau.S = rep(0, numIterations),
                            tau.G = rep(0, numIterations),
                            p = rep(0, numIterations),
                            l = rep(0, numIterations),
                            #scases = matrix(0, nrow = numWeeks, ncol = numIterations),
                            #ecases = matrix(0, nrow = numWeeks, ncol = numIterations),
                            cutHeight = rep(0, numIterations), # 27.02.2020
                            sBlockSize = matrix(0, nrow = 2, ncol = numIterations)) # block size, starting point 17.03.2020
  if(dimBeta == 123){
    accept$B <- array(0, dim = c(numWeeksGroups, numRegionsGroups, numSequenceGroups))
    reject$B <- array(0, dim = c(numWeeksGroups, numRegionsGroups, numSequenceGroups))
    parameters.stored$B <- array(0, dim = c(numWeeksGroups, numRegionsGroups, numSequenceGroups))
  }
  if(autocorrPrior == 1){ # autocorrelation X update 13.04.2020
    accept$p <- NULL
    reject$p <- NULL
    parameters.stored$p <- NULL
    accept$p10 <- 0
    reject$p10 <- 0
    parameters.stored$p10 <- rep(0, numIterations)
    accept$p01 <- 0
    reject$p01 <- 0
    parameters.stored$p01 <- rep(0, numIterations)
  }
  
  storage <- list(accept = accept,
                  reject = reject,
                  parameters = parameters.stored)
  
  return(storage)
}

adaptiveConfig <- function(ifAdaptive = FALSE){
  # TODO I am not using this for now... organise it
  if(ifAdaptive){
    # Auxiliar storage for adaptive updates (from 05_MCMC/RCode_201212/05.R/02.)
    cSdProp <- 3 # 0.05 TODO what's that???
    if(adaptiveMCMC){
      # Based on info found in 05_MCMC/AdMCMC_Lecture.pdf or http://www.mafy.lut.fi/study/sam/statmode09.pdf
      # I kept the notation (nPrev, nminus1, n, ...)
      sdim <- 2.4^2/numSequences
      eps <- 0.01
      cumulativeAlphaSum_n <- rep(0, numSequences)
      cumulativeAlphaSum_nminus1 <- rep(0, numSequences)
      cumulativeCovariance <- matrix(0, numSequences, numSequences)
      
      config <- list(cSdProp = cSdProp,
                     sdim = sdim,
                     eps = eps,
                     cumulativeAlphaSum_n = cumulativeAlphaSum_n,
                     cumulativeAlphaSum_nminus1 = cumulativeAlphaSum_nminus1,
                     cumulativeCovariance = cumulativeCovariance)
    }
  }else{
    config <- NULL
  }
  return(config)
}

