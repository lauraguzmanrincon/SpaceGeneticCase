# Requires: config, sigmaJumps, storage, parameters
# I literally ignored II. and III.

# 1: single updates.  2: block update

adaptTypeG <- 1

if(config$ifGUpdate){
  if(adaptTypeG == 1){
    # I. (one-by-one - normal jumps)
    for (k in 1:numSequences) {
    #for (k in c(28)) {
      proposalG <- rnorm(1, parameters$G[k], sd = sigmaJumps$G[k])
      # CHOOSE:
      # Gaussian process for G prior:
      larG <- lpriorG(k, proposalG) + llG(k, proposalG) - lpriorG(k, parameters$G[k]) - llG(k, parameters$G[k])
      if(k == 1){
        cat(sprintf("%.*f", 2, larG))
      }
      # 02042019 Independent priors for G:
      #larG <- lpriorGIndependent(k, proposalG) + llG(k, proposalG) - lpriorGIndependent(k, parameters$G[k]) - llG(k, parameters$G[k])
      u <- runif(1)
      if (log(u) < larG) {
        parameters$G[k] <- proposalG
        storage$accept$G[k] <- storage$accept$G[k] + 1
        sigmaJumps$G[k] <- sigmaJumps$G[k]*x
        matrixG[,,k] <- proposalG
        matrixGexp[,,k] <- exp(proposalG)
      } else {
        storage$reject$G[k] <- storage$reject$G[k] + 1
        sigmaJumps$G[k] <- sigmaJumps$G[k]*x^(-0.7857)
      }
    }
    cat("G updated single ")
  }else if(adaptTypeG == 2){
    # copied from code below on 18.07.2019
    #proposalG <- rmvnorm.spam(n = 1, Sigma = deltaMatrix/parameters$tau.G)[1,]
    proposalG <- rmvnorm.spam(n = 1, Q = invDeltaMatrix*parameters$tau.G)[1,] # 19.02.2020
  }
  
  
  
  
  
  
}else if(config$ifGUpdate && 0 && (it%%freqGUpdate == 1 | freqGUpdate == 1)){
  # II. Block update G and l (G: using Ruew,2001 algorithm described in Knorr-Held,2002 and using SPAM package)
  # Proposal l
  repeat{
    proposJump <- rnorm(1, parameters$l, sd = l.sigma)
    if(proposJump > 0 & proposJump < M){
      proposalL <- proposJump
      break
    }
  }
  
  # Proposal G
  Q <- chol(as.spam(parameters$tau.G*invDeltaMatrix))
  cat("chol done ")
  z <- rnorm(numSequences, mean = 0, sd = 1)
  proposalG <- backsolve.spam(Q, z)
  
  # Accept/reject
  deltaMatrixProposal <- as.spam(deltaMatrixFn(proposalL, maternParameter) + 0.01*diag(numSequences))
  invDeltaMatrixProposal <- as.spam(solve(deltaMatrixProposal)) # slow if very large l
  matrixS <- matrix(parameters$S, nrow = numWeeks, ncol = numSequences, byrow = FALSE)
  matrixB <- matrix(parameters$B, nrow = numWeeks, ncol = numSequences, byrow = TRUE)
  normConstantRatioLJump <- log(pnorm(M - parameters$l) - pnorm(- parameters$l)) - log(pnorm(M - proposalL) - pnorm(- proposalL))
  cat("proposal inversion done ")
  lar <- llGBlock(proposalG, matrixS, matrixB, parameters$X) - llGBlock(parameters$G, matrixS, matrixB, parameters$X) +
    lpriorLGBlock(proposalG, invDeltaMatrixProposal) - lpriorLGBlock(parameters$G, as.spam(invDeltaMatrix)) + normConstantRatioLJump
  cat("lar for G and l done ", exp(lar), " ")
  u <- runif(1)
  if (log(u) < lar) {
    parameters$l <- proposalL
    parameters$G <- proposalG
    accept$l <- accept$l + 1
    accept$G <- accept$G + 1
    l.sigma <- l.sigma*x # TODO adaptation of l jump if block updates with g??
    cat("G and rho accepted ")
  } else {
    reject$l <- reject$l + 1
    reject$G <- reject$G + 1
    l.sigma <- l.sigma*x^(-0.7857)
    cat("G and rho rejected ")
  }
  
}else if(config$ifGUpdate && 0 && (it%%freqGUpdate == 1 | freqGUpdate == 1)){
  # III. Block update G and l (G: normal jumps)
  # Proposal l
  repeat{
    proposJump <- rnorm(1, parameters$l, sd = l.sigma)
    if(proposJump > 0 & proposJump < M){
      proposalL <- proposJump
      break
    }
  }
  proposalL <- parameters$l
  
  # Proposal G (adaptive)
  if(!adaptiveMCMC){
    #proposalG <- rep(0, numSequences)
    #for (k in 1:numSequences) {
    #  proposalG[k] <- rnorm(1, parameters$G[k], sd = G.sigma[k])
    #}
    proposalG <- rmvnorm.spam(n = 1, Sigma = deltaMatrix/parameters$tau.G)[1,]
  }else{
    nPrev <- it - 1 # how many we are looking at
    cumulativeAlphaSum_nminus1 <- cumulativeAlphaSum_n
    cumulativeAlphaSum_n <- cumulativeAlphaSum_n + parameters$G
    if(it > 2){
      cumulativeMean_n <- cumulativeAlphaSum_n/nPrev
      cumulativeMean_nminus1 <- cumulativeAlphaSum_nminus1/(nPrev - 1)
      cumulativeCovariance <- (nPrev - 1)/nPrev*cumulativeCovariance +
        sdim/nPrev*(nPrev*cumulativeMean_nminus1%*%t(cumulativeMean_nminus1) -
                      (nPrev + 1)*cumulativeMean_n%*%t(cumulativeMean_n) +
                      parameters$G%*%t(parameters$G) +
                      eps*diag(1, numSequences, numSequences))
      cat("cumulating G! ")
    }
    if(it < 50 | accept$G[1] == 0){
      proposalG <- rnorm(numSequences, mean = parameters$G, sd = cSdProp/parameters$tau.G)
      cat("first G! ")
    }else{
      cumulativeCovarianceSpam <- as.spam(cumulativeCovariance)
      proposalG <- mvrnorm(1, parameters$G, cumulativeCovariance)
      cat("adaptive G! ")
    }
  }
  
  # Accept/reject TODOOO
  #deltaMatrixProposal <- as.spam(deltaMatrixFn(proposalL, maternParameter) + 0.01*diag(numSequences))
  cat("n")
  #invDeltaMatrixProposal <- as.spam(solve(deltaMatrixProposal)) # slow if very large l
  deltaMatrixProposal <- as.spam(deltaMatrix)
  invDeltaMatrixProposal <- as.spam(invDeltaMatrix)
  cat("n")
  matrixS <- matrix(parameters$S, nrow = numWeeks, ncol = numSequences, byrow = FALSE)
  matrixB <- matrix(parameters$B, nrow = numWeeks, ncol = numSequences, byrow = TRUE)
  normConstantRatioLJump <- log(pnorm(M - parameters$l) - pnorm(- parameters$l)) - log(pnorm(M - proposalL) - pnorm(- proposalL))
  cat("proposal inversion done ")
  lar <- llGBlock(proposalG, matrixS, matrixB, parameters$X) - llGBlock(parameters$G, matrixS, matrixB, parameters$X) +
    lpriorLGBlock(proposalG, invDeltaMatrixProposal) - lpriorLGBlock(parameters$G, as.spam(invDeltaMatrix)) + normConstantRatioLJump
  cat(formatC(lar, format = "e", digits = 2))
  #cat(formatC(lpriorLGBlock(proposalG, invDeltaMatrixProposal), format = "e", digits = 2))
  #cat(formatC(lpriorLGBlock(parameters$G, as.spam(invDeltaMatrix)), format = "e", digits = 2))
  if(it == 11) Sys.sleep(10)
  cat("lar for G and l done ", exp(lar), " ")
  u <- runif(1)
  if (log(u) < lar) {
    parameters$l <- proposalL
    parameters$G <- proposalG
    accept$l <- accept$l + 1
    accept$G <- accept$G + 1
    l.sigma <- l.sigma*x
    G.sigma <- G.sigma*x
    cat("G and rho ad. accepted ")
  } else {
    reject$l <- reject$l + 1
    reject$G <- reject$G + 1
    l.sigma <- l.sigma*x^(-0.7857)
    G.sigma <- G.sigma*x^(-0.7857)
    cat("G and rho ad. rejected ")
  }
}
