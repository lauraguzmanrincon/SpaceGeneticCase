# Requires: config, sigmaJumps, storage, parameters
# I literally ignored II. and III.

# 1: single updates.  2: 'global' update.  3: block update.

#adaptTypeG <- 3
adaptTypeG <- ifelse(it%%2, 3, 1)
parameters$cutHeight <- 0 # only for adaptTypeG = 3 # 27.02.2020

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
  }else if(adaptTypeG == 3){
    # IV. Block update as in Knorr-Held,1999 and Fahrmeir,2001
    # This code was constructed in sandbox 45.R/2.
    
    # Random cut
    # treeCuts has three levels: original groups (groupDupId), low group (dim3Cases), and cut clusters (cutClust)
    #cutHeight <- sample(10:130, size = 1, replace = FALSE)
    cutHeight <- sample(10:1000, size = 1, replace = FALSE)
    treeCuts <- data.table(groupDupId = 1:length(hClustOut$order), cutClust = cutree(hClustOut, h = cutHeight))
    # Merge cut with dimension (as in 1.b.)
    temp <- casesForModels[, .N, .(groupDupId, dim3Cases)]
    setkey(treeCuts, groupDupId)
    setkey(temp, groupDupId)
    treeCuts[temp, dim3Cases := dim3Cases]
    #treeCuts[, .N, .(cutClust, dim3Cases)]
    
    listBlocks <- treeCuts[, .N, cutClust][order(cutClust)]
    listBlocks[, id := .I]
    cutToDim3 <- treeCuts[, .N, .(cutClust, dim3Cases)]
    
    parameters$cutHeight <- cutHeight # 27.02.2020
    
    # Blocks loop
    # K = invDeltaMatrix = AMatrix%*%solve(deltaMatrixProposal)%*%AMatrix that is already computed at this point
    # Extract indexes of the chosen block (clustIndex): dim3ClusterIndexes
    # Compute solve(K[indexes,indexes]), K[indexes,-indexes], G[-indexes]
    # Compute mean and covariance of prior conditional
    # Propose a jump from prior conditional
    # Compute acceptance ratio using the likelihood only
    # Draw random number to accept/reject
    # Accept/reject
    for(clustIndex in 1:nrow(listBlocks)){
      # Note that for blocks of size 1, this is not equivalent as single RW jumps, since conditional updates depend on the neighbours state
      cutCluster <- listBlocks[id == clustIndex, cutClust]
      dim3ClusterIndexes <- cutToDim3[cutClust == cutCluster, dim3Cases] # TODO should be ordered? don't think so...
      dim3ClusterComplement <- setdiff(1:numDims[3], dim3ClusterIndexes)
      
      muA <- solve(invDeltaMatrix[dim3ClusterIndexes, dim3ClusterIndexes]) # TODO invertible? Prove!
      muB <- invDeltaMatrix[dim3ClusterIndexes, dim3ClusterComplement]
      muC <- parameters$G[dim3ClusterComplement]
      
      muSample <- muA%*%muB%*%muC
      sigmaSample <- muA/parameters$tau.G
      proposalGKnorr <- mvrnorm(n = 1, mu = muSample, Sigma = sigmaSample)
      
      larGKnorr <- llGKnorr(dim3ClusterIndexes, proposalGKnorr) - llGKnorr(dim3ClusterIndexes, parameters$G[dim3ClusterIndexes])
      
      u <- runif(1)
      if (log(u) < larGKnorr) {
        parameters$G[dim3ClusterIndexes] <- proposalGKnorr
        storage$accept$GcondUp[dim3ClusterIndexes] <- storage$accept$GcondUp[dim3ClusterIndexes] + 1
        #NO?sigmaJumps$G[dim3ClusterIndexes] <- sigmaJumps$G[dim3ClusterIndexes]*x
        matrixG[,,dim3ClusterIndexes] <- proposalGKnorr
        matrixGexp[,,dim3ClusterIndexes] <- exp(proposalGKnorr)
      } else {
        storage$reject$GcondUp[dim3ClusterIndexes] <- storage$reject$GcondUp[dim3ClusterIndexes] + 1
        #NO?sigmaJumps$G[dim3ClusterIndexes] <- sigmaJumps$G[dim3ClusterIndexes]*x^(-0.7857)
      }
    }
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
