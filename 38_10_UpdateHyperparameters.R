# Requires: config, sigmaJumps, storage, parameters
# Hyperparameters: Tr, Ts, Tk, l, p
# TODO

# Update l and tauG (block) ----

if(config$ifLTauGUpdate && it > 20){ # ifLTauGUpdate    ifLTauGUpdate && it%%20 == 0
  # L update (positive proposal)
  repeat{
    proposJump <- rnorm(1, parameters$l, sd = sigmaJumps$l)
    if(proposJump > 0){
      proposalL <- proposJump
      break
    }
  }
  #deltaMatrixProposal <- deltaMatrixFn(proposalL, maternParameter) # + 0.5*diag(numSequences)
  # SPAM deltaMatrixProposal <- as.spam(deltaMatrixFn(proposalL, distanceMatrixMCMC, sqrDistanceMatrix, config$maternParameter) + 0.01*diag(numSequences))
  deltaMatrixProposal <- deltaMatrixFn(proposalL, distanceMatrixMCMC, sqrDistanceMatrix, config$maternParameter) + 0.01*diag(numSequences) # ??? 03.03.2020 17.45
  #invDeltaMatrixProposal <- as.spam(solve(deltaMatrixProposal)) # slow if very large l
  invDeltaMatrixProposal <- AMatrix%*%solve(deltaMatrixProposal)%*%AMatrix # 22.01.2020 new precision matrix, no spam (see 43.R)
  
  # Tau update
  extraTerm <- 0.5*(t(parameters$G)%*%invDeltaMatrixProposal%*%parameters$G) # MUST be greater than 0
  proposalTau <- rgamma(1, shape = constants$aG + (numSequences/2), rate = constants$bG + extraTerm)
  
  cat(sprintf("%.*f", 2, proposalL), " ", sprintf("%.*f", 2, proposalTau), " ")
  
  # Accept/reject
  normConstantRatioLJump <- log(pnorm(parameters$l)) - log(pnorm(proposalL))
  larLT <- lpriorLTauGBlock(invDeltaMatrixProposal, proposalTau, proposalL) -
    lpriorLTauGBlock(invDeltaMatrix, parameters$tau.G, parameters$l) + normConstantRatioLJump # 22.01.2020 no as.spam(invDeltaMatrix) (see 43.R)
  #temp <- as.spam(invDeltaMatrix)
  #0.5*determinant.spam(temp)$modulus - 0.5*proposalTau*(t(parameters$G)%*%temp%*%parameters$G)
  #0.5*determinant.spam(invDeltaMatrixProposal)$modulus - 0.5*proposalTau*(t(parameters$G)%*%invDeltaMatrixProposal%*%parameters$G)
  # shall we assume det is always positive? TODO
  cat(sprintf("%.*f", 2, larLT), " ")
  uL <- runif(1)
  if (log(uL) < larLT) {
    parameters$l <- proposalL
    parameters$tau.G <- proposalTau
    storage$accept$l <- storage$accept$l + 1
    sigmaJumps$l <- sigmaJumps$l*x # TODO ???
    deltaMatrix <- deltaMatrixProposal # debug 18.12.2019 :)
    invDeltaMatrix <- invDeltaMatrixProposal # debug 18.12.2019 :)
    cat("LT accepted ")
  } else {
    storage$reject$l <- storage$reject$l + 1
    #sigmaJumps$l <- ifelse(sigmaJumps$l > 0.05, sigmaJumps$l*x^(-0.7857), sigmaJumps$l)
    sigmaJumps$l <- sigmaJumps$l*x^(-0.7857) # adapted jump added on 02.03.2020
    cat("LT rejected ")
  }
}

# Update tauS (gamma) ----
if(config$ifTauSUpdate){
  extraTermS <- 0.5*(t(parameters$S)%*%seasonalCoefficientMatrixSqr%*%parameters$S)
  proposalTauS <- rgamma(1, shape = constants$aS + ((numWeeks - 2)/2), rate = constants$bS + extraTermS)
  parameters$tau.S <- proposalTauS
}

# Update tauR (gamma) ----
if(config$ifTauRUpdate){
  # check notes 26.06.2019
  extraTermR <- 0.25*sum((matrixRR - t(matrixRR))^2)
  proposalTauR <- rgamma(1, shape = constants$aR + 0.5*numRegions, rate = constants$bR + extraTermR)
  parameters$tau.R <- proposalTauR
}

# Update p (beta) ----
if(config$ifPUpdate){
  # Discussion on whether or not we should use Xij or xsigmatheta... conclusion: the second
  #str(new$mod[[1]][[1]]$X)
  parameters$p <- rbeta(1, shape1 = sum(parameters$X) + constants$aP, shape2 = sum(1 - parameters$X) + constants$bP)
}






# OUTDATED Update l (truncated normal jumps) ----
# OUTDATED: has a truncated jump with M
# TODO update
if(config$ifLUpdate){
  repeat{
    proposJump <- rnorm(1, parameters$l, sd = 1)
    if(proposJump > 0 & proposJump < M){
      proposalL <- proposJump
      break
    }
  }
  
  #deltaMatrixProposal <- deltaMatrixFn(proposalL, maternParameter) # + 0.5*diag(numSequences)
  deltaMatrixProposal <- as.spam(deltaMatrixFn(proposalL, maternParameter) + 0.01*diag(numSequences))
  invDeltaMatrixProposal <- as.spam(solve(deltaMatrixProposal)) # slow if very large l
  normConstantRatioLJump <- log(pnorm(M - parameters$l) - pnorm(- parameters$l)) - log(pnorm(M - proposalL) - pnorm(- proposalL))
  lar <- lpriorl(invDeltaMatrixProposal) - lpriorl(as.spam(invDeltaMatrix)) + normConstantRatioLJump
  uL <- runif(1)
  if (log(uL) < lar) {
    parameters$l <- proposalL
    accept$l <- accept$l + 1
    l.sigma <- l.sigma*x
  } else {
    reject$l <- reject$l + 1
    l.sigma <- l.sigma*x^(-0.7857)
  }
}

# UNUSED Update tauG (gamma) ----
# TODO update
if(config$ifTauGUpdate){
  extraTerm <- 0.5*(t(parameters$G)%*%invDeltaMatrix%*%parameters$G) # MUST be greater than 0
  proposalTau <- rgamma(1, shape = a_tau + (numSequences/2), rate = b_tau + extraTerm)
  parameters$tau.G <- proposalTau
}
