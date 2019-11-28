# Requires: config, sigmaJumps, storage, parameters

#if(config$ifSUpdate && (it%%2)){
if(config$ifSUpdate){
  # I. (normal jumps)
  for (i in 1:numWeeks){
    proposalSN <- rnorm(1, parameters$S[i], sd = sigmaJumps$S[i])
    larSN <- lpriorS(i, proposalSN) + llS(i, proposalSN) - lpriorS(i, parameters$S[i]) - llS(i, parameters$S[i])
    uSN <- runif(1)
    if (log(uSN) < larSN) {
      parameters$S[i] <- proposalSN
      storage$accept$S[i] <- storage$accept$S[i] + 1
      sigmaJumps$S[i] <- sigmaJumps$S[i]*x
      matrixS[i,,] <- proposalSN
      matrixSexp[i,,] <- exp(proposalSN)
    } else {
      storage$reject$S[i] <- storage$reject$S[i] + 1
      sigmaJumps$S[i] <- sigmaJumps$S[i]*x^(-0.7857)
    }
  }
  cat("S normal updated ")
  # II. (conditional proposal)
  #}else if(config$ifSUpdate && (!it %%2)){
}else if(config$ifSUpdate && 0){
  for (i in 1:numWeeks){
    tempS <- parameters$S
    tempS[i] <- 0
    proposalSC <- rnorm(1, -(tempS%*%seasonalCoefficientMatrixSqr[i,])/(1*seasonalCoefficientMatrixSqr[i,i]), # 2 cancels
                        sd = sqrt(1/(parameters$tau.S*seasonalCoefficientMatrixSqr[i,i])))
    larSC <- llS(i, proposalSC) - llS(i, parameters$S[i])
    uSC <- runif(1)
    if (log(uSC) < larSC) {
      parameters$S[i] <- proposalSC
      storage$accept$ScondUp[i] <- storage$accept$ScondUp[i] + 1
      matrixS[i,,] <- proposalSC
      matrixSexp[i,,] <- exp(proposalSC)
    } else {
      storage$reject$ScondUp[i] <- storage$reject$ScondUp[i] + 1
    }
  }
  cat("S cond. updated ")
}
