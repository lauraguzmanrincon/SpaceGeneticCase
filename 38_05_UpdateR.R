# Requires: config, sigmaJumps, storage, parameters

if(config$ifRUpdate && (it%%2)){
#if(config$ifRUpdate){
  # I. (normal jumps)
  for (j in 1:numRegions){
  #for (j in c(257)){
    proposalRN <- rnorm(1, mean = parameters$R[j], sd = sigmaJumps$R[j])
    larRN <- lpriorR(j, proposalRN) + llR(j, proposalRN) - lpriorR(j, parameters$R[j]) - llR(j, parameters$R[j])
    uRN <- runif(1)
    if (log(uRN) < larRN) {
      parameters$R[j] <- proposalRN
      storage$accept$R[j] <- storage$accept$R[j] + 1
      sigmaJumps$R[j] <- sigmaJumps$R[j]*x
      matrixR[,j,] <- proposalRN
      matrixRexp[,j,] <- exp(proposalRN)
      matrixRR[j,] <- proposalRN*matrixForGMRF[j,]
    } else {
      storage$reject$R[j] <- storage$reject$R[j] + 1
      sigmaJumps$R[j] <- sigmaJumps$R[j]*x^(-0.7857)
    }
  }
  cat("R normal updated ")
  # II. (conditional proposal)
  #}else if(config$ifRUpdate && (!it%%2)){
}else if(config$ifRUpdate && (!it%%2)){
  for (j in 1:numRegions){
  #for (j in c(257)){
    proposalRC <- rnorm(n = 1, mean = sum(matrixRR[,j])/neighLength[j], sd = sqrt(1/(parameters$tau.R*neighLength[j])))
    larRC <- llR(j, proposalRC) - llR(j, parameters$R[j])
    uRC <- runif(1)
    if (log(uRC) < larRC) {
      parameters$R[j] <- proposalRC
      storage$accept$RcondUp[j] <- storage$accept$RcondUp[j] + 1
      matrixR[,j,] <- proposalRC
      matrixRexp[,j,] <- exp(proposalRC)
      matrixRR[j,] <- proposalRC*matrixForGMRF[j,]
    } else {
      storage$reject$RcondUp[j] <- storage$reject$RcondUp[j] + 1
    }
  }
  cat("R cond. updated ")
}
