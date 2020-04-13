# Requires: config, sigmaJumps, storage, parameters
# I only modified the auxTerm. Check notes on 30.03.2019 to recall
# NOO Also, note that X is now an array, and so the auxiliar matrices

# Update X block (at the same time since conditionally independent) ----
if(config$ifXUpdate){
  #matrixS <- array(parameters$S, dim = c(numWeeks, numRegions, numSequences))
  #matrixR <- aperm(array(parameters$R, dim = c(numRegions, numWeeks, numSequences)), perm = c(2,1,3))
  #matrixB <- aperm(array(parameters$B[iToGroups], dim = c(numRegions, numWeeks, numSequences)), perm = c(2,1,3))
  #matrixG <- aperm(array(parameters$G, dim = c(numSequences, numWeeks, numRegions)), perm = c(2,3,1))
  
  # Large one
  ##termXis1 <- parameters$p*exp((y*matrixB) - exp(matrixG + matrixS + matrixB))
  ##termXis0 <- (1 - parameters$p)*exp(- exp(matrixG + matrixS))
  ##termXis1[is.infinite(termXis1) & !is.infinite(termXis0)] <- 1
  ##termXis0[is.infinite(termXis1) & !is.infinite(termXis0)] <- 0
  ##termXis1[!is.infinite(termXis1) & is.infinite(termXis0)] <- 0
  ##termXis0[!is.infinite(termXis1) & is.infinite(termXis0)] <- 1
  ##termXis1[is.infinite(termXis1) & is.infinite(termXis0)] <- 0.5
  ##termXis0[is.infinite(termXis1) & is.infinite(termXis0)] <- 0.5
  ##probaXis1 <- termXis1/(termXis1 + termXis0)
  ##uX <- matrix(runif(numSequences*numWeeks), nrow = numWeeks, ncol = numSequences)
  ##proposalX <- 1*(uX < probaXis1)
  # Avoid infinite problems 28.03.2019 WRONG:
  ##Q0_1 <- (1-parameters$p)*exp(- (y*matrixB) + exp(matrixG + matrixS + matrixB) - exp(matrixG + matrixS))/parameters$p
  ##logProbaXis1 <- - log(1 + Q0_1)
  ##uX <- matrix(runif(numSequences*numWeeks), nrow = numWeeks, ncol = numSequences)
  ##proposalX <- 1*(log(uX) < logProbaXis1)
  
  # Aproximation using notes on 30.03.2019NO 25.06.2019 (better using rbinom thatn runif? TODO):
  # (Note that auxTerm is the collapsed sum of a cube (auxTermCube))
  # TODONah move this block to a function since it contains likelihoods ?? Nah
  auxTermCube <- - (matrixPop*exp(parameters$a)*matrixGexp*matrixSexp*matrixRexp) - y*matrixB + (matrixPop*exp(parameters$a)*matrixGexp*matrixSexp*matrixRexp*exp(matrixB))
  
  #tempCollapsedByK <- apply(auxTermCube, 1:2, sum)
  # here I'm assuming dimToInclude has 2 dimensions and it's ordered
  tempCollapsedByK <- drop(auxTermCube)
  auxTerm <- outer(1:numBlockDims[dimToInclude[1]], 1:numBlockDims[dimToInclude[2]], FUN = auxFnUpdateX)
  
  if(0){
    # Approximation avoiding infinites
    # TODO check if bug
    expAuxTerm <- exp(auxTerm)
    expAuxTermIsInf <- is.infinite(expAuxTerm)
    expAuxTermNoInf <- expAuxTerm
    expAuxTermNoInf[expAuxTermIsInf] <- 0
    logProbaXis1 <- expAuxTermIsInf*(-log(1 - parameters$p) + log(parameters$p) - auxTerm) + # log(A/(A+B)) ~ -log(B/A) ???
      (!expAuxTermIsInf)*(-log(1 + (1-parameters$p)*expAuxTermNoInf/parameters$p)) # log(A/(A+B)) = -log(1 + B/A)
    uX <- matrix(runif(numBlockDims[dimToInclude[1]]*numBlockDims[dimToInclude[2]]), nrow = numBlockDims[dimToInclude[1]], ncol = numBlockDims[dimToInclude[2]])
    proposalX <- 1*(log(uX) < logProbaXis1)
  }else{
    expAuxTerm <- exp(auxTerm)
    probaXis1 <- parameters$p/(parameters$p + (1-parameters$p)*expAuxTerm)
    cat(mean(probaXis1, na.rm = T), " ")
    uX <- matrix(runif(numBlockDims[dimToInclude[1]]*numBlockDims[dimToInclude[2]]), nrow = numBlockDims[dimToInclude[1]], ncol = numBlockDims[dimToInclude[2]])
    proposalX <- 1*(uX < probaXis1)
    cat(sum(proposalX, na.rm = T), " ")
  }
  
  parameters$X <- proposalX
  
  matrixX <- array(parameters$X[allToGroups], dim = numDims)
  matrixXBexp <- exp(matrixX*matrixB)
  cat("X updated ")
}