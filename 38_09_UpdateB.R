# Requires: config, sigmaJumps, storage, parameters
# WARNING Assumes dimBeta is only 2 or 3

if(config$ifBUpdate){
  for(l in 1:numBeta){
  #for(l in which.max(colSums(parameters$X[1,,]))){
    #if(sum(parameters$X[,l]) > 0){
      repeat{
        proposJump <- rnorm(1, parameters$B[l], sd = sigmaJumps$B[l])
        if(proposJump > 0){
          proposalB <- proposJump
          break
        }
      }
      larB <- lpriorB(l, proposalB) + llB(l, proposalB) -
        lpriorB(l, parameters$B[l]) - llB(l, parameters$B[l]) + log(pnorm(parameters$B[l])) - log(pnorm(proposalB))
      uB <- runif(1)
      if (log(uB) < larB) {
        parameters$B[l] <- proposalB
        storage$accept$B[l] <- storage$accept$B[l] + 1
        sigmaJumps$B[l] <- sigmaJumps$B[l]*x
        if(dimBeta == 2){
          matrixB[,jToGroups == l,] <- proposalB
          matrixXBexp[,jToGroups == l,] <- exp(matrixX[,jToGroups == l,]*proposalB)
        }else{
          matrixB[,,kToGroups == l] <- proposalB
          matrixXBexp[,,kToGroups == l] <- exp(matrixX[,,kToGroups == l]*proposalB)
        }
      } else {
        storage$reject$B[l] <- storage$reject$B[l] + 1
        sigmaJumps$B[l] <- sigmaJumps$B[l]*x^(-0.7857)
      }
    #}
  }
  cat("B updated ")
}

# This code changed considerably from previous MCMC versions
# To be used if B is a matrix
if(0){
  if(config$ifBUpdate){
    for(j in 1:numRegionsGroups){
      for (k in 1:numSequencesGroups) {
        if(sum(parameters$X[,j,k]) > 0){
          repeat{
            proposJump <- rnorm(1, parameters$B[j,k], sd = sigmaJumps$B[j,k])
            if(proposJump > 0){
              proposalB <- proposJump
              break
            }
          }
          larB <- lpriorB(j, k, proposalB) + llB(j, k, proposalB) -
            lpriorB(j, k, parameters$B[j,k]) - llB(j, k, parameters$B[j, k]) + log(pnorm(parameters$B[j,k])) - log(pnorm(proposalB))
          uB <- runif(1)
          if (log(uB) < larB) {
            parameters$B[j,k] <- proposalB
            storage$accept$B[j,k] <- storage$accept$B[j,k] + 1
            sigmaJumps$B[j,k] <- sigmaJumps$B[j,k]*x
            #matrixB[,j,k] <- proposalB
            matrixXBexp[,j,k] <- exp(matrixX[,j,k]*proposalBkansbdkfa)
          } else {
            storage$reject$B[j,k] <- storage$reject$B[j,k] + 1
            sigmaJumps$B[j,k] <- sigmaJumps$B[j,k]*x^(-0.7857)
          }
        }
      }
    }
    cat("B updated ")
  }
}