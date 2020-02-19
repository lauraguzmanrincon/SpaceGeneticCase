# Requires: config, sigmaJumps, storage, parameters
# WARNING Assumes dimBeta is only 2 or 3 or 123

if(dimBeta != 123){
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
}

# This code changed considerably from previous MCMC versions
# To be used if B is a matrix DONE # 07.01.2020 :D
if(dimBeta == 123){
  if(config$ifBUpdate){
    for(i in 1:numWeeksGroups){
      for(j in 1:numRegionsGroups){
        for (k in 1:numSequenceGroups) {
          #if(sum(parameters$X[i,j,k]) > 0){
          repeat{
            proposJump <- rnorm(1, parameters$B[i,j,k], sd = sigmaJumps$B[i,j,k])
            if(proposJump > 0){
              proposalB <- proposJump
              break
            }
          }
          larB <- lpriorBijk(i, j, k, proposalB) + llBijk(i, j, k, proposalB) -
            lpriorBijk(i, j, k, parameters$B[i,j,k]) - llBijk(i, j, k, parameters$B[i,j, k]) + log(pnorm(parameters$B[i,j,k])) - log(pnorm(proposalB))
          uB <- runif(1)
          if (log(uB) < larB) {
            parameters$B[i,j,k] <- proposalB
            storage$accept$B[i,j,k] <- storage$accept$B[i,j,k] + 1
            sigmaJumps$B[i,j,k] <- sigmaJumps$B[i,j,k]*x
            matrixB[iToGroups == i,jToGroups == j,kToGroups == k] <- proposalB
            matrixXBexp[iToGroups == i,jToGroups == j,kToGroups == k] <- exp(matrixX[iToGroups == i,jToGroups == j,kToGroups == k]*proposalB)
          } else {
            storage$reject$B[i,j,k] <- storage$reject$B[i,j,k] + 1
            sigmaJumps$B[i,j,k] <- sigmaJumps$B[i,j,k]*x^(-0.7857)
          }
          #}
        }
      }
    }
    cat("B updated ")
  }
}