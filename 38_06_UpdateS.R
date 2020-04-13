# Requires: config, sigmaJumps, storage, parameters

# Choose update method
adaptTypeS <- ifelse(it%%2, 3, 1) # do not update 2 and 3 at the same time (storage issues to be fixed)
#adaptTypeS <- ifelse(floor(it/5)%%2, 3, 1)
parameters$sBlockSize <- 0 # only for adaptTypeS = 3

if(config$ifSUpdate){
  if(adaptTypeS == 1){
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
  }else if(adaptTypeS == 2){
    # II. (conditional proposal)
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
  }else if(adaptTypeS == 3){
    # II. (conditional block proposal)
    # UPDATE 30.03.2020: update all ST groups independently (**)
    randomSize <- sample(x = 1:11, size = 1, replace = TRUE)
    #randomStart <- sample(x = randomSize, size = 1, replace = TRUE)
    randomStart <- 1
    parameters$sBlockSize <- c(randomSize, randomStart)
    itNumBlocksS <- (randomStart != 1) + ceiling((numWeeks/numSTGroups - randomStart + 1)/randomSize) # **
    
    for(scm in 1:numSTGroups){ # **
      scmIndexes <- (scm - 1)*(numWeeks/numSTGroups) + 1:(numWeeks/numSTGroups) # **
      for(sBlockIndex in 1:itNumBlocksS){
        #if(sBlockIndex == 1 & randomStart != 1){
        #  sBlockIndexes <- 1:(randomStart - 1)
        #}else{
        #sBlockIndexes <- intersect(randomStart - 1 + (sBlockIndex - 1 - (randomStart != 1))*randomSize + (1:randomSize), 1:numWeeks)
        #}
        sBlockIndexes <- intersect((scm - 1)*(numWeeks/numSTGroups) + randomStart - 1 + (sBlockIndex - 1 - (randomStart != 1))*randomSize + (1:randomSize), 1:numWeeks) # **
        sBlockComplement <- setdiff(scmIndexes, sBlockIndexes) # **
        
        # Fahrmeir: construct mu and sigma for proposal
        muA <- solve(seasonalCoefficientMatrixSqr[sBlockIndexes, sBlockIndexes])
        muB <- seasonalCoefficientMatrixSqr[sBlockIndexes, sBlockComplement]
        muC <- parameters$S[sBlockComplement]
        
        muSample <- -muA%*%muB%*%muC
        sigmaSample <- muA/parameters$tau.S
        
        # Create proposal
        proposalSKnorr <- mvrnorm(n = 1, mu = muSample, Sigma = sigmaSample)
        larSKnorr <- llSKnorr(sBlockIndexes, proposalSKnorr) - llSKnorr(sBlockIndexes, parameters$S[sBlockIndexes])
        
        u <- runif(1)
        if(log(u) < larSKnorr){
          parameters$S[sBlockIndexes] <- proposalSKnorr
          storage$accept$ScondUp[sBlockIndexes] <- storage$accept$ScondUp[sBlockIndexes] + 1
          matrixS[sBlockIndexes,,] <- proposalSKnorr
          matrixSexp[sBlockIndexes,,] <- exp(proposalSKnorr)
          #cat("Accepted:")
        }else{
          storage$reject$ScondUp[sBlockIndexes] <- storage$reject$ScondUp[sBlockIndexes] + 1
          #cat("Rejected:")
        }
      }
    }
    
    cat("S cond. block updated ")
  }
}
