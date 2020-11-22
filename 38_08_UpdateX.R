# Requires: config, sigmaJumps, storage, parameters
# Independent prior notes:  I only modified the auxTerm. Check notes on 30.03.2019 to recall
#                           NOO Also, note that X is now an array, and so the auxiliar matrices

# Update X block (at the same time since conditionally independent) ----
if(config$ifXUpdate & it > 20){
  if(autocorrPrior == 0){
    # Independent prior
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
  }else{
    # Autocorrelated prior (13.04.2020) - Check notes 13.04.2020
    # Here I assume dimToInclude = c(1,#). That is, the first dim to include is 1 (tested in 38_00.R).
    
    proposalX <- array(0, numBlockDims[dimToInclude])
    uX <- matrix(runif(numBlockDims[dimToInclude[1]]*numBlockDims[dimToInclude[2]]), nrow = numBlockDims[dimToInclude[1]], ncol = numBlockDims[dimToInclude[2]])
    
    # Stationary distribution
    I0 <- parameters$p10/(parameters$p10 + parameters$p01)
    I1 <- 1- I0
    
    # Likelihood matrices (idea as in 'Intependent prior'*)
    # (*) "here I'm assuming dimToInclude has 2 dimensions and it's ordered"
    # Seems to cost ~2min per 1000 iterataions
    auxTermCube <- y*matrixB - (matrixPop*exp(parameters$a)*matrixGexp*matrixSexp*matrixRexp*exp(matrixB))
    tempCollapsedByK <- drop(auxTermCube)
    logL1 <- outer(1:numBlockDims[dimToInclude[1]], 1:numBlockDims[dimToInclude[2]], FUN = auxFnUpdateX)
    auxTermCube <- -matrixPop*exp(parameters$a)*matrixGexp*matrixSexp*matrixRexp
    tempCollapsedByK <- drop(auxTermCube)
    logL0 <- outer(1:numBlockDims[dimToInclude[1]], 1:numBlockDims[dimToInclude[2]], FUN = auxFnUpdateX)
    
    dimsSTSet <- c(numBlockDims[dimToInclude[1]]/numSTGroups, numBlockDims[dimToInclude[2]])
    numWeeksGroupsSTSet <- dimsSTSet[1]
    for(stIndex in 1:numSTGroups){ # loop added to adapt FFBS to ST cases (only change)
      # Storage
      X0 <- array(0, dimsSTSet)
      X1 <- array(0, dimsSTSet)
      W0 <- array(0, dimsSTSet)
      W1 <- array(0, dimsSTSet)
      #logW0 <- array(0, numBlockDims[dimToInclude]) # for approximate(?) version
      #logW1 <- array(0, numBlockDims[dimToInclude]) # for approximate(?) version
      
      # Forward Filtering
      # 1.1
      X0[1,] <- rep(I0, dimsSTSet[2])
      X1[1,] <- rep(I1, dimsSTSet[2])
      
      for(tFF in 1:numWeeksGroupsSTSet){
        # Complete version
        if(0){
          # t.1
          if(tFF != 1){
            X0[tFF,] <- (1 - parameters$p01)*W0[tFF - 1,] + parameters$p10*W1[tFF - 1,]
            X1[tFF,] <- parameters$p01*W0[tFF - 1,] + (1 - parameters$p10)*W1[tFF - 1,]
          }
          # t.2
          Y0t <- L0[tFF,]*X0[tFF,]
          Y1t <- L1[tFF,]*X1[tFF,]
          YYt <- Y0t + Y1t
          # t.3
          W0[tFF,] <- Y0t/YYt
          W1[tFF,] <- Y1t/YYt
        }
        # Approximate(?) version
        else{
          # t.1
          if(tFF != 1){
            X0[tFF,] <- (1 - parameters$p01)*W0[tFF - 1,] + parameters$p10*W1[tFF - 1,]
            X1[tFF,] <- parameters$p01*W0[tFF - 1,] + (1 - parameters$p10)*W1[tFF - 1,]
          }
          # t.2
          logY0t <- logL0[tFF,] + log(X0[tFF,])
          logY1t <- logL1[tFF,] + log(X1[tFF,])
          # t.3
          W0[tFF,] <- exp(logY0t) # here W now contains Y
          W1[tFF,] <- exp(logY1t)
        }
      }
      
      # Backward Sampling
      # T.4
      #proposalX[stIndex*numWeeksGroupsSTSet,] <- 1*(uX[numWeeksGroupsSTSet,] < W1[numWeeksGroupsSTSet,]) # Complete version
      proposalX[stIndex*numWeeksGroupsSTSet,] <- 1*(uX[numWeeksGroupsSTSet,] < W1[numWeeksGroupsSTSet,]/(W0[numWeeksGroupsSTSet,] + W1[numWeeksGroupsSTSet,])) #Approximate(?) vers.
      
      for(tBS in (numWeeksGroupsSTSet - 1):1){
        # t.4
        tBSinProposal <- (stIndex - 1)*numWeeksGroupsSTSet + tBS
        Z0t <- W0[tBS,]*(((1 - parameters$p01)/X0[tBS + 1,])*(1 - proposalX[tBSinProposal + 1,]) + (parameters$p01/X1[tBS + 1,])*proposalX[tBSinProposal + 1,])
        Z1t <- W1[tBS,]*((parameters$p10/X0[tBS + 1,])*(1 - proposalX[tBSinProposal + 1,]) + ((1 - parameters$p10)/X1[tBS + 1,])*proposalX[tBSinProposal + 1,])
        
        # Draw sample
        proposalX[tBSinProposal,] <- 1*(uX[tBSinProposal,] < Z1t/(Z0t + Z1t))
      }
    }
    
    # To fix ST_TW_20_#
    proposalX[is.na(proposalX)] <- 0
    
    # Store
    # TODO Note initial state of X is 3D while parameters$X is 2D
    # Note that a dropped parameters$X is a consistent format
    parameters$X <- proposalX
    matrixX <- array(parameters$X[allToGroups], dim = numDims)
    matrixXBexp <- exp(matrixX*matrixB)
    cat("X autocorr updated ")
  }
}