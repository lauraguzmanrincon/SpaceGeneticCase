# Priors and likelihood of parameters and hyperparameters
# Parameters: a, Gk, Si, Rj, Xij, Bl
# Hyperparameters: Tr, Ts, Tk, l, p
# These functions compute the priors and likelihoods of a, Gk, Rj, Si, Bl
# 30.03.2020: added , na.rm = T to all likelihoods

# Parameter a ----
lpriorA <- function(a){
  value <- -(a - constants$ma)^2/(2*constants$va)
  return(value)
}
llA <- function(a){
  value <- sum(y*a - matrixPop*exp(a)*matrixGexp*matrixSexp*matrixRexp*matrixXBexp, na.rm = T)
  return(value)
}

# Parameter Gk ----
lpriorG <- function(k, Gk){
  #value <- - parameters$tau.G*sum(Gk*parameters$G*invDeltaMatrix[k,]) + invDeltaMatrix[k,k]*(parameters$G[k]*Gk - 0.5*Gk^2) # correction 02042019 WRONG
  value <- - parameters$tau.G*sum(Gk*parameters$G*invDeltaMatrix[k,]) + parameters$tau.G*invDeltaMatrix[k,k]*(parameters$G[k]*Gk - 0.5*Gk^2) #correction 16072019
  return(value)
}
llG <- function(k, Gk){
  value <- sum(y[,,k]*Gk - matrixPop[,,k]*exp(parameters$a)*exp(Gk)*matrixSexp[,,k]*matrixRexp[,,k]*matrixXBexp[,,k], na.rm = T)
  return(value)
}

# Parameter Si ----
lpriorS <- function(i, Si){
  vectorS <- parameters$S
  vectorS[i] <- Si
  value <- - 0.5*parameters$tau.S*Si*(2*(seasonalCoefficientMatrixSqr[i,]%*%vectorS) - seasonalCoefficientMatrixSqr[i,i]*Si)
  return(value)
}
llS <- function(i, Si){
  value <- sum(y[i,,]*Si - matrixPop[i,,]*exp(parameters$a)*matrixGexp[i,,]*exp(Si)*matrixRexp[i,,]*matrixXBexp[i,,], na.rm = T)
  return(value)
}

# Parameter Bsigmaj ----
lpriorB <- function(j, Bl){
  value <- (constants$aB - 1)*log(Bl) - constants$bB*Bl # 07.01.2020 (constants$aB + 1) was a BUG!
  return(value)
}
llB <- function(l, Bl){
  if(dimBeta == 2){
    value <- sum(y[,jToGroups == l,]*matrixX[,jToGroups == l,]*Bl - matrixPop[,jToGroups == l,]*exp(parameters$a)*matrixGexp[,jToGroups == l,]*
                   matrixSexp[,jToGroups == l,]*matrixRexp[,jToGroups == l,]*exp(matrixX[,jToGroups == l,]*Bl), na.rm = T)
  }else{
    value <- sum(y[,,kToGroups == l]*matrixX[,,kToGroups == l]*Bl - matrixPop[,,kToGroups == l]*exp(parameters$a)*matrixGexp[,,kToGroups == l]*
                   matrixSexp[,,kToGroups == l]*matrixRexp[,,kToGroups == l]*exp(matrixX[,,kToGroups == l]*Bl), na.rm = T)
  }
  return(value)
}

# Parameter Bsigmaijk block ----
lpriorBijk <- function(i, j, k, Bijk){
  value <- (constants$aB - 1)*log(Bijk) - constants$bB*Bijk
  return(value)
}
llBijk <- function(i, j, k, Bijk){
  value <- sum(y[iToGroups == i,jToGroups == j,kToGroups == k]*matrixX[iToGroups == i,jToGroups == j,kToGroups == k]*Bijk - matrixPop[iToGroups == i,jToGroups == j,kToGroups == k]*
                 exp(parameters$a)*matrixGexp[iToGroups == i,jToGroups == j,kToGroups == k]*matrixSexp[iToGroups == i,jToGroups == j,kToGroups == k]*
                 matrixRexp[iToGroups == i,jToGroups == j,kToGroups == k]*exp(matrixX[iToGroups == i,jToGroups == j,kToGroups == k]*Bijk), na.rm = T)
  return(value)
}

# Parameter Rj ----
lpriorR <- function(j, Rj){
  value <- -0.5*parameters$tau.R*neighLength[j]*(Rj - sum(matrixRR[,j])/neighLength[j])^2
  return(value)
}
llR <- function(j, Rj){
  value <- sum(y[,j,]*Rj - matrixPop[,j,]*exp(parameters$a)*matrixGexp[,j,]*matrixSexp[,j,]*exp(Rj)*matrixXBexp[,j,], na.rm = T)
  return(value)
}

# Update block l and tauG ----
#' Prior for an informative l (gamma) and tauG
lpriorLTauGBlock <- function(invMatrix, tau, l){
  #value <- 0.5*numSequences*log(tau) + 0.5*log(det(invMatrix)) - 0.5*tau*(t(parameters$G)%*%invMatrix%*%parameters$G) +
  value <- 0.5*numSequences*log(tau) + 0.5*determinant(invMatrix, logarithm = T)$modulus - 0.5*tau*(t(parameters$G)%*%invMatrix%*%parameters$G) +
    (constants$aG - 1)*log(tau) - constants$bG*tau + (constants$aL - 1)*log(l) - constants$bL*l # 22.01.2020 no determinant.spam (see 43.R)
  return(value)
}

# Build covariance matrix ----
#' maternParam: if 0 the function returns the squared exponential. If not 0, it is the parameter defining the Matern function
#' inputMatrix: if maternParam is 0, it is the squared of distances. If maternParam is not 0, ot is the matrix of distances
#' sqrInputMatrix (19.01.2020) squared of distance matrix
deltaMatrixFn <- function(l, inputMatrix, sqrInputMatrix, maternParam){
  if(maternParam == 0){
    value <- exp(-sqrInputMatrix/(2*l^2))
  }
  else if(maternParam == 1/2){
    value <- exp(-inputMatrix/l)
  }else if(maternParam == 3/2){
    value <- (1 + (sqrt(3)*inputMatrix/l))*exp(-sqrt(3)*inputMatrix/l)
  }else{
    warning("Wrong l choice")
  }
  return(value)
}
#deltaMatrixFn <- function(l, maternParam = 0){
#  if(maternParam == 0){
#    value <- exp(-sqrDistanceMatrix/(2*l^2))
#  }
#  else if(maternParam == 1/2){
#    value <- exp(-distanceMatrixMCMC/l)
#  }
#  return(value)
#}

# 18.12.2019 :)
# Overall likelihood ----
lll <- function(){
  value <- sum(y*(parameters$a + matrixG + matrixS + matrixR + matrixB*matrixX) - matrixPop*exp(parameters$a)*matrixGexp*matrixSexp*matrixRexp*matrixXBexp, na.rm = T)
  return(value)
}

# 26.02.2020 (from sandbox 45.R/2.)
# Parameter Gk for Knorr block update ----
llGKnorr <- function(kVector, GkVector){
  GkMatrix <- drop(aperm(array(GkVector, dim = c(length(GkVector), numWeeks, numRegions)), perm = c(2,3,1)))
  # Note there is a drop we didn't use before... it seems it was not relevant until we use a block of indices???
  value <- sum(y[,,kVector]*GkMatrix - matrixPop[,,kVector]*exp(parameters$a)*exp(GkMatrix)*matrixSexp[,,kVector]*matrixRexp[,,kVector]*matrixXBexp[,,kVector], na.rm = T)
  return(value)
}

# 17.03.2020 (Similar to llGKnorr)
# Parameter Si for Knorr block update ----
llSKnorr <- function(iVector, SiVector){
  SiMatrix <- drop(array(SiVector, dim = c(length(SiVector), numRegions, numSequences)))
  # Note there is a drop we didn't use before... it seems it was not relevant until we use a block of indices???
  value <- sum(y[iVector,,]*SiMatrix - matrixPop[iVector,,]*exp(parameters$a)*exp(SiMatrix)*matrixGexp[iVector,,]*matrixRexp[iVector,,]*matrixXBexp[iVector,,], na.rm = T)
  return(value)
}

# Old function to be removed ----
#' Prior for an uninformative l
lpriorl <- function(invMatrix){
  value <- 0.5*log(det(invMatrix)) - 0.5*parameters$tau.G*(t(parameters$G)%*%invMatrix%*%parameters$G)
  return(value)
}

lpriorLGBlock <- function(Gvector, invMatrix){
  value <- 0.5*log(det(invMatrix)) - 0.5*parameters$tau.G*(t(Gvector)%*%invMatrix%*%Gvector)
  return(value)
}

llGBlock <- function(Gvector, S, B, X){
  # S,B,X matrices
  G <- matrix(Gvector, nrow = numWeeks, ncol = numSequences, byrow = TRUE)
  value <- sum((y*G) - exp(parameters$a + G + S + X*B))
  return(value)
}

#' 02042019 Independent prior for Gk
lpriorGIndependent <- function(k, Gk){
  value <- - 0.5*parameters$tau.G*Gk^2
  return(value)
}
