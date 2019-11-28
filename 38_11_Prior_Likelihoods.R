# Priors and likelihood of parameters and hyperparameters
# Parameters: a, Gk, Si, Rj, Xij, Bl
# Hyperparameters: Tr, Ts, Tk, l, p
# These functions compute the priors and likelihoods of a, Gk, Rj, Si, Bl

# Parameter a ----
lpriorA <- function(a){
  value <- -(a - constants$ma)^2/(2*constants$va)
  return(value)
}
llA <- function(a){
  value <- sum(y*a - matrixPop*exp(a)*matrixGexp*matrixSexp*matrixRexp*matrixXBexp)
  return(value)
}

# Parameter Gk ----
lpriorG <- function(k, Gk){
  #value <- - parameters$tau.G*sum(Gk*parameters$G*invDeltaMatrix[k,]) + invDeltaMatrix[k,k]*(parameters$G[k]*Gk - 0.5*Gk^2) # correction 02042019 WRONG
  value <- - parameters$tau.G*sum(Gk*parameters$G*invDeltaMatrix[k,]) + parameters$tau.G*invDeltaMatrix[k,k]*(parameters$G[k]*Gk - 0.5*Gk^2) #correction 16072019
  return(value)
}
llG <- function(k, Gk){
  value <- sum(y[,,k]*Gk - matrixPop[,,k]*exp(parameters$a)*exp(Gk)*matrixSexp[,,k]*matrixRexp[,,k]*matrixXBexp[,,k])
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
  value <- sum(y[i,,]*Si - matrixPop[i,,]*exp(parameters$a)*matrixGexp[i,,]*exp(Si)*matrixRexp[i,,]*matrixXBexp[i,,])
  return(value)
}

# Parameter Bsigmaj ----
lpriorB <- function(j, Bl){
  value <- (constants$aB + 1)*log(Bl) - constants$bB*Bl
  return(value)
}
llB <- function(l, Bl){
  if(dimBeta == 2){
    value <- sum(y[,jToGroups == l,]*matrixX[,jToGroups == l,]*Bl - matrixPop[,jToGroups == l,]*exp(parameters$a)*matrixGexp[,jToGroups == l,]*
                   matrixSexp[,jToGroups == l,]*matrixRexp[,jToGroups == l,]*exp(matrixX[,jToGroups == l,]*Bl))
  }else{
    value <- sum(y[,,kToGroups == l]*matrixX[,,kToGroups == l]*Bl - matrixPop[,,kToGroups == l]*exp(parameters$a)*matrixGexp[,,kToGroups == l]*
                   matrixSexp[,,kToGroups == l]*matrixRexp[,,kToGroups == l]*exp(matrixX[,,kToGroups == l]*Bl))
  }
  return(value)
}

# Parameter Rj ----
lpriorR <- function(j, Rj){
  value <- -0.5*parameters$tau.R*neighLength[j]*(Rj - sum(matrixRR[,j])/neighLength[j])^2
  return(value)
}
llR <- function(j, Rj){
  value <- sum(y[,j,]*Rj - matrixPop[,j,]*exp(parameters$a)*matrixGexp[,j,]*matrixSexp[,j,]*exp(Rj)*matrixXBexp[,j,])
  return(value)
}

# Update block l and tauG ----
#' Prior for an informative l (gamma) and tauG
lpriorLTauGBlock <- function(invMatrix, tau, l){
  #value <- 0.5*numSequences*log(tau) + 0.5*log(det(invMatrix)) - 0.5*tau*(t(parameters$G)%*%invMatrix%*%parameters$G) +
  value <- 0.5*numSequences*log(tau) + 0.5*determinant.spam(invMatrix)$modulus - 0.5*tau*(t(parameters$G)%*%invMatrix%*%parameters$G) +
    (constants$aG - 1)*log(tau) - constants$bG*tau + (constants$aL - 1)*log(l) - constants$bL*l
  return(value)
}

# Build covariance matrix ----
#' maternParam: if 0 the function returns the squared exponential. If not 0, it is the parameter defining the Matern function
#' inputMatrix: if maternParam is 0, it is the squared of distances. If maternParam is not 0, ot is the matrix of distances
deltaMatrixFn <- function(l, inputMatrix, maternParam){
  if(maternParam == 0){
    value <- exp(-inputMatrix/(2*l^2))
  }
  else if(maternParam == 1/2){
    value <- exp(-inputMatrix/l)
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
