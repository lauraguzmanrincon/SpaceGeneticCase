# Plot distributions ----
plotGamma <- function(shape, rate, xlab = "", returnPlot = 1){
  val <- seq(0, qgamma(0.95, shape = shape, rate = rate), length.out = 51)
  if(returnPlot == 0){
    return(list(x = val, y = dgamma(val, shape = shape, rate = rate)))
  }else{
    pp <- ggplot(data.table(x = val, y = dgamma(val, shape = shape, rate = rate)), aes(x = x, y = y)) + geom_path() + theme_laura() +
      labs(x = xlab, y = "probability distribution")
    return(pp)
  }
}

plotNorm <- function(mean, sd, xlab = "", returnPlot = 1){
  val <- seq(qnorm(0.05, mean = mean, sd = sd), qnorm(0.95, mean = mean, sd = sd), length.out = 51)
  if(returnPlot == 0){
    return(list(x = val, y = dnorm(val, mean = mean, sd = sd)))
  }else{
    pp <- ggplot(data.table(x = val, y = dnorm(val, mean = mean, sd = sd)), aes(x = x, y = y)) + geom_path() + theme_laura() +
      labs(x = xlab, y = "probability distribution")
    return(pp)
  }
}

plotBeta <- function(shape1, shape2, xlab = "", returnPlot = 1){
  val <- seq(0, 1, length.out = 51)
  if(returnPlot == 0){
    return(list(x = val, y = dbeta(val, shape1 = shape1, shape2 = shape2)))
  }else{
    pp <- ggplot(data.table(x = val, y = dbeta(val, shape1 = shape1, shape2 = shape2)), aes(x = x, y = y)) + geom_path() + theme_laura() +
      labs(x = xlab, y = "probability distribution")
    return(pp)
  }
}

plotLogNorm <- function(mean, sd, xlab = "", returnPlot = 1){
  val <- seq(qnorm(0.05, mean = mean, sd = sd), qnorm(0.95, mean = mean, sd = sd), length.out = 51)
  if(returnPlot == 0){
    return(list(x = val, y = dnorm(val, mean = mean, sd = sd)))
  }else{
    pp <- ggplot(data.table(x = exp(val), y = dnorm(val, mean = mean, sd = sd)), aes(x = x, y = y)) + geom_path() + theme_laura() +
      labs(x = xlab, y = "probability distribution")
    return(pp)
  }
}

plotGammaTruncated <- function(shape, rate, upp, xlab = ""){
  val <- seq(0, upp, length.out = 51)
  normalConst <- pgamma(upp, shape = shape, rate = rate)
  pp <- ggplot(data.table(x = val, y = dgamma(val, shape = shape, rate = rate)/normalConst), aes(x = x, y = y)) + geom_path() + theme_laura() +
    labs(x = xlab, y = "probability distribution")
  return(pp)
}

#ttt <- matrix(1:6, ncol = 2)
#ttt1 <- c(1,2,2) # matrix(c(1,2,2,1,2,2), ncol = 2)
#ttt2 <- c(1,2) # matrix(c(1,1,1,2,2,2), ncol = 2)
#holahola <- Vectorize(function(x, y) sum(ttt[ttt1 == x, ttt2 == y]))
#outer(1:2,1:2, FUN = holahola)
# Auxiliar function only for the X update, Requires the existence of tempCollapsedByK
#OLDauxFnUpdateX <- Vectorize(function(iInd, jInd, zInd) sum(tempCollapsedByK[iToGroups == iInd, jToGroups == jInd, kToGroups == kInd]))
auxFnUpdateX <- Vectorize(function(iInd, jInd) sum(tempCollapsedByK[listDims[[dimToInclude[1]]] == iInd, listDims[[dimToInclude[2]]] == jInd]))

collapseAuxFn <- Vectorize(function(iInd, jInd) sum(matToCollapse[listDims[[dimToInclude[1]]] == iInd, listDims[[dimToInclude[2]]] == jInd]))
collapseMatrixFn <- function() outer(1:numBlockDims[dimToInclude[1]], 1:numBlockDims[dimToInclude[2]], FUN = collapseAuxFn)

# 24.03.2020
# Taken from epiclustR TEMPORARLY
# TODO
#' block_precision_order_2
#' 
#' Computes the precision matrix for a random walk of order 2
#' of the given block size
#'
#' @param n the block size
#' @return The precision matrix K corresponding to a random walk of order 2.
block_precision_order_2 <- function(n) {
  # Generates the structure matrix for a random walk of order 2
  K<-6*diag(2 + n + 2) # need 2 on either side for ends
  for (j in 1:nrow(K)) {
    if (j>2) {K[j,j-2]<-1}
    if (j>1) {K[j,j-1]<--4}
    if (j+1<=ncol(K)) {K[j,j+1]<--4}
    if (j+2<=ncol(K)) {K[j,j+2]<-1}
  }
  K[1,1:2] <- K[n+4,n+4:3] <- c(1, -2)
  K[2,1:2] <- K[n+3,n+4:3] <- c(-2, 5)
  K
}

