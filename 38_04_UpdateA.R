# Requires: config, sigmaJumps, storage, parameters

# Update
if(config$ifAUpdate){
  proposalA <- rnorm(1, parameters$a, sd = sigmaJumps$a)
  larA <- lpriorA(proposalA) + llA(proposalA) - lpriorA(parameters$a) - llA(parameters$a)
  u <- runif(1)
  if (log(u) < larA) {
    parameters$a <- proposalA
    storage$accept$a <- storage$accept$a + 1
    sigmaJumps$a <- sigmaJumps$a*x
  } else {
    storage$reject$a <- storage$reject$a + 1
    sigmaJumps$a <- sigmaJumps$a*x^(-0.7857)
  }
  #cat(sprintf("%.*f", 2, larA), " ")
  #cat(sprintf("%.*f", 2, log(u) < larA), " ")
  cat("A updated ")
}
