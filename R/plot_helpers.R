.FindPredictedDensityNeal <- function(
    x, c, phi, tau, alpha, mu_phi, lambda, nu1, nu2) {
  nci <- .NumberOfObservationsInEachCluster(c)

  # Find probability that new theta[i+1] lies in a particular cluster or a new one
  pci <- c(nci, alpha) # Could form new cluster
  pci <- pci / sum(pci)

  dens <- .MixtureDens(x, w = pci[1:nc], mu = phi, sd = 1 / sqrt(tau)) +
    .PredNewDens(x, w = pci[nc + 1], mu_phi, lambda, nu1, nu2)
}


.FindPredictedDensityWalker <- function(
    x, w, phi, tau, mu_phi, lambda, nu1, nu2) {
  prob_new_clust <- 1 - sum(w)

  dens <- .MixtureDens(x, w = w, mu = phi, sd = 1 / sqrt(tau)) +
    .PredNewDens(x, w = prob_new_clust, mu_phi, lambda, nu1, nu2)
  return(dens)
}


# Find mixture density
# Function which where you pass it a set of means, sds and weights and it
# returns the density of the corresponding mixture of normals
# Arguments:
# x - vector of values at which to evaluate mixture density
# w - vector of weights
# mu - the means
# sd - the sds
.MixtureDens <- function(x, w, mu, sd) {
  DTemp <- mapply(
    function(mu, sig, w, x) w * dnorm(x, mean = mu, sd = sig),
    mu,
    sd,
    w,
    MoreArgs = list(x = x))
  apply(DTemp, 1, sum) # Sum up the various mixtures
}


# The predictive for a new observation is a scaled t-distribution
.PredNewDens <- function(x, w, mu_phi, lambda, nu1, nu2) {
  w * exp(.LogMarginalNormalGamma(x, mu_phi, lambda, nu1, nu2))
}

