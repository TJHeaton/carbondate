################################################
# Function which updates the ith element of c given current theta and phi
# Arguments:
# i - the element to update in c
# c - the vector of cluster identifiers
# phi - means of current clusters
# tau - precision of current clusters
# theta - the given calendar ages where in a specific cluster
# theta ~ N(phi, sigma^2)
# lambda, nu1, nu2 - parameters in NormalGamma prior on (phi, tau)
# mu_phi - current estimate of mu_phi
# alpha - parameter in stick breaking
# Returns:
# c - updated vector of cluster identifiers
# phi - updated vector of current cluster means
# tau - updated vector of current cluster precisions
# Note: phi and tau will need to stored as a list on return as variable length

.BivarUpdateClusterIdentifier <- function(
    i, c, phi, tau, theta, lambda, nu1, nu2, mu_phi, alpha, use_cpp) {
  nc <- length(phi)
  ci <- c[i] # Cluster of element to update
  cminus <- c[-i] # the other cluster elements
  if (sum(cminus == ci) == 0) {
    phi <- phi[-ci] # Remove phi for ci if no other elements
    tau <- tau[-ci] # Remove tau for ci
    cminus[cminus > ci] <- cminus[cminus > ci] - 1 # Adjust labelling
    nc <- nc - 1 # Adjust n levels
  }
  nci <- .NumberOfObservationsInEachCluster(as.integer(cminus))

  # Likelihood of theta given phi and tau
  cprob <- stats::dnorm(theta, mean = phi, sd = 1 / sqrt(tau))
  logmarg <- .LogMarginalNormalGamma(theta, mu_phi, lambda, nu1, nu2)
  # Concatenate marginal of theta for a new cluster
  cprob <- c(cprob, exp(logmarg))

  # weight by number in class (or alpha for new cluster)
  cprob <- cprob * c(nci, alpha)
  class <- sample.int(nc + 1, 1, prob = cprob)
  if (class == (nc + 1)) {
    # We have sampled a new state
    # - create new phi and tau from posterior given theta
    nu1new <- nu1 + 0.5
    nu2new <- nu2 + (lambda * (theta - mu_phi)^2) / (2 * (lambda + 1))
    lambdanew <- lambda + 1
    muphinew <- (lambda * mu_phi + theta) / (lambda + 1)
    taunew <- stats::rgamma(1, shape = nu1new, rate = nu2new)
    phinew <- stats::rnorm(
      1, mean = muphinew, sd = 1 / sqrt(lambdanew * taunew))
    phi <- c(phi, phinew)
    tau <- c(tau, taunew)
  }
  # Now update the return class variables
  c[-i] <- cminus
  c[i] <- class
  retlist <- list(c = c, phi = phi, tau = tau)
  return(retlist)
}


# Function as above but using a gamma prior on the value of alpha
.UpdateAlphaGammaPrior <- function(n, alpha, prshape = 0.5, prrate = 1, nci) {
  propsd = 1
  # Sample new alpha from truncated normal distribution
  repeat {
    alphanew <- stats::rnorm(1, alpha, propsd)
    if (alphanew > 0) {
      break
    }
  }
  logprrat <- stats::dgamma(
      alphanew, shape = prshape, rate = prrate, log = TRUE) -
    stats::dgamma(alpha, shape = prshape, rate = prrate, log = TRUE)
  loglikrat <- .AlphaLogLiklihood(n, alphanew, nci) - .AlphaLogLiklihood(n, alpha, nci)
  # Adjust for non-symmetric truncated normal proposal
  logproprat <- stats::pnorm(alpha / propsd, log.p = TRUE) -
    stats::pnorm(alphanew / propsd, log.p = TRUE)
  HR <- exp(logprrat + loglikrat + logproprat)
  if (stats::runif(1) < HR) {
    return(alphanew)
  }
  return(alpha)
}

# This function will work out the likelihood of a particular Dir(alpha)
# given a partition c
# Arguments are:
# c - vector of the classes of each observation
# alpha - parameter in Dir(alpha)
# Return the likelihood
.AlphaLogLiklihood <- function(n, alpha, nci) {
  nc <- length(nci)
  # Note we have to use pmax(nci-1, 1) here to account for clusters of
  # size 1 have 0! = 1
  loglik <- nc*log(alpha) + sum(
    sapply(
      pmax(nci-1, 1),
      function(x) sum(log(1:x)))) - sum(log(alpha:(alpha+n-1)))
  return(loglik)
}


.NumberOfObservationsInEachCluster <- function(cluster_identifiers) {
  n_clust <- max(cluster_identifiers)
  nci <- apply(
    t(as.matrix(1:n_clust)), # A row vector
    2,
    function(x, c) sum(c == x),
    c = cluster_identifiers)
}
