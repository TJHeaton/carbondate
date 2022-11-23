# Function which performs one iteration of Walker DP update on (w, phi, tau)
# Input arguments:
# theta- the current observations theta_i ~ N(phi_delta[i], prec = tau_delta[i])
# w - the current weights
# v - the current vs
# delta - the current cluster allocations
# n_clust - current n_clust
# alpha, mu_phi, lambda, nu1, nu2 - the current DP parameters (not updated here)
.DPWalkerUpdate_cpp <- function(
    theta, w, v, delta, n_clust, alpha, mu_phi, lambda, nu1, nu2) {

  retlist <- WalkerUpdateWeights_cpp(w, v, delta, n_clust, alpha)

  u = retlist$u
  v = retlist$v
  w = retlist$weight
  n_clust = retlist$n_clust

  retlist <- WalkerUpdateClusterPhiTau_cpp(n_clust, as.double(theta), delta, mu_phi, lambda, nu1, nu2)

  phi <- retlist$phi;
  tau <- retlist$tau;

  delta = WalkerUpdateClusterIdentifiers_cpp(as.double(theta), u, w, phi, tau)

  # Return w, v, delta, phi and tau
  return_list <- list(
    w = w, v = v, delta = delta, phi = phi, tau = tau, n_clust = n_clust)
  return(return_list)
}


# Function which updates the Dirichlet process parameter alpha
# Does this via Metropolis-Hastings with truncated proposal distribution
# Arguments are:
# delta - vector of the classes of each observation
# alpha - current alpha parameter in Dir(alpha)
# prshape and prrate - prior on alpha ~ Gamma(shape, rate)
# propsd - proposal sd for new alpha
.WalkerUpdateAlpha <- function(
    delta, alpha, prshape = 0.5, prrate = 1, propsd = 1) {
  # Sample new alpha from truncated normal distribution
  repeat {
    alpha_new <- stats::rnorm(1, alpha, propsd)
    if (alpha_new > 0) {
      break
    }
  }
  # Find distinct number of populated clusters
  d <- length(unique(delta))
  n <- length(delta)
  # Now find the likelihood and prior
  logprrat <- stats::dgamma(alpha_new, shape = prshape, rate = prrate, log = TRUE) -
    stats::dgamma(alpha, shape = prshape, rate = prrate, log = TRUE)
  loglikrat <- .WalkerAlphaLogLiklihood(d = d, alpha = alpha_new, n = n) -
    .WalkerAlphaLogLiklihood(d = d, alpha = alpha, n = n)
  # Adjust for non-symmetric truncated normal proposal
  logproprat <- stats::pnorm(alpha / propsd, log.p = TRUE) -
    stats::pnorm(alpha_new / propsd, log.p = TRUE)
  HR <- exp(logprrat + loglikrat + logproprat)
  if (stats::runif(1) < HR) {
    return(alpha_new) # Accept alpha_new
  }
  return(alpha) # Or reject and keep alpha
}


.SampleNewPhi <- function(mu, sigma, range = c(0, pkg.globals$MAX_YEAR_BP)) {
  inrange <- FALSE
  while (!inrange) {
    newphi <- stats::rnorm(1, mu, sigma)
    inrange <- (newphi < max(range)) & (newphi > min(range))
  }
  return(newphi)
}


.WalkerAlphaLogLiklihood <- function(d, alpha, n) {
  return(d * log(alpha) + lgamma(alpha) - lgamma(alpha + n))
}
