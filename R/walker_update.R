# Function which performs one iteration of Walker DP update on (w, phi, tau)
# Input arguments:
# theta- the current observations theta_i ~ N(phi_delta[i], prec = tau_delta[i])
# w - the current weights
# v - the current vs
# delta - the current cluster allocations
# phi - the current cluster means
# tau - the current cluster precisions
# n_clust - current n_clust
# c, mu_phi, lambda, nu1, nu2 - the current DP parameters (not updated here)
.DPWalkerUpdate <- function(
    theta, w, v, delta, phi, tau, n_clust, c, mu_phi, lambda, nu1, nu2) {
  # Create relevant variables
  n <- length(theta)

  # Create auxiliary u variables
  u <- stats::runif(n, min = 0, max = w[delta])
  ustar <- min(u)


  # Now update the weights
  wnew <- c()
  j <- 1
  brprod <- 1 # This is the current product of (1-v[1])...(1-v[j-1])
  # Need to work out alpha and beta internally

  # Iteratively update the weights until we have all we need
  # To update v[j] we do the following
  while (sum(wnew) < 1 - ustar) {

    if (j <= n_clust) {
      # We have to work out alpha and beta and sample a new v[j] by inverse cdf

      # Find the indices we need to search over i.e. which are in cluster j and
      # which are above cluster j
      clustj <- which(delta == j)
      aboveclj <- which(delta > j)

      # Find alpha_j
      if (length(clustj) == 0) {
        alpha <- 0
      } else {
        alpha <- max(u[clustj] / brprod)
      }

      # Find beta_j
      if (length(aboveclj) == 0) {
        beta <- 1
      } else {
        # Vector to aid denominator: ith entry is Prod_{l < i, l != j} (1- v_l)
        prodvtemp <- cumprod(1 - v) / (1 - v[j])
        beta <- 1 -
          max(
            u[aboveclj] / (v[delta[aboveclj]] * prodvtemp[delta[aboveclj] - 1]))
      }

      # Now sample from the correct posterior by inversion
      A <- (1 - alpha)^c
      B <- A - (1 - beta)^c
      utemp <- stats::runif(1)
      v[j] <- 1 - (A - B * utemp)^(1 / c)
    } else {
      # We are in the case that j > k^star (for current n_clust) and v[j] is just
      # from the prior beta
      v[j] <- stats::rbeta(1, 1, c)
    }

    # We now have v so we need to extend w (and update brprod for next iteration)
    wnew <- c(wnew, brprod * v[j])
    brprod <- brprod * (1 - v[j])

    if (sum(is.na(wnew)) != 0) {
      cat("c = ", c, "\n")
      cat("j < k^star is", j <= n_clust, "\n")
      cat("v[j] = ", v[j], "\n")
      browser()
      stop("Weird weights")
    }
    j <- j + 1
  }

  # Now update n_clust (the number of weights we have) and truncate w and v at the
  # correct values
  n_clust <- length(wnew)
  w <- wnew
  v <- v[1:n_clust]

  # Now update the cluster means and precisions (note we have to introduce new
  # ones for the new states without observations)
  for (i in 1:n_clust) {
    # Find which observations belong to this cluster
    clusti <- which(delta == i)
    if (length(clusti) == 0) {
      # No observations in this cluster so sample from the prior
      tau[i] <- stats::rgamma(1, shape = nu1, rate = nu2)
      phi[i] <- stats::rnorm(1, mean = mu_phi, sd = 1 / sqrt(lambda * tau[i]))
    } else {
      # There are some observations we need to update the phi and tau for this
      # cluster (conjugate according to NormalGamma prior)
      gibbs_parameters <- .UpdatePhiTau(
        theta[clusti], mu_phi = mu_phi, lambda = lambda, nu1 = nu1, nu2 = nu2)
      phi[i] <- gibbs_parameters$phi
      tau[i] <- gibbs_parameters$tau
    }
  }
  # Only store the values we need
  phi <- phi[1:n_clust]
  tau <- tau[1:n_clust]

  # Now update the allocations for each observation by sampling from them
  for (i in 1:n) {
    possid <- which(w > u[i])
    dens <- stats::dnorm(
      phi[possid], mean = theta[i], sd = 1 / sqrt(tau[possid]))
    dens[is.na(dens)] <- 0 # Fudge to remove erroneous problems
    delta[i] <- possid[sample.int(length(possid), 1, prob = dens)]
  }

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
