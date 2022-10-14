################################################
# Function which updates the ith element of c given current theta and phi
# Arguments:
# i - the element to update in c
# c - the vector of cluster identifiers
# phi - means of current clusters
# tau - precision of current clusters
# theta - the given calendar ages where in a spcific cluster theta ~ N(phi, sigma^2)
# lambda, nu1, nu2 - parameters in NormalGamma prior on (phi, tau)
# muphi - current estimate of muphi
# alpha - parameter in stick breaking
# Returns:
# c - updated vector of cluster identifiers
# phi - updated vector of current cluster means
# tau - updated vector of current cluster precisions
# Note: phi and tau will need to stored as a list on return as variable length

.BivarUpdatec <- function(i, c, phi, tau, theta, lambda, nu1, nu2, muphi, alpha) {
  nc <- length(phi)
  ci <- c[i] # Cluster of element to update
  cminus <- c[-i] # the other cluster elements
  if (sum(cminus == ci) == 0) {
    phi <- phi[-ci] # Remove phi for ci if no other elements
    tau <- tau[-ci] # Remove tau for ci
    cminus[cminus > ci] <- cminus[cminus > ci] - 1 # Adjust labelling
    nc <- nc - 1 # Adjust n levels
  }
  nci <- apply(as.row(1:nc), 2, function(x, cminus) sum(cminus == x), cminus = cminus)

  cprob <- stats::dnorm(theta, mean = phi, sd = 1 / sqrt(tau)) # Likelihood of theta given phi and tau
  logmarg <- .Logmargnormgamma(theta, muphi, lambda, nu1, nu2)
  cprob <- c(cprob, exp(logmarg)) # Concatenate marginal of theta for a new cluster

  cprob <- cprob * c(nci, alpha) # weight by number in class (or alpha for new cluster)
  class <- sample(1:(nc + 1), 1, prob = cprob)
  if (class == (nc + 1)) { # We have sampled a new state - create new phi and tau from posterior given theta
    nu1new <- nu1 + 0.5
    nu2new <- nu2 + (lambda * (theta - muphi)^2) / (2 * (lambda + 1))
    lambdanew <- lambda + 1
    muphinew <- (lambda * muphi + theta) / (lambda + 1)
    taunew <- stats::rgamma(1, shape = nu1new, rate = nu2new) # Sample new tau
    phinew <- stats::rnorm(1, mean = muphinew, sd = 1 / sqrt(lambdanew * taunew)) # and then phi
    phi <- c(phi, phinew) # Create updated phi
    tau <- c(tau, taunew) # and tau
  }
  # Now update the return class variables
  c[-i] <- cminus
  c[i] <- class
  retlist <- list(c = c, phi = phi, tau = tau)
  return(retlist)
}

# Function which works out the marginal of theta when theta ~ N(phi, sd = sqrt(1/tau)) and (phi,tau) are NormalGamma
.Logmargnormgamma <- function(theta, muphi, lambda, nu1, nu2) {
  margprec <- (nu1 * lambda) / (nu2 * (lambda + 1))
  margdf <- 2 * nu1

  A <- lgamma((margdf + 1) / 2) - lgamma(margdf / 2)
  B <- 0.5 * (log(margprec) - log(margdf) - log(pi))
  C <- -((margdf + 1) / 2) * log(1 + (margprec * (theta - muphi)^2) / margdf)
  logden <- A + B + C
}


# This function will update DP process parameter alpha by MH
# Arguments are:
# c - vector of the classes of each observation
# alpha - parameter in Dir(alpha)
# mualpha - mean of log-normal
# sigalpha - sd of lognormal
.Updatealphalognormpr <- function(c, alpha, mualpha = -3, sigalpha = 1, propsd = 1) {
  uold <- log(alpha)
  unew <- stats::rnorm(1, uold, propsd)
  alphanew <- exp(unew)
  logprrat <- stats::dnorm(unew, mean = mualpha, sd = sigalpha, log = TRUE) - stats::dnorm(uold, mean = mualpha, sd = sigalpha, log = TRUE)
  loglikrat <- .LogLikalpha(c, alphanew) - .LogLikalpha(c, alpha)
  HR <- exp(logprrat + loglikrat)
  if (is.na(HR)) {
    cat(logprrat, "    ", loglikrat, "\n")
    cat(alpha, "     ", alphanew, "\n")
    cat(c, "\n")
    stop("Not a number")
  }

  if (stats::runif(1) < HR) {
    return(alphanew) # Accept alphanew
  }
  return(alpha) # Or reject and keep alpha
}


# Function as above but using a gamma prior on the value of alpha
.Updatealphagammapr <- function(c, alpha, prshape = 0.5, prrate = 1, propsd = 1) {
  # Sample new alpha from truncated normal distribution
  repeat {
    alphanew <- stats::rnorm(1, alpha, propsd)
    if (alphanew > 0) {
      break
    }
  }
  logprrat <- stats::dgamma(alphanew, shape = prshape, rate = prrate, log = TRUE) - stats::dgamma(alpha, shape = prshape, rate = prrate, log = TRUE)
  loglikrat <- .LogLikalpha(c, alphanew) - .LogLikalpha(c, alpha)
  logproprat <- stats::pnorm(alpha / propsd, log.p = TRUE) - stats::pnorm(alphanew / propsd, log.p = TRUE) # Adjust for non-symmetric truncated normal proposal
  HR <- exp(logprrat + loglikrat + logproprat)
  if (stats::runif(1) < HR) {
    return(alphanew) # Accept alphanew
  }
  return(alpha) # Or reject and keep alpha
}

# This function will work out the likelihood of a particular Dir(alpha)
# given a partition c
# Arguments are:
# c - vector of the classes of each observation
# alpha - parameter in Dir(alpha)
# Return the likelihood
.LogLikalpha <- function(c, alpha) {
  n <- length(c)
  nc <- max(c)
  nci <- apply(as.row(1:nc), 2, function(x, c) sum(c == x), c = c)
  #  cat("nci = ", nci, "\n")
  #  lik <- ((alpha^nc)*prod(factorial(nci-1))) / prod(alpha:(alpha+n-1)) # Need to work with log-likelihood
  loglik <- nc*log(alpha) + sum(sapply(pmax(nci-1, 1), function(x) sum(log(1:x)))) - sum(log(alpha:(alpha+n-1)))
  # Note we have to use pmax(nci-1, 1) here to account for clusters of size 1 have 0! = 1
  return(loglik)
}
