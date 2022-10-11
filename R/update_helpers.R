############################################################
# This function will sample from the posterior of (phi, tau) given a set of calendar ages from that cluster
# For each cluster in turn we have
# theta ~ N(phi, sd = 1/sqrt(tau))
# Arguments:
# theta - vector of object calendar ages (all belonging to same cluster)
# muphi, lambda, nu1, nu2 - parameters in NormalGamma for (phi, tau)
.Updatephitau <- function(theta, muphi, lambda, nu1, nu2) {
  nclus <- length(theta)
  thetabar <- mean(theta)
  s <- mean((theta - thetabar)^2)
  # Update parameters according to conjugate prior
  nu1new <- nu1 + nclus / 2
  nu2new <- nu2 + 0.5 * ((nclus * s) + (lambda * nclus * ((thetabar - muphi)^2) / (lambda + nclus)))
  lambdanew <- lambda + nclus
  muphinew <- (lambda * muphi + nclus * thetabar) / (lambda + nclus)
  # Now sample new values for precision tau and mean phi
  taunew <- stats::rgamma(1, shape = nu1new, rate = nu2new) # Sample new tau
  phinew <- stats::rnorm(1, mean = muphinew, sd = 1 / sqrt(lambdanew * taunew)) # and then phi

  retlist <- list(tau = taunew, phi = phinew)
  return(retlist)
}


# Function which updates muphi via Gibbs sampling based upon current (phi, tau) values
# Arguments:
# phi - current vector of phi values (means of clusters)
# tau - current vector of precisions
# lambda - current value of lambda
# A - prior mean of muphi
# B - prior precision for muphi
.Updatemuphi <- function(phi, tau, lambda, A, B) {
  locprec <- tau * lambda # We then observe phi_i ~ N(muphi, prec = lambda*tau_i)

  postmean <- (A * B + sum(locprec * phi)) / (B + sum(locprec))
  postprec <- B + sum(locprec)
  stats::rnorm(1, mean = postmean, sd = 1 / sqrt(postprec))
}


# Create a function which works out the log-likelihood
# (i.e. log-target for any z (theta) given other parameters
# A fast version which assume IntCal is constant on a year scale
.thetaloglikfast <- function(z, prmean, prsig, c14obs, c14sig, mucalallyr, sigcalallyr) {
  yr <- floor(z) # Assume the cal likelihood is constant on year scale
  # Return -Inf if outside range
  if (yr < 1 | yr > 50000) {
    return(-Inf)
  }

  mucal <- mucalallyr[yr]
  sigcal <- sigcalallyr[yr]
  loglik <- stats::dnorm(z, mean = prmean, sd = prsig, log = TRUE) +
    stats::dnorm(c14obs, mean = mucal, sd = sqrt(sigcal^2 + c14sig^2), log = TRUE)
  return(loglik)
}
