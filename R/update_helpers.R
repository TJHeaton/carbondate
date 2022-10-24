pkg.globals <- new.env()

pkg.globals$MAX_YEAR_BP <- 50000


############################################################
# This function will sample from the posterior of (phi, tau) given a set of
# calendar ages from that cluster
# For each cluster in turn we have
# theta ~ N(phi, sd = 1/sqrt(tau))
# Arguments:
# theta - vector of object calendar ages (all belonging to same cluster)
# mu_phi, lambda, nu1, nu2 - parameters in NormalGamma for (phi, tau)
.UpdatePhiTau <- function(theta, mu_phi, lambda, nu1, nu2) {
  nclus <- length(theta)
  thetabar <- mean(theta)
  s <- mean((theta - thetabar)^2)
  # Update parameters according to conjugate prior
  nu1new <- nu1 + nclus / 2
  nu2new <- nu2 + 0.5 *
    (
      (nclus * s)
      + (lambda * nclus * ((thetabar - mu_phi)^2) / (lambda + nclus)))
  lambdanew <- lambda + nclus
  muphinew <- (lambda * mu_phi + nclus * thetabar) / (lambda + nclus)
  # Now sample new values for precision tau and mean phi
  taunew <- stats::rgamma(1, shape = nu1new, rate = nu2new)
  phinew <- stats::rnorm(1, mean = muphinew, sd = 1 / sqrt(lambdanew * taunew))

  return_list <- list(tau = taunew, phi = phinew)
  return(return_list)
}


# Function which updates mu_phi via Gibbs sampling based upon current (phi, tau)
# values
# Arguments:
# phi - current vector of phi values (means of clusters)
# tau - current vector of precisions
# lambda - current value of lambda
# A - prior mean of mu_phi
# B - prior precision for mu_phi
.UpdateMuPhi <- function(phi, tau, lambda, A, B) {
  locprec <- tau * lambda
  # We then observe phi_i ~ N(mu_phi, prec = lambda*tau_i)

  postmean <- (A * B + sum(locprec * phi)) / (B + sum(locprec))
  postprec <- B + sum(locprec)
  stats::rnorm(1, mean = postmean, sd = 1 / sqrt(postprec))
}


# Calculates the log-target for any calendar year (theta) given other parameters
# A fast version which assume the calibration curve is constant on a year scale
.ThetaLogLikelihood <- function(
    theta, prmean, prsig, c14obs, c14sig, mucalallyr, sigcalallyr) {
  # Assume the cal likelihood is constant on year scale
  yr <- floor(theta)

  if (yr < 1 | yr > pkg.globals$MAX_YEAR_BP) {
    return(-Inf)
  }

  mucal <- mucalallyr[yr]
  sigcal <- sigcalallyr[yr]
  loglik <- stats::dnorm(theta, mean = prmean, sd = prsig, log = TRUE) +
    stats::dnorm(
      c14obs, mean = mucal, sd = sqrt(sigcal^2 + c14sig^2), log = TRUE)
  return(loglik)
}
