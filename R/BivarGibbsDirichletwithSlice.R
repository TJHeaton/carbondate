#' Performs the Gibbs sampler to estimate the density of
#' a set of observations using a Dirichlet Mixture Density
#'
#' THIS CONSIDERS BOTH MEAN AND VARIANCE OF THE CLUSTERS TO BE UNKNOWN
#'
#' @param x TODO
#' @param xsig TODO
#' @param lambda,nu1,nu2  Hyperparameters for the prior on the means
#' \eqn{\phi_j} and precision \eqn{\tau_j} of each individual calendar age
#' cluster \eqn{j}.
#' \deqn{(\phi_j, \tau_j)|\mu_{\phi} \sim
#' \textrm{NormalGamma}(\mu_{\phi}, \lambda, \nu_1, \nu_2)} where
#' \eqn{\mu_{\phi}} is the overall cluster centering.
#' @param A,B  Prior on \eqn{\mu_{\phi}} giving the mean and precision of the
#' overall centering \eqn{\mu_{\phi} \sim N(A, B^{-1})} i.e.
#' B small is uninformative.
#' @param mualpha TODO
#' @param sigalpha TODO
#' @param alphaprshape TODO
#' @param alphaprrate TODO
#' @param niter  The number of MCMC iterations (optional). Default is 100.
#' @param nthin  How much to thin the output (optional). 1 is no thinning,
#' a larger number is more thinning. Default is 10. Must choose an integer more
#' than 1 and not too close to `n_iter`, since after burn-in there are
#' \eqn{(n_{\textrm{iter}}/n_{\textrm{thin}})/2} samples from posterior to
#' potentially use.
#' @param theta The initial estimate for the underlying calendar ages
#' (optional). If supplied it must be a vector with the same length as
#' `c14_determinations`. Will be overridden if `sensible_initialisation` is
#' `TRUE`.
#' @param alphapr TODO
#' @param w TODO
#' @param m TODO
#' @param calibration_curve TODO
#' @param nclusinit TODO
#' @param sensibleinit TODO
#' @param showprogress Whether to show a progress bar in the console during
#' execution. Default is `TRUE`.
#'
#' @return A list containing 6 items
#' @export
#'
BivarGibbsDirichletwithSlice <- function(x, xsig,
                                         lambda, nu1, nu2, A, B, mualpha, sigalpha, alphaprshape, alphaprrate,
                                         niter = 100, nthin = 10, theta = NA, alphapr = "gamma",
                                         w = 200, m = 50, calibration_curve, nclusinit = 10,
                                         sensibleinit = TRUE, showprogress = TRUE) {

  # Find the number of observations
  nobs <- length(x)

  # Initialise variables c, muphi, alpha and (phi, tau)
  if (is.na(nclusinit)) nclusinit <- 1 # Single cluster

  # Sample classes (but make sure all are represented)
  allclass <- FALSE
  while (!allclass) { # Edited to allow different weights for mixture distribution
    c <- sample(1:nclusinit, nobs, replace = TRUE) # Sample some classes (make sure all are represented)
    allclass <- (length(unique(c)) == nclusinit)
  }

  # Initalise theta by choosing maximum posterior from independent coarse version
  if (sensibleinit) {
    initial_probabilities <- mapply(
      CalibrateSingleDetermination,
      x,
      xsig,
      MoreArgs = list(calibration_curve=calibration_curve))
    indices_of_max_probability = apply(initial_probabilities, 2, which.max)
    theta <- calibration_curve$calendar_age[indices_of_max_probability]
    muphi <- stats::median(theta)
    # Over-ride A and B in prior for muphi from range of theta
    A <- stats::median(theta)
    B <- 1 / (max(theta) - min(theta))^2
  } else {
    muphi <- mean(x) * 8267 / 8033 # Start off with muphi set at mean(x)*8267/8033
    if (is.na(theta[1])) theta <- x * 8267 / 8033 # Initial assumption that age = 14C determination * 8267/8033
  }

  # Set a value for alpha from prior
  alpha <- switch(alphapr,
                  lognorm = exp(stats::rnorm(1, mualpha, sd = sigalpha)),
                  gamma = 0.0001, # stats::rgamma(1, shape = alphaprshape, rate = alphaprrate),
                  stop("Unknown form for prior on gamma")
  )

  nclus <- length(unique(c))
  tau <- rep(nclus, 1 / (diff(range(x)) / 4)^2)
  phi <- stats::rnorm(nclus, mean = muphi, sd = diff(range(x)) / 2) # min = 0.8*min(x), max = 1.2*max(x) ) * 8267/8033

  phiout <- list(phi) # Needs to be a ragged array i.e. list
  tauout <- list(tau)
  cout <- matrix(NA, nrow = floor(niter / nthin) + 1, ncol = nobs)
  thetaout <- matrix(NA, nrow = floor(niter / nthin) + 1, ncol = nobs)
  alphaout <- rep(NA, length = floor(niter / nthin) + 1)
  muphiout <- rep(NA, length = floor(niter / nthin) + 1)
  outi <- 1

  cout[outi, ] <- c
  thetaout[outi, ] <- theta
  alphaout[outi] <- alpha
  muphiout[outi] <- muphi

  ## Interpolate onto single year grid to speed up updating thetas
  integer_cal_year_curve <- InterpolateCalibrationCurve(
    1:pkg.globals$MAX_YEAR_BP, calibration_curve)
  mucalallyr <- integer_cal_year_curve$c14_age
  sigcalallyr <- integer_cal_year_curve$c14_sig


  # Create a progress bar if we want to show progress
  if (showprogress) {
    pb <- utils::txtProgressBar(min = 0, max = niter, style = 3)
  }

  for (iter in 1:niter) {
    if (showprogress) {
      # Update progress bar every 1000 iterations
      if (iter %% 100 == 0) {
        utils::setTxtProgressBar(pb, iter)
      }
    }

    # Update c
    for (i in 1:nobs) {
      newclusters <- .BivarUpdatec(i, c = c, phi = phi, tau = tau, theta = theta[i], lambda = lambda, nu1 = nu1, nu2 = nu2, muphi = muphi, alpha = alpha)
      c <- newclusters$c
      phi <- newclusters$phi
      tau <- newclusters$tau
      if (max(c) != length(phi)) stop("Lengths do not match")
    }
    # Update phi and tau values
    for (j in 1:length(phi)) {
      GibbsParams <- .UpdatePhiTau(theta[c == j], mu_phi = muphi, lambda = lambda, nu1 = nu1, nu2 = nu2)
      phi[j] <- GibbsParams$phi
      tau[j] <- GibbsParams$tau
    }
    # Update muphi based on current phi and tau values
    muphi <- .UpdateMuPhi(phi = phi, tau = tau, lambda = lambda, A = A, B = B)


    # Update y values (the calendar ages of the objects) by slice sampling
    for (k in 1:nobs) {
      theta[k] <- .SliceSample(
        TARGET = .ThetaLogLikelihood,
        x0 = theta[k],
        slice_width = w,
        slice_multiplier = m,
        type = "log",
        prmean = phi[c[k]],
        prsig = 1 / sqrt(tau[c[k]]),
        c14obs = x[k],
        c14sig = xsig[k],
        mucalallyr = mucalallyr,
        sigcalallyr = sigcalallyr)
    }

    # Update alpha
    alpha <- switch(alphapr,
                    lognorm = .Updatealphalognormpr(c, alpha, mualpha = mualpha, sigalpha = sigalpha),
                    gamma = .Updatealphagammapr(c, alpha, prshape = alphaprshape, prrate = alphaprrate),
                    stop("Unknown form for prior on gamma")
    )

    # Store thinned output
    if (iter %% nthin == 0) {
      outi <- outi + 1
      phiout[[outi]] <- phi
      tauout[[outi]] <- tau
      cout[outi, ] <- c
      thetaout[outi, ] <- theta
      alphaout[outi] <- alpha
      muphiout[outi] <- muphi
    }
  }
  retlist <- list(c = cout, phi = phiout, tau = tauout, theta = thetaout, alpha = alphaout, muphi = muphiout)

  if (showprogress) close(pb)
  return(retlist)
}
