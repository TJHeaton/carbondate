#' Implements calibration using Walker updating of the DPMM
#'
#' @param c14_determinations A vector containing the radiocarbon determinations
#' @param c14_uncertainties A vector containing the radiocarbon determination uncertainties.
#' Must be the same length as `c14_determinations`.
#' @param calibration_curve A dataframe which should contain one column entitled
#' c14_age and one column entitled c14_sig.
#' This format matches [carbondate::intcal20].
#' @param lambda  TODO
#' @param nu1  TODO
#' @param nu2  TODO
#' @param A  TODO
#' @param B  TODO
#' @param cprshape  TODO
#' @param cprrate  TODO
#' @param niter  TODO
#' @param nthin  TODO
#' @param theta  TODO
#' @param slicew  TODO
#' @param m  TODO
#' @param kstar TODO
#' @param sensibleinit TODO
#' @param showprogress  TODO
#'
#' @return a list TODO add more
#' @export
#'
#' @examples WalkerBivarDirichlet(
#'   c14_determinations=c(602, 805, 1554),
#'   c14_uncertainties=c(35, 34, 45),
#'   calibration_curve=intcal20,
#'   lambda=1,
#'   nu1=1,
#'   nu2=1,
#'   A=1,
#'   B=1,
#'   cprshape=1,
#'   cprrate=1)
WalkerBivarDirichlet <- function(
    c14_determinations,
    c14_uncertainties,
    calibration_curve,
    lambda,
    nu1,
    nu2,
    A,
    B,
    cprshape,
    cprrate,
    niter = 100,
    nthin = 10,
    theta = NA,
    slicew = 1000,
    m = 10,
    kstar = 10,
    sensibleinit = TRUE,
    showprogress = TRUE) {

  # Initialise parameters
  nobs <- length(c14_determinations)

  # Initalise theta by choosing maximum from independent coarse version
  if (sensibleinit) {
    initprobs <- mapply(
      CalibrateSingleDetermination,
      c14_determinations,
      c14_uncertainties,
      MoreArgs = list(calibration_curve=calibration_curve))
    theta <- calibration_curve$calendar_age[apply(initprobs, 2, which.max)]
    muphi <- stats::median(theta)
    # Over-ride A and B from range of theta
    A <- stats::median(theta)
    B <- 1 / (max(theta) - min(theta))^2
  } else {
    # TODO: Where do these numbers come from?
    muphi <- mean(c14_determinations) * 8267 / 8033 # Start off with muphi set at mean(x)*8267/8033
    if (is.na(theta[1])) theta <- c14_determinations * 8267 / 8033 # Initial assumption that age = 14C determination * 8267/8033
  }

  # Choose 4 initial cluster
  if (is.na(kstar)) kstar <- 4

  # Sample initial value of Dirichlet parameter from gamma prior (do not allow v small values of c as this causes crashes)
  c <- 2 # stats::rgamma(1, shape = cprshape, rate = cprrate)

  # Choose initial tau and phi
  tau <- stats::rgamma(kstar, shape = nu1, rate = nu2)
  phi <- stats::rnorm(kstar, mean = muphi, sd = 1 / sqrt(lambda * tau))

  # Initialise v, w and delta
  v <- stats::rbeta(kstar, 1, c)
  w <- v * c(1, cumprod(1 - v)[-kstar]) # These are the weights(w1, ..., w_k*)
  delta <- sample(1:kstar, nobs, replace = TRUE) # Initial allocations

  ##############################################################################
  #### Create storage for output

  # Create an output list (ragged array or matrix)
  phiout <- list(phi)
  tauout <- list(tau)
  wout <- list(w)
  deltaout <- matrix(NA, nrow = floor(niter / nthin) + 1, ncol = nobs)
  cout <- rep(NA, length = floor(niter / nthin) + 1)
  ndclustout <- rep(NA, length = floor(niter / nthin) + 1)
  muphiout <- rep(NA, length = floor(niter / nthin) + 1)
  thetaout <- matrix(NA, nrow = floor(niter / nthin) + 1, ncol = nobs)

  outi <- 1
  deltaout[outi, ] <- delta
  cout[outi] <- c
  ndclustout[outi] <- length(unique(delta))
  muphiout[outi] <- muphi
  thetaout[outi, ] <- theta

  ################################################################################

  ## Interpolate onto single year grid to speed up updating thetas
  IntCalyrgrid <- InterpolateCalibrationCurve(1:50000, calibration_curve)
  mucalallyr <- IntCalyrgrid$c14_age
  sigcalallyr <- IntCalyrgrid$c14_sig

  if (showprogress) {
    pb <- utils::txtProgressBar(min = 0, max = niter, style = 3)
  }

  #####################################
  # Now the calibration and DPMM
  for (MHiter in 1:niter) {
    if (showprogress) {
      # Update progress bar every 1000 iterations
      if (MHiter %% 100 == 0) {
        utils::setTxtProgressBar(pb, MHiter)
      }
    }
    DPMMUpdate <- .DPWalkerUpdate(
      theta = theta, w = w, v = v, delta = delta, phi = phi, tau = tau, kstar = kstar, c = c, muphi = muphi,
      lambda = lambda, nu1 = nu1, nu2 = nu2
    )
    w <- DPMMUpdate$w
    delta <- DPMMUpdate$delta
    phi <- DPMMUpdate$phi
    tau <- DPMMUpdate$tau
    v <- DPMMUpdate$v
    kstar <- DPMMUpdate$kstar

    c <- .WalkerUpdatec(delta = delta, c = c, prshape = cprshape, prrate = cprrate)
    muphi <- .Updatemuphi(phi = phi, tau = tau, lambda = lambda, A = A, B = B)

    # Update theta values (the calendar ages of the objects) by slice sampling
    for (k in 1:nobs) {
      theta[k] <- .SliceSample(
        TARGET = .thetaloglikfast, x0 = theta[k],
        w = slicew, m = m, type = "log",
        prmean = phi[delta[k]], prsig = 1 / sqrt(tau[delta[k]]),
        c14obs = c14_determinations[k], c14sig = c14_uncertainties[k],
        mucalallyr = mucalallyr, sigcalallyr = sigcalallyr
      )
    }

    # Save the output on every nthin iteration
    if (MHiter %% nthin == 0) {
      outi <- outi + 1
      deltaout[outi, ] <- delta
      cout[outi] <- c
      ndclustout[outi] <- length(unique(delta))
      phiout[[outi]] <- phi
      tauout[[outi]] <- tau
      thetaout[outi, ] <- theta
      wout[[outi]] <- w
      muphiout[outi] <- muphi
    }
  }
  retlist <- list(delta = deltaout, c = cout, nclust = ndclustout, phi = phiout, tau = tauout, theta = thetaout, w = wout, muphi = muphiout)
  if (showprogress) close(pb)
  return(retlist)
}
