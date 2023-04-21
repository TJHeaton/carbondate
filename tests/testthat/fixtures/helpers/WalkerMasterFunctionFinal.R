# This file contains a function which implements calibration using Walker updating of the DPMM

# Arguments are:
# x, xsig - the radiocarbon determinations and their uncertainties
# lambda, nu1, nu2 - the hyperparameters for the prior on the (phi, tau) ~ NormalGamma(muphi, lambda, nu1, nu2)
# A, B - prior on muphi ~ N(A, B^-1) i.e. B small is uninformative
# cprshape, cprrate - the hyperparameters on the DP parameters c ~ Gamma(shape, rate)
# niter, nthin - the number of MCMC iterations and how much to thin
# theta - initial estimate of the underlying calibrated dates
# wslice, m - parameters for slice sampling
# calcurve - the calibration curve we are using
# kstar - intiial number of clusters
# sensibleinit - logical as to wheteher we choose sensible starting values and adaptive prior on muphi (A, B)
WalkerBivarDirichlet <- function(x, xsig,
                                 lambda, nu1, nu2, A, B, cprshape, cprrate,
                                 niter = 100, nthin = 10, theta = NA,
                                 slicew = 1000, m = 10, calcurve, kstar = 10,
                                 sensibleinit = TRUE, showprogress = TRUE) {

  # Initialise parameters
  nobs <- length(x)

  # Initalise theta by choosing maximum from independent coarse version
  if(sensibleinit) {
    IntCalyrgrid <- FindCal(1:50000, calcurve$c14ag, calcurve$calage, calcurve$c14sig)
    initprobs <- mapply(calibind, x, xsig, MoreArgs = list(calmu = IntCalyrgrid$mu, calsig = IntCalyrgrid$sigma))
    theta <- apply(initprobs, 2, which.max)
    muphi <- median(theta)
    # Over-ride A and B from range of theta
    A <- median(theta)
    B <- 1 / (max(theta) - min(theta))^2
  } else {
    muphi <- mean(x) * 8267/8033 # Start off with muphi set at mean(x)*8267/8033
    if(is.na(theta[1])) theta <- x*8267/8033 # Initial assumption that age = 14C determination * 8267/8033
  }

  # Choose 4 initial cluster
  if(is.na(kstar)) kstar <- 4

  # Sample initial value of Dirichlet parameter from gamma prior (do not allow v small vallues of c as this causes crashes)
  c <- 2 # rgamma(1, shape = cprshape, rate = cprrate)

  # Choose initial tau and phi
  tau <- rgamma(kstar, shape = nu1, rate = nu2)
  phi <- rnorm(kstar, mean = muphi, sd = 1/sqrt(lambda*tau))

  # Initialise v, w and delta
  v <- rbeta(kstar, 1, c)
  w <- v * c(1, cumprod(1-v)[-kstar]) # These are the weights(w1, ..., w_k*)
  delta <- sample(1:kstar, nobs, replace = TRUE) # Initial allocations

  ##############################################################################
  #### Create storage for output

  # Create an output list (ragged array or matrix)
  phiout <- list(phi)
  tauout <- list(tau)
  wout <- list(w)
  deltaout <- matrix(NA, nrow = floor(niter/nthin) + 1, ncol = nobs)
  cout <- rep(NA, length = floor(niter/nthin) + 1)
  ndclustout <- rep(NA, length = floor(niter/nthin) + 1)
  muphiout <- rep(NA, length = floor(niter/nthin) + 1)
  thetaout <- matrix(NA, nrow = floor(niter/nthin) + 1, ncol = nobs)

  outi <- 1
  deltaout[outi,] <- delta
  cout[outi] <- c
  ndclustout[outi] <- length(unique(delta))
  muphiout[outi] <- muphi
  thetaout[outi,] <- theta

  ################################################################################
  # Interpolate calibration curve onto single year to speed up code
  # Read in the calibration curve and store output on annual grid
  caltheta <- calcurve$calage
  calmu <- calcurve$c14age
  calsig <- calcurve$c14sig

  ## Interpolate onto single year grid to speed up updating thetas
  IntCalyrgrid <- FindCal(1:50000, calmu, caltheta, calsig)
  mucalallyr <- IntCalyrgrid$mu
  sigcalallyr <- IntCalyrgrid$sigma

  # Create a progress bar if we want to show progress
  if(showprogress) {
    pb <- txtProgressBar(min = 0, max = niter, style = 3)
  }

  #####################################
  # Now the calibration and DPMM
  for(MHiter in 1:niter) {
    if(showprogress) {
      # Update progress bar every 1000 iterations
      if(MHiter %% 100 == 0){
        setTxtProgressBar(pb, MHiter)
      }
    }
    # if(MHiter == 1307) browser()
    # Perform a Walker update of the DPMM
    DPMMUpdate <- DPWalkerUpdate(theta = theta, w = w, v = v, delta = delta, phi = phi, tau = tau, kstar = kstar, c = c, muphi = muphi,
                                 lambda = lambda, nu1 = nu1, nu2 = nu2)
    w <- DPMMUpdate$w
    delta <- DPMMUpdate$delta
    phi <- DPMMUpdate$phi
    tau <- DPMMUpdate$tau
    v <- DPMMUpdate$v
    kstar <- DPMMUpdate$kstar

    # Now update the other parameters (c and muphi)
    # Update the Dirichlet parameter c (stick breaking part )
    c <- WalkerUpdatec(delta = delta, c = c, prshape = cprshape, prrate = cprrate)

    # Update muphi based on current phi and tau values
    muphi <- Updatemuphi(phi = phi, tau = tau, lambda = lambda, A = A, B = B)

    # Update theta values (the calendar ages of the objects) by slice sampling
    for(k in 1:nobs) {
      theta[k] <- SliceSample(TARGET = thetaloglikfast, x0 = theta[k],
                              w = slicew, m = m, type = "log",
                              prmean = phi[delta[k]], prsig = 1/sqrt(tau[delta[k]]),
                              c14obs = x[k], c14sig = xsig[k],
                              mucalallyr = mucalallyr, sigcalallyr = sigcalallyr)
    }

    # Save the output on every nthin iteration
    if(MHiter %% nthin == 0) {
      outi <- outi+1
      deltaout[outi,] <- delta
      cout[outi] <- c
      ndclustout[outi] <- length(unique(delta))
      phiout[[outi]] <- phi
      tauout[[outi]] <- tau
      thetaout[outi,] <- theta
      wout[[outi]] <- w
      muphiout[outi] <- muphi
    }
  }
  retlist <- list(delta = deltaout, c = cout, nclust = ndclustout, phi = phiout, tau = tauout, theta = thetaout, w = wout, muphi = muphiout)
  if(showprogress) close(pb)
  return(retlist)
}
