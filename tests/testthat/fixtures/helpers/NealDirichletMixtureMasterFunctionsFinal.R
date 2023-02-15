#### IN THIS VERSION WE HAVE TWO FUNCTIONS:
## UniGibbsDirichlet --- ALLOWS JUST THE MEAN OF THE CLUSTERS TO VARY (with variance identical, fixed and known for all clusters)
## BivarGibbsDirichlet --- ALLOWS BOTH THE MEAN AND VARIANCE OF THE CLUSTERS TO VARY

# This file will code up Dirichlet Process density estimation
# We have variables:
# Observed
# X - the radicarbon determinations of the observations
# Unknown
# c - the classification/cluster parameter
# phi - the means of the normal distribution for each cluster
# theta - calendar ages of the objects - prior theta[i] ~ N(phi[c[i]], sigma^2)


###########################################################################################
# Function which will perform the Gibbs sampler to estimate the density of
# a set of observations using a Dirichlet Mixture Density
## THIS CONSIDERS JUST THE MEAN OF THE CLUSTERS TO BE UNKNOWN
# Arguments:
# x - the set of radicoarbon determinations
# xsig - sd on radiocarbon determinations
# sigma - the sd on underlying calendar ages in each cluster i.e. theta ~ N(phi, sd = sigma)
# prmean and prsd - the mean and sd of G0 i.e. phi ~ N(prmean, sd = prsd)
# alphainit - initial guess for parameter in stick breaking - updated within sampler
# propsd - the standard deviation for proposing a new calendar age in MH
# theta - an initial guess for calendar ages if we have one
# alphapr - the type of prior on alpha in DP - "lognorm" or "gamma"
UniGibbsDirichlet <- function(x, xsig, sigma = 1,
                              prmean = 0, alpha = 1, prsd = 1,
                              niter = 100, nthin = 10, propsd = 1,
                              theta = NA, alphapr = "lognorm") {
  nobs <- length(x)

  # Create an initial guess for theta (the underlying calendar ages)
  ## Might want to change this to make a better initial guess (e.g. individual posterior means given X)
  if(is.na(theta[1])) theta <- rnorm(nobs, mean = prmean, sd = sigma)

  # Initialise c and phi
  c <- rep(1, nobs) # All in the same class
  phi <- rnorm(length(unique(c)), mean = prmean, sd = prsd) # Initialise from G0


  phiout <- list(phi) # Needs to be a ragged array i.e. list
  cout <- matrix(NA, nrow = floor(niter/nthin) + 1, ncol = nobs)
  thetaout <- matrix(NA, nrow = floor(niter/nthin) + 1, ncol = nobs)
  alphaout <- rep(NA, length = floor(niter/nthin) + 1)
  outi <- 1
  cout[outi,] <- c
  thetaout[outi,] <- theta
  alphaout[outi] <- alpha

  # Read in calibration curve
  calcurve <- read.csv("intcal13.csv", sep = ",", header = TRUE)
  caltheta <- calcurve$calage
  calmu <- calcurve$c14age
  calsig <- calcurve$c14sig

  # Create a progress bar
  pb <- txtProgressBar(min = 0, max = niter, style = 3)

  for(iter in 1:niter) {
    # Update progress bar every 1000 iterations
    if(iter %% 1000 == 0){
      setTxtProgressBar(pb, iter)
    }

    # Update c
    for(i in 1:nobs) {
      newclusters <- Updatec(i, c = c, phi= phi, theta = theta[i], sigma = sigma,
                             prmean = prmean, prsd = prsd, alpha = alpha)
      c <- newclusters$c
      phi <- newclusters$phi
      if(max(c) != length(phi)) stop("Lengths do not match")
    }
    # Update phi values
    for(j in 1:length(phi)) {
      phi[j] <- Updatephi(theta[c == j], sigma = sigma, prmean = prmean, prsd = prsd)
    }
    # Update y values (the calendar ages of the objects) by MH
    for(k in 1:nobs) {
      theta[k] <- Updatecalage(theta = theta[k], phi = phi[c[k]], sigma = sigma, x = x[k], xsig = xsig[k],
                               calmu = calmu, caltheta = caltheta, calsig = calsig, propsd = propsd)
    }
    # Update alpha
    alpha <- switch(alphapr,
                    lognorm = Updatealphalognormpr(c, alpha),
                    gamma = Updatealphagammapr(c, alpha),
                    stop("Unknown form for prior on gamma"))
    if(iter %% nthin == 0) {
      cat("Iteration = ", iter, "\n") # Progress just to show it is still working
      outi <- outi+1
      phiout[[outi]] <- phi
      cout[outi,] <- c
      thetaout[outi,] <- theta
      alphaout[outi] <- alpha
    }
  }
  close(pb)
  retlist <- list(c = cout, phi = phiout, theta = thetaout, alpha = alphaout)
  return(retlist)
}



###########################################################################################
# Function which will perform the Gibbs sampler to estimate the density of
# a set of observations using a Dirichlet Mixture Density
## THIS CONSIDERS BOTH MEAN AND VARIANCE OF THE CLUSTERS TO BE UNKNOWN
# We use a normal-gamma prior for the clusters which has the following 5 hyperparameters
# phi|muphi, tau ~ N(muphi, precision = (lambda * tau))
# and
# tau ~ gamma(shape = nu1, rate = nu2)
# muphi ~ N(A, precision = B)

##### Arguments:
## Observed:
# x - the set of radicoarbon determinations
# xsig - sd on radiocarbon determinations
## Hyperparameters:
# lambda - prior on precision on phi
# nu1, nu2 - shap and rate of tau in normal-gamma for phi
# A, B - prior parameters for muphi i.e. muphi ~ N(A, precision = B)
# mualpha, musigma - prior parameters for alpha (Dirichlet process)
## Tuning and user selected
# niter - number of iterations to perform
# nthin - output every nthin iteration
# theta - an initial guess for calendar ages if we have one
# alphapr - the type of prior on alpha in DP - "lognorm" or "gamma"
## Slice sampling parameters
# w - estimate of typical slice size
# m - integer limiting slice size to mw
# sensibleinit - logical as to whether we choose sensible starting values and adaptive prior on muphi (A, B)
BivarGibbsDirichletwithSlice <- function(x, xsig,
                                         lambda, nu1, nu2, A, B, mualpha, sigalpha, alphaprshape, alphaprrate,
                                         niter = 100, nthin = 10, theta = NA, alphapr = "gamma",
                                         w = 200, m = 50, calcurve, nclusinit = 10,
                                         sensibleinit = TRUE, showprogress = TRUE) {

  # Find the number of observations
  nobs <- length(x)

  # Initialise variables c, muphi, alpha and (phi, tau)
  if(is.na(nclusinit)) nclusinit <- 1 # Single cluster

  # Sample classes (but make sure all are represented)
  allclass <- FALSE
  while(!allclass) { # Edited to allow different weights for mixture distribution
   c <- sample(1:nclusinit, nobs, replace = TRUE) # Sample some classes (make sure all are represented)
   allclass <- (length(unique(c)) == nclusinit)
  }

  # Initalise theta by choosing maximum posterior from independent coarse version
  if(sensibleinit) {
    initprobs <- mapply(calibind, x, xsig, MoreArgs = list(calmu = calcurve$c14age, calsig = calcurve$c14sig))
    theta <- calcurve$calage[apply(initprobs, 2, which.max)]
    muphi <- median(theta)
    # Over-ride A and B in prior for muphi from range of theta
    A <- median(theta)
    B <- 1 / (max(theta) - min(theta))^2
  } else {
    muphi <- mean(x) * 8267/8033 # Start off with muphi set at mean(x)*8267/8033
    if(is.na(theta[1])) theta <- x*8267/8033 # Initial assumption that age = 14C determination * 8267/8033
  }

  # Set a value for alpha from prior
  alpha <- switch(alphapr,
                  lognorm = exp(rnorm(1, mualpha, sd = sigalpha)),
                  gamma = 0.0001, # rgamma(1, shape = alphaprshape, rate = alphaprrate),
                  stop("Unknown form for prior on gamma"))

  nclus <- length(unique(c))
  tau <- rep( 1 / (diff(range(x))/4)^2 , nclus)
  phi <- rnorm(nclus, mean = muphi, sd = diff(range(x))/2)  # min = 0.8*min(x), max = 1.2*max(x) ) * 8267/8033

  phiout <- list(phi) # Needs to be a ragged array i.e. list
  tauout <- list(tau)
  cout <- matrix(NA, nrow = floor(niter/nthin) + 1, ncol = nobs)
  thetaout <- matrix(NA, nrow = floor(niter/nthin) + 1, ncol = nobs)
  alphaout <- rep(NA, length = floor(niter/nthin) + 1)
  muphiout <- rep(NA, length = floor(niter/nthin) + 1)
  outi <- 1

  cout[outi,] <- c
  thetaout[outi,] <- theta
  alphaout[outi] <- alpha
  muphiout[outi] <- muphi

  ## Read in calibration curve
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

  for(iter in 1:niter) {
    if(showprogress) {
      # Update progress bar every 1000 iterations
      if(iter %% 100 == 0){
        setTxtProgressBar(pb, iter)
      }
    }

    # Update c
    for(i in 1:nobs) {
      newclusters <- BivarUpdatec(i, c = c, phi = phi, tau = tau, theta = theta[i], lambda = lambda, nu1 = nu1, nu2 = nu2, muphi = muphi, alpha = alpha)
      c <- newclusters$c
      phi <- newclusters$phi
      tau <- newclusters$tau
      if(max(c) != length(phi)) stop("Lengths do not match")
    }
    # Update phi and tau values
    for(j in 1:length(phi)) {
      GibbsParams <- BivarUpdatephitau(theta[c == j], muphi = muphi, lambda = lambda, nu1 = nu1, nu2 = nu2)
      phi[j] <- GibbsParams$phi
      tau[j] <- GibbsParams$tau
    }
    # Update muphi based on current phi and tau values
    muphi <- Updatemuphi(phi = phi, tau = tau, lambda = lambda, A = A, B = B)


    # Update y values (the calendar ages of the objects) by slice sampling
    for(k in 1:nobs) {
      theta[k] <- SliceSample(TARGET = thetaloglikfast, x0 = theta[k],
                              w = w, m = m, type = "log",
                              prmean = phi[c[k]], prsig = 1/sqrt(tau[c[k]]),
                              c14obs = x[k], c14sig = xsig[k],
                              mucalallyr = mucalallyr, sigcalallyr = sigcalallyr)
    }

    # Update alpha
    alpha <- switch(alphapr,
                    lognorm = Updatealphalognormpr(c, alpha, mualpha = mualpha, sigalpha = sigalpha),
                    gamma = Updatealphagammapr(c, alpha, prshape = alphaprshape, prrate = alphaprrate),
                    stop("Unknown form for prior on gamma"))

    # Store thinned output
    if(iter %% nthin == 0) {
      outi <- outi+1
      phiout[[outi]] <- phi
      tauout[[outi]] <- tau
      cout[outi,] <- c
      thetaout[outi,] <- theta
      alphaout[outi] <- alpha
      muphiout[outi] <- muphi
    }
  }
  retlist <- list(c = cout, phi = phiout, tau = tauout, theta = thetaout, alpha = alphaout, muphi = muphiout)

  if(showprogress) close(pb)
  return(retlist)
}

# Find mixture density
# Function which where you pass it a set of means, sds and weights and it returns the density of the corresponding mixture of normals
# Arguments:
# x - vector of values at which to evaluate mixture density
# w - vector of weights
# mu - the means
# sd - the sds
MixtureDens <- function(x, w, mu, sd) {
  DTemp <- mapply(function(mu, sig, w, x) w * dnorm(x, mean = mu, sd = sig), mu, sd, w, MoreArgs = list(x = x))
  apply(DTemp, 1, sum) # Sum up the various mixtures
}

# The predictive for a new observation is a scaled t-distribution
PredNewDens <- function(x, w, muphi, lambda, nu1, nu2) {
  w *exp(Logmargnormgamma(x, muphi, lambda, nu1, nu2))
}

NealFindpred <- function(x, c, phi, tau, alpha, muphi, lambda, nu1, nu2) {
  # Find the nci - number of current observations in each cluster
  nc <- length(phi)
  nci <- apply(as.row(1:nc), 2, function(x, c) sum(c == x), c = c)
  # Find probability that new theta[i+1] lies in a particular cluster or a new one
  pci <- c(nci, alpha) # Could form new cluster
  pci <- pci/sum(pci) # Normalise
  # Find density from existing clusters
  dens <- MixtureDens(x, w = pci[1:nc], mu = phi, sd = 1/sqrt(tau)) + PredNewDens(x, w = pci[nc+1], muphi, lambda, nu1, nu2)
}

