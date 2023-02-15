######################################################
# Functions are called by the master Gibbs functions
######################################################



# The as.row function which we will use to count number of diferent clusters
as.row <- function(x) {
  y <- NULL
  x <- as.matrix(x)
  if (nrow(x) > ncol(x)) {
    y <- t(x)
  }
  else {
    y <- x
  }
  y
}

######################################################
### Univariate Versions of Functions
######################################################

######## Updatec ###################
# Function which updates the ith element of c given current theta and phi
# Arguments:
# i - the element to update in c
# c - the vector of cluster identifiers
# phi - vector of current parameters
# theta - the given calendar ages where in a spcific cluster theta ~ N(phi, sigma^2)
# sigma - the sd on y ~ N(phi, sd = sigma)
# prmean and prsd - the mean and sd of G0 i.e. phi ~ N(prmean, prsd^2)
# alpha - parameter in stick breaking (considered fixed here although can give hyperprior)
# Returns:
# c - updated vector of cluster identifiers
# phi - updated vector of current cluster parameters
# Note: phi will need to stored as a list on return as variable length
Updatec <- function(i, c, phi, theta, sigma = 1, prmean = 0, prsd = 1, alpha = 1) {
  nc <- length(phi)
  ci <- c[i] # Cluster of element to update
  cminus <- c[-i] # the other cluster elements
  if(sum(cminus == ci) == 0) {
    phi <- phi[-ci] # Remove phi for ci if no other elements
    cminus[cminus > ci] <- cminus[cminus > ci] - 1 # Adjust labelling
    nc <- nc - 1 # Adjust n levels
  }
  nci <- apply(as.row(1:nc), 2, function(x, cminus) sum(cminus == x), cminus = cminus)
  cprob <- c(dnorm(theta, mean = phi, sd = sigma), dnorm(theta, mean = prmean, sd = sqrt(sigma^2 + prsd^2)))
  cprob <- cprob * c(nci, alpha) # weight by number in class
  class <- sample(1:(nc+1), 1, prob = cprob)
  if(class == (nc+1)) { # We have sampled a new state
    hsd <- sqrt( 1 / (1/prsd ^2 + 1/sigma^2) )
    hmean <- hsd^2 * ( prmean/(prsd^2) + theta/(sigma^2) )
    phi <- c(phi, rnorm(1, mean = hmean, sd = hsd)) # Sample new phi
  }
  # Now update the return class variables
  c[-i] <- cminus
  c[i] <- class
  retlist <- list(c = c, phi = phi)
  return(retlist)
}

############################################################
# This function will sample from the posterior of phi given a set of calendar ages from that cluster
# For each cluster in turn we have
# theta ~ N(phi, sigma^2)
# Arguments:
# theta - vector of object calendar ages (all belonging to same cluster)
# sigma - the sd on theta ~ N(phi, sd = sigma)
# prmean and prsd - the mean and sd of G0
Updatephi <- function(theta, sigma, prmean, prsd) {
  nclus <- length(theta)
  postsd <- sqrt( 1 / (1/prsd ^2 + nclus/sigma^2) )
  postmean <- postsd^2 * ( prmean/(prsd^2) + sum(theta)/(sigma^2) )
  rnorm(1, mean = postmean, sd = postsd)
}


######################################################
### Bivariate Versions of Functions - ie updates unknown cluster means and variances
######################################################



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

BivarUpdatec <- function(i, c, phi, tau, theta, lambda, nu1, nu2, muphi, alpha) {
  nc <- length(phi)
  ci <- c[i] # Cluster of element to update
  cminus <- c[-i] # the other cluster elements
  if(sum(cminus == ci) == 0) {
    phi <- phi[-ci] # Remove phi for ci if no other elements
    tau <- tau[-ci] # Remove tau for ci
    cminus[cminus > ci] <- cminus[cminus > ci] - 1 # Adjust labelling
    nc <- nc - 1 # Adjust n levels
  }
  nci <- apply(as.row(1:nc), 2, function(x, cminus) sum(cminus == x), cminus = cminus)

  cprob <- dnorm(theta, mean = phi, sd = 1/sqrt(tau)) # Likelihood of theta given phi and tau
  logmarg <- Logmargnormgamma(theta, muphi, lambda, nu1, nu2)
  cprob <- c(cprob, exp(logmarg)) # Concatenate marginal of theta for a new cluster

  cprob <- cprob * c(nci, alpha) # weight by number in class (or alpha for new cluster)
  class <- sample(1:(nc+1), 1, prob = cprob)
  if(class == (nc+1)) { # We have sampled a new state - create new phi and tau from posterior given theta
    nu1new <- nu1 + 0.5
    nu2new <- nu2 + (lambda * (theta-muphi)^2)/(2*(lambda + 1))
    lambdanew <- lambda + 1
    muphinew <- (lambda*muphi + theta)/(lambda+1)
    taunew <- rgamma(1, shape = nu1new, rate = nu2new) # Sample new tau
    phinew <- rnorm(1, mean = muphinew, sd = 1/sqrt(lambdanew*taunew)) # and then phi
    phi <- c(phi, phinew) # Create updated phi
    tau <- c(tau, taunew) # and tau
  }
  # Now update the return class variables
  c[-i] <- cminus
  c[i] <- class
  retlist <- list(c = c, phi = phi, tau = tau)
  return(retlist)
}

#####

# Function which works out the marginal of theta when theta ~ N(phi, sd = sqrt(1/tau)) and (phi,tau) are NormalGamma
Logmargnormgamma <- function(theta, muphi, lambda, nu1, nu2) {
  margprec <- (nu1 * lambda)/(nu2 * (lambda +1))
  margdf <- 2 * nu1

  A <- lgamma((margdf + 1)/2) - lgamma(margdf/2)
  B <- 0.5 * (log(margprec) - log(margdf) - log(pi))
  C <- -((margdf + 1)/2) * log(1 + (margprec*(theta-muphi)^2)/margdf)
  logden <- A + B + C
}


############################################################
# This function will sample from the posterior of (phi, tau) given a set of calendar ages from that cluster
# For each cluster in turn we have
# theta ~ N(phi, sd = 1/sqrt(tau))
# Arguments:
# theta - vector of object calendar ages (all belonging to same cluster)
# muphi, lambda, nu1, nu2 - parameters in NormalGamma for (phi, tau)
BivarUpdatephitau <- function(theta, muphi, lambda, nu1, nu2) {
  nclus <- length(theta)
  thetabar <- mean(theta)
  s <- mean((theta-thetabar)^2)
  # Update parameters according to conjugate prior
  nu1new <- nu1 + nclus/2
  nu2new <- nu2 + 0.5 * ((nclus*s) + (lambda * nclus * ((thetabar-muphi)^2)/(lambda + nclus)))
  lambdanew <- lambda + nclus
  muphinew <- (lambda*muphi + nclus*thetabar)/(lambda+nclus)
  # Now sample new values for precision tau and mean phi
  taunew <- rgamma(1, shape = nu1new, rate = nu2new) # Sample new tau
  phinew <- rnorm(1, mean = muphinew, sd = 1/sqrt(lambdanew*taunew)) # and then phi

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
Updatemuphi <- function(phi, tau, lambda, A, B) {

  locprec <- tau * lambda # We then observe phi_i ~ N(muphi, prec = lambda*tau_i)

  # Manual adjustment to get rid of tiny taus which seem to be causing a rounding problem
  # Effectively remove these clusters
  # locprec[1/sqrt(locprec) > 50000] <- 0

  postmean <- (A*B + sum(locprec*phi)) / (B + sum(locprec))
  postprec <- B + sum(locprec)
  rnorm(1, mean = postmean, sd = 1/sqrt(postprec))
}
### May need to write a univariate function if we want that to also allow updating of prior mean


#####################################################################
### Potentially shared functions used by both univariate and bivariate
#####################################################################


############################################################
# This function will perform MH to update calendar age of the objects
# given radiocarbon determination and cluster mean
# Arguments:
# theta - current cal age estimate
# phi - mean of cluster theta is thought to come from
# sigma - sd of cluster
# x - radiocarbon determination
# xsig - sd of radicoarbon determination
# calmu - posterior mean of calibration curve
# caltheta - corresponding calendar ages
# calsig - posterior sd of calibration curve
# propsd - the proposal sd for new calendar age
Updatecalage <- function(theta, phi, sigma, x, xsig, calmu, caltheta, calsig, propsd) {
  # Propose a new calendar age
  thetanew <- rnorm(1, mean = theta, sd = propsd)

  # Bit of code which finds mu(theta) and sdcal(theta) by interpolation
  calold <- FindCal(theta, calmu, caltheta, calsig)
  calnew <- FindCal(thetanew, calmu, caltheta, calsig)

  # Find log-likelihoods
  loglikold <- dnorm(theta, mean = phi, sd = sigma, log = TRUE) + dnorm(x, mean = calold$mu, sd = sqrt(calold$sigma^2 + xsig^2), log = TRUE)
  logliknew <- dnorm(thetanew, mean = phi, sd = sigma, log = TRUE) + dnorm(x, mean = calnew$mu, sd = sqrt(calnew$sigma^2 + xsig^2), log = TRUE)

  # Find HR
  HR <- exp(logliknew - loglikold)
  # Decide whether to accept or reject
  if(runif(1) < HR) return(thetanew) # Accept new thetanew
  return(theta) # Else reject and return theta
}


############################################################
# This function will perform a SLICE SAMPLER to
# to update calendar age of the objects
# given radiocarbon determination and cluster mean
# Arguments:
# theta0 - current cal age estimate
# phi - mean of cluster theta is thought to come from
# sigma - sd of cluster
# x - radiocarbon determination
# xsig - sd of radicoarbon determination
# calmu - posterior mean of calibration curve
# caltheta - corresponding calendar ages
# calsig - posterior sd of calibration curve
## Slice sampling parameters
# w - estimate of typical slice size
# m - integer limiting slice size to mw
SliceUpdatecalage <- function(theta0, phi, sigma, x, xsig,
                              calmu, caltheta, calsig,
                              w, m) {

  # Now perform slice update using this function
  SliceSample(TARGET = thetaloglik, x0 = theta0,
              w = w, m = m, type = "log",
              phi = phi, sigma = sigma, x = x, xsig = xsig,
              calmu = calmu, caltheta = caltheta, calsig = calsig)
}


# Create a function which works out the log-likelihood
# (i.e. log-target for any z (theta) given other parameters

thetaloglik <- function(z, prmean, prsig, c14obs, c14sig, calmu, caltheta, calsig) {
  # Return -Inf if outside range
  if(z < 0 | z > 50000) return(-Inf)

  # Find mu(theta) and sdcal(theta) by interpolation
  calvals <- FindCal(z, calmu, caltheta, calsig)
  loglik <- dnorm(z, mean = prmean, sd = prsig, log = TRUE) +
    dnorm(c14obs, mean = calvals$mu, sd = sqrt(calvals$sigma^2 + c14sig^2), log = TRUE)
  return(loglik)
}

# A fast version which assume IntCal is constant on a year scale
thetaloglikfast <- function(z, prmean, prsig, c14obs, c14sig, mucalallyr, sigcalallyr) {
  yr <- floor(z) # Assume the cal likelihood is constant on year scale
  # Return -Inf if outside range
  if(yr < 1 | yr > 50000) return(-Inf)

  mucal <- mucalallyr[yr]
  sigcal <- sigcalallyr[yr]
  loglik <- dnorm(z, mean = prmean, sd = prsig, log = TRUE) +
    dnorm(c14obs, mean = mucal, sd = sqrt(sigcal^2 + c14sig^2), log = TRUE)
  return(loglik)
}



# This function will interpolate a value for the calibration curve
# at age theta based upon given grid
# Arguments:
# theta - the age at which you want to interpolate
# caltheta - the grid of ages for the calibration curve
# calmu - the posterior mean for the calibration curve
# calsig - the posterior sd for the calibration curve
FindCal <- function(theta, calmu, caltheta, calsig) {
  mutheta <- approx(caltheta, calmu, xout = theta, rule = 2)$y
  musigma <- approx(caltheta, calsig, xout = theta, rule = 2)$y
  retlist <- list(mu = mutheta, sigma = musigma)
  return(retlist)
}




# This function will work out the likelihood of a particular Dir(alpha)
# given a partition c
# Arguments are:
# c - vector of the classes of each observation
# alpha - parameter in Dir(alpha)
# Return the likelihood
LogLikalpha <- function(c, alpha) {
  n <- length(c)
  nc <- max(c)
  nci <- apply(as.row(1:nc), 2, function(x, c) sum(c == x), c = c)
#  cat("nci = ", nci, "\n")
#  lik <- ((alpha^nc)*prod(factorial(nci-1))) / prod(alpha:(alpha+n-1)) # Need to work with log-likelihood
  loglik <- nc*log(alpha) + sum(sapply(pmax(nci-1, 1), function(x) sum(log(1:x)))) - sum(log(alpha:(alpha+n-1)))
  # Note we have to use pmax(nci-1, 1) here to account for clusters of size 1 have 0! = 1
  return(loglik)
}




# This function will update DP process parameter alpha by MH
# Arguments are:
# c - vector of the classes of each observation
# alpha - parameter in Dir(alpha)
# mualpha - mean of log-normal
# sigalpha - sd of lognormal
Updatealphalognormpr <- function(c, alpha, mualpha = -3, sigalpha = 1, propsd = 1) {
  uold <- log(alpha)
  unew <- rnorm(1, uold, propsd)
  alphanew <- exp(unew)
  logprrat <- dnorm(unew, mean = mualpha, sd = sigalpha, log = TRUE) - dnorm(uold, mean = mualpha, sd = sigalpha, log = TRUE)
  loglikrat <- LogLikalpha(c, alphanew) - LogLikalpha(c, alpha)
  HR <- exp(logprrat + loglikrat)
  if(is.na(HR)) {
    cat(logprrat, "    ", loglikrat, "\n")
    cat(alpha, "     ", alphanew, "\n")
    cat(c, "\n")
    stop("Not a number")
  }

  if(runif(1) < HR) {
    return(alphanew) # Accept alphanew
  }
  return(alpha) # Or reject and keep alpha
}

# Function as above but using a gamma prior on the value of alpha
Updatealphagammapr <- function(c, alpha, prshape = 0.5, prrate = 1, propsd = 1) {
# Sample new alpha from truncated normal distribution
  repeat {
    alphanew <- rnorm(1, alpha, propsd)
    if (alphanew>0)
      break
  }
  logprrat <- dgamma(alphanew, shape = prshape, rate = prrate, log = TRUE) - dgamma(alpha, shape = prshape, rate = prrate, log = TRUE)
  loglikrat <- LogLikalpha(c, alphanew) - LogLikalpha(c, alpha)
  logproprat <- pnorm(alpha/propsd, log.p = TRUE) - pnorm(alphanew/propsd, log.p = TRUE) # Adjust for non-symmetric truncated normal proposal
  HR <- exp(logprrat + loglikrat + logproprat)
  if(runif(1) < HR) {
    return(alphanew) # Accept alphanew
  }
  return(alpha) # Or reject and keep alpha
}

# Function which updates the Dirichlet process parameter c
# Does this via MH with truncated proposal distribution
# Arguments are:
# delta - vector of the classes of each observation
# c - current c parameter in Dir(c)
# prshape and prrate - prior on c ~ Gamma(shape, rate)
# propsd - proposal sd for new c
WalkerUpdatec <- function(delta, c, prshape = 0.5, prrate = 1, propsd = 1) {
  # Sample new c from truncated normal distribution
  repeat {
    cnew <- rnorm(1, c, propsd)
    if (cnew>0)
      break
  }
  # Find distinct number of populated clusters
  d <- length(unique(delta))
  n <- length(delta)
  # Now find the likelihood and prior
  logprrat <- dgamma(cnew, shape = prshape, rate = prrate, log = TRUE) - dgamma(c, shape = prshape, rate = prrate, log = TRUE)
  loglikrat <- WalkerLogLikc(d = d, c = cnew, n = n) - WalkerLogLikc(d = d, c = c, n = n)
  logproprat <- pnorm(c/propsd, log.p = TRUE) - pnorm(cnew/propsd, log.p = TRUE) # Adjust for non-symmetric truncated normal proposal
  HR <- exp(logprrat + loglikrat + logproprat)
  if(runif(1) < HR) {
    return(cnew) # Accept cnew
  }
  return(c) # Or reject and keep c
}

WalkerLogLikc <- function(d, c, n) {
  d*log(c) + lgamma(c) - lgamma(c+n)
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

WalkerFindpred <- function(x, w, phi, tau, muphi, lambda, nu1, nu2) {
  pnewclust <- 1 - sum(w)
  # Find density from existing clusters
  MixtureDens(x, w = w, mu = phi, sd = 1/sqrt(tau)) + PredNewDens(x, w = pnewclust, muphi, lambda, nu1, nu2)
}


# Function which performs one iteration of Walker DP update on (w, phi, tau)
# Input arguments:
# theta- the current observations theta_i ~ N(phi_delta[i], prec = tau_delta[i])
# w - the current weights
# v - the current vs
# delta - the current cluster allocations
# phi - the current cluster means
# tau - the current cluster precisions
# kstar - current kstar
# c, muphi, lambda, nu1, nu2 - the current DP parameters (not updated here)
DPWalkerUpdate <- function(theta, w, v, delta, phi, tau, kstar, c, muphi, lambda, nu1, nu2) {
  # Create relevant variables
  n <- length(theta)

  # Create auxiliary u variables
  u <- runif(n, min = 0, max = w[delta])
  ustar <- min(u)


  # Now update the weights
  wnew <- c()
  j <- 1
  brprod <- 1 # This is the current product of (1-v[1])...(1-v[j-1])
  # Need to work out alpha and beta internally

  # Iteratively update the weights until we have all we need
  # To update v[j] we do the following
  while(sum(wnew) < 1 - ustar) {
    if(j <= kstar) { # We have to work out alpha and beta and sample a new v[j] by inverse cdf
      # Find the indices we need to search over i.e. which are in cluster j and which are above cluster j
      clustj <- which(delta == j)
      aboveclj <- which(delta > j)

      # Find alpha_j
      if(length(clustj) == 0) {
        alpha <- 0
      } else {
        alpha <- max(u[clustj]/brprod)
      }

      # Find beta_j
      if(length(aboveclj) == 0) {
        beta <- 1
      } else {
        # Find beta_j
        prodvtemp <- cumprod(1-v) / (1-v[j]) # Vector to aid denominator. The ith entry is Prod_{l < i, l != j} (1- v_l)
        beta <- 1 - max(u[aboveclj] / (v[delta[aboveclj]] * prodvtemp[delta[aboveclj]-1]))
      }

      # Now sample from the correct posterior by inversion
      A <- (1 - alpha)^c
      B <- A - (1 - beta)^c
      utemp <- runif(1)
      v[j] <- 1 - (A - B * utemp)^(1/c)
    } else { # We are in the case that j > k^star (for current kstar) and v[j] is just from the prior beta
      v[j] <- rbeta(1, 1, c)
    }

    # We now have v so we need to extend w (and update brprod for next iteration)
    wnew <- c(wnew, brprod * v[j] )
    brprod <- brprod * (1- v[j])

    if(sum(is.na(wnew)) != 0) {
      cat("c = ",c, "\n")
      cat("j < k^star is", j <= kstar, "\n")
      cat("v[j] = ", v[j], "\n")
      browser()
      stop("Weird weights")
    }
    j <- j+1
  }

  # Now update kstar (the number of weights we have) and truncate w and v at the correct values
  kstar <- length(wnew)
  w <- wnew
  v <- v[1:kstar]

  # Now update the cluster means and precisions (note we have to introduce new ones for the new states without observations)
  for(i in 1:kstar) {
    # Find which observations belong to this cluster
    clusti <- which(delta == i)
    if(length(clusti)  == 0) { # No observations in this cluster so sample from the prior
      tau[i] <- rgamma(1, shape = nu1, rate = nu2)
      phi[i] <- rnorm(1, mean = muphi, sd = 1/sqrt(lambda*tau[i])) # sampnewphi(mu = muphi, sigma = 1/sqrt(lambda*tau)) # rnorm(1, mean = muphi, sd = 1/sqrt(lambda*tau))
    } else { # There are some observations we need to update the phi and tau for this cluster (conjugate according to NormalGamma prior)
      GibbsParams <- BivarUpdatephitau(theta[clusti], muphi = muphi, lambda = lambda, nu1 = nu1, nu2 = nu2)
      phi[i] <- GibbsParams$phi
      tau[i] <- GibbsParams$tau
    }
  }
  # Only store the values we need
  phi <- phi[1:kstar]
  tau <- tau[1:kstar]

  # Now update the allocations for each observation by sampling from them
  for(i in 1:n) {
    possid <- which(w > u[i])
    dens <- dnorm(phi[possid], mean = theta[i], sd = 1/sqrt(tau[possid]))
    dens[is.na(dens)] <- 0 # Fudge to remove erroneous problems
    delta[i] <- possid[sample.int(length(possid), 1, prob = dens)]
  }

  # Return w, v, delta, phi and tau
  retlist  <- list(w = w, v = v, delta = delta, phi = phi, tau = tau, kstar = kstar)
  return(retlist)
}

# Fudge function which makes sure when we sample a new phi it has to lie within [0,50000]
sampnewphi <- function(mu, sigma, range = c(0,50000)) {
  inrange <- FALSE
  while(!inrange) {
    newphi <- rnorm(1, mu, sigma)
    inrange <- (newphi < max(range)) & (newphi > min(range))
  }
  newphi
}


# Function which plots the individual posterior calendar ages for the Walker/Neal DP
# Arguments:
# MCMCoutput - the Walker/Neal output
# ident - the determination you want to show the individual posrterior cal age for
# y, er - vector of c14ages and c14sigs for all observations
# calcurve - the calibraiton curve you are using
# nburn and nbreaks - choices of the burn in length and the resolution of the histogram
# Output:
# A plot of the histogram of the posterior calendar age
plotindpost <- function(MCMCoutput, ident, y, er, calcurve, nburn = NA, resolution = 5) {

  resolution = ceiling(max(resolution, 1))
  par(mfrow = c(1,1))
  theta <- MCMCoutput$theta[, ident]
  c14age <- y[ident]
  c14sig <- er[ident]

  if(is.na(nburn)) {
    cat("Setting burn in as half MCMC chain\n")
    nburn <- floor(length(theta)/2)
  }
  if(nburn > length(theta)) stop("Burn in is longer than MCMC chain")
  theta <- theta[-(1:nburn)]

  # Choose the number of breaks
  nbreaks <- min(100, floor(length(theta)/10))

  # Find the calendar age range to plot
  xrange <- range(theta)
  xrange[1] = floor(xrange[1])
  while (xrange[1] %% resolution != 0) xrange[1] = xrange[1] - 1

  xrange[2] = ceiling(xrange[2])
  while (xrange[2] %% resolution != 0) xrange[2] = xrange[2] + 1

  temp <- c(which.min(abs(calcurve$calage - xrange[1])),which.min(abs(calcurve$calage - xrange[2])))
  idrange <- temp[1]:temp[2]
  calcurve$ub<- calcurve$c14age + 1.96*calcurve$c14sig
  calcurve$lb <- calcurve$c14age - 1.96*calcurve$c14sig
  yrange<- range(calcurve$ub[idrange], calcurve$lb[idrange])

  plot(calcurve$calage, calcurve$c14age, col = "blue",
       ylim = yrange, xlim = rev(xrange),
       xlab = "Calendar Age (cal yr BP)", ylab = expression(paste(""^14, "C", " age (", ""^14, "C yr BP)")),
       type = "l", main = substitute(paste("Posterior of ", i^th, " determination ", c14age, "\u00B1", c14sig, ""^14,"C yr BP"),
                                     list(i = ident, c14age = c14age,c14sig = c14sig)),
       xaxs = "i", yaxs = "i")
  lines(calcurve$calage, calcurve$ub, lty = 2, col = "blue")
  lines(calcurve$calage, calcurve$lb, lty = 2, col = "blue")

  # Plot the 14C determination on the y-axis
  yfromto <- seq(c14age - 3 * c14sig, c14age + 3 * c14sig, by = 1)
  radpol <- cbind( c(0, dnorm(yfromto, mean = c14age, sd = c14sig), 0),
                   c(min(yfromto), yfromto, max(yfromto)))
  radpol[,1] <- radpol[,1] * 0.1 * (xrange[2] - xrange[1]) / max(radpol[,1])
  radpol[,1] <-  xrange[2] - radpol[,1]
  polygon(radpol, col = rgb(1,0,0,.5))

  # Plot the posterior cal age on the x-axis
  par(new = TRUE, las = 1)
  # Create hist but do not plot - works out senssible ylim and allows to work out nbreaks for desired resolution
  breaks <-seq(xrange[1], xrange[2], by=resolution)
  temphist <- hist(theta, breaks = breaks, plot = FALSE)
  finalhist <-hist(theta, prob = TRUE, breaks = breaks, xlim = rev(xrange),
       axes = FALSE, xlab = NA, ylab = NA, main = "", xaxs = "i", ylim = c(0,2.5 * max(temphist$density)))
  return(finalhist)
}



