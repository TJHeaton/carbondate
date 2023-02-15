# Function which creates raw data
CreateDataVaryWeight <- function(nobs, nclus, DPparams, nu1, nu2, muphi, lambda, calcurve, xsig) {
  # Create some nclus clusters with means and variances
  tautrue <- rgamma(nclus, nu1, nu2)
  phitrue <- rnorm(nclus, mean = muphi, sd = 1/sqrt(lambda*tautrue))
  
  # Now create mixture weights
  if(length(DPparams) != nclus) stop("Incorrect parameters to sample from Dirichlet distribution")
  weights <- rgamma(nclus, shape = DPparams, rate = 1)
  weights <- weights/sum(weights)
  
  # Create some observed data from the clusters according to probabilities
  ctrue <- sample(1:nclus, nobs, replace = TRUE, prob = weights)
  theta <- rnorm(nobs, mean = phitrue[ctrue], sd = 1/sqrt(tautrue[ctrue]))
  
  #### Now create some radiocarbon determinations x
  
  # Interpolate calobration curve mean and sd at theta values
  calinterp <- FindCal(theta, calmu = calcurve$c14age, caltheta = calcurve$calage, calsig = calcurve$c14sig)
  
  # Sample some calibration curve values
  xcalcurve <- rnorm(nobs, calinterp$mu, calinterp$sigma)
  
  x <- rnorm(nobs, mean = xcalcurve, sd = xsig)  
  retlist <- list(x = x, thetatrue = theta, tautrue = tautrue, phitrue = phitrue, weights = weights, ctrue = ctrue)
  return(retlist)
}

# Single cluster example
# Function which creates raw data
#CreateDataSingleCluster <- function(nobs, murange, calcurve, a, g, h) {
CreateDataSingleCluster <- function(nobs, murange, calcurve, masternu1 = 2, masternu2 = 4000, xsig) {
  # Choose a value for the mean of the cluster
  phitrue <- runif(1, murange[1], murange[2])
  
  # Now estimate a precision
  #tempbeta <- rgamma(1, g, h)
  #tautrue <- rgamma(1, a, tempbeta)
  tautrue <- rgamma(1, masternu1, masternu2)
  
  # Create some observed data from the clusters according to probabilities
  weights <- 1
  ctrue <- rep(1, nobs)
  theta <- rnorm(nobs, mean = phitrue[ctrue], sd = 1/sqrt(tautrue[ctrue]))
  
  #### Now create some radiocarbon determinations x
  
  # Interpolate calobration curve mean and sd at theta values
  calinterp <- FindCal(theta, calmu = calcurve$c14age, caltheta = calcurve$calage, calsig = calcurve$c14sig)
  
  # Sample some calibration curve values
  xcalcurve <- rnorm(nobs, calinterp$mu, calinterp$sigma)
  
  x <- rnorm(nobs, mean = xcalcurve, sd = xsig)  
  retlist <- list(x = x, thetatrue = theta, tautrue = tautrue, phitrue = phitrue, weights = weights, ctrue = ctrue)
  return(retlist)
}


# Uniform distribution example
CreateSingleUniformPhaseData <- function(nobs, startrange, calcurve, phaserange, xsig) {
  # Choose a value for the start of the phase and length of the pahse 
  phasestart <- runif(1, min = startrange[1], max = startrange[2]) # Start ~ Unif[50, 14000]
  phaselength <- runif(1, min = phaserange[1], max = phaserange[2]) # Range ~ Unif[50,1000]
  
  theta <- runif(nobs, min = phasestart, max = phasestart + phaselength) # Sample nobs from the uniform  phase

  #### Now create some radiocarbon determinations x
  
  # Interpolate calobration curve mean and sd at theta values
  calinterp <- FindCal(theta, calmu = calcurve$c14age, caltheta = calcurve$calage, calsig = calcurve$c14sig)
  
  # Sample some calibration curve values
  xcalcurve <- rnorm(nobs, calinterp$mu, calinterp$sigma)
  
  x <- rnorm(nobs, mean = xcalcurve, sd = xsig)  
  retlist <- list(x = x, thetatrue = theta, phasestart = phasestart, phaselength = phaselength)
  return(retlist)
}




# A function which works out the loss E[(thetahat -  theta)^2]
# We assume we have nobs things we are calibrating and npost posterior estimates
# Arguments:
# thetaout - npost * nobs matrix of posterior estimates (each row shows a single set of estimated thetas for all objects)
# true - vector (length nobs) of the true values of theta
FindLoss <- function(thetaout, true, nburn = NA) {
  if(is.na(nburn)) nburn <- floor(dim(thetaout)[1]/2)
  thetaout <- thetaout[-(1:nburn),]
  dist <- sweep(thetaout, MARGIN = 2, true,`-`)
  l2loss <- mean(dist^2)
  l1loss <- mean(abs(dist))
  retlist <- list(l1loss = l1loss, l2loss = l2loss)
  return(retlist)
}


# To calibrate all the determinations individually then we need to:
# 1 - Find calibration mu and sigma for each calendar age in yrange
# 2 - For each radicoarbon determination find the density of a normal 

# Function which finds the calibrtion curve at any y
FindCalCurve <- function(theta, calcurve) { 
  calcurvemean <- approx(calcurve$calage, calcurve$c14age, theta, rule = 2)$y
  calcurvesd <- approx(calcurve$calage, calcurve$c14sig, theta, rule = 2)$y   
  retlist <- list(curvemean = calcurvemean, curvesd = calcurvesd)
  return(retlist)
}

# Function which calibrates a single determination x ~ N(cal(theta), xsig^2)
# Arguments:
# x and xsig the observed radiocarbon determination and uncertainty
# calmu - a vector of posterior mean of the calibration curve over range of theta
# calsig - corresponding vector of posterior sd of calibration curve over range of theta 
calibind <- function(x, xsig, calmu, calsig) {
  indprob <- dnorm(x, mean = calmu, sd = sqrt(calsig^2 + xsig^2))
  indprob <- indprob/sum(indprob)
  return(indprob)
} 

# Find the expected loss of individual calibration
FindIndLoss <- function(probs, true, yfromto)
{
  dist <- (yfromto - true)
  l2loss <- sum(probs * dist^2)
  l1loss <- sum(probs * abs(dist))
  retvec <- c(l1loss, l2loss)
  return(retvec)
}
