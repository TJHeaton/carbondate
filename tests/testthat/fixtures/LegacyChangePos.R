## Proposal 1: Alter the position of a randomly chosen internal changepoint
# Arguments:
# th - the observed thetas
# s - the current changepoints in the rate
# h - the heights
# intrate - integral_0^L  nu(t) dt
LegacyChangePos <- function(th, s, h, intrate)
{
  source(test_path("fixtures", "LegacyHelpers.R"))

  ns <- length(s)
  nth <- length(th)

  if(ns < 3) stop("Have proposed to move an internal changepoint when there are none")
  # Select internal changepoint at random
  j <- resample(2:(ns-1), 1) # Care as may pass single integer j and with sample() will then pick from 1:j which we do not want

  # Propose new changepoint
  sjold <- s[j]
  sjnew <- runif(1, min = s[j-1], max = s[j+1])
  # Store new set of changepoints
  snew <- s
  snew[j] <- sjnew

  # Find prior s ratio
  psratio <- ((s[j+1] - sjnew)*(sjnew - s[j-1])) / ((s[j+1] - sjold)*(sjold - s[j-1]))

  # Find the likelihood of the theta data given both sets of changepoints (only changes in period )

  # Find intratenew by adjusting previous intrate version
  adjint <- (h[j] - h[j-1]) * (s[j] - sjnew) # Diff in changepoint location * Diff in height
  intratenew <- intrate + adjint

  # Find which ths will contribute to the likelihood ratio and find their log-likelihood
  if(sjnew < sjold) { # Have shifted sjnew towards t = 0
    minsj <- sjnew
    maxsj <- sjold
    nuseth <- sum(th > minsj & th < maxsj) # Number of thetas that have been affected
    loglikold <- (nuseth * log(h[j-1])) - intrate # In old these have rate h[j-1] as before s[j]
    logliknew <- (nuseth * log(h[j])) - intratenew # In new they have rate h[j]as after sjnew
  } else { # Have shifted sjnew away from t = 0
    minsj <- sjold
    maxsj <- sjnew
    nuseth <- sum(th > minsj & th < maxsj) # Number of thetas that have been affected
    loglikold <- (nuseth * log(h[j])) - intrate # In old these have rate h[j] as after s[j]
    logliknew <- (nuseth * log(h[j-1])) - intratenew # In new they have rate h[j-1] as before sjnew
  }

  thlikratio <- exp(logliknew - loglikold)

  # Find acceptance probability
  HR <- thlikratio * psratio

  # Determine acceptance and return result
  if(runif(1) < HR) {
    retlist <- list(s = snew, intrate = intratenew)	# Accept and return modified changepoints and intrate
  } else {
    retlist <- list(s = s, intrate = intrate) # Else reject and return old changepoints and intrate
  }
  return(retlist)
}

