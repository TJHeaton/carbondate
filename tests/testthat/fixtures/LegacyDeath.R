## Proposal 4: Kill a current changepoint
LegacyDeath <- function(th, s, h, intrate, alpha, beta, lambda, propratio)
{
  source(test_path("fixtures", "LegacyHelpers.R"))

  ns <- length(s)
  k <- ns - 2
  # Select a changepoint to remove
  j <- resample(2:(ns-1), 1) # Need to take care as 2:(ns-1) can be a single integer and will so not work as desired if we pass it to sample

  # Create new height for interval s[j-1], s[j+1] via inverse of birth step
  hjnew <- exp( 1/(s[j+1] - s[j-1]) * ( (s[j] - s[j-1]) * log(h[j-1]) + (s[j+1]-s[j]) * log(h[j]) ) )

  # Store proposed changepoint and height vectors
  snew <- s[-j]
  hnew <- append(h[-c(j-1, j)], hjnew, after = j-2) # Care as could change first/last element

  # Find the prior ratio for dimension
  logpkratio <- dpois(k-1 , lambda, log = TRUE) - dpois(k, lambda, log = TRUE)

  # Find the prior ratio for changepoint spacing
  psratio <- (s[ns] - s[1])^2 / (2 * k * (2*k+1) ) * (s[j+1] - s[j-1]) / ((s[j+1] - s[j]) * (s[j] - 	s[j-1]))

  # Find the prior ratio for the heights NEED CARE WITH ROUNDING
  phratio <- (gamma(alpha)/(beta^alpha)) / exp( (alpha-1)*(log(h[j])+log(h[j-1])-log(hjnew)) - beta*(h[j] + h[j-1] - hjnew))

  # Jacobian
  Jac <- hjnew / ((h[j-1] + h[j])^2)

  # Find the likelihood of the data
  intrateadj <- ((hjnew - h[j-1])*(s[j] - s[j-1])) + ((hjnew - h[j])*(s[j+1] - s[j]))
  intratenew <- intrate + intrateadj

  nuseth1 <- sum(th <= s[j] & th > s[j-1]) # in old these will have rate h[j-1]
  nuseth2 <- sum(th < s[j+1] & th > s[j])  # in old these will have rate h[j]

  loglikold <- (nuseth1 * log(h[j-1])) + (nuseth2 * log(h[j]))  - intrate # In old they have rate h[j-1] or h[j] dependent upon if after s[j]
  logliknew <- ((nuseth1 + nuseth2) * log(hjnew)) - intratenew # In new these all have rate hnew

  logthlikratio <- logliknew - loglikold

  # Find acceptance probability
  HR <- exp(logthlikratio + logpkratio) * psratio * phratio * Jac * propratio

  # Determine acceptance and return result
  if(runif(1) < HR)	{
    retlist <- list(
      s = snew,
      h = hnew,
      intrate = intratenew,
      hastings_ratio = HR)	# Accept and return modified changepoints + heights
  } else {
    retlist <- list(
      s = s,
      h = h,
      intrate = intrate,
      hastings_ratio = HR) # Else reject and return old heights
  }
  return(retlist)
}
