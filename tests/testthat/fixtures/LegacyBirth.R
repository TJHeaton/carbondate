## Proposal 3: Give birth to a new changepoint
# Arguments:
# th - the observed thetas
# s - the current changepoints in the rate
# h - the heights
# intrate - integral_0^L  nu(t) dt
# alpha, beta - parameters on heights
LegacyBirth <- function(th, s, h, intrate, alpha, beta, lambda, propratio)
{
  ns <- length(s)
  k <- ns - 2 # dimension of space i.e. number of internal changepoints

  # Propose new changepoint
  sstar <- runif(1, min = s[1], max = s[ns])

  # Find which interval it falls into
  j <- max(which(s < sstar)) # Don't need to worry about boundary

  # Current height between s_j and s_{j+1}
  hjold <- h[j]

  # Sample u = U[0,1] and adjust heights
  u <- runif(1)
  h1new <- hjold * (u/(1-u))^((s[j+1]-sstar)/(s[j+1] - s[j]))
  h2new <- h1new * (1-u) / u

  # Create new changepoint and height vector in sorted order
  snew <- c(s[1:j], sstar, s[(j+1):ns]) # No care as no change to first/last element i.e. j = 1, ns -1
  hnew <- append(h[-j], c(h1new, h2new), after = j-1) # Care as could change first or last elements

  # Find the prior ratio for dimension
  logpkratio <- dpois(k+1, lambda, log = TRUE) - dpois(k, lambda, log = TRUE)

  # Find the prior ratio for changepoint spacing
  psratio <- (2 * (k+1) * (2*k + 3) / (s[ns] - s[1])^2 )  * (sstar - s[j]) * (s[j+1] - sstar)/(s[j+1] - s[j])

  # Find the prior ratio for the heights NEED CARE WITH ROUNDING
  phratio <- ((beta^alpha)/gamma(alpha)) * exp( (alpha-1)*(log(h1new)+log(h2new)-log(hjold)) - beta*(h1new + h2new - hjold))

  # Jacobian
  Jac <- (h1new + h2new)^2 / hjold

  # Find the likelihood of the data
  intrateadj <- ((h1new - hjold)*(sstar - s[j])) + ((h2new - hjold)*(s[j+1] - sstar))
  intratenew <- intrate + intrateadj

  nuseth1 <- sum(th <= sstar & th > s[j])
  nuseth2 <- sum(th < s[j+1] & th > sstar)

  loglikold <- ((nuseth1 + nuseth2) * log(hjold)) - intrate # In old these all have rate hold
  logliknew <- (nuseth1 * log(h1new)) + (nuseth2 * log(h2new)) - intratenew # In new they have rate h1new or h2new dependent upon if after sstar

  logthlikratio <- logliknew - loglikold

  # Find acceptance probability
  HR <- exp(logpkratio + logthlikratio)*psratio*phratio*Jac*propratio
  #	print(c(HR, logpkratio, logmulikratio, psratio, phratio, Jac, propratio))
  # Determine acceptance and return result
  if(runif(1) < HR)	{
    retlist <- list(s = snew, h = hnew, intrate = intratenew, hastings_ratio = HR)	# Accept and return modified changepoints + heights
  } else {
    retlist <- list(s = s, h = h, intrate = intrate, hastings_ratio = HR) # Else reject and return old heights
  }
  return(retlist)
}
