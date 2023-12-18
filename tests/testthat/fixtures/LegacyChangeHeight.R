## Proposal 2: Alter the height of a randomly chosen step
# Arguments:
# th - the observed thetas
# s - the current changepoints in the rate
# h - the heights
# intrate - integral_0^L  nu(t) dt
# alpha, beta - parameters on heights
LegacyChangeHe <- function(th, s, h, intrate, alpha, beta)
{
  nh <- length(h)
  nth <- length(th)

  # Select step height to alter - will change height between s[j] and s[j+1]
  j <- sample(nh, 1) # Can use sample as just pass a integer so will pick from 1:nh

  # Propose new height of step so that log(hjnew/hjold) ~ Unif[-0.5, 0.5]
  hjold <- h[j]
  u <- runif(1, min = -0.5, max = 0.5)
  hjnew <- hjold * exp(u)

  # Store new set of heights
  hnew <- h
  hnew[j] <- hjnew

  # Find prior h ratio
  logphratio <- alpha*u - beta*(hjnew -hjold)

  # Find the likelihood of the data given the new height in the step [s_j , s_j+1]
  # Adjust the ingetrated rate to account for new height
  intratenew <- intrate + (hjnew - hjold)*(s[j+1] - s[j])

  # Find how many thetas are affected
  nuseth <- sum(th < s[j+1] & th > s[j])

  loglikold <- (nuseth * log(hjold)) - intrate # In old these have rate h[j-1] as before s[j]
  logliknew <- (nuseth * log(hjnew)) - intratenew # In new they have rate h[j]as after sjnew

  logthlikratio <- logliknew - loglikold

  # Find acceptance probability (better to use log results rather than multiply as logphratio is simple format)
  HR <- exp(logphratio + logthlikratio)

  # Determine acceptance and return result
  if(runif(1) < HR)	{
    retlist <- list(h = hnew, intrate = intratenew)	# Accept and return modified heights and intrate
  }	else {
    retlist <- list(h = h, intrate = intrate)	# Reject and return old heights and intrate
  }
  return(retlist)
}
