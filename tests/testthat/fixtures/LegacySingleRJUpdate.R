# This file will contain the code to perform the 4 proposal steps associated with the RJMCMC method which models the random walk variance as a series of step functions. The 4 moves are:

# 1) Alter the height of a randomly chosen step
# 2) Alter the position of a randomly chosen step
# 3) Create a new step at a randomly chosen location - Birth step
# 4) Kill one of the current steps at random - Death step

# This function takes as variables
#  th        	- calendar ages of objects (DOES NOT NEED TO BE SORTED)
#  s          - Positions of the changepoints (MUST BE SORTED)
#  h          - Values of the PP rate in each step
# intrate     - The integrated rate of the PP
# alpha, beta - shape and rate of Gamma distr for h
# lambda      - expected number of changepoints
# pMove       - calculated prob of 4 different types of MH move

LegacyRJMCMCVar <- function(th, s, h, intrate, alpha, beta, lambda, pMove)
{
  # Read in legacy individual updates
  source(test_path("fixtures", "LegacyChangePos.R"))
  source(test_path("fixtures", "LegacyChangeHeight.R"))
  source(test_path("fixtures", "LegacyBirth.R"))
  source(test_path("fixtures", "LegacyDeath.R"))
  source(test_path("fixtures", "LegacyHelpers.R"))

  ns <- length(s)
  nh <- length(h)
  k <- ns - 2
  if(nh != k + 1) {
    stop("Error in matching dimension of s and h")
  }
  u <- runif(1)
  if(u < pMove$Pos[nh]) # Propose moving the position of a changepoint
  {
    Res <- LegacyChangePos(th = th, s = s, h = h, intrate = intrate)
    s <- Res$s
    intrate <- Res$intrate
  }
  else if(u < pMove$Pos[nh] + pMove$He[nh]) # Propose moving the height of a step
  {
    Res <- LegacyChangeHe(th = th, s = s, h = h, intrate = intrate, alpha = alpha, beta = beta)
    h <- Res$h
    intrate <- Res$intrate
  }
  else if(u < pMove$Pos[nh] + pMove$He[nh] + pMove$Birth[nh]) # Propose a birth step
  {
    Res <- LegacyBirth(th = th, s = s, h = h, intrate = intrate, alpha = alpha, beta = beta,
                 lambda = lambda, propratio = (pMove$Death[nh+1] * (s[ns] - s[1])) / (pMove$Birth[nh] * (k+1)) )
    s <- Res$s
    h <- Res$h
    intrate <- Res$intrate
  }
  else # Propose a death step
  {
    Res <- LegacyDeath(th = th, s = s, h = h, intrate = intrate, alpha = alpha, beta = beta,
                 lambda = lambda, propratio = (pMove$Birth[nh-1]*k)/(pMove$Death[nh]*(s[ns] - s[1])) )
    s <- Res$s
    h <- Res$h
    intrate <- Res$intrate
  }
  list(s = s, h = h, intrate = intrate)
}
