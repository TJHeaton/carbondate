## This function will work out the probabilities for each move in the RJ MCMC sampler
# Return 4 vectors in list with each move probability for number of heights h = 1,..., kmax+1
# Note we use the number of heights to avoid confusion as to the fact that sometimes nk = 0 and S does not naturally create vectors with index p[0]
# Also note that nk = nh - 1

LegacyFindMoveProb <- function(lambda, kmax)
{
  pMovePos <- rep(NA, length = kmax + 1)
  pMoveHe <- rep(NA, length = kmax + 1)
  pMoveBirth <- rep(NA, length = kmax + 1)
  pMoveDeath <- rep(NA, length = kmax + 1)

  # Fix constraints
  pMoveBirth[kmax+1] <- 0
  pMoveDeath[1] <- 0

  # Now find other probabilities of birth and death ignoring constant c
  pMoveBirth[1:kmax] <- pmin(1, dpois(1:kmax, lambda = lambda)/dpois(0:(kmax -1), lambda = lambda))
  pMoveDeath[2:(kmax+1)] <- pmin(1, dpois(0:(kmax-1), lambda = lambda)/dpois(1:kmax, lambda = lambda))

  c <- 0.9/ max(pMoveBirth+pMoveDeath)

  # Rescale to allow other moves a reasonable probability
  pMoveBirth <- c * pMoveBirth
  pMoveDeath <- c * pMoveDeath

  # Find other move probabilities
  pMovePos <- (1 - pMoveBirth -pMoveDeath)/2
  pMoveHe <- pMovePos

  # Finally deal with case that if nh = 1 (i.e. nk = 0) then cannot move position of a changepoint
  pMovePos[1] <- 0
  pMoveHe[1] <- 2 * pMoveHe[1]

  # Check that all add up to 1
  #	if(any(pMoveBirth + pMoveDeath + pMoveHe + pMovePos != 1)) print("Error")

  # Return a list with the probabilities
  list(Pos = pMovePos, He = pMoveHe, Birth = pMoveBirth, Death = pMoveDeath, sum = pMoveBirth + pMoveDeath + pMoveHe + pMovePos)
}
