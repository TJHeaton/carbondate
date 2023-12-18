# Find prior on theta for given Poisson Process
LegacyFindCalPrior <- function(s, h, t) {
  NumChange <- length(s) - 2
  # Prior is proportional to rate at given calendar age
  priort <- approx(x = s[1:(NumChange+2)], y = c(h[1:(NumChange+1)], 0), xout = t, method = "constant")$y
  return(priort)
}

