# A wrapper which performs all of the slice sampler
# Arguments:
# TARGET - function proportional to what we want to sample from
# x0 - current x value
# type - either "raw" or "log" dependent upon if TARGET is log(f(x))
# Returns:
# x1 - updated x value
.SliceSample <- function(
    x0,
    w,
    m,
    prmean,
    prsig,
    c14obs,
    c14sig,
    mucalallyr,
    sigcalallyr) {
  # Find height of slice
  y <- .ThetaLogLikelihood(
    x0, prmean, prsig, c14obs, c14sig, mucalallyr, sigcalallyr) - stats::rexp(1)
  I <- .FindSliceInterval(
    x0, y, w, m, prmean, prsig, c14obs, c14sig, mucalallyr, sigcalallyr)
  sampled_value = .SampleFromSliceInterval(
    x0, y, I, prmean, prsig, c14obs, c14sig, mucalallyr, sigcalallyr)
  return(sampled_value)
}

# Find the slice interval [L,R] from which we will sample
# Arguments:
# x0 - the current value of x
# y - current vertical level of slice
# w - estimate of typical slice size
# m - integer limiting slice size to mw
# Returns:
# [L,R] - the slice interval
.FindSliceInterval <- function(
    x0, y, w, m, prmean, prsig, c14obs, c14sig, mucalallyr, sigcalallyr) {
  # Initalise slice
  L <- x0 - w * stats::runif(1)
  R <- L + w

  # Decide how far we can extend on both sides
  J <- floor(m * stats::runif(1)) # Max steps on LHS
  K <- (m - 1) - J # Max steps on RHS

  # Now perform the stepping out
  while (J > 0 & y < .ThetaLogLikelihood(
    L, prmean, prsig, c14obs, c14sig, mucalallyr, sigcalallyr)) {
    L <- L - w
    J <- J - 1
  } # LHS stepping out

  while (K > 0 & y < .ThetaLogLikelihood(
    R, prmean, prsig, c14obs, c14sig, mucalallyr, sigcalallyr)) {
    R <- R + w
    K <- K - 1
  } # RHS stepping out

  return(slice_interval = c(L, R))
}

# Arguments:
# x0 - the current value of x
# y - current vertical level of slice
# slice_interval - the slice interval in form (L,R)
.SampleFromSliceInterval <- function(
    x0, y, slice_interval, prmean, prsig, c14obs, c14sig, mucalallyr, sigcalallyr) {
  Lbar <- slice_interval[1]
  Rbar <- slice_interval[2]

  repeat{
    x1 <- Lbar + stats::runif(1) * (Rbar - Lbar)

    # Break loop if we have sampled satisfactorily
    if (y < .ThetaLogLikelihood(
      x1, prmean, prsig, c14obs, c14sig, mucalallyr, sigcalallyr)) {
      return(x1)
    }

    # Else shrink the interval
    if (x1 < x0) {
      Lbar <- x1
    } else {
      Rbar <- x1
    }
  }
}
