# A wrapper which performs all of the slice sampler
# Arguments:
# TARGET - function proportional to what we want to sample from
# x0 - current x value
# type - either "raw" or "log" dependent upon if TARGET is log(f(x))
# Returns:
# x1 - updated x value
.SliceSample <- function(TARGET, x0, w, m, type = "raw", ...) {
  # Find height of slice (depends on if target is on raw or log scale)
  y <- switch(type,
              raw = stats::runif(1, max = TARGET(x0, ...)),
              log = TARGET(x0, ...) - stats::rexp(1)
  )
  # Find interval
  I <- .FindInterval(TARGET, x0 = x0, y = y, w = w, m = m, ...)
  # Sample with shrinkage and return
  .SamplefromInt(TARGET, x0 = x0, y = y, sliceint = I, ...)
}

# Find the slice interval [L,R] from which we will sample
# Arguments:
# TARGET - the function proportional to the density
# x0 - the current value of x
# y - current vertical level of slice
# w - estimate of typical slice size
# m - integer limiting slice size to mw
# Returns:
# [L,R] - the slice interval
.FindInterval <- function(TARGET, x0, y, w, m, ...) {
  # Initalise slice
  L <- x0 - w * stats::runif(1)
  R <- L + w

  # Decide how far we can extend on both sides
  J <- floor(m * stats::runif(1)) # Max steps on LHS
  K <- (m - 1) - J # Max steps on RHS

  # Now perform the stepping out
  while (J > 0 & y < TARGET(L, ...)) {
    L <- L - w
    J <- J - 1
  } # LHS stepping out

  while (K > 0 & y < TARGET(R, ...)) {
    R <- R + w
    K <- K - 1
  } # RHS stepping out

  return(sliceint = c(L, R))
}

# A function which now samples from the slice interval
# Arguments:
# TARGET - the function proportional to the density
# x0 - the current value of x
# y - current vertical level of slice
# sliceint - the slice interval in form (L,R)
.SamplefromInt <- function(TARGET, x0, y, sliceint, ...) {
  Lbar <- sliceint[1]
  Rbar <- sliceint[2]

  repeat{
    x1 <- Lbar + stats::runif(1) * (Rbar - Lbar)

    # Break loop if we have sampled satisfactorily
    if (y < TARGET(x1, ...)) {
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
