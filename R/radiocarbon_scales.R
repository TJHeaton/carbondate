.AddF14cColumns <- function(data) {
  if (any(names(data) == "f14c") && any(names(data) == "f14c_sig")) {
    return(data)
  }
  data$f14c <- exp(-data$c14_age / 8033)
  data$f14c_sig <- data$f14c * data$c14_sig / 8033
  return(data)
}


.Convert14CageToF14c <- function(c14_determinations, c14_sigmas) {

  f14c <- exp(-c14_determinations / 8033)
  f14c_sig <- f14c * c14_sigmas / 8033

  return(data.frame(f14c = f14c, f14c_sig = f14c_sig))
}


.AddC14ageColumns <- function(data) {
  if (any(names(data) == "c14_age") && any(names(data) == "c14_sig")) {
    return(data)
  }
  data$c14_age <- -8033 * log(data$f14c)
  data$c14_sig <- 8033 * data$f14c_sig / data$f14c
  return(data)
}


.ConvertF14cTo14Cage <- function(f14c_determinations, f14c_sigmas) {

  c14_age <- -8033 * log(f14c_determinations)
  c14_sig <- 8033 * f14c_sigmas / f14c_determinations

  return(data.frame(c14_age = c14_age, c14_sig = c14_sig))
}

