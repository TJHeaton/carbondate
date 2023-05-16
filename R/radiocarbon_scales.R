.AddF14cColumns <- function(data) {
  if (any(names(data) == "f14c") && any(names(data) == "f14c_sig")) {
    return(data)
  }
  data$f14c = exp(-data$c14_age / 8033)
  data$f14c_sig = data$f14c * data$c14_sig / 8033
  return(data)
}


.AddC14ageColumns <- function(data) {
  if (any(names(data) == "c14_age") && any(names(data) == "c14_sig")) {
    return(data)
  }
  data$c14_age =  -8033 * log(data$f14c)
  data$c14_sig = 8033 * data$f14c_sig / data$f14c
  return(data)
}


.ConvertF14cTo14Cage <- function(c14_determinations, c14_sigmas) {

  f14c <- -8033 * log(c14_determinations)
  f14c_sig <- 8033 * c14_sigmas / c14_determinations

  return(data.frame(f14c = f14c, f14c_sig = f14c_sig))
}


.Convert14CageToF14c <- function(f14c_determinations, f14c_sigmas) {

  c14_age <- exp(-f14c_determinations / 8033)
  c14_sig <- c14_age * f14c_sigmas / 8033

  return(data.frame(c14_age = c14_age, c14_sig = c14_sig))
}
