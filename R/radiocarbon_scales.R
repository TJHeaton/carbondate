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
