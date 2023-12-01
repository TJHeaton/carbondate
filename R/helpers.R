.SetNOut <- function (n_iter, n_thin) {
  n_out <- floor(n_iter / n_thin) + 1
}


.SetNBurn <- function(n_burn, n_iter, n_thin) {
  if (is.na(n_burn)) {
    n_burn <- ceiling(n_iter / (2 * n_thin))
  } else {
    n_burn <- floor(n_burn / n_thin)
  }
  return(n_burn)
}


.SetNEnd <- function(n_end, n_iter, n_thin) {
  if (is.na(n_end)) {
    n_end <- .SetNOut(n_iter, n_thin)
  } else {
    n_end <- floor(n_end / n_thin)
  }
  return(n_end)
}