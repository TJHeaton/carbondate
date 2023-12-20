# Small function which performs sampling as we want it to
# We need to take care with using the sample command as sometimes we pass a single integer j.
# If we use sample() then we will draw from 1:j which is not what we want
# This resample function will stop this happening
resample <- function(x, size, ...)
{
  if(length(x) <= 1) {
    if(!missing(size) && size == 0) x[FALSE]
    else x
  }
  else sample(x, size, ...)
}
