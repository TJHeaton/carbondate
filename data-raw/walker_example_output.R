set.seed(8)

library(carbondate)


walker_example_output <- WalkerBivarDirichlet(
  c14_determinations = kerr$c14_ages,
  c14_sigmas = kerr$c14_sig,
  calibration_curve=intcal20,
  n_iter = 100000,
  n_thin = 100)


## Cut off the first half of the data to reduce file size
## thin_id <- seq(nburn, n_iter, length = 100)

usethis::use_data(walker_example_output, overwrite = TRUE, compress = "bzip2")
