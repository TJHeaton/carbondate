set.seed(8)

library(carbondate)

walker_example_output <- WalkerBivarDirichlet(
  rc_determinations = kerr$c14_age,
  rc_sigmas = kerr$c14_sig,
  F14C_inputs = FALSE,
  calibration_curve=intcal20,
  n_iter = 1e5,
  n_thin = 100)


## Cut off the first half of the data to reduce file size
n_out = length(walker_example_output$alpha)
thin_id <- seq((n_out+1)/2, n_out, by = 1)

walker_example_output$cluster_identifiers = walker_example_output$cluster_identifiers[thin_id,]
walker_example_output$alpha = walker_example_output$alpha[thin_id]
walker_example_output$n_clust = walker_example_output$n_clust[thin_id]
walker_example_output$phi = walker_example_output$phi[thin_id]
walker_example_output$tau = walker_example_output$tau[thin_id]
walker_example_output$weight = walker_example_output$weight[thin_id]
walker_example_output$calendar_ages = walker_example_output$calendar_ages[thin_id,]
walker_example_output$mu_phi = walker_example_output$mu_phi[thin_id]

usethis::use_data(walker_example_output, overwrite = TRUE, compress = "bzip2")
