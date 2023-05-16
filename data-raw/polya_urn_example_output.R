set.seed(8)

library(carbondate)

polya_urn_example_output <- PolyaUrnBivarDirichlet(
  rc_determinations = kerr$c14_age,
  rc_sigmas = kerr$c14_sig,
  F14C_inputs = FALSE,
  calibration_curve=intcal20,
  n_iter = 1e5,
  n_thin = 100)


## Cut off the first half of the data to reduce file size
n_out = length(polya_urn_example_output$alpha)
thin_id <- seq((n_out+1)/2, n_out, by = 1)

polya_urn_example_output$cluster_identifiers = polya_urn_example_output$cluster_identifiers[thin_id]
polya_urn_example_output$observations_per_cluster = polya_urn_example_output$observations_per_cluster[thin_id]
polya_urn_example_output$alpha = polya_urn_example_output$alpha[thin_id]
polya_urn_example_output$n_clust = polya_urn_example_output$n_clust[thin_id]
polya_urn_example_output$phi = polya_urn_example_output$phi[thin_id]
polya_urn_example_output$tau = polya_urn_example_output$tau[thin_id]
polya_urn_example_output$calendar_ages = polya_urn_example_output$calendar_ages[thin_id,]
polya_urn_example_output$mu_phi = polya_urn_example_output$mu_phi[thin_id]

usethis::use_data(polya_urn_example_output, overwrite = TRUE, compress = "bzip2")
