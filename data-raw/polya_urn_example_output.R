set.seed(8)

library(carbondate)

polya_urn_example_output <- PolyaUrnBivarDirichlet(
  rc_determinations = two_normals$c14_age,
  rc_sigmas = two_normals$c14_sig,
  F14C_inputs = FALSE,
  calibration_curve=intcal20,
  n_iter = 1e5,
  n_thin = 100)

usethis::use_data(polya_urn_example_output, overwrite = TRUE, compress = "bzip2")
