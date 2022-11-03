test_that("gives expected result for walker output", {
  load(test_path("testdata", "NPWalker_output.rda"))
  load(test_path("testdata", "NPWalker_postprocessing.rda"))

  set.seed(seednum)
  walker_output = list(
    cluster_identifiers = WalkerTemp$delta,
    alpha = WalkerTemp$c,
    n_clust = WalkerTemp$nclust,
    phi = WalkerTemp$phi,
    tau = WalkerTemp$tau,
    calendar_ages = WalkerTemp$theta,
    weight = WalkerTemp$w,
    mu_phi = WalkerTemp$muphi,
    update_type = "walker",
    lambda = lambda,
    nu1 = nu1,
    nu2 = nu2)

  density_list = PlotPredictiveCalendarAgeDensity(
    c14_determinations = x,
    c14_uncertainties = xsig,
    calibration_curve = intcal20,
    output_data = walker_output,
    n_posterior_samples = npostsum)

  expect_equal(density_list[[1]]$density_mean, postden)
  expect_equal(density_list[[1]]$density_ci_lower, postdenCI[1,])
  expect_equal(density_list[[1]]$density_ci_upper, postdenCI[2,])

})


test_that("gives expected result for neal output", {
  load(test_path("testdata", "NPNeal_output.rda"))
  load(test_path("testdata", "NPNeal_postprocessing.rda"))

  set.seed(seednum)
  neal_output = list(
    cluster_identifiers = NealTemp$c,
    alpha = NealTemp$alpha,
    n_clust = c(),
    phi = NealTemp$phi,
    tau = NealTemp$tau,
    calendar_ages = NealTemp$theta,
    mu_phi = NealTemp$muphi,
    update_type = "neal",
    lambda = lambda,
    nu1 = nu1,
    nu2 = nu2)
  n_out = length(neal_output$alpha)
  for (i in 1:n_out) {
    cluster_identifiers = neal_output$cluster_identifiers[i,]
    neal_output$n_clust[i] = length(unique(cluster_identifiers))
  }

  density_list = PlotPredictiveCalendarAgeDensity(
    c14_determinations = x,
    c14_uncertainties = xsig,
    calibration_curve = intcal20,
    output_data = neal_output,
    n_posterior_samples = npostsum)

  expect_equal(density_list[[1]]$density_mean, postden)
  expect_equal(density_list[[1]]$density_ci_lower, postdenCI[1,])
  expect_equal(density_list[[1]]$density_ci_upper, postdenCI[2,])

})
