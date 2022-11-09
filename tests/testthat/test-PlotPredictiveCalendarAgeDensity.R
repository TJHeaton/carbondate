test_that("gives expected result for walker output", {
  load(test_path("testdata", "NPWalker_input.rda"))
  set.seed(seednum)
  # NOTE this test will fail if the test for the function below which we use to
  # set up the data fails
  walker_output = WalkerBivarDirichlet(
    c14_determinations=x,
    c14_sigmas=xsig,
    calibration_curve=intcal20,
    n_iter=niter,
    n_thin=nthin)
  load(test_path("testdata", "NPWalker_postprocessing.rda"))

  set.seed(seednum)

  density_list = PlotPredictiveCalendarAgeDensity(
    output_data = walker_output,
    n_posterior_samples = npostsum,
    interval_width = "bespoke",
    bespoke_probability = 0.95)

  expect_equal(density_list[[1]]$density_mean, postden)
  expect_equal(density_list[[1]]$density_ci_lower, postdenCI[1,])
  expect_equal(density_list[[1]]$density_ci_upper, postdenCI[2,])

})


test_that("gives expected result for neal output", {
  load(test_path("testdata", "NPNeal_input.rda"))
  set.seed(seednum)
  # NOTE this test will fail if the test for the function below which we use to
  # set up the data fails
  neal_output = PolyaUrnBivarDirichlet(
    c14_determinations = x,
    c14_sigmas = xsig,
    calibration_curve = intcal20,
    n_iter = niter,
    n_thin = nthin)
  load(test_path("testdata", "NPNeal_postprocessing.rda"))

  set.seed(seednum)

  density_list = PlotPredictiveCalendarAgeDensity(
    output_data = neal_output,
    n_posterior_samples = npostsum,
    interval_width = "bespoke",
    bespoke_probability = 0.95)

  expect_equal(density_list[[1]]$density_mean, postden)
  expect_equal(density_list[[1]]$density_ci_lower, postdenCI[1,])
  expect_equal(density_list[[1]]$density_ci_upper, postdenCI[2,])

})
