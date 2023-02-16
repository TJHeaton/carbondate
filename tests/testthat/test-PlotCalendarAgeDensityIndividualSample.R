test_that("gives expected result for walker output", {
  load(test_path("fixtures", "walker_input.rda"))
  set.seed(seednum)
  # NOTE this test will fail if the test for the function below which we use to
  # set up the data fails
  walker_output = WalkerBivarDirichlet(
    c14_determinations=x,
    c14_sigmas=xsig,
    calibration_curve=intcal20,
    n_iter=niter,
    n_thin=nthin)
  load(test_path("fixtures", "walker_independent_posterior.rda"))

  density_hist = PlotCalendarAgeDensityIndividualSample(
    ident = ident,
    output_data = walker_output,
    resolution = resolution)

  expect_equal(density_hist$breaks, indpost$breaks)
  expect_equal(density_hist$counts, indpost$counts)
  expect_equal(density_hist$density, indpost$density)
})


test_that("gives expected result for neal output", {
  load(test_path("fixtures", "polya_urn_input.rda"))
  set.seed(seednum)
  # NOTE this test will fail if the test for the function below which we use to
  # set up the data fails
  neal_output = PolyaUrnBivarDirichlet(
    c14_determinations = x,
    c14_sigmas = xsig,
    calibration_curve = intcal20,
    n_iter = niter,
    n_thin = nthin)
  load(test_path("fixtures", "polya_urn_independent_posterior.rda"))

  set.seed(seednum)

  density_hist = PlotCalendarAgeDensityIndividualSample(
    ident = ident,
    output_data = neal_output,
    resolution = resolution)

  expect_equal(density_hist$breaks, indpost$breaks)
  expect_equal(density_hist$counts, indpost$counts)
  expect_equal(density_hist$density, indpost$density)
})
