test_that("BivarGibbsDirichletwithSlice gives expected result", {
  # Here we load input and output data from the functions in the original
  # nonparametric calibration code (that we know is working correctly). We run
  # the function again here and compare results to test that the function is
  # working correctly

  set.seed(14)
  load(test_path("testdata", "NPNeal_input.rda"))
  calculated_neal_temp = BivarGibbsDirichletwithSlice(
    c14_determinations = x,
    c14_uncertainties = xsig,
    calibration_curve = intcal20,
    lambda = lambda,
    nu1 = nu1,
    nu2 = nu2,
    A = A,
    B = B,
    alpha_shape = alphaprshape,
    alpha_rate = alphaprrate,
    n_iter = niter,
    n_thin = nthin,
    calendar_ages = inittheta,
    slice_width = w,
    slice_multiplier = m,
    n_clust = nclusinit)

  load(test_path("testdata", "NPNeal_output.rda"))
  expect_equal(calculated_neal_temp$cluster_identifiers, NealTemp$c)
  expect_equal(calculated_neal_temp$phi, NealTemp$phi)
  expect_equal(calculated_neal_temp$tau, NealTemp$tau)
  expect_equal(calculated_neal_temp$calendar_ages, NealTemp$theta)
  expect_equal(calculated_neal_temp$alpha, NealTemp$alpha)
  expect_equal(calculated_neal_temp$mu_phi, NealTemp$muphi)
  n_out = length(calculated_neal_temp$alpha)
  for (i in 1:n_out) {
    cluster_identifiers = calculated_neal_temp$cluster_identifiers[i,]
    expect_equal(
      calculated_neal_temp$n_clust[i],
      length(unique(cluster_identifiers)))
  }
})


# TODO: test with a lognorm type for alpha. This is not currently tested.
