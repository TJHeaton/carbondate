test_that("WalkerBivarDirichlet gives expected result", {
  # Here we load input and output data from the functions in the original
  # nonparametric calibration code (that we know is working correctly). We run
  # the function again here and compare results to test that the function is
  # working correctly

  set.seed(14)
  load(test_path("testdata", "NPWalker_input.rda"))
  calculated_walker_temp = WalkerBivarDirichlet(
    c14_determinations=x,
    c14_uncertainties=xsig,
    calibration_curve=intcal20,
    lambda=lambda,
    nu1=nu1,
    nu2=nu2,
    A=A,
    B=B,
    alpha_shape=cprshape,
    alpha_rate=cprrate,
    n_iter=niter,
    n_thin=nthin,
    calendar_ages=inittheta,
    slice_width=slicew,
    slice_multiplier=m,
    n_clust=kstar)

  load(test_path("testdata", "NPWalker_output.rda"))
  expect_equal(calculated_walker_temp$cluster_identifiers, WalkerTemp$delta)
  expect_equal(calculated_walker_temp$alpha, WalkerTemp$c)
  expect_equal(calculated_walker_temp$n_clust, WalkerTemp$nclust)
  expect_equal(calculated_walker_temp$phi, WalkerTemp$phi)
  expect_equal(calculated_walker_temp$tau, WalkerTemp$tau)
  expect_equal(calculated_walker_temp$calendar_ages, WalkerTemp$theta)
  expect_equal(calculated_walker_temp$weight, WalkerTemp$w)
  expect_equal(calculated_walker_temp$mu_phi, WalkerTemp$muphi)

})
