test_that("BivarGibbsDirichletwithSlice gives expected result", {
  # Here we load input and output data from the functions in the original
  # nonparametric calibration code (that we know is working correctly). We run
  # the function again here and compare results to test that the function is
  # working correctly

  set.seed(14)
  load(test_path("testdata", "NPNeal_input.rda"))
  calculated_neal_temp = BivarGibbsDirichletwithSlice(
    x = x,
    xsig = xsig,
    lambda = lambda,
    nu1 = nu1,
    nu2 = nu2,
    A = A,
    B = B,
    mualpha = mualpha,
    sigalpha = sigalpha,
    alphaprshape = alphaprshape,
    alphaprrate = alphaprrate,
    niter = niter,
    nthin = nthin,
    theta = inittheta,
    w = w,
    m = m,
    calibration_curve = intcal20,
    nclusinit = nclusinit)

  load(test_path("testdata", "NPNeal_output.rda"))
  expect_equal(calculated_neal_temp$c, NealTemp$c)
  expect_equal(calculated_neal_temp$phi, NealTemp$phi)
  expect_equal(calculated_neal_temp$tau, NealTemp$tau)
  expect_equal(calculated_neal_temp$theta, NealTemp$theta)
  expect_equal(calculated_neal_temp$alpha, NealTemp$alpha)
  expect_equal(calculated_neal_temp$muphi, NealTemp$muphi)

})
