test_that("WalkerBivarDirichlet gives expected result", {
  # Here we load input and output data from the functions provided in the original
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
    cprshape=cprshape,
    cprrate=cprrate,
    niter=niter,
    nthin=nthin,
    theta=inittheta,
    slicew=slicew,
    m=m,
    kstar=kstar,
    showprogress = TRUE)

  load(test_path("testdata", "NPWalker_output.rda"))
  expect_equal(calculated_walker_temp, WalkerTemp)
})
