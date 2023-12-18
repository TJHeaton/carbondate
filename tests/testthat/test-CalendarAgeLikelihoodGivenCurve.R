test_that("CalendarAgeLikelihoodGivenCurve gives same as legacy code", {

  rc_determination <- 3000
  rc_sigma <- 40

  set.seed(15)
  theta <- runif(n = 100, min = 0, max = 50000)

  revised_likelihood <- .CalendarAgeLikelihoodGivenCurve(
    rc_determination = rc_determination,
    rc_sigma = rc_sigma,
    theta = theta,
    F14C_inputs = FALSE,
    calibration_curve = intcal20)

  source(test_path("fixtures", "LegacyCalibrate.R"))

  legacy_cal_curve <- intcal20
  names(legacy_cal_curve)[1:3] <- c("calage", "c14age", "c14sig")

  legacy_likelihood <- LegacyCalibrate(
    x = rc_determination,
    xsig = rc_sigma,
    t = theta,
    calcurve = legacy_cal_curve
  )

  expect_identical(revised_likelihood, legacy_likelihood)

})
