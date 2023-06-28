test_that("SPD estimate works with intcal20 data, c14 inputs", {
  expect_silent(FindSummedProbabilityDistribution(
    calendar_age_range_BP=c(600, 1600),
    rc_determinations=c(602, 805, 1554),
    rc_sigmas=c(35, 34, 45),
    F14C_inputs=FALSE,
    calibration_curve=intcal20))
})


test_that("SPD probabilities sum to one, c14 inputs", {
  SPD = FindSummedProbabilityDistribution(
    calendar_age_range_BP=c(600, 1600),
    rc_determinations=c(602, 805, 1554),
    rc_sigmas=c(35, 34, 45),
    F14C_inputs=FALSE,
    calibration_curve=intcal20)

  expect_equal(sum(SPD$probability), 1)
})


test_that("SPD estimate works with intcal20 data, f14c inputs", {
  expect_silent(FindSummedProbabilityDistribution(
    calendar_age_range_BP=c(600, 16000),
    rc_determinations=c(0.602, 0.805, 0.1554),
    rc_sigmas=c(35, 34, 45),
    F14C_inputs=TRUE,
    calibration_curve=intcal20))
})


test_that("SPD probabilities sum to one, f14c inputs", {
  SPD = FindSummedProbabilityDistribution(
    calendar_age_range_BP=c(600, 16000),
    rc_determinations=c(0.602, 0.805, 0.1554),
    rc_sigmas=c(0.0011, 0.0013, 0.0012),
    F14C_inputs=TRUE,
    calibration_curve=intcal20)

  expect_equal(sum(SPD$probability), 1)
})
