test_that("SPD estimate works with intcal20 data", {
  expect_silent(FindSummedProbabilityDistribution(
    calendar_age_range=c(600, 1600),
    c14_determinations=c(602, 805, 1554),
    c14_sigmas=c(35, 34, 45),
    calibration_curve=intcal20))
})


test_that("SPD probabilities sum to one", {
  SPD = FindSummedProbabilityDistribution(
    calendar_age_range=c(600, 1600),
    c14_determinations=c(602, 805, 1554),
    c14_sigmas=c(35, 34, 45),
    calibration_curve=intcal20)

  expect_equal(sum(SPD$probability), 1)
})
