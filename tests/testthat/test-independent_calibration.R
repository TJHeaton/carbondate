test_that("single determination probabilities same length as calibration data", {
  test_calibration_data = data.frame(
    c14_age = c(2, 3, 4, 5, 6, 7, 8),
    c14_sig = c(0.2, 0.1, 0.3, 0.3, 0.4, 0.5, 0.3)
  )
  probabilities = CalibrateSingleDetermination(4, 0.1, test_calibration_data)
  expect_equal(length(probabilities), length(test_calibration_data$c14_age))
})

test_that("single determination probabilities sum to one", {
  test_calibration_data = data.frame(
    c14_age = c(2, 3, 4, 5, 6, 7, 8),
    c14_sig = c(0.2, 0.1, 0.3, 0.3, 0.4, 0.5, 0.3)
  )
  probabilities = CalibrateSingleDetermination(4, 0.1, test_calibration_data)
  expect_equal(sum(probabilities), 1)
})

test_that("single calibration probabilities work with intcal20 data", {
  expect_silent(CalibrateSingleDetermination(51020, 35, intcal20))
})


test_that("SPD estimate works with intcal20 data", {
  expect_silent(FindSPD(
    calendar_age_range=c(600, 1600),
    c14_determinations=c(602, 805, 1554),
    c14_uncertainties=c(35, 34, 45),
    calibration_data=intcal20))
})


test_that("SPD probabilities sum to one", {
  SPD = FindSPD(
    calendar_age_range=c(600, 1600),
    c14_determinations=c(602, 805, 1554),
    c14_uncertainties=c(35, 34, 45),
    calibration_data=intcal20)

  expect_equal(sum(SPD$probability), 1)
})
