test_that("single determination probabilities same length as calibration data, c14 inputs", {
  test_calibration_data <- data.frame(
    calendar_age_BP = c(3, 4, 6, 5, 9, 8, 7),
    c14_age = c(2, 3, 4, 5, 6, 7, 8),
    c14_sig = c(0.2, 0.1, 0.3, 0.3, 0.4, 0.5, 0.3)
  )
  cali_result <- CalibrateSingleDetermination(4, 0.1, test_calibration_data, FALSE)
  probabilities <- cali_result$probability
  expect_equal(length(probabilities), length(test_calibration_data$c14_age))
})

test_that("single determination probabilities sum to one, c14 inputs", {
  test_calibration_data <- data.frame(
    calendar_age_BP = c(3, 4, 6, 5, 9, 8, 7),
    c14_age = c(2, 3, 4, 5, 6, 7, 8),
    c14_sig = c(0.2, 0.1, 0.3, 0.3, 0.4, 0.5, 0.3)
  )
  cali_result <- CalibrateSingleDetermination(4, 0.1, test_calibration_data, FALSE)
  probabilities <- cali_result$probability
  expect_equal(sum(probabilities), 1)
})

test_that("single calibration probabilities work with intcal20 data, c14 inputs", {
  expect_silent(CalibrateSingleDetermination(51020, 35, intcal20, FALSE))
})

test_that("single determination probabilities same length as calibration data, f14c inputs", {
  test_calibration_data <- data.frame(
    calendar_age_BP = c(3, 4, 6, 5, 9, 8, 7),
    f14c = c(2, 3, 4, 5, 6, 7, 8),
    f14c_sig = c(0.2, 0.1, 0.3, 0.3, 0.4, 0.5, 0.3)
  )
  cali_result <- CalibrateSingleDetermination(4, 0.1, test_calibration_data, TRUE)
  probabilities <- cali_result$probability
  expect_equal(length(probabilities), length(test_calibration_data$f14c))
})

test_that("single determination probabilities sum to one, f14c inputs", {
  test_calibration_data <- data.frame(
    calendar_age_BP = c(3, 4, 6, 5, 9, 8, 7),
    f14c = c(2, 3, 4, 5, 6, 7, 8),
    f14c_sig = c(0.2, 0.1, 0.3, 0.3, 0.4, 0.5, 0.3)
  )
  cali_result <- CalibrateSingleDetermination(4, 0.1, test_calibration_data, TRUE)
  probabilities <- cali_result$probability
  expect_equal(sum(probabilities), 1)
})

test_that("single calibration probabilities work with intcal20 data, f14c inputs", {
  expect_silent(CalibrateSingleDetermination(0.54, 0.002, intcal20, TRUE))
})
