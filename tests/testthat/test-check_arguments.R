test_that("test check calibration curve has correct headings - F14C inputs is NA - fails", {
  calibration_curve <- data.frame(calendar_age_BP = 1:3, f14c = 1:3)
  arg_check <- .makeAssertCollection()

  .CheckCalibrationCurve(arg_check, calibration_curve, NA)

  expect_equal(
    arg_check$getMessages(),
    "The calibration curve must have the columns ('f14c', 'f14c_sig') and/or the columns ('c14_age', 'c14_sig')"
  )
})


test_that("test check calibration curve has correct headings - F14C inputs is NA - passes", {
  calibration_curve <- data.frame(calendar_age_BP = 1:3, f14c = 1:3, f14c_sig = 1:3, c14_age = 1:3, c14_sig = 1:3)
  arg_check <- .makeAssertCollection()

  .CheckCalibrationCurve(arg_check, calibration_curve, NA)

  expect_true(arg_check$isEmpty())
})


test_that("test check calibration curve has correct headings - F14C inputs is TRUE - fails", {
  calibration_curve <- data.frame(calendar_age_BP = 1:3, f14c = 1:3)
  arg_check <- .makeAssertCollection()

  .CheckCalibrationCurve(arg_check, calibration_curve, TRUE)

  expect_equal(
    arg_check$getMessages(),
    "The calibration curve must have the columns ('f14c', 'f14c_sig')"
  )
})


test_that("test check calibration curve has correct headings - F14C inputs is TRUE - passes", {
  calibration_curve <- data.frame(calendar_age_BP = 1:3, f14c = 1:3, f14c_sig = 1:3)
  arg_check <- .makeAssertCollection()

  .CheckCalibrationCurve(arg_check, calibration_curve, NA)

  expect_true(arg_check$isEmpty())
})


test_that("test check calibration curve has correct headings - F14C inputs is FALSE - fails", {
  calibration_curve <- data.frame(calendar_age_BP = 1:3, c14_age = 1:3)
  arg_check <- .makeAssertCollection()

  .CheckCalibrationCurve(arg_check, calibration_curve, FALSE)

  expect_equal(
    arg_check$getMessages(),
    "The calibration curve must have the columns ('c14_age', 'c14_sig')"
  )
})


test_that("test check calibration curve has correct headings - F14C inputs is FALSE - passes", {
  calibration_curve <- data.frame(calendar_age_BP = 1:3, c14_age = 1:3, c14_sig = 1:3)
  arg_check <- .makeAssertCollection()

  .CheckCalibrationCurve(arg_check, calibration_curve, FALSE)

  expect_true(arg_check$isEmpty())
})


test_that("test check integer - fails", {
  my_char <- "a"
  my_vec <- 1:10L
  my_double <- 2.4
  my_int <- 9L
  arg_check <- .makeAssertCollection()

  .CheckInteger(arg_check, my_char)
  .CheckInteger(arg_check, my_vec)
  .CheckInteger(arg_check, my_double)
  .CheckInteger(arg_check, my_int, lower = 10)
  .CheckInteger(arg_check, my_int, upper = 5)

  expect_equal(
    arg_check$getMessages(),
    c(
      "my_char must be an integer",
      "my_vec must be an integer",
      "my_double must be an integer",
      "my_int must be more than or equal to 10",
      "my_int must be less than or equal to 5"
    )
  )
})


test_that("test check integer - passes", {
  my_integerish <- 10
  my_integer <- 5L
  arg_check <- .makeAssertCollection()

  .CheckInteger(arg_check, my_integerish)
  .CheckInteger(arg_check, my_integer)
  .CheckInteger(arg_check, my_integer, upper = 5)
  .CheckInteger(arg_check, my_integer, upper = 10)
  .CheckInteger(arg_check, my_integer, lower = 5)
  .CheckInteger(arg_check, my_integer, lower = 2)

  expect_true(arg_check$isEmpty())
})


test_that("test check number - fails", {
  my_char <- "a"
  my_vec <- 1:10
  arg_check <- .makeAssertCollection()

  .CheckNumber(arg_check, my_char)
  .CheckNumber(arg_check, my_vec)

  expect_equal(
    arg_check$getMessages(), c("my_char must be a number", "my_vec must be a number")
  )
})


test_that("test check number - passes", {
  my_num <- 3
  arg_check <- .makeAssertCollection()

  .CheckNumber(arg_check, my_num)

  expect_true(arg_check$isEmpty())
})


test_that("test check flag - fails", {
  my_char <- "a"
  my_num <- 3
  my_vec <- c(TRUE, FALSE)
  arg_check <- .makeAssertCollection()

  .CheckFlag(arg_check, my_char)
  .CheckFlag(arg_check, my_num)
  .CheckFlag(arg_check, my_vec)

  expect_equal(
    arg_check$getMessages(),
    c(
        "my_char must be a single logical value (TRUE, FALSE OR NA)",
        "my_num must be a single logical value (TRUE, FALSE OR NA)",
        "my_vec must be a single logical value (TRUE, FALSE OR NA)"
    )
  )
})


test_that("test check flag - passes", {
  arg_check <- .makeAssertCollection()

  .CheckFlag(arg_check, TRUE)
  .CheckFlag(arg_check, FALSE)
  .CheckFlag(arg_check, NA)

  expect_true(arg_check$isEmpty())
})


test_that("test number array - fails", {
  my_char <- c("a", "b")
  my_vec <- c(1, 2, 3)
  my_vec_with_na <- c(my_vec, NA)
  arg_check <- .makeAssertCollection()

  .CheckNumberVector(arg_check, my_char)
  .CheckNumberVector(arg_check, my_vec, min_length = 10)
  .CheckNumberVector(arg_check, my_vec, len = 4)
  .CheckNumberVector(arg_check, my_vec_with_na)

  expect_equal(
    arg_check$getMessages(),
    c(
        "my_char must have numeric entries (and not be NA)",
        "my_vec must have at least 10 elements",
        "my_vec must have exactly 4 elements",
        "my_vec_with_na must have numeric entries (and not be NA)"
    )
  )
})


test_that("test number array - passes", {
  my_vec <- c(1, 2.5, 3)
  arg_check <- .makeAssertCollection()

  .CheckNumberVector(arg_check, 3)
  .CheckNumberVector(arg_check, my_vec)
  .CheckNumberVector(arg_check, my_vec, min_length = 2)
  .CheckNumberVector(arg_check, my_vec, min_length = 3)
  .CheckNumberVector(arg_check, my_vec, len = 3)

  expect_true(arg_check$isEmpty())
})


test_that("test check n_burn - fails", {
  arg_check <- .makeAssertCollection()

  n_burn <- "a"
  .CheckNBurnAndNEnd(arg_check, n_burn, NA, n_iter = 1200, n_thin = 10)
  n_burn <- 500
  .CheckNBurnAndNEnd(arg_check, n_burn, NA, n_iter = 1200, n_thin = 10)
  n_burn <- -1
  .CheckNBurnAndNEnd(arg_check, n_burn, NA, n_iter = 1200, n_thin = 10)
  n_end <- 1300
  .CheckNBurnAndNEnd(arg_check, NA, n_end, n_iter = 1200, n_thin = 10)
  n_end <- 700
  .CheckNBurnAndNEnd(arg_check, NA, n_end, n_iter = 1200, n_thin = 15)

  expect_equal(
    arg_check$getMessages(),
    c(
        "n_burn must be an integer",
        "n_burn must be less than or equal to 200",
        "n_burn must be more than or equal to 0",
        "n_end must be less than or equal to 1200",
        "n_end must be more than or equal to 750"
    )
  )
})


test_that("test check n_burn - passes", {
  arg_check <- .makeAssertCollection()

  .CheckNBurnAndNEnd(arg_check, NA, NA, n_iter = 1200, n_thin = 10)
  .CheckNBurnAndNEnd(arg_check, 5000, NA, n_iter = 10000, n_thin = 10)
  .CheckNBurnAndNEnd(arg_check, NA, 9000, n_iter = 10000, n_thin = 10)
  .CheckNBurnAndNEnd(arg_check, 4000, 8000, n_iter = 10000, n_thin = 10)

  expect_true(arg_check$isEmpty())
})