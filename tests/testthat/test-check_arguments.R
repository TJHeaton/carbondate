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

  expect_equal(
    arg_check$getMessages(),
    c(
      "my_char must be an integer",
      "my_vec must be an integer",
      "my_double must be an integer",
      "my_int must be more than 10")
  )
})

test_that("test check integer - passes", {
  my_integerish <- 10
  my_integer <- 5L
  arg_check <- .makeAssertCollection()

  .CheckInteger(arg_check, my_integerish)
  .CheckInteger(arg_check, my_integer)

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
        "my_vec must be a single logical value (TRUE, FALSE OR NA)")
  )
})

test_that("test check flag - passes", {
  arg_check <- .makeAssertCollection()

  .CheckFlag(arg_check, TRUE)
  .CheckFlag(arg_check, FALSE)
  .CheckFlag(arg_check, NA)

  expect_true(arg_check$isEmpty())
})