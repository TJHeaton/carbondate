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
  my_num <- 2.5
  arg_check <- .makeAssertCollection()

  .CheckNumber(arg_check, my_char)
  .CheckNumber(arg_check, my_vec)
  .CheckNumber(arg_check, my_num, lower = 3.1)
  .CheckNumber(arg_check, my_num, upper = 2.1)

  expect_equal(
    arg_check$getMessages(),
    c(
      "my_char must be a number",
      "my_vec must be a number",
      "my_num must be more than or equal to 3.1",
      "my_num must be less than or equal to 2.1"
    )
  )
})


test_that("test check number - passes", {
  my_num <- 3
  arg_check <- .makeAssertCollection()

  .CheckNumber(arg_check, my_num)
  .CheckNumber(arg_check, my_num, lower = 2.9, upper = 3.1)

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


test_that("test vector - fails", {
  my_vec <- c(1, 2, 3)
  arg_check <- .makeAssertCollection()

  .CheckVector(arg_check, my_vec, min_length = 10)
  .CheckVector(arg_check, my_vec, len = 4)

  expect_equal(
    arg_check$getMessages(),
    c("my_vec must have at least 10 elements", "my_vec must have exactly 4 elements")
  )
})


test_that("test vector - passes", {
  my_vec <- c(1, 2.5, 3)
  my_char <- c("a", "b", "c", "d")
  arg_check <- .makeAssertCollection()

  .CheckVector(arg_check, 3)
  .CheckVector(arg_check, "a")
  .CheckVector(arg_check, my_vec, min_length = 3)
  .CheckVector(arg_check, my_char, min_length = 3)
  .CheckVector(arg_check, my_char, len = 4)

  expect_true(arg_check$isEmpty())
})


test_that("test number vector - fails", {
  my_char <- c("a", "b")
  my_vec <- c(1, 2, 3)
  my_vec_with_na <- c(my_vec, NA)
  arg_check <- .makeAssertCollection()

  .CheckNumberVector(arg_check, my_char)
  .CheckNumberVector(arg_check, my_vec, min_length = 10)
  .CheckNumberVector(arg_check, my_vec, len = 4)
  .CheckNumberVector(arg_check, my_vec, lower = 2)
  .CheckNumberVector(arg_check, my_vec_with_na)

  expect_equal(
    arg_check$getMessages(),
    c(
        "my_char must have numeric entries (and not be NA)",
        "my_vec must have at least 10 elements",
        "my_vec must have exactly 4 elements",
        "all entries of my_vec must be more than or equal to 2",
        "my_vec_with_na must have numeric entries (and not be NA)"
    )
  )
})


test_that("test number vector - passes", {
  my_vec <- c(1, 2.5, 3)
  arg_check <- .makeAssertCollection()

  .CheckNumberVector(arg_check, 3)
  .CheckNumberVector(arg_check, my_vec)
  .CheckNumberVector(arg_check, my_vec, min_length = 2)
  .CheckNumberVector(arg_check, my_vec, min_length = 3)
  .CheckNumberVector(arg_check, my_vec, len = 3)
  .CheckNumberVector(arg_check, my_vec, lower = 0)

  expect_true(arg_check$isEmpty())
})


test_that("test check choice", {
  allowed_choices <- c("apple", "banana", "kiwi")
  output_data <- list(fruit = "melon")

  arg_check <- .makeAssertCollection()
  .CheckChoice(arg_check, output_data$fruit, allowed_choices)
  expect_equal(arg_check$getMessages(), "output_data$fruit must be one of: apple, banana, kiwi")

  arg_check <- .makeAssertCollection()
  .CheckChoice(arg_check, "banana", allowed_choices)
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


test_that("test check slice parameters", {
  arg_check <- .makeAssertCollection()

  .CheckSliceParameters(arg_check, slice_width = 10, slice_multiplier = NA, sensible_initialisation = FALSE)
  .CheckSliceParameters(arg_check, slice_width = 10, slice_multiplier = 0.1, sensible_initialisation = TRUE)

  .CheckSliceParameters(arg_check, slice_width = NA, slice_multiplier = 10, sensible_initialisation = FALSE)
  .CheckSliceParameters(arg_check, slice_width = 0.1, slice_multiplier = 10, sensible_initialisation = FALSE)
  .CheckSliceParameters(arg_check, slice_width = 0.1, slice_multiplier = 10, sensible_initialisation = TRUE)

  expect_equal(
    arg_check$getMessages(),
    c(
        "slice_multiplier must be a number",
        "slice_multiplier must be more than or equal to 1",
        "slice_width must be a number",
        "slice_width must be more than or equal to 1",
        "slice_width must be more than or equal to 1"
    )
  )

  arg_check <- .makeAssertCollection()
  .CheckSliceParameters(arg_check, slice_width = 100, slice_multiplier = 10, sensible_initialisation = FALSE)
  .CheckSliceParameters(arg_check, slice_width = 100, slice_multiplier = 10, sensible_initialisation = TRUE)
  .CheckSliceParameters(arg_check, slice_width = NA, slice_multiplier = 10, sensible_initialisation = TRUE)
  expect_true(arg_check$isEmpty())

})
