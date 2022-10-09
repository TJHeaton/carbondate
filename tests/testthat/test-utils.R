test_that("interpolating calibration curve works", {
  # Choose ages not in calibration curve
  new_calendar_ages = c(54770, 54290)

  new_curve = InterpolateCalibrationCurve(new_calendar_ages, intcal20)

  # Manually checked these results
  expected_c14_ages = c(49905, 49660.5)
  expected_c14_sigs = c(963, 823.5)

  expect_equal(new_curve, data.frame(
    calendar_age=new_calendar_ages,
    c14_age=expected_c14_ages,
    c14_sig=expected_c14_sigs))
})
