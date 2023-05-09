test_that("interpolating calibration curve works, c14 ages", {
  # Choose ages not in calibration curve
  new_calendar_ages = c(54770, 54290)

  new_curve = InterpolateCalibrationCurve(new_calendar_ages, intcal20, F14C_outputs = FALSE)

  # Manually checked these results
  expected_c14_ages = c(49905, 49660.5)
  expected_c14_sigs = c(963, 823.5)

  expect_equal(new_curve, data.frame(
    calendar_age_BP=new_calendar_ages,
    c14_age=expected_c14_ages,
    c14_sig=expected_c14_sigs))
})

test_that("interpolating calibration curve works, f14c ages", {
  # Choose ages not in calibration curve
  new_calendar_ages = c(54770, 54290)

  new_curve = InterpolateCalibrationCurve(new_calendar_ages, intcal20, F14C_outputs = TRUE)

  # Manually checked these results
  expected_f14c_ages = c(
    (intcal20$f14c[12] + intcal20$f14c[13])/2,
    (intcal20$f14c[36] + intcal20$f14c[37])/2)
  expected_f14c_sigs = c(
    (intcal20$f14c_sig[12] + intcal20$f14c_sig[13])/2,
    (intcal20$f14c_sig[36] + intcal20$f14c_sig[37])/2)

  expect_equal(new_curve, data.frame(
    calendar_age_BP=new_calendar_ages,
    f14c=expected_f14c_ages,
    f14c_sig=expected_f14c_sigs))
})


test_that("interpolating calibration curve works, both age scales", {
  # Choose ages not in calibration curve
  new_calendar_ages = c(54770, 54290)

  new_curve = InterpolateCalibrationCurve(new_calendar_ages, intcal20, F14C_outputs = NA)

  # Manually checked these results
  expected_c14_ages = c(49905, 49660.5)
  expected_c14_sigs = c(963, 823.5)
  expected_f14c_ages = c(
    (intcal20$f14c[12] + intcal20$f14c[13])/2,
    (intcal20$f14c[36] + intcal20$f14c[37])/2)
  expected_f14c_sigs = c(
    (intcal20$f14c_sig[12] + intcal20$f14c_sig[13])/2,
    (intcal20$f14c_sig[36] + intcal20$f14c_sig[37])/2)

  expect_equal(new_curve, data.frame(
    calendar_age_BP=new_calendar_ages,
    c14_age=expected_c14_ages,
    c14_sig=expected_c14_sigs,
    f14c=expected_f14c_ages,
    f14c_sig=expected_f14c_sigs))
})
