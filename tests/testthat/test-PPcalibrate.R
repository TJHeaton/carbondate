test_that("PPcalibrate gives expected outcomes", {
  skip_on_cran()
  # Checking
  set.seed(15)
  n_observed <- 20
  n_iter <- 100
  n_thin <- 10
  calibration_curve <- intcal20
  calendar_grid_resolution <- 10
  test_rc_sigma <- 15
  rc_sigmas <- rep(test_rc_sigma, n_observed)

  # Create artificial rc_determinations
  calendar_age_range <- c(1800,2150)
  observed_age_range <- c(1900,2000)
  true_theta <- seq(from = observed_age_range[1],
                    to = observed_age_range[2],
                    length = n_observed)

  # Value of curve itself
  rc_determinations <- approx(x = calibration_curve$calendar_age_BP,
                              y = calibration_curve$c14_age,
                              xout = true_theta)$y

  PP_fit_output <- PPcalibrate(
    rc_determinations = rc_determinations,
    rc_sigmas = rc_sigmas,
    calibration_curve = calibration_curve,
    calendar_age_range = calendar_age_range,
    calendar_grid_resolution = calendar_grid_resolution,
    n_iter = 1000,
    show_progress = FALSE
  )

  post_n_internal_changes <- PP_fit_output$n_internal_changes # Double
  post_rate_h <- PP_fit_output$rate_h # Length is integer
  post_rate_s <- PP_fit_output$rate_s
  lower_rate_s <- sapply(post_rate_s, function(x) x[1])
  upper_rate_s <- sapply(post_rate_s, function(x) x[length(x)])
  max_theta <- max(PP_fit_output$calendar_ages)
  min_theta <- min(PP_fit_output$calendar_ages)


  expect_identical(post_n_internal_changes,
                   sapply(post_rate_h, length) - 1) # Converts int to double
  expect_identical(sapply(post_rate_s, length) - 1L,
                   sapply(post_rate_h, length))
  expect_length(unique(lower_rate_s), 1L)
  expect_length(unique(upper_rate_s), 1L)
  expect_lt(lower_rate_s[1], min_theta)
  expect_gt(upper_rate_s[1], max_theta)
})

test_that("PPcalibrate supports one calibration curve per sample", {
  skip_on_cran()

  set.seed(15)
  n_observed <- 20
  calendar_grid_resolution <- 10
  test_rc_sigma <- 15
  rc_sigmas <- rep(test_rc_sigma, n_observed)
  calendar_age_range <- c(1800, 2150)
  true_theta <- seq(from = 1900, to = 2000, length = n_observed)

  sample_curves <- c(rep(list(intcal20), n_observed / 2), rep(list(shcal20), n_observed / 2))
  rc_determinations <- mapply(
    function(curve, theta) {
      approx(
        x = curve$calendar_age_BP,
        y = curve$c14_age,
        xout = theta
      )$y
    },
    sample_curves,
    true_theta
  )

  PP_fit_output <- PPcalibrate(
    rc_determinations = rc_determinations,
    rc_sigmas = rc_sigmas,
    calibration_curve = sample_curves,
    calendar_age_range = calendar_age_range,
    calendar_grid_resolution = calendar_grid_resolution,
    n_iter = 1000,
    show_progress = FALSE
  )

  expect_identical(PP_fit_output$update_type, "RJPP")
  expect_equal(ncol(PP_fit_output$calendar_ages), n_observed)
  expect_true(all(PP_fit_output$calendar_ages >= min(calendar_age_range)))
  expect_true(all(PP_fit_output$calendar_ages <= max(calendar_age_range)))
})

test_that("PPcalibrate validates per-sample calibration_curve length", {
  skip_on_cran()

  rc_determinations <- c(3000, 3050)
  rc_sigmas <- c(20, 20)

  expect_error(
    PPcalibrate(
      rc_determinations = rc_determinations,
      rc_sigmas = rc_sigmas,
      calibration_curve = list(intcal20),
      n_iter = 100,
      show_progress = FALSE
    ),
    "If calibration_curve is a list, it must have exactly 2 elements"
  )
})
