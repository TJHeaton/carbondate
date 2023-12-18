test_that("ChangeHeight gives expected outcomes", {

  set.seed(14)

  rate_s <- initial_rate_s <- c(0, 2, 5, 7, 10)
  rate_h <- initial_rate_h <- c(5, 2, 4, 6)

  # Choose sensible prior on h
  # Prior mean matches mean of rate_h
  prior_h_shape <- 0.1 * mean(rate_h)
  prior_h_rate <- 0.1 # Bit disperse

  n_heights <- length(rate_h)
  integrated_rate <- .FindIntegral(
    initial_rate_s,
    initial_rate_h)

  n_theta <- 2 * integrated_rate # Double expected number

  calendar_ages <- stats::runif(
    n_theta,
    min = min(initial_rate_s),
    max = max(initial_rate_s))

  for(i in 1:1000) {
    return_val <- .ChangeHeight(
      theta = calendar_ages,
      rate_s = rate_s,
      rate_h = rate_h,
      integrated_rate = integrated_rate,
      prior_h_shape = prior_h_shape,
      prior_h_rate = prior_h_rate)

    rate_h <- return_rate_h <- return_val$rate_h
    integrated_rate <- return_integrated_rate <- return_val$integrated_rate
  }

  # Tests that heights are all strictly positive
  expect_true( all(return_rate_h > 0) )

  # Tests that number of heights hasn't varied
  expect_equal(
    length(return_rate_h),
    n_heights)


  # Test that have updated integrated rate correctly
  expect_equal(
    return_integrated_rate,
    .FindIntegral(rate_s, return_rate_h)
  )

  # Test as to whether it has updated the heights
  expect_false(identical(return_rate_h, initial_rate_h))

})


### Second test
test_that("ChangeHeight gives same as legacy code", {

  set.seed(14)
  rate_s <- initial_rate_s <- c(0, 2, 5, 7, 10)
  rate_h <- initial_rate_h <- c(5, 2, 4, 6)

  # Choose sensible prior on h
  # Prior mean matches mean of rate_h
  prior_h_shape <- 0.1 * mean(rate_h)
  prior_h_rate <- 0.1 # Var = 10 * mean(h) disperse

  n_heights <- length(rate_h)
  integrated_rate <- initial_integrated_rate <- .FindIntegral(
    initial_rate_s,
    initial_rate_h)

  n_theta <- 2 * integrated_rate # Double expected number

  calendar_ages <- stats::runif(
    n_theta,
    min = min(initial_rate_s),
    max = max(initial_rate_s))

  set.seed(11)
  for(i in 1:1000) {
    return_val <- .ChangeHeight(
      theta = calendar_ages,
      rate_s = rate_s,
      rate_h = rate_h,
      integrated_rate = integrated_rate,
      prior_h_shape = prior_h_shape,
      prior_h_rate = prior_h_rate)

    rate_h <- revised_rate_h <- return_val$rate_h
    integrated_rate <- revised_integrated_rate <- return_val$integrated_rate
  }

  set.seed(11)
  source(test_path("fixtures", "LegacyChangeHeight.R"))
  # Reset the starting values
  rate_h <- initial_rate_h
  rate_s <- initial_rate_s
  integrated_rate <- initial_integrated_rate

  for(i in 1:1000) {
    return_val <- LegacyChangeHe(
      th = calendar_ages,
      s = rate_s,
      h = rate_h,
      intrate = integrated_rate,
      alpha = prior_h_shape,
      beta = prior_h_rate)

    rate_h <- legacy_rate_h <- return_val$h
    integrated_rate <- legacy_integrated_rate <- return_val$intrate
  }

  # Test revised code gives same rate_h after running multiple time
  expect_identical(revised_rate_h, legacy_rate_h)

  # Test revised code gives same integrated rate after running multiple time
  expect_identical(revised_integrated_rate, legacy_integrated_rate)

})



