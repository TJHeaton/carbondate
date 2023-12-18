test_that("Death gives expected outcomes", {

  set.seed(14)

  # Set some true values for simulation of theta
  time_start <- 0
  time_end <- 10
  true_rate_s <- c(time_start , time_end)
  true_rate_h <- 20

  # Choose sensible prior on h
  # Prior mean matches mean of rate_h
  prior_h_shape <- 0.1 * mean(true_rate_h)
  prior_h_rate <- 0.1 # Bit disperse

  prior_n_internal_changepoints_lambda <- 60
  proposal_ratio <- 0.1

  true_integrated_rate <- .FindIntegral(
    true_rate_s,
    true_rate_h)

  n_theta <- 2 * true_integrated_rate # Double expected number

  calendar_ages <- stats::runif(
    n = n_theta,
    min = time_start,
    max = time_end
  )

  # Create too many initial changepoints
  initial_n_internal_changepoints <- 300
  initial_rate_s <- sort(
    c(
      time_start,
      runif(initial_n_internal_changepoints,
            min = time_start,
            max = time_end),
      time_end
      )
    )
  initial_rate_h <- jitter(
    rep(true_rate_h, initial_n_internal_changepoints + 1),
    amount = 2
  )

  # Test death step with initial values
  rate_s <- initial_rate_s
  rate_h <- initial_rate_h
  n_heights <- length(rate_h)
  integrated_rate <- .FindIntegral(
    rate_s,
    rate_h)

  n_iters <- initial_n_internal_changepoints - 1

  # Create variable to store progress
  n_changes <- rep(NA, n_iters)
  n_heights <- rep(NA, n_iters)
  are_heights_positive <- rep(NA, n_iters)
  are_changepoints_increasing <- rep(NA, n_iters)
  are_changepoints_bounds_correct <- rep(NA, n_iters)
  are_rate_lengths_compatible <- rep(NA, n_iters)

  set.seed(11)
  for(i in 1:n_iters) {
    return_val <- .Death(
      theta = calendar_ages,
      rate_s = rate_s,
      rate_h = rate_h,
      integrated_rate = integrated_rate,
      prior_h_shape = prior_h_shape,
      prior_h_rate = prior_h_rate,
      prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
      proposal_ratio = proposal_ratio)

    rate_h <- return_rate_h <- return_val$rate_h
    rate_s <- return_rate_s <- return_val$rate_s
    integrated_rate <- return_integrated_rate <- return_val$integrated_rate

    n_changes[i] <- length(rate_s)
    n_heights[i] <- length(rate_h)

    are_heights_positive[i] <- all(return_rate_h >= 0)
    are_changepoints_increasing[i] <- all(diff(return_rate_s) > 0)
    are_changepoints_bounds_correct[i] <- (
      (min(return_rate_s) == min(initial_rate_s)) &&
        (max(return_rate_s) == max(initial_rate_s))
    )
  }

  # Tests that heights are all strictly positive
  expect_true( all(are_heights_positive) )

  # Tests that changepoint locations are strictly increasing
  expect_true( all(are_changepoints_increasing) )

  # Tests that keeps same initial and final changepoints
  expect_true( all(are_changepoints_bounds_correct) )

  # Tests that number of height is always one less than number of changepoints
  expect_identical(n_changes - 1L, n_heights)

  # Test that have updated integrated rate correctly
  expect_equal(
    return_integrated_rate,
    .FindIntegral(return_rate_s, return_rate_h)
  )

  # Test as to whether it has updated the heights and changepoints
  # Whether is passes or fails will depend upon seed and initialisation point
  # This version should pass (as accepts some changes)
  expect_false(identical(return_rate_h, initial_rate_h))
  expect_false(identical(return_rate_s, initial_rate_s))

  # That each running either increases or keeps number changepoints the same
  expect_true(all(diff(n_changes) == 0 | diff(n_changes) == (-1)))

})

### Second test
test_that("Death gives same as legacy code", {

  set.seed(14)

  # Set some true values for simulation of theta
  time_start <- 0
  time_end <- 10
  true_rate_s <- c(time_start , time_end)
  true_rate_h <- 20

  # Choose sensible prior on h
  # Prior mean matches mean of rate_h
  prior_h_shape <- 0.1 * mean(true_rate_h)
  prior_h_rate <- 0.1 # Bit disperse

  prior_n_internal_changepoints_lambda <- 60
  proposal_ratio <- 0.1

  true_integrated_rate <- .FindIntegral(
    true_rate_s,
    true_rate_h)

  n_theta <- 2 * true_integrated_rate # Double expected number

  calendar_ages <- stats::runif(
    n = n_theta,
    min = time_start,
    max = time_end
  )

  # Create too many initial changepoints
  initial_n_internal_changepoints <- 300
  initial_rate_s <- sort(
    c(
      time_start,
      runif(initial_n_internal_changepoints,
            min = time_start,
            max = time_end),
      time_end
    )
  )
  initial_rate_h <- jitter(
    rep(true_rate_h, initial_n_internal_changepoints + 1),
    amount = 2
  )
  initial_integrated_rate <- .FindIntegral(
    initial_rate_s,
    initial_rate_h)

  set.seed(11)
  rate_s <- initial_rate_s
  rate_h <- initial_rate_h
  integrated_rate <- initial_integrated_rate

  for(i in 1:1000) {
    return_val <- .Death(
      theta = calendar_ages,
      rate_s = rate_s,
      rate_h = rate_h,
      integrated_rate = integrated_rate,
      prior_h_shape = prior_h_shape,
      prior_h_rate = prior_h_rate,
      prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
      proposal_ratio = proposal_ratio)

    rate_h <- revised_rate_h <- return_val$rate_h
    rate_s <- revised_rate_s <- return_val$rate_s
    integrated_rate <- revised_integrated_rate <- return_val$integrated_rate
  }

  set.seed(11)
  source(test_path("fixtures", "LegacyDeath.R"))
  # Reset the starting values
  rate_s <- initial_rate_s
  rate_h <- initial_rate_h
  integrated_rate <- initial_integrated_rate

  for(i in 1:1000) {
    return_val <- LegacyDeath(
      th = calendar_ages,
      s = rate_s,
      h = rate_h,
      intrate = integrated_rate,
      alpha = prior_h_shape,
      beta = prior_h_rate,
      lambda = prior_n_internal_changepoints_lambda,
      propratio = proposal_ratio)

    rate_h <- legacy_rate_h <- return_val$h
    rate_s <- legacy_rate_s <- return_val$s
    integrated_rate <- legacy_integrated_rate <- return_val$intrate
  }

  # Check has changed things
  expect_false(identical(revised_rate_h, initial_rate_h))
  expect_false(identical(revised_rate_s, initial_rate_s))

  # Test revised code gives same rate_s after running multiple time
  expect_identical(revised_rate_s, legacy_rate_s)

  # Test revised code gives same rate_h after running multiple time
  expect_identical(revised_rate_h, legacy_rate_h)

  # Test revised code gives same integrated rate after running multiple time
  expect_identical(revised_integrated_rate, legacy_integrated_rate)

})

