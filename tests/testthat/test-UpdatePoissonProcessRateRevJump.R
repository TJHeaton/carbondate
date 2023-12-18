test_that("UpdatePoissonProcessRateRevJump gives expected outcomes", {

  # Set up "true" values
  set.seed(14)
  n_iters <- 1000

  true_rate_s <- 100 * c(0, 2, 5, 7, 10)
  true_rate_h <- c(5, 2, 4, 6)

  time_start <- min(true_rate_s)
  time_end <- max(true_rate_s)

  # Prior mean matches mean of rate_h
  prior_h_shape <- 0.1 * mean(true_rate_h)
  prior_h_rate <- 0.1 # Bit disperse

  prior_n_internal_changepoints_lambda <- 6
  k_max_internal_changepoints <- 20

  n_sections <- length(true_rate_h)
  n_obs_per_section <- rpois(
    n = n_sections,
    lambda = true_rate_h * diff(true_rate_s))

  calendar_ages <- c()
  for(i in 1:n_sections) {
    calendar_ages <- c(
      calendar_ages,
      runif(n_obs_per_section[i],
            min = true_rate_s[i],
            max = true_rate_s[i+1])
    )
  }

  # Now initialise (as overdispersed) and run
  prob_move <- .FindMoveProbability(
    prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
    k_max_internal_changepoints = k_max_internal_changepoints,
    rescale_factor = 0.9)

  # Create too many initial changepoints
  initial_n_internal_changepoints <- pmax(0, k_max_internal_changepoints - 5)
  initial_rate_s <- sort(
    c(
      time_start,
      runif(initial_n_internal_changepoints,
            min = time_start,
            max = time_end),
      time_end
    )
  )
  initial_rate_h <- rgamma(n = initial_n_internal_changepoints + 1,
                           shape = prior_h_shape,
                           rate = prior_h_rate)

  initial_integrated_rate <- .FindIntegral(
    initial_rate_s,
    initial_rate_h)

  # Create variable to store progress
  n_changes <- rep(NA, n_iters)
  n_heights <- rep(NA, n_iters)
  are_heights_positive <- rep(NA, n_iters)
  are_changepoints_increasing <- rep(NA, n_iters)
  are_changepoints_bounds_correct <- rep(NA, n_iters)
  are_rate_lengths_compatible <- rep(NA, n_iters)

  set.seed(11)
  rate_s <- initial_rate_s
  rate_h <- initial_rate_h
  integrated_rate <- initial_integrated_rate

  for(i in 1:n_iters) {
    return_val <- UpdatePoissonProcessRateRevJump(
      theta = calendar_ages,
      rate_s = rate_s,
      rate_h = rate_h,
      integrated_rate = integrated_rate,
      prior_h_shape = prior_h_shape,
      prior_h_rate = prior_h_rate,
      prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
      prob_move = prob_move
    )

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

  expect_true( all(are_heights_positive) )
  expect_true( all(are_changepoints_increasing) )
  expect_true( all(are_changepoints_bounds_correct) )
  expect_identical(n_changes - 1L, n_heights)
  expect_equal(
    return_integrated_rate,
    .FindIntegral(return_rate_s, return_rate_h)
  )

  # Test as to whether it has updated the heights and changepoints
  # Whether is passes or fails will depend upon seed and initialisation point
  # This version should pass (as accepts some changes)
  expect_false(identical(return_rate_h, initial_rate_h))
  expect_false(identical(return_rate_s, initial_rate_s))

  expect_true(all(diff(n_changes) == 0
                  | diff(n_changes) == (-1)
                  | diff(n_changes) == 1)
  )

})


### Second test
test_that("UpdatePoissonProcessRateRevJump gives same as legacy code", {

  # Set up "true" values
  set.seed(14)
  n_iters <- 1000

  true_rate_s <- 100 * c(0, 2, 5, 7, 10)
  true_rate_h <- c(5, 2, 4, 6)

  time_start <- min(true_rate_s)
  time_end <- max(true_rate_s)

  # Prior mean matches mean of rate_h
  prior_h_shape <- 0.1 * mean(true_rate_h)
  prior_h_rate <- 0.1 # Bit disperse

  prior_n_internal_changepoints_lambda <- 6
  k_max_internal_changepoints <- 20

  n_sections <- length(true_rate_h)
  n_obs_per_section <- rpois(
    n = n_sections,
    lambda = true_rate_h * diff(true_rate_s))

  calendar_ages <- c()
  for(i in 1:n_sections) {
    calendar_ages <- c(
      calendar_ages,
      runif(n_obs_per_section[i],
            min = true_rate_s[i],
            max = true_rate_s[i+1])
    )
  }

  # Create too many initial changepoints
  initial_n_internal_changepoints <- pmax(0, k_max_internal_changepoints - 5)
  initial_rate_s <- sort(
    c(
      time_start,
      runif(initial_n_internal_changepoints,
            min = time_start,
            max = time_end),
      time_end
    )
  )
  initial_rate_h <- rgamma(n = initial_n_internal_changepoints + 1,
                           shape = prior_h_shape,
                           rate = prior_h_rate)

  initial_integrated_rate <- .FindIntegral(
    initial_rate_s,
    initial_rate_h)

  # Run revised
  prob_move <- .FindMoveProbability(
    prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
    k_max_internal_changepoints = k_max_internal_changepoints,
    rescale_factor = 0.9)

  set.seed(11)
  rate_s <- initial_rate_s
  rate_h <- initial_rate_h
  integrated_rate <- initial_integrated_rate

  for(i in 1:n_iters) {
    return_val <- UpdatePoissonProcessRateRevJump(
      theta = calendar_ages,
      rate_s = rate_s,
      rate_h = rate_h,
      integrated_rate = integrated_rate,
      prior_h_shape = prior_h_shape,
      prior_h_rate = prior_h_rate,
      prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
      prob_move = prob_move
    )

    rate_h <- revised_rate_h <- return_val$rate_h
    rate_s <- revised_rate_s <- return_val$rate_s
    integrated_rate <- revised_integrated_rate <- return_val$integrated_rate
  }

  # Run legacy
  source(test_path("fixtures", "LegacyFindMoveProb.R"))
  source(test_path("fixtures", "LegacySingleRJUpdate.R"))
  pMove <- LegacyFindMoveProb(
    lambda = prior_n_internal_changepoints_lambda,
    kmax = k_max_internal_changepoints
  )

  set.seed(11)
  # Reset starting values
  rate_s <- initial_rate_s
  rate_h <- initial_rate_h
  integrated_rate <- initial_integrated_rate

  for(i in 1:n_iters) {
    return_val <- LegacyRJMCMCVar(
      th = calendar_ages,
      s = rate_s,
      h = rate_h,
      intrate = integrated_rate,
      alpha = prior_h_shape,
      beta = prior_h_rate,
      lambda = prior_n_internal_changepoints_lambda,
      pMove = pMove
    )

    rate_h <- legacy_rate_h <- return_val$h
    rate_s <- legacy_rate_s <- return_val$s
    integrated_rate <- legacy_integrated_rate <- return_val$intrate
  }

  # Check has changed things
  expect_false(identical(revised_rate_h, initial_rate_h))
  expect_false(identical(revised_rate_s, initial_rate_s))

  expect_identical(revised_rate_s, legacy_rate_s)
  expect_identical(revised_rate_h, legacy_rate_h)
  expect_identical(revised_integrated_rate, legacy_integrated_rate)

})
