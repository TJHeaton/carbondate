test_that("FindCalendarAgePriorGivenPoissonProcess gives expected outcomes", {

  rate_s <- c(0,5, 10)
  rate_h <- c(4, 7)
  n_theta <- 100

  n_heights <- length(rate_h)
  theta <- seq(min(rate_s), max(rate_s), length = n_theta)

  prior_for_theta <- .FindCalendarAgePriorGivenPoissonProcess(
    rate_s = rate_s,
    rate_h = rate_h,
    theta = theta
  )

  expect_length(prior_for_theta, n_theta)
  expect_identical(prior_for_theta[1], rate_h[1])
  expect_identical(prior_for_theta[n_theta - 1], rate_h[n_heights])
  expect_identical(prior_for_theta[n_theta], 0)

})

### Second test
# The legacy code is slightly different as will give h[nh] if theta is s[nh+1]
test_that("FindCalendarAgePriorGivenPoissonProcess gives same as legacy code", {

  rate_s <- c(0, 5, 10)
  rate_h <- c(4, 7)
  n_theta <- 100

  n_heights <- length(rate_h)
  theta <- seq(min(rate_s), max(rate_s), length = n_theta)

  revised_prior_for_theta <- .FindCalendarAgePriorGivenPoissonProcess(
    rate_s = rate_s,
    rate_h = rate_h,
    theta = theta
  )

  source(test_path("fixtures", "LegacyFindCalPrior.R"))
  legacy_prior_for_theta <- LegacyFindCalPrior(
    s = rate_s,
    h = rate_h,
    t = theta
  )

  expect_identical(revised_prior_for_theta,
                   legacy_prior_for_theta)

})
