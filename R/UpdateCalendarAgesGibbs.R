.UpdateCalendarAgesGibbs <- function(
    likelihood_calendar_ages_from_calibration_curve,
    calendar_age_grid,
    rate_s,
    rate_h)
{
  prior_calendar_ages <- .FindCalendarAgePriorGivenPoissonProcess(
    rate_s = rate_s,
    rate_h = rate_h,
    theta = calendar_age_grid)

  posterior_calendar_ages <- sweep(
    likelihood_calendar_ages_from_calibration_curve,
    MARGIN=1,
    prior_calendar_ages,
    FUN = `*`)

  updated_calendar_ages <- apply(
    posterior_calendar_ages,
    MARGIN = 2,
    FUN = function(probs, theta) {
      sample(theta, 1, prob = probs)
    },
    theta = calendar_age_grid)

  return(updated_calendar_ages)
}


.TrimmedUpdateCalendarAgesGibbs <- function(
    calendar_age_grid,
    rate_s,
    rate_h,
    likelihood_values,
    likelihood_offsets)
{
  prior_calendar_ages <- .FindCalendarAgePriorGivenPoissonProcess(
    rate_s = rate_s,
    rate_h = rate_h,
    theta = calendar_age_grid)

  updated_calendar_ages <- UpdateCalendarAgesGibbsCpp(
    prior_calendar_ages, calendar_age_grid, likelihood_values, likelihood_offsets)

  return(updated_calendar_ages)
}


# Find (un-normalised) prior on theta (on a specific calendar grid) for a given Poisson Process
.FindCalendarAgePriorGivenPoissonProcess <- function(
    rate_s,
    rate_h,
    theta)
{
  # Prior is proportional to rate at given calendar age
  prior_theta <- stats::approx(
    x = rate_s,
    y = c(rate_h, 0), # approx() interpolates from left
    xout = theta,
    method = "constant")$y

  return(prior_theta)
}

