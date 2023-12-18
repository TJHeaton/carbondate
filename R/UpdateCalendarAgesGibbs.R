#' Title
#'
#' @param likelihood_calendar_ages_from_calibration_curve A matrix with column i being vector of
#' likelihood of theta (for each value in calendar_age_grid) of ith rc_determination
#'
#' @param calendar_age_grid The possible set of calendar ages the observations can take
#'
#' @param rate_s Current set of (calendar age) changepoints in rate of
#' piecewise constant Poisson process
#' @param rate_h Current set of heights/rates in each section of
#' piecewise constant Poisson process
#'
#' @return TODO
#' @export
#'
#' @examples # TODO
UpdateCalendarAgesGibbs <- function(
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
    trimmed_likelihood_calendar_ages_from_calibration_curve,
    calendar_age_grid,
    rate_s,
    rate_h)
{
  prior_calendar_ages <- .FindCalendarAgePriorGivenPoissonProcess(
    rate_s = rate_s,
    rate_h = rate_h,
    theta = calendar_age_grid)

  trimmed_posterior_calendar_ages <- lapply(
    trimmed_likelihood_calendar_ages_from_calibration_curve,
    FUN = function(trimmed_likelihood, prior_calendar_ages) {
      trimmed_prior_calendar_ages = prior_calendar_ages[
        trimmed_likelihood$start_index:trimmed_likelihood$end_index]
      list(
        values = trimmed_likelihood$values * trimmed_prior_calendar_ages,
        start_index = trimmed_likelihood$start_index,
        end_index = trimmed_likelihood$end_index)
    },
    prior_calendar_ages = prior_calendar_ages)

  trimmed_updated_calendar_ages <- sapply(
    trimmed_posterior_calendar_ages,
    FUN = function(trimmed_probs, theta) {
      sample(theta[trimmed_probs$start_index:trimmed_probs$end_index], 1, prob=trimmed_probs$values)
    },
    theta = calendar_age_grid,
    simplify = TRUE)

  return(trimmed_updated_calendar_ages)
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

