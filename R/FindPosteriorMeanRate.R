#' Find Posterior Mean Rate of Sample Occurrence for Poisson Process Model
#'
#' @description
#' Given output from the Poisson process fitting function [carbondate::PPcalibrate] calculate
#' the posterior mean rate of sample occurrence (i.e., the underlying Poisson process
#' rate \eqn{\lambda(t)}) together with specified probability intervals, on a given calendar age
#' grid (provided in cal yr BP).
#'
#' \strong{Note:} If you want to calculate and plot the result, use
#' [carbondate::PlotPosteriorMeanRate] instead.
#'
#' For more information read the vignette: \cr
#' \code{vignette("Poisson-process-modelling", package = "carbondate")}
#'
#' @inheritParams PlotPosteriorMeanRate
#' @param calendar_age_sequence A vector containing the calendar age grid on which to
#' calculate the posterior mean rate.
#'
#' @return A list, each item containing a data frame of the `calendar_age_BP`, the `rate_mean`
#' and the confidence intervals for the rate - `rate_ci_lower` and `rate_ci_upper`.
#'
#' @export
#'
#' @seealso [carbondate::PlotPosteriorMeanRate]
#'
#' @examples
#' # NOTE: All these examples are shown with a small n_iter and n_posterior_samples
#' # to speed up execution.
#' # Try n_iter and n_posterior_samples as the function defaults.
#'
#' pp_output <- PPcalibrate(
#'     pp_uniform_phase$c14_age,
#'     pp_uniform_phase$c14_sig,
#'     intcal20,
#'     n_iter = 1000,
#'     show_progress = FALSE)
#'
#' # Default plot with 2 sigma interval
#' FindPosteriorMeanRate(pp_output, seq(450, 640, length=10), n_posterior_samples = 100)
FindPosteriorMeanRate <- function(
    output_data,
    calendar_age_sequence,
    n_posterior_samples = 5000,
    interval_width = "2sigma",
    bespoke_probability = NA,
    n_burn = NA,
    n_end = NA) {

  n_iter <- output_data$input_parameters$n_iter
  n_thin <- output_data$input_parameters$n_thin

  arg_check <- .InitializeErrorList()
  .CheckOutputData(arg_check, output_data, "RJPP")
  .CheckIntervalWidth(arg_check, interval_width, bespoke_probability)
  .CheckCalendarAgeSequence(arg_check, calendar_age_sequence)
  .CheckNumber(arg_check, min(calendar_age_sequence), lower = min(output_data$rate_s[[1]]))
  .CheckNumber(arg_check, max(calendar_age_sequence), upper = max(output_data$rate_s[[1]]))
  .CheckInteger(arg_check, n_posterior_samples, lower = 10)
  .CheckNBurnAndNEnd(arg_check, n_burn, n_end, n_iter, n_thin)
  .ReportErrors(arg_check)

  n_burn <- .SetNBurn(n_burn, n_iter, n_thin)
  n_end <- .SetNEnd(n_end, n_iter, n_thin)

  indices <- sample((n_burn + 1):n_end, n_posterior_samples, replace = ((n_end - n_burn) < n_posterior_samples))
  rate <- matrix(NA, nrow = n_posterior_samples, ncol = length(calendar_age_sequence))
  for (i in 1:n_posterior_samples) {
    ind <- indices[i]
    rate[i,] <- stats::approx(
      x = output_data$rate_s[[ind]],
      y = c(output_data$rate_h[[ind]], 0),
      xout = calendar_age_sequence,
      method = "constant")$y
  }
  edge_width <- switch(
    interval_width,
    "1sigma" = 1 - stats::pnorm(1),
    "2sigma"  = 1 - stats::pnorm(2),
    "bespoke" = (1 - bespoke_probability)/2
  )
  posterior_rate <- data.frame(
    calendar_age_BP = calendar_age_sequence,
    rate_mean = apply(rate, 2, mean),
    rate_ci_lower = apply(rate, 2, stats::quantile, probs = edge_width),
    rate_ci_upper = apply(rate, 2, stats::quantile, probs = 1 - edge_width)
  )

  return(posterior_rate)
}
