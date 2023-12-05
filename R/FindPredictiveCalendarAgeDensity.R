#' Find the predictive calendar age density of from the output data
#'
#' Given a calendar age sequence and output data from one of the package upate
#' functions, calculates the predicitive calendar age density and confidence
#' intervals. Note if you want to calculate and plot the result, use
#' [carbondate::PlotPredictiveCalendarAgeDensity] instead.
#'
#' @inheritParams PlotPredictiveCalendarAgeDensity
#' @param calendar_age_sequence A vector containing the calendar ages to
#' calculate the predictive density for.
#'
#' @return A data frame of the `calendar_age`, the
#' `density_mean` and the confidence intervals for the density
#' `density_ci_lower` and `density_ci_upper`.
#'
#' @export
#'
#' @examples
#' # Find results for example output, 2-sigma confidence interval (default)
#' FindPredictiveCalendarAgeDensity(
#'   polya_urn_example_output, seq(600, 1700, length=12), 500)
#'
#' # Find results for example output, 1-sigma confidence interval (default)
#' FindPredictiveCalendarAgeDensity(
#'   polya_urn_example_output, seq(600, 1700, length=12), 500, "1sigma")
#'
#' # Find results for example output, 95% confidence interval (default)
#' FindPredictiveCalendarAgeDensity(
#'   polya_urn_example_output, seq(600, 1700, length=12), 500, "bespoke", 0.95)
FindPredictiveCalendarAgeDensity <- function(
    output_data,
    calendar_age_sequence,
    n_posterior_samples,
    interval_width = "2sigma",
    bespoke_probability = NA,
    n_burn = NA,
    n_end = NA) {

  arg_check <- .InitializeErrorList()

  .CheckOutputData(arg_check, output_data)
  n_iter <- output_data$input_parameters$n_iter
  n_thin <- output_data$input_parameters$n_thin

  .CheckCalendarAgeSequence(arg_check, calendar_age_sequence)
  .CheckIntervalWidth(arg_check, interval_width, bespoke_probability)
  .CheckNBurnAndNEnd(arg_check, n_burn, n_end, n_iter, n_thin)
  .CheckInteger(arg_check, n_posterior_samples, lower = 10)
  .ReportErrors(arg_check)

  edge_width <- switch(
    interval_width,
    "1sigma" = 1 - stats::pnorm(1),
    "2sigma"  = 1 - stats::pnorm(2),
    "bespoke" = (1 - bespoke_probability)/2)

  n_burn <- .SetNBurn(n_burn, n_iter, n_thin)
  n_end <- .SetNEnd(n_end, n_iter, n_thin)

  if (output_data$update_type == "Walker") {
    return(
      FindPredictiveDensityAndCIWalker(
        calendar_ages = calendar_age_sequence,
        weights = output_data$weight,
        phis = output_data$phi,
        taus = output_data$tau,
        mu_phis = output_data$mu_phi,
        lambda = output_data$input_parameters$lambda,
        nu1 = output_data$input_parameters$nu1,
        nu2 = output_data$input_parameters$nu2,
        n_posterior_samples = n_posterior_samples,
        quantile_edge_width = edge_width,
        n_burn = n_burn,
        n_end = n_end
      )
    )
  } else {
    return(
      FindPredictiveDensityandCIPolyaUrn(
        calendar_ages = calendar_age_sequence,
        observations_per_clusters = output_data$observations_per_cluster,
        phis = output_data$phi,
        taus = output_data$tau,
        alphas = output_data$alpha,
        mu_phis = output_data$mu_phi,
        n_obs = length(output_data$cluster_identifiers[[1]]),
        lambda = output_data$input_parameters$lambda,
        nu1 = output_data$input_parameters$nu1,
        nu2 = output_data$input_parameters$nu2,
        n_posterior_samples = n_posterior_samples,
        quantile_edge_width = edge_width,
        n_burn = n_burn,
        n_end = n_end
      )
    )
  }
}
