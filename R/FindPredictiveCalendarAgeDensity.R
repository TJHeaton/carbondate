#' Find Predictive Estimate of Shared Calendar Age Density from Bayesian Non-Parametric
#' DPMM Output
#'
#' @description
#' Given output from one of the Bayesian non-parametric summarisation functions (either
#' [carbondate::PolyaUrnBivarDirichlet] or [carbondate::WalkerBivarDirichlet]) calculate the
#' predictive (summarised/shared) calendar age density and probability intervals
#' on a given calendar age grid.
#'
#' \strong{Note:} If you want to calculate and plot the result, use
#' [carbondate::PlotPredictiveCalendarAgeDensity] instead.
#'
#' @inheritParams PlotPredictiveCalendarAgeDensity
#' @param calendar_age_sequence A vector containing the calendar age grid on which to
#' calculate the predictive (summarised/shared) density.
#'
#' @return A data frame of the `calendar_age`, the
#' `density_mean` and the confidence intervals for the density
#' `density_ci_lower` and `density_ci_upper`.
#'
#' @export
#'
#' @seealso [carbondate::PlotPredictiveCalendarAgeDensity]
#'
#' @examples
#' # NOTE: All these examples are shown with a small n_iter and n_posterior_samples
#' # to speed up execution.
#' # Try n_iter and n_posterior_samples as the function defaults.
#'
#' # First generate output data
#' polya_urn_output <- PolyaUrnBivarDirichlet(
#'     two_normals$c14_age,
#'     two_normals$c14_sig,
#'     intcal20,
#'     n_iter = 100,
#'     show_progress = FALSE)
#'
#' # Find results for example output, 2-sigma confidence interval (default)
#' FindPredictiveCalendarAgeDensity(
#'     polya_urn_output, seq(3600, 4700, length=12), n_posterior_samples = 500)
FindPredictiveCalendarAgeDensity <- function(
    output_data,
    calendar_age_sequence,
    n_posterior_samples = 5000,
    interval_width = "2sigma",
    bespoke_probability = NA,
    n_burn = NA,
    n_end = NA) {

  arg_check <- .InitializeErrorList()

  .CheckOutputData(arg_check, output_data,  c("Polya Urn", "Walker"))
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
