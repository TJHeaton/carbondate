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
#'   walker_example_output, seq(600, 1700, length=12), 500)
#'
#' # Find results for example output, 1-sigma confidence interval (default)
#' FindPredictiveCalendarAgeDensity(
#'   walker_example_output, seq(600, 1700, length=12), 500, "1sigma")
#'
#' # Find results for example output, 95% confidence interval (default)
#' FindPredictiveCalendarAgeDensity(
#'   walker_example_output, seq(600, 1700, length=12), 500, "bespoke", 0.95)
FindPredictiveCalendarAgeDensity <- function(
    output_data,
    calendar_age_sequence,
    n_posterior_samples,
    interval_width = "2sigma",
    bespoke_probability = NA) {

  arg_check <- checkmate::makeAssertCollection()

  .CheckOutputData(arg_check, output_data)
  .CheckCalendarAgeSequence(arg_check, calendar_age_sequence)
  .CheckIntervalWidth(arg_check, interval_width, bespoke_probability)
  checkmate::assertInt(n_posterior_samples, lower = 10, add = arg_check)
  checkmate::reportAssertions(arg_check)

  density_matrix <- .FindDensityPerSampleID(
    output_data, calendar_age_sequence, n_posterior_samples)

  edge_width = switch(
    interval_width,
    "1sigma" = 1 - stats::pnorm(1),
    "2sigma"  = 1 - stats::pnorm(2),
    "bespoke" = (1 - bespoke_probability)/2)
  density_confidence_intervals <- apply(
    density_matrix,
    1,
    stats::quantile,
    probs = c(edge_width, 1 - edge_width))
  density_mean <- apply(density_matrix, 1, mean)

  return(data.frame(
    calendar_age=calendar_age_sequence,
    density_mean=density_mean,
    density_ci_lower=density_confidence_intervals[1, ],
    density_ci_upper=density_confidence_intervals[2, ]))
}


# Creates a matrix where each column is the density for a particular sample id
.FindDensityPerSampleID <- function(
    output_data, calendar_age_sequence, n_posterior_samples) {
  n_out <- length(output_data$alpha)
  n_burn <- floor(n_out / 2)

  sample_ids <- sample(
    x = n_burn:n_out,
    size = n_posterior_samples,
    replace = n_posterior_samples > (n_out - n_burn))

  density_matrix <- apply(
    matrix(sample_ids, 1, n_posterior_samples),
    2,
    function(i, output_data, x) {
      if (output_data$update_type == "Walker") {
        .FindPredictiveDensityWalker(
          x,
          weight = output_data$weight[[i]],
          phi = output_data$phi[[i]],
          tau = output_data$tau[[i]],
          mu_phi = output_data$mu_phi[i],
          lambda = output_data$input_parameters$lambda,
          nu1 = output_data$input_parameters$nu1,
          nu2 = output_data$input_parameters$nu2)
      } else {
        .FindPredictiveDensityPolyaUrn(
          x,
          cluster_identifiers = output_data$cluster_identifiers[i, ],
          phi = output_data$phi[[i]],
          tau = output_data$tau[[i]],
          alpha = output_data$alpha[i],
          mu_phi = output_data$mu_phi[i],
          lambda = output_data$input_parameters$lambda,
          nu1 = output_data$input_parameters$nu1,
          nu2 = output_data$input_parameters$nu2)
      }
    },
    output_data = output_data, x = calendar_age_sequence)
  return(density_matrix)
}


.FindPredictiveDensityPolyaUrn <- function(
    x, cluster_identifiers, phi, tau, alpha, mu_phi, lambda, nu1, nu2) {
  n_clust <- length(phi)
  nci <- .NumberOfObservationsInEachCluster(cluster_identifiers)

  # Find probability that new theta[i+1] is in a particular cluster or a new one
  pci <- c(nci, alpha) # Could form new cluster
  pci <- pci / sum(pci)

  return(FindPredictiveDensityPolyaUrn_cpp(
    calendar_ages = x,
    cluster_identifiers = as.integer(cluster_identifiers),
    phi = phi,
    tau = tau,
    alpha = alpha,
    mu_phi = mu_phi,
    lambda = lambda,
    nu1 = nu1,
    nu2 = nu2))
}


.FindPredictiveDensityWalker <- function(
    x, weight, phi, tau, mu_phi, lambda, nu1, nu2) {
  prob_new_clust <- 1 - sum(weight)
  return(FindPredictiveDensityWalker_cpp(
    calendar_ages = x,
    weight = weight,
    phi = phi,
    tau = tau,
    mu_phi = mu_phi,
    lambda = lambda,
    nu1 = nu1,
    nu2 = nu2))
}


.NumberOfObservationsInEachCluster <- function(cluster_identifiers) {
  n_clust <- max(cluster_identifiers)
  nci <- apply(
    t(as.matrix(1:n_clust)), # A row vector
    2,
    function(x, c) sum(c == x),
    c = cluster_identifiers)
}
