#' Plots the number of clusters
#'
#' Once a function has been run to calibrate a set of radiocarbon
#' determinations, the estimated number of clusters can be plotted using this
#' function.
#'
#' @inheritParams PlotPredictiveCalendarAgeDensity
#'
#' @return None
#'
#' @export
#'
#' @examples
#' # Plot results for the number of clusters
#' polya_urn_output <- PolyaUrnBivarDirichlet(
#'     two_normals$c14_age, two_normals$c14_sig, intcal20, n_iter = 1e4, show_progress = FALSE)
#' PlotNumberOfClusters(polya_urn_output)
PlotNumberOfClusters <- function(output_data, n_burn = NA, n_end = NA) {

  arg_check <- .InitializeErrorList()
  .CheckOutputData(arg_check, output_data, c("Polya Urn", "Walker"))
  n_iter <- output_data$input_parameters$n_iter
  n_thin <- output_data$input_parameters$n_thin
  .CheckNBurnAndNEnd(arg_check, n_burn, n_end, n_iter, n_thin)
  .ReportErrors(arg_check)

  n_burn <- .SetNBurn(n_burn, n_iter, n_thin)
  n_end <- .SetNEnd(n_end, n_iter, n_thin)

  n_clusters <- output_data$n_clust[(n_burn + 1):n_end]
  graphics::hist(
    n_clusters,
    xlab = "Number of Clusters",
    probability = TRUE,
    breaks = seq(0.5, max(n_clusters) + 1, by = 1),
    main = NA)
}
