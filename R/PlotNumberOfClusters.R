#' Plots the number of clusters
#'
#' Once a function has been run to calibrate a set of radiocarbon
#' determinations, the estimated number of clusters can be plotted using this
#' function.
#'
#' @inheritParams PlotPredictiveCalendarAgeDensity
#' @param title The title to use for the plot. If nothing is given it will use
#' whatever is assigned to the `label` field of the `output_data`. If there is
#' no label, and nothing is provided as an argument, then no title will be
#' shown.
#'
#' @return No return value
#' @export
#'
#' @examples
#' # Plot results for the number of clusters
#' PlotNumberOfClusters(output_data = walker_example_output)
PlotNumberOfClusters <- function(output_data, title = NULL, n_burn = NA, n_end = NA) {

  .CheckOutputData(NULL, output_data)
  n_iter = output_data$input_parameters$n_iter
  n_thin = output_data$input_parameters$n_thin

  checkmate::assertInt(n_burn, lower = 0, upper = n_iter - 100 * n_thin, na.ok = TRUE)

  if (is.null(title)) {
    title = output_data$label
  }

  n_out <- length(output_data$n_clust)
  if (is.na(n_burn)) {
    n_burn = floor(n_out / 2)
  } else {
    n_burn = floor(n_burn / n_thin)
  }
  if (is.na(n_end)) {
    n_end = n_out
  } else {
    n_end = floor(n_end / n_thin)
  }

  n_clusters <- output_data$n_clust[(n_burn + 1):n_end]
  graphics::hist(n_clusters,
       xlab = "Number of Clusters",
       main = title,
       probability = TRUE,
       breaks = seq(0.5, max(n_clusters) + 1, by = 1))
}
