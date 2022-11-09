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
PlotNumberOfClusters <- function(output_data, title = NULL) {

  .CheckOutputData(NULL, output_data)

  if (is.null(title)) {
    title = output_data$label
  }

  n_out <- length(output_data$n_clust)
  n_burn <- floor(n_out / 2)

  n_clusters <- output_data$n_clust[n_burn:n_out]
  graphics::hist(n_clusters,
       xlab = "Number of Clusters",
       main = title,
       probability = TRUE,
       breaks = seq(0.5, max(n_clusters) + 1, by = 1))
}
