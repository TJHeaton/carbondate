#' Plots the number of clusters
#'
#' Once a function has been run to calibrate a set of radiocarbon
#' determinations, the estimated number of clusters can be plotted using this
#' function.
#'
#' @param output_data Data returned from one of the updating functions e.g.
#' [carbondate::WalkerBivarDirichlet] or
#' [carbondate::BivarGibbsDirichletwithSlice].
#'
#' @return No return value
#' @export
#'
#' @examples
#' # Plot results for the number of clusters
#' PlotNumberOfClusters(output_data = walker_example_output)
PlotNumberOfClusters <- function(output_data) {
  n_out <- length(output_data$n_clust)
  n_burn <- floor(n_out / 2)

  n_clusters <- output_data$n_clust[n_burn:n_out]
  graphics::hist(n_clusters,
       xlab = "Number of Clusters",
       main = NA,
       probability = TRUE,
       breaks = seq(0.5, max(n_clusters) + 1, by = 1))
}
