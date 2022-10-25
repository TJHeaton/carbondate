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
#' # First generate some output data
#' walker_temp = WalkerBivarDirichlet(
#'   c14_determinations = c(602, 805, 954),
#'   c14_uncertainties = c(35, 34, 45),
#'   calibration_curve = intcal20,
#'   lambda = 0.1,
#'   nu1 = 0.25,
#'   nu2 = 10,
#'   alpha_shape = 1,
#'   alpha_rate = 1)
#'
#' # Plot results for the number of clusters
#' PlotNumberOfClusters(output_data = walker_temp)
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
