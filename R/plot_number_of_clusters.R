#' Plots the number of clusters
#'
#' @param output_data Data returned from one of the updating functions e.g.
#' [carbondate::WalkerBivarDirichlet] or
#' [carbondate::BivarGibbsDirichletwithSlice].
#'
#' @return No return value
#' @export
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
