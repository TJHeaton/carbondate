#' Plot Number of Calendar Age Clusters Estimated
#' in Bayesian Non-Parametric DPMM Output
#'
#' @description
#' Given output from one of the Bayesian non-parametric summarisation functions (either
#' [carbondate::PolyaUrnBivarDirichlet] or [carbondate::WalkerBivarDirichlet]) plot the
#' estimated number of calendar age clusters represented by the \eqn{{}^{14}}C samples.
#'
#' For more information read the vignette: \cr
#' \code{vignette("Non-parametric-summed-density", package = "carbondate")}
#'
#' @inheritParams PlotPredictiveCalendarAgeDensity
#'
#' @return None
#'
#' @export
#'
#' @seealso [carbondate::PlotPredictiveCalendarAgeDensity] and
#' [carbondate::PlotCalendarAgeDensityIndividualSample] for more plotting functions using DPMM output.
#'
#' @examples
#' # NOTE: these examples are shown with a small n_iter to speed up execution.
#' # When you run ensure n_iter gives convergence (try function default).
#'
#' polya_urn_output <- PolyaUrnBivarDirichlet(
#'     two_normals$c14_age,
#'     two_normals$c14_sig,
#'     intcal20,
#'     n_iter = 500,
#'     n_thin = 2,
#'     show_progress = FALSE)
#'
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
