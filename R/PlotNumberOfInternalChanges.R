#' Plot Number of Changepoints in Rate of Sample Occurrence for Poisson Process Model
#'
#' @description
#' Given output from the Poisson process fitting function [carbondate::PPcalibrate], plot
#' the posterior distribution for the number of internal changepoints in the underlying rate of
#' sample occurrence (i.e., in \eqn{\lambda(t)}) over the period under study.
#'
#' For more information read the vignette: \cr
#' \code{vignette("Poisson-process-modelling", package = "carbondate")}
#'
#' @inheritParams PlotPosteriorMeanRate
#'
#' @return None
#'
#' @export
#'
#' @examples
#' # NOTE: This example is shown with a small n_iter to speed up execution.
#' # Try n_iter and n_posterior_samples as the function defaults.
#'
#' pp_output <- PPcalibrate(
#'     pp_uniform_phase$c14_age,
#'     pp_uniform_phase$c14_sig,
#'     intcal20,
#'     n_iter = 1000,
#'     show_progress = FALSE)
#'
#' PlotNumberOfInternalChanges(pp_output)
PlotNumberOfInternalChanges <- function(output_data, n_burn = NA, n_end = NA) {

  arg_check <- .InitializeErrorList()
  .CheckOutputData(arg_check, output_data, "RJPP")
  n_iter <- output_data$input_parameters$n_iter
  n_thin <- output_data$input_parameters$n_thin
  .CheckNBurnAndNEnd(arg_check, n_burn, n_end, n_iter, n_thin)
  .ReportErrors(arg_check)

  n_burn <- .SetNBurn(n_burn, n_iter, n_thin)
  n_end <- .SetNEnd(n_end, n_iter, n_thin)

  n_changes <- output_data$n_internal_changes[(n_burn + 1):n_end]
  graphics::hist(n_changes,
       xlab = "Number of Internal Changes",
       probability = TRUE,
       breaks = seq(0.5, max(n_changes) + 1, by = 1),
       main = "")
}
