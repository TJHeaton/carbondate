#' Plots the number of internal changes
#'
#' Once a function has been run to calibrate a set of radiocarbon determinations,
#' the estimated number of internal changes can be plotted using this function.
#'
#' @inheritParams PlotPosteriorChangePoints
#'
#' @return No return value
#' @export
#'
#' @examples
#' # TODO
PlotNumberOfInternalChanges <- function(output_data, n_burn = NA, n_end = NA) {

  arg_check <- .InitializeErrorList()
  .CheckRJPPOutputData(arg_check, output_data)
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
