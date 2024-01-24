#' Plot Calendar Ages of Changes in Rate of Sample Occurrence for Poisson Process Model
#'
#' @description
#' Given output from the Poisson process fitting function [carbondate::PPcalibrate], plot the
#' posterior density estimates for the calendar ages at which there are internal changepoints in the rate
#' of sample occurrence \eqn{\lambda(t)}. These density estimates are calculated \strong{conditional}
#' upon the number of internal changepoints within the period under study (which is specified as an input
#' to the function).
#'
#' Having conditioned on the number of changes, `n_change`, the code will extract all realisations
#' from the the posterior of the MCMC sampler which have that number of internal changepoints in the
#' estimate of \eqn{\lambda(t)}. It will then provide density estimates for the (ordered) calendar ages
#' of those internal changepoints. These density estimates are obtained using a Gaussian kernel.
#'
#' \strong{Note: These graphs will become harder to interpret as the specified number of changepoints
#' increases}
#'
#' For more information read the vignette: \cr
#' \code{vignette("Poisson-process-modelling", package = "carbondate")}
#'
#' @inheritParams PlotPosteriorMeanRate
#'
#' @param n_changes Number of internal changepoints to condition on, and plot for. A vector
#' which can contain at most 4 elements, with values in the range 1 to 6. If not given, then
#' `c(1, 2, 3)` will be used.
#' @param kernel_bandwidth (Optional) The bandwidth used for the (Gaussian) kernel smoothing of
#' the calendar age densities. If not given, then 1/50th of the overall calendar age range will be used.
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
#' # Plot the posterior change points for only 2 or 3 internal changes
#' PlotPosteriorChangePoints(pp_output, n_changes = c(2, 3))
PlotPosteriorChangePoints <- function(
    output_data,
    n_changes = c(1, 2, 3),
    n_burn = NA,
    n_end = NA,
    kernel_bandwidth = NA) {

  n_iter <- output_data$input_parameters$n_iter
  n_thin <- output_data$input_parameters$n_thin

  arg_check <- .InitializeErrorList()
  .CheckOutputData(arg_check, output_data, "RJPP")
  .CheckIntegerVector(arg_check, n_changes, lower = 1, upper = 6, max_length = 4)
  .CheckNBurnAndNEnd(arg_check, n_burn, n_end, n_iter, n_thin)
  if (!is.na(kernel_bandwidth)) .CheckNumber(arg_check, kernel_bandwidth, lower = 0)
  .ReportErrors(arg_check)

  n_burn <- .SetNBurn(n_burn, n_iter, n_thin)
  n_end <- .SetNEnd(n_end, n_iter, n_thin)

  max_density <- 0
  all_densities <- list()

  posterior_n_internal_changes <- output_data$n_internal_changes[n_burn + 1:n_end]
  posterior_rate_s <- output_data$rate_s[n_burn + 1:n_end]

  cal_age_range <- sort(output_data$input_parameters$pp_cal_age_range)
  if (is.na(kernel_bandwidth)) kernel_bandwidth <- diff(cal_age_range) / 50

  colors <- c("blue", "darkgreen", "red", "purple", "cyan3", "darkgrey")
  legend <- NULL

  for (n_change in n_changes) {
    legend <- c(legend, paste(n_change, "internal changes"))

    index <- which(posterior_n_internal_changes == n_change)
    if (length(index) == 0) {
      warning(paste("No posterior samples with", n_change, "internal changes"))
      next
    }

    extracted_posteriors <- do.call(rbind, posterior_rate_s[index])
    for (j in 2:(n_change + 1)) {
      tryCatch({
        smoothed_density <- stats::density(
          extracted_posteriors[, j], bw = kernel_bandwidth, from = cal_age_range[1], to = cal_age_range[2])
        if (max(smoothed_density$y) > max_density) max_density <- max(smoothed_density$y)
        this_line <- list(x = smoothed_density$x, y = smoothed_density$y, n_change = n_change, j = j)
        all_densities  <- append(all_densities, list(this_line))
      },
      error = function(cond) {
        message(paste("Could not calculate density for n_change =", n_change, ", j =",j, ": ", conditionMessage(cond)))
      })
    }
  }

  graphics::plot(
    x = NA,
    y = NA,
    xlim = rev(cal_age_range),
    ylim = c(0, max_density * 1.2),
    xlab = "Calendar Age (cal yr BP)",
    ylab = "Density",
    type = "n",
  )
  for (line in all_densities) {
    graphics::lines(line$x, line$y, lty = line$n_change, col = colors[line$n_change], lwd = 2)
  }
  graphics::legend("topright", legend = legend, lty = n_changes, col = colors[n_changes])
}
