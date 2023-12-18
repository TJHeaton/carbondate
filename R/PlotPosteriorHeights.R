#' Plots the posterior heights
#'
#' Description here
#'
#' @inheritParams PlotPosteriorChangePoints
#' @param kernel_bandwidth The bandwidth used for the KDE of the density (optional). If not give 1/50th of the
#' maximum height will be used.
#'
#' @export
#'
#' @examples
#' # Plot results for a single calibration
#' # TODO
PlotPosteriorHeights <- function(
    output_data,
    n_changes = c(1, 2, 3),
    n_burn = NA,
    n_end = NA,
    kernel_bandwidth = NA) {

  n_iter <- output_data$input_parameters$n_iter
  n_thin <- output_data$input_parameters$n_thin

  arg_check <- .InitializeErrorList()
  .CheckRJPPOutputData(arg_check, output_data)
  .CheckIntegerVector(arg_check, n_changes, lower = 1, upper = 6, max_length = 4)
  .CheckNBurnAndNEnd(arg_check, n_burn, n_end, n_iter, n_thin)
  if (!is.na(kernel_bandwidth)) .CheckNumber(arg_check, kernel_bandwidth, lower = 0)
  .ReportErrors(arg_check)

  n_burn <- .SetNBurn(n_burn, n_iter, n_thin)
  n_end <- .SetNEnd(n_end, n_iter, n_thin)

  max_density <- 0
  all_densities <- list()

  posterior_n_internal_changes <- output_data$n_internal_changes[n_burn + 1:n_end]
  posterior_rate_h <- output_data$rate_h[n_burn + 1:n_end]

  max_height <- 0
  for (n_change in n_changes) {
    index <- which(posterior_n_internal_changes == n_change)
    if (length(index) != 0) {
      extracted_posteriors <- do.call(rbind, posterior_rate_h[index])
      max_n_height <- max(extracted_posteriors)
      if (max(max_n_height) > max_height) max_height <- max(max_n_height)
    }
  }
  if (is.na(kernel_bandwidth)) kernel_bandwidth <- max_height / 50

  colors <- c("blue", "darkgreen", "red", "purple", "cyan3", "darkgrey")
  legend <- NULL

  for (n_change in n_changes) {
    legend <- c(legend, paste(n_change, "internal changes"))

    index <- which(posterior_n_internal_changes == n_change)
    if (length(index) == 0) {
      warning(paste("No posterior samples with", n_change, "internal changes"))
      next
    }

    extracted_posteriors <- do.call(rbind, posterior_rate_h[index])
    for (j in 1:(n_change + 1)) {
      tryCatch({
        smoothed_density <- stats::density(extracted_posteriors[, j], bw = kernel_bandwidth, from = 0, to = max_height)
        if (max(smoothed_density$y) > max_density) max_density <- max(smoothed_density$y)
        this_line <- list(x = smoothed_density$x, y = smoothed_density$y, n_change = n_change)
        all_densities  <- append(all_densities, list(this_line))
      },
      error = function(cond) {
        message(paste("Could not calculate density for n_change =", n_change, ", j =",j, ": ", conditionMessage(cond)))
      })
    }
  }

  plot(
    x = NA,
    y = NA,
    xlim = c(0, max_height),
    ylim = c(0, max_density * 1.2),
    xlab = "Rate of process (events per year)",
    ylab = "Density",
    type = "n"
  )
  for (line in all_densities) {
    graphics::lines(line$x, line$y, lty = line$n_change, col = colors[line$n_change], lwd = 2)
  }
  legend("topright", legend = legend, lty = n_changes, col = colors[n_changes])
}
