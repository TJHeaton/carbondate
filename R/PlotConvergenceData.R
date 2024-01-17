#' Plots data that can be used to estimate convergence of the MCMC sampler
#'
#' This plots the Kullback-Leibler divergence between the initial predictive density
#' and the density calculated during a run of a calibration function. This can be used to
#' estimate whether the calibration function has converged.
#'
#' @param output_data The return value from one of the updating functions e.g.
#' [carbondate::PolyaUrnBivarDirichlet] or
#' [carbondate::WalkerBivarDirichlet].
#' @param n_initial The number of iterations to use for the 'initial' predictive density,
#' which is then compared with all subsequent data points. If not given the minimum of 1000
#' iterations, or 1 / 10 of the the total number of data points, will be used.
#'
#' @return None
#'
#' @export
#'
#' @examples
#' # Plot results for the example data
#' polya_urn_output <- PolyaUrnBivarDirichlet(
#'     two_normals$c14_age, two_normals$c14_sig, intcal20, n_iter = 1e4, show_progress = FALSE)
#' PlotConvergenceData(polya_urn_output)
PlotConvergenceData <- function(output_data, n_initial = NA) {

  arg_check <- .InitializeErrorList()
  .CheckOutputData(arg_check, output_data, c("Polya Urn", "Walker"))

  n_iter <- output_data$input_parameters$n_iter
  n_thin <- output_data$input_parameters$n_thin
  n_out <- length(output_data$mu_phi)

  if (is.na(n_initial)) {
    n_initial <- min(floor(n_out / 10), floor(1000 / n_thin))
  } else {
    .CheckInteger(arg_check, n_initial, lower = 10 * n_thin, upper = floor(n_out / 2))
    n_initial <- floor(n_initial / n_thin)
  }

  .ReportErrors(arg_check)

  densities <- output_data$density_data$densities
  iters <- c(1, seq(n_thin, n_iter, by = n_thin))

  final_initial_index <- min(which(iters > n_initial))
  mean_initial_density <- apply(densities[1:final_initial_index, ], 2, sum) / final_initial_index

  kld_instant <- rep(NA, n_out - final_initial_index)
  kld_instant_avg <- rep(0, ceiling((n_out - final_initial_index) / 100))
  iters_avg <- rep(NA, ceiling((n_out - final_initial_index) / 100))
  avg_index <- 1
  for (i in (final_initial_index + 1):n_out) {
    kld_instant[i - final_initial_index] <- .KLD(mean_initial_density, densities[i, ])
    if (avg_index <= length(iters_avg)) {
      kld_instant_avg[avg_index] <- kld_instant_avg[avg_index] + kld_instant[i - final_initial_index]
    }
    if (i %% 100 == 0) {
      kld_instant_avg[avg_index] <- kld_instant_avg[avg_index] / 100
      iters_avg[avg_index] <- iters[i]
      avg_index <- avg_index + 1
    }
  }

  graphics::plot(
    iters[(final_initial_index + 1):n_out]/1000,
    (kld_instant),
    ylab = "Kullback Leibler divergence",
    xlab = "k iterations", type="p", pch=".",
    main = paste(
      "Kullback Leibler divergence between \n initial and a later iteration"))
  graphics::lines(iters_avg/1000, (kld_instant_avg), col="red")
  graphics::legend(
    "topright", legend = c("Single KLD", "Average of 100 KLDs"),
    lty = c(-1, 1), pch = c(".", NA), col = c("black", "red"))

}

.KLD <- function(P, Q) {
  P <- P/sum(P)
  Q <- Q/sum(Q)
  return (sum(P * log(P / Q)))
}
