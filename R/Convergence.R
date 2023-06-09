#' Plots a histogram of the Gelman–Rubin convergence diagnostic for a single chain
#'
#' This plots a histogram of the potential scale reduction factor (PSRF) for each of
#' the calendar age observations for a single MCMC chain, by splitting the chain into 3 parts after
#' n_burn and comparing the within-chain variance with the between-chains variance.
#' If the chain have converged to the target posterior distribution, then PSRF should be close to 1.
#'
#' @param output_data The return value from one of the updating functions e.g.
#' [carbondate::WalkerBivarDirichlet] or
#' [carbondate::PolyaUrnBivarDirichlet].
#' @param n_burn The number of samples required for burn-in.
#'
#' @return No value
#'
#' @export
#'
#' @examples
#' # Plot results for the example data
#' PlotConvergenceData(walker_example_output, 500)
PlotGelmanRubinDiagnostic <- function(output_data, n_burn) {
  n_obs = length(output$input_data$rc_determinations)
  n_thin = output$input_parameters$n_thin
  R = rep(0, n_obs)

  for (i in 1:n_obs) {
    R[i] = single_chain_psrf(output$calendar_ages[, i], n_burn, n_thin)
  }

  hist(R, xlab = "")
}


single_chain_psrf <- function(theta, n_burn, n_thin) {
  n_out = length(theta)
  n_burn = floor(n_burn / n_thin) + 1
  n_chain = floor((n_out - n_burn) / 3)

  chain_one = theta[(n_burn + 1):(n_burn + n_chain)]
  chain_two = theta[(n_burn + n_chain + 1):(n_burn + 2 * n_chain)]
  chain_three = theta[(n_burn + 2 * n_chain + 1):(n_burn + 3 * n_chain)]

  overall_mean = mean(theta[n_burn:(n_burn + 3 * n_chain)])

  within_chain_variance = 0.5 * (var(chain_one) + var(chain_two) + var(chain_three))
  between_chain_variance = 0.5 * n_chain * (
    (mean(chain_one) - overall_mean)^2
    + (mean(chain_two) - overall_mean)^2
    + (mean(chain_three) - overall_mean)^2
  )

  pooled_variance = (n_chain - 1) / n_chain * within_chain_variance
  pooled_variance = pooled_variance + 4 / (3 * n_chain) * between_chain_variance;

  return(pooled_variance / within_chain_variance)
}
