#' Plots a histogram of the Gelman–Rubin convergence diagnostic for a single chain
#'
#' This plots a histogram of the potential scale reduction factor (PSRF) for each of
#' the calendar age observations for a single MCMC chain, by splitting the chain into segments after
#' `n_burn` and comparing the within-chain variance with the between-chains variance of the segments.
#' If the chain have converged to the target posterior distribution, then PSRF should be close to 1.
#'
#' @param output_data The return value from one of the updating functions e.g.
#' [carbondate::PolyaUrnBivarDirichlet] or
#' [carbondate::WalkerBivarDirichlet].
#' @param n_burn The number of iterations required for burn-in  - any iterations before this
#' are not used to calculate the PSRF. If not given, the first half of the
#' MCMC chain is discarded. Note that the maximum
#' value that can be chosen is `n_iter - 100 * n_thin` (where `n_iter` and `n_thin` are the
#' arguments given to [carbondate::PolyaUrnBivarDirichlet] or c[carbondate::WalkerBivarDirichlet]).
#' @param n_segments The number of segments to split the chain into. Default is 3, must be a
#' number between 2 and 10.
#'
#' @return No value
#'
#' @export
#'
#' @examples
#' # Plot results for the example data
#' PlotGelmanRubinDiagnosticSingleChain(polya_urn_example_output)
PlotGelmanRubinDiagnosticSingleChain <- function(output_data, n_burn = NA, n_segments = 3) {

  arg_check <- .InitializeErrorList()

  .CheckOutputData(arg_check, output_data, c("Polya Urn", "Walker"))
  n_iter <- output_data$input_parameters$n_iter
  n_thin <- output_data$input_parameters$n_thin
  n_obs <- length(output_data$input_data$rc_determinations)

  .CheckNBurnAndNEnd(arg_check, n_burn, NA, n_iter, n_thin)
  .CheckInteger(arg_check, n_segments, lower = 2, upper = 10)
  .ReportErrors(arg_check)

  n_burn <- .SetNBurn(n_burn, n_iter, n_thin)

  R <- rep(0, n_obs)
  for (i in 1:n_obs) {
    R[i] <- .SingleChainPsrf(output_data$calendar_ages[, i], n_burn, n_segments)
  }

  graphics::hist(R, xlab = "", main = "Histogram of PSRF for posterior calendar ages")
}


#' Plots a histogram of the Gelman–Rubin convergence diagnostic for a multiple chains
#'
#' This plots a histogram of the potential scale reduction factor (PSRF) for each of
#' the calendar age observations for a multiple chains by comparing the within-chain variance with
#' the between-chains variance after `n_burn` iterations.
#' If the chain have converged to the target posterior distribution, then PSRF should be close to 1
#' for all calendar ages.
#'
#' @param output_data_list A list, each item containing the return value from one of the updating
#' functions e.g. [carbondate::PolyaUrnBivarDirichlet] or [carbondate::WalkerBivarDirichlet].
#' The minimum number of elements in the list is 2.
#' @param n_burn The number of iterations required for burn-in  - any iterations before this
#' are not used to calculate the PSRF. If not given, the first half of the
#' MCMC chain is discarded. Note that the maximum
#' value that can be chosen is `n_iter - 100 * n_thin` (where `n_iter` and `n_thin` are the
#' arguments given to [carbondate::PolyaUrnBivarDirichlet] or c[carbondate::WalkerBivarDirichlet]).
#' @param n_burn The number of samples required for burn-in.
#'
#' @return No value
#'
#' @export
#'
#' @examples
#' # Plot results for the many chains
#' po = list()
#' for (i in 1:3) po[[i]] = PolyaUrnBivarDirichlet(kerr$c14_age, kerr$c14_sig, intcal20, n_iter=1e4)
#' PlotGelmanRubinDiagnosticMultiChain(po)
PlotGelmanRubinDiagnosticMultiChain <- function(output_data_list, n_burn = NA) {

  arg_check <- .InitializeErrorList()
  number_of_output_data <- length(output_data_list)
  .CheckInteger(arg_check, number_of_output_data, lower = 2)
  .CheckMultipleOutputDataConsistent(arg_check, output_data_list)
  for (i in 1:number_of_output_data) {
    .CheckOutputData(arg_check, output_data_list[[i]], c("Polya Urn", "Walker"))
  }
  n_iter <- output_data_list[[1]]$input_parameters$n_iter
  n_thin <- output_data_list[[1]]$input_parameters$n_thin
  n_obs <- length(output_data_list[[1]]$input_data$rc_determinations)

  .CheckNBurnAndNEnd(arg_check, n_burn, NA, n_iter, n_thin)
  .ReportErrors(arg_check)

  n_burn <- .SetNBurn(n_burn, n_iter, n_thin)
  R <- rep(0, n_obs)
  get_calendar_ages <- function(output_data, i) return(output_data$calendar_ages[, i])
  for (i in 1:n_obs) {
    theta_list <- lapply(output_data_list, get_calendar_ages, i)
    R[i] <- .MultiChainPsrf(theta_list, n_burn)
  }

  graphics::hist(R, xlab = "", main = "Histogram of PSRF for posterior calendar ages")
}


.SingleChainPsrf <- function(theta, n_burn, n_segments) {
  n_out <- length(theta)
  n_chain <- floor((n_out - n_burn) / n_segments)

  chains <- list()
  for (i in 1:n_segments) {
    chains[[i]] <- theta[(n_burn + (i - 1) * n_chain + 1):(n_burn + i * n_chain)]
  }

  return(.FindPsrf(chains, n_segments, n_chain))
}


.MultiChainPsrf <- function(theta_list, n_burn) {
  M <- length(theta_list)
  n_out <- length(theta_list[[1]])
  n_chain <- n_out - n_burn

  chains <- list()
  for (i in 1:M) {
    chains[[i]] <- theta_list[[i]][(n_burn + 1):n_out]
  }

  return(.FindPsrf(chains, M, n_chain))
}


.FindPsrf <- function(chains, M, N) {
  overall_mean <- mean(unlist(chains))
  diff_from_mean <- function(chain, overall_mean) return((mean(chain) - overall_mean)^2)

  within_chain_variance <- sum(unlist(lapply(chains, stats::var))) / M
  between_chain_variance <- N * sum(unlist(lapply(chains, diff_from_mean, overall_mean))) / (M - 1)

  pooled_variance <- (N - 1) / N * within_chain_variance
  pooled_variance <- pooled_variance + (M + 1) / (M * N) * between_chain_variance

  return(pooled_variance / within_chain_variance)
}
