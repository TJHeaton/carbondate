#' Plot Histogram of the Gelman-Rubin Convergence Diagnostic for a Single MCMC Chain
#'
#' @description
#' This plots a histogram of the potential scale reduction factors (PSRF) for each of the individual posterior
#' calendar age estimates within a single MCMC chain. Achieved by splitting the chain into segments
#' after `n_burn` and comparing the \emph{within-chain} variance with the \emph{between-chains}
#' variance of the segments. The PSRF of each sample's posterior calendar age is calculated.
#' If the chain have converged to the target posterior distribution, then PSRF should be close to 1
#' for all of the samples (a stringent condition is that all values are less than 1.1).
#'
#' For more information read the vignette: \cr
#' \code{vignette("determining-convergence", package = "carbondate")}
#'
#' @param output_data The return value from one of the updating functions, e.g.,
#' [carbondate::PolyaUrnBivarDirichlet], [carbondate::WalkerBivarDirichlet] or [carbondate::PPcalibrate].
#' @param n_burn The number of MCMC iterations that should be discarded for burn-in. This relates to
#' the total number of iterations `n_iter` when running the original update functions (not the
#' thinned `output_data`). Any MCMC iterations before this are not used in the calculations of the PSRF.
#' If not given, the first half of the MCMC chain is discarded. Note: The maximum value that the function
#' will allow is `n_iter - 100 * n_thin` (where `n_iter` and `n_thin` are the arguments that were given to
#' [carbondate::PPcalibrate]) which would leave only 100 of the (thinned) values in `output_data`.
#' @param n_segments The number of segments to split the chain into. Default is 3, must be a
#' number between 2 and 10.
#'
#' @return None
#'
#' @export
#'
#' @examples
#' # Plot results for the example data - n_iter is too small for convergence
#' # Try increasing n_iter to see the values of the PSRF decrease
#' polya_urn_output <- PolyaUrnBivarDirichlet(
#'     two_normals$c14_age,
#'     two_normals$c14_sig,
#'     intcal20,
#'     n_iter = 500,
#'     show_progress = FALSE)
#' PlotGelmanRubinDiagnosticSingleChain(polya_urn_output)
PlotGelmanRubinDiagnosticSingleChain <- function(output_data, n_burn = NA, n_segments = 3) {

  arg_check <- .InitializeErrorList()

  .CheckOutputData(arg_check, output_data, c("Polya Urn", "Walker", "RJPP"))
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


#' Plot Histogram of the Gelman-Rubin Convergence Diagnostic for Multiple Independent
#' MCMC Chains
#'
#' @description
#' This plots a histogram of the potential scale reduction factors (PSRF) for each of the individual
#' posterior calendar age estimates for multiple independent MCMC chains. Achieved by comparing the
#' \emph{within-chain} variance with the \emph{between-chains} variance after `n_burn` iterations.
#' The PSRF of each sample's posterior calendar age is calculated.
#' If the chain have converged to the target posterior distribution, then PSRF should
#' be close to 1 for all of the samples (a stringent condition is that all values are less than 1.1).
#'
#' For more information read the vignette: \cr
#' \code{vignette("determining-convergence", package = "carbondate")}
#'
#' @param output_data_list A list, each item containing the return value from one of the updating
#' functions e.g. [carbondate::PolyaUrnBivarDirichlet], [carbondate::WalkerBivarDirichlet] or
#' [carbondate::PPcalibrate]. The minimum number of elements in the list is 2.
#' @param n_burn The number of MCMC iterations that should be discarded for burn-in. This relates to
#' the total number of iterations `n_iter` when running the original update functions (not the
#' thinned `output_data`). Any MCMC iterations before this are not used in the calculations of the PSRF.
#' If not given, the first half of the MCMC chain is discarded. Note: The maximum value that the function
#' will allow is `n_iter - 100 * n_thin` (where `n_iter` and `n_thin` are the arguments that were given to
#' [carbondate::PPcalibrate]) which would leave only 100 of the (thinned) values in `output_data`.
#'
#' @return None
#'
#' @export
#'
#'
#' @examples
#' # Plot results for the example data - n_iter is too small for convergence
#' # Try increasing n_iter to see the values of the PSRF decrease
#' po = list()
#' for (i in 1:3) {
#'     set.seed(i)
#'     po[[i]] <- PolyaUrnBivarDirichlet(
#'         two_normals$c14_age,
#'         two_normals$c14_sig,
#'         intcal20,
#'         n_iter=400,
#'         show_progress = FALSE)
#' }
#' PlotGelmanRubinDiagnosticMultiChain(po)
PlotGelmanRubinDiagnosticMultiChain <- function(output_data_list, n_burn = NA) {

  arg_check <- .InitializeErrorList()
  number_of_output_data <- length(output_data_list)
  .CheckInteger(arg_check, number_of_output_data, lower = 2)
  .CheckMultipleOutputDataConsistent(arg_check, output_data_list)
  for (i in 1:number_of_output_data) {
    .CheckOutputData(arg_check, output_data_list[[i]], c("Polya Urn", "Walker", "RJPP"))
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
