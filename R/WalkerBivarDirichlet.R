#' Calibration of a set of individual radiocarbon samples using Walker updating
#' of the DPMM
#'
#' @description
#' This function takes as an input a set of radiocarbon determinations and
#' associated 1-sigma uncertainties, as well as the calibration curve which
#' should be used, and returns output data that can be sampled to estimate the
#' joint calendar age density and cluster.
#'
#' For more information read the vignette: \cr
#' \code{vignette("Non-parametric-summed-density", package = "carbondate")}
#'
#' @inheritParams FindSummedProbabilityDistribution
#' @param n_iter  The number of MCMC iterations (optional). Default is 100,000.
#' @param n_thin  How much to thin the output (optional). 1 is no thinning,
#' a larger number is more thinning. Default is 10. Must choose an integer more
#' than 1 and not too close to `n_iter`, to ensure there are enought samples from
#' posterior to potentially use.
#' @param use_F14C_space Whether the calculations are carried out in F14C space (default is TRUE).
#' If FALSE, calculations are carried out in \eqn{{}^{14}}C yr BP space.
#' @param slice_width Parameter for slice sampling (optional). If not given a value
#' is chosen intelligently based on the spread of the initial calendar ages.
#' Must be given if `sensible_initialisation` is `FALSE`.
#' @param slice_multiplier  Integer parameter for slice sampling (optional).
#' Default is 10. Limits the slice size to `slice_multiplier * slice_width`.
#' @param n_clust The initial number of clusters (optional). Must
#' be less than the length of `rc_determinations`. Default is 10 or the length
#' of `rc_determinations` if it is less than 10.
#' @param show_progress Whether to show a progress bar in the console during
#' execution. Default is `TRUE`.
#' @param sensible_initialisation Whether to use sensible start values and
#' adaptive prior on \eqn{\mu_{\phi}} and  (A, B).
#' If this is `TRUE` (the default), then all the remaining arguments below are
#' ignored.
#' @param lambda,nu1,nu2  Hyperparameters for the prior on the means
#' \eqn{\phi_j} and precision \eqn{\tau_j} of each individual calendar age
#' cluster \eqn{j}.
#' \deqn{(\phi_j, \tau_j)|\mu_{\phi} \sim
#' \textrm{NormalGamma}(\mu_{\phi}, \lambda, \nu_1, \nu_2)} where
#' \eqn{\mu_{\phi}} is the overall cluster centering. Required if
#' `sensible_initialisation` is `FALSE`.
#' @param A,B  Prior on \eqn{\mu_{\phi}} giving the mean and precision of the
#' overall centering \eqn{\mu_{\phi} \sim N(A, B^{-1})} i.e.
#' B small is uninformative. Required if `sensible_initialisation` is `FALSE`.
#' @param mu_phi Initial value of the overall cluster centering \eqn{\mu_{\phi}}.
#' Required if `sensible_initialisation` is `FALSE`.
#' @param alpha_shape,alpha_rate Hyperparameters for the shape and rate on prior
#' for DP concentration, \eqn{\alpha}, determining the number of clusters we
#' expect to observe among our n sampled objects.
#' \eqn{\alpha \sim \Gamma(\eta_1, \eta_2)} where \eqn{\eta_1, \eta_2} are
#' the `alpha_shape` and `alpha_rate`. A small alpha means more concentrated
#' (i.e. few clusters) while a large alpha means not concentrated (i.e. many
#' clusters).  Required if `sensible_initialisation` is `FALSE`.
#' @param calendar_ages  The initial estimate for the underlying calendar ages
#' (optional). If supplied it must be a vector with the same length as
#' `rc_determinations`.  Required if `sensible_initialisation` is `FALSE`.
#'
#' @return A list with 11 items. The first 8 items contain output data, each of
#' which have one dimension of size \eqn{n_{\textrm{out}} =
#' \textrm{floor}( n_{\textrm{iter}}/n_{\textrm{thin}}) + 1}, each row storing
#' the result from every \eqn{n_{\textrm{thin}}}th iteration:
#'
#' \describe{
#'  \item{`cluster_identifiers`}{An \eqn{n_{\textrm{out}}} by
#'     \eqn{n_{\textrm{obs}}} integer matrix. Gives the cluster allocation
#'      - an integer between 1 and n_clust - for each observation.}
#'  \item{`alpha`}{A double vector of length \eqn{n_{\textrm{out}}} giving the DP
#'     concentration parameter.}
#'  \item{`n_clust`}{An integer vector of length \eqn{n_{\textrm{out}}} giving
#'      the number of clusters.}
#'  \item{`phi`}{A list of length \eqn{n_{\textrm{out}}} each entry giving
#'      a vector of the cluster means \eqn{\phi_j}.}
#'  \item{`tau`}{A list of length \eqn{n_{\textrm{out}}} each entry giving
#'      a vector of the cluster uncertainties \eqn{\tau_j}.}
#'  \item{`weight`}{A list of length \eqn{n_{\textrm{out}}} each entry giving
#'      the mixing weights of each cluster.}
#'  \item{`calendar_ages`}{An \eqn{n_{\textrm{out}}} by \eqn{n_{\textrm{obs}}}
#'     integer matrix. Gives the calendar age for each observation.}
#'  \item{`mu_phi`}{A vector of length \eqn{n_{\textrm{out}}} giving the overall
#'      centering \eqn{\mu_{\phi}} of the clusters.}
#' }
#' where \eqn{n_{\textrm{obs}}} is the number of radiocarbon observations i.e.
#' the length of `rc_determinations`.
#'
#' The remaining items give information about input data, input parameters (or
#' those calculated using `sensible_initialisation`) and update_type
#'
#' \describe{
#'  \item{`update_type`}{A string that always has the value "Walker".}
#'  \item{`input_data`}{a list containing the C14 data used and the name of
#'  the calibration curve used.}
#'  \item{`input_parameters`}{A list containing the values of the fixed
#'  hyperparameters `lambda`, `nu1`, `nu2`, `A`, `B`, `alpha_shape`, and
#'  `alpha_rate`, and the slice parameters `slice_width` and
#'  `slice_multiplier`.}
#' }
#'
#' @export
#'
#' @examples
#' # Note these examples are shown with a small n_iter to speed up execution.
#' # When you run ensure n_iter gives convergence (try function default).
#'
#' walker_output <- WalkerBivarDirichlet(
#'     two_normals$c14_age,
#'     two_normals$c14_sig,
#'     intcal20,
#'     n_iter = 100,
#'     show_progress = FALSE)
#'
#' # The radiocarbon determinations can be given as F14C concentrations
#' walker_output <- WalkerBivarDirichlet(
#'     two_normals$f14c,
#'     two_normals$f14c_sig,
#'     intcal20,
#'     F14C_inputs = TRUE,
#'     n_iter = 100,
#'     show_progress = FALSE)
WalkerBivarDirichlet <- function(
    rc_determinations,
    rc_sigmas,
    calibration_curve,
    F14C_inputs = FALSE,
    n_iter = 1e5,
    n_thin = 10,
    use_F14C_space = TRUE,
    slice_width = NA,
    slice_multiplier = 10,
    show_progress = TRUE,
    sensible_initialisation = TRUE,
    lambda = NA,
    nu1 = NA,
    nu2 = NA,
    A = NA,
    B = NA,
    alpha_shape = NA,
    alpha_rate = NA,
    mu_phi = NA,
    calendar_ages = NA,
    n_clust = min(10, length(rc_determinations))) {

  ##############################################################################
  # Check input parameters
  num_observations <- length(rc_determinations)

  arg_check <- .InitializeErrorList()

  .CheckInputData(arg_check, rc_determinations, rc_sigmas, F14C_inputs)
  .CheckCalibrationCurve(arg_check, calibration_curve, NA)
  .CheckDPMMParameters(
    arg_check,
    sensible_initialisation,
    num_observations,
    lambda,
    nu1,
    nu2,
    A,
    B,
    alpha_shape,
    alpha_rate,
    mu_phi,
    calendar_ages,
    n_clust)
  .CheckIterationParameters(arg_check, n_iter, n_thin)
  .CheckSliceParameters(arg_check, slice_width, slice_multiplier, sensible_initialisation)
  .CheckFlag(arg_check, F14C_inputs)
  .CheckFlag(arg_check, use_F14C_space)

  .ReportErrors(arg_check)

  ##############################################################################
  ## Interpolate cal curve onto single year grid to speed up updating thetas
  integer_cal_year_curve <- InterpolateCalibrationCurve(NA, calibration_curve, use_F14C_space)
  interpolated_calendar_age_start <- integer_cal_year_curve$calendar_age_BP[1]
  if (use_F14C_space) {
    interpolated_rc_age <- integer_cal_year_curve$f14c
    interpolated_rc_sig <- integer_cal_year_curve$f14c_sig
  } else {
    interpolated_rc_age <- integer_cal_year_curve$c14_age
    interpolated_rc_sig <- integer_cal_year_curve$c14_sig
  }

  ##############################################################################
  # Save input data
  input_data <- list(
    rc_determinations = rc_determinations,
    rc_sigmas = rc_sigmas,
    F14C_inputs = F14C_inputs,
    calibration_curve_name = deparse(substitute(calibration_curve)))

  ##############################################################################
  # Convert the scale of the initial determinations to F14C or C14_age as appropriate
  # if they aren't already

  if (F14C_inputs == use_F14C_space) {
    rc_determinations <- as.double(rc_determinations)
    rc_sigmas <- as.double(rc_sigmas)
  } else if (F14C_inputs == FALSE) {
    converted <- .Convert14CageToF14c(rc_determinations, rc_sigmas)
    rc_determinations <- converted$f14c
    rc_sigmas <- converted$f14c_sig
  } else {
    converted <- .ConvertF14cTo14Cage(rc_determinations, rc_sigmas)
    rc_determinations <- converted$c14_age
    rc_sigmas <- converted$c14_sig
  }

  ##############################################################################
  # Initialise parameters
  initial_probabilities <- mapply(
    .ProbabilitiesForSingleDetermination,
    rc_determinations,
    rc_sigmas,
    MoreArgs = list(F14C_inputs=use_F14C_space, calibration_curve=integer_cal_year_curve))

  spd <- apply(initial_probabilities, 1, sum)
  cumulative_spd <- cumsum(spd) / sum(spd)
  spd_range_1_sigma <- c(
    integer_cal_year_curve$calendar_age_BP[min(which(cumulative_spd > (1 - 0.683)/2))],
    integer_cal_year_curve$calendar_age_BP[max(which(cumulative_spd < (1 + 0.683)/2))])
  spd_range_2_sigma <- c(
    integer_cal_year_curve$calendar_age_BP[min(which(cumulative_spd > (1 - 0.954)/2))],
    integer_cal_year_curve$calendar_age_BP[max(which(cumulative_spd < (1 + 0.954)/2))])
  spd_range_3_sigma <- c(
    integer_cal_year_curve$calendar_age_BP[min(which(cumulative_spd > (1 - 0.997)/2))],
    integer_cal_year_curve$calendar_age_BP[max(which(cumulative_spd < (1 + 0.997)/2))])

  plot_range <- spd_range_3_sigma + c(-1, 1) * diff(spd_range_3_sigma) * 0.1
  plot_range <- c(
    max(plot_range[1], min(calibration_curve$calendar_age_BP)),
    min(plot_range[2], max(calibration_curve$calendar_age_BP)))

  densities_cal_age_sequence <- seq(plot_range[1], plot_range[2], length.out=100)

  if (sensible_initialisation) {
    indices_of_max_probability <- apply(initial_probabilities, 2, which.max)
    calendar_ages <- integer_cal_year_curve$calendar_age_BP[indices_of_max_probability]

    mu_phi <- mean(spd_range_2_sigma)
    A <- mean(spd_range_2_sigma)
    B <- 1 / (diff(spd_range_2_sigma))^2

    tempspread <- 0.05 * diff(spd_range_1_sigma)
    tempprec <- 1/(tempspread)^2

    lambda <- (100 / diff(spd_range_3_sigma))^2
    nu1 <- 0.25
    nu2 <- nu1 / tempprec

    alpha_shape <- 1
    alpha_rate <- 1

    if (is.na(slice_width)) slice_width <- diff(spd_range_3_sigma)
  }

  # do not allow very small values of alpha as this causes crashes
  alpha <- 2

  tau <- stats::rgamma(n_clust, shape = nu1, rate = nu2)
  phi <- stats::rnorm(n_clust, mean = mu_phi, sd = 1 / sqrt(lambda * tau))

  v <- stats::rbeta(n_clust, 1, alpha)
  weight <- v * c(1, cumprod(1 - v)[-n_clust])
  cluster_identifiers <- as.integer(sample(1:n_clust, num_observations, replace = TRUE))
  calendar_ages <- as.double(calendar_ages)

  ##############################################################################
  # Save input parameters
  input_parameters <- list(
    lambda = lambda,
    nu1 = nu1,
    nu2 = nu2,
    A = A,
    B = B,
    alpha_shape = alpha_shape,
    alpha_rate = alpha_rate,
    slice_width = slice_width,
    slice_multiplier = slice_multiplier,
    n_iter = n_iter,
    n_thin = n_thin)

  ##############################################################################
  # Create storage for output
  n_out <- .SetNOut(n_iter, n_thin)

  phi_out <- list(phi)
  tau_out <- list(tau)
  w_out <- list(weight)
  cluster_identifiers_out <- matrix(NA, nrow = n_out, ncol = num_observations)
  alpha_out <- rep(NA, length = n_out)
  n_clust_out <- rep(NA, length = n_out)
  mu_phi_out <- rep(NA, length = n_out)
  theta_out <- matrix(NA, nrow = n_out, ncol = num_observations)
  densities_out <- matrix(NA, nrow = n_out, ncol = length(densities_cal_age_sequence))

  output_index <- 1
  cluster_identifiers_out[output_index, ] <- cluster_identifiers
  alpha_out[output_index] <- alpha
  n_clust_out[output_index] <- length(unique(cluster_identifiers))
  mu_phi_out[output_index] <- mu_phi
  theta_out[output_index, ] <- calendar_ages
  densities_out[output_index, ] <- FindInstantPredictiveDensityWalker(
    densities_cal_age_sequence, weight, phi, tau, mu_phi, lambda, nu1, nu2)
  ##############################################################################
  # Now the calibration and DPMM
  if (show_progress) {
    progress_bar <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
  }
  for (iter in 1:n_iter) {
    if (show_progress) {
      if (iter %% 100 == 0) {
        utils::setTxtProgressBar(progress_bar, iter)
      }
    }
    DPMM_update <- WalkerUpdateStep(
      as.double(calendar_ages),
      weight,
      v,
      cluster_identifiers,
      alpha,
      mu_phi,
      alpha_shape,
      alpha_rate,
      lambda,
      nu1,
      nu2,
      A,
      B,
      slice_width,
      slice_multiplier,
      rc_determinations,
      rc_sigmas,
      interpolated_calendar_age_start,
      interpolated_rc_age,
      interpolated_rc_sig)
    weight <- DPMM_update$weight
    cluster_identifiers <- DPMM_update$cluster_ids
    phi <- DPMM_update$phi
    tau <- DPMM_update$tau
    v <- DPMM_update$v
    alpha <- DPMM_update$alpha
    mu_phi <- DPMM_update$mu_phi
    calendar_ages <- DPMM_update$calendar_ages


    if (iter %% n_thin == 0) {
      output_index <- output_index + 1
      cluster_identifiers_out[output_index, ] <- cluster_identifiers
      alpha_out[output_index] <- alpha
      n_clust_out[output_index] <- length(unique(cluster_identifiers))
      phi_out[[output_index]] <- phi
      tau_out[[output_index]] <- tau
      theta_out[output_index, ] <- calendar_ages
      w_out[[output_index]] <- weight
      mu_phi_out[output_index] <- mu_phi
      densities_out[output_index, ] <- FindInstantPredictiveDensityWalker(
        densities_cal_age_sequence, weight, phi, tau, mu_phi, lambda, nu1, nu2)
    }
  }

  density_data <- list(densities = densities_out, calendar_ages = densities_cal_age_sequence)

  return_list <- list(
    cluster_identifiers = cluster_identifiers_out,
    alpha = alpha_out,
    n_clust = n_clust_out,
    phi = phi_out,
    tau = tau_out,
    weight = w_out,
    calendar_ages = theta_out,
    mu_phi = mu_phi_out,
    update_type="Walker",
    input_data = input_data,
    input_parameters = input_parameters,
    density_data = density_data)
  if (show_progress) close(progress_bar)
  return(return_list)
}
