#' Calibration of a set of individual radiocarbon samples by performing the
#' Gibbs sampler using a DPMM
#'
#'
#' @description This function takes as an input a set of radiocarbon
#' determinations and associated 1-sigma uncertainties, as well as the
#' calibration curve which should be used, and returns output data that can be
#' sampled to estimate the joint calendar age density and clusters. This method
#' considers both the mean and the variance of the clusters to be unknown.
#'
#' @inheritParams WalkerBivarDirichlet
#'
#' @return A list with 10 items. The first 7 items contain output data, each of
#' which have one dimension of size \eqn{n_{\textrm{out}} =
#' \textrm{floor}( n_{\textrm{iter}}/n_{\textrm{thin}}) + 1}, each row storing
#' the result from every \eqn{n_{\textrm{thin}}}th iteration:
#'
#' \describe{
#'  \item{`cluster_identifiers`}{A list of length \eqn{n_{\textrm{out}}} each entry
#'      giving the cluster allocation - an integer between 1 and n_clust - for each
#'      observation.}
#'  \item{`alpha`}{A double vector of length \eqn{n_{\textrm{out}}} giving the
#'      DP concentration parameter.}
#'  \item{`n_clust`}{An integer vector of length \eqn{n_{\textrm{out}}} giving
#'      the number of clusters.}
#'  \item{`phi`}{A list of length \eqn{n_{\textrm{out}}} each entry giving
#'      a vector of length n_clust of the cluster means \eqn{\phi_j}.}
#'  \item{`tau`}{A list of length \eqn{n_{\textrm{out}}} each entry giving
#'      a vector of length n_clust of the cluster uncertainties \eqn{\tau_j}.}
#'  \item{`observations_per_cluster`}{A list of length \eqn{n_{\textrm{out}}} each entry giving
#'      a vector of length n_clust of the number of observations for that cluster.}
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
#'  \item{`update_type`}{A string that always has the value "Polya Urn".}
#'  \item{`input_data`}{a list containing the C14 data used and the name of
#'  the calibration curve used.}
#'  \item{`input_parameters`}{A list containing the values of the fixed
#'  hyperparameters `lambda`, `nu1`, `nu2`, `A`, `B`, `alpha_shape`,
#'  `alpha_rate` and `mu_phi`, and the slice parameters `slice_width` and
#'  `slice_multiplier`.}
#' }
#'
#' @export
#'
#' @examples
#' # Basic usage making use of sensible initialisation to set most values and
#' # using a saved example data set and the IntCal20 curve. See the result by plotting the predictive
#' # density
#' output = PolyaUrnBivarDirichlet(two_normals$c14_age, two_normals$c14_sig, intcal20, n_iter = 1e4)
#' PlotPredictiveCalendarAgeDensity(output)
#'
#' # The radiocarbon determinations can be given as F14C concentrations
#' output = PolyaUrnBivarDirichlet(
#'     two_normals$f14c,
#'     two_normals$f14c_sig,
#'     intcal20,
#'     F14C_inputs = TRUE,
#'     n_iter = 1e4)
#' PlotPredictiveCalendarAgeDensity(output)
PolyaUrnBivarDirichlet <- function(
    rc_determinations,
    rc_sigmas,
    calibration_curve,
    F14C_inputs = FALSE,
    n_iter = 1e5,
    n_thin = 10,
    use_F14C_space = TRUE,
    slice_width = NA,
    slice_multiplier = 10,
    n_clust = min(10, length(rc_determinations)),
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
    calendar_ages = NA) {

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
  .CheckFlag(arg_check, show_progress)
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
  all_clusters_represented <- FALSE
  while (!all_clusters_represented) {
    cluster_identifiers <- sample(1:n_clust, num_observations, replace = TRUE)
    all_clusters_represented <- (length(unique(cluster_identifiers)) == n_clust)
  }

  observations_per_cluster <- rep(0, n_clust)
  for (obs in cluster_identifiers)
    observations_per_cluster[obs] <- observations_per_cluster[obs] + 1

  initial_probabilities <- mapply(
    .ProbabilitiesForSingleDetermination,
    rc_determinations,
    rc_sigmas,
    MoreArgs = list(F14C_inputs=use_F14C_space, calibration_curve=integer_cal_year_curve))

  spd <- apply(initial_probabilities, 1, sum)
  cumulative_spd <- cumsum(spd) / sum(spd)
  spd_range_1_sigma <- c(
    integer_cal_year_curve$calendar_age_BP[min(which(cumulative_spd > 0.16))],
    integer_cal_year_curve$calendar_age_BP[max(which(cumulative_spd < 0.84))])
  spd_range_2_sigma <- c(
    integer_cal_year_curve$calendar_age_BP[min(which(cumulative_spd > 0.025))],
    integer_cal_year_curve$calendar_age_BP[max(which(cumulative_spd < 0.975))])
  spd_range_3_sigma <- c(
    integer_cal_year_curve$calendar_age_BP[min(which(cumulative_spd > 0.001))],
    integer_cal_year_curve$calendar_age_BP[max(which(cumulative_spd < 0.999))])

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

  alpha <- 0.0001

  tau <- rep(1 / (diff(spd_range_3_sigma) / 4)^2, n_clust)
  phi <- stats::rnorm(n_clust, mean = mu_phi, sd = diff(spd_range_3_sigma) / 2)

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
  cluster_identifiers_out <- list(as.integer(cluster_identifiers))
  observations_per_cluster_out <- list(as.integer(observations_per_cluster))
  calendar_ages_out <- matrix(NA, nrow = n_out, ncol = num_observations)
  alpha_out <- rep(NA, length = n_out)
  mu_phi_out <- rep(NA, length = n_out)
  n_clust_out <- rep(NA, length = n_out)
  densities_out <- matrix(NA, nrow = n_out, ncol = length(densities_cal_age_sequence))

  output_index <- 1
  calendar_ages_out[output_index, ] <- calendar_ages
  alpha_out[output_index] <- alpha
  mu_phi_out[output_index] <- mu_phi
  n_clust_out[output_index] <- length(unique(cluster_identifiers))
  densities_out[output_index, ] <- FindInstantPredictiveDensityPolyaUrn(
    densities_cal_age_sequence,
    as.integer(observations_per_cluster),
    phi,
    tau,
    alpha,
    mu_phi,
    num_observations,
    lambda,
    nu1,
    nu2)
  ##############################################################################
  # Now the calibration
  if (show_progress) {
    progress_bar <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
  }
  for (iter in 1:n_iter) {
    if (show_progress) {
      if (iter %% 100 == 0) {
        utils::setTxtProgressBar(progress_bar, iter)
      }
    }
    DPMM_update <- PolyaUrnUpdateStep(
      as.double(calendar_ages),
      as.integer(cluster_identifiers),
      phi,
      tau,
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
      as.double(rc_determinations),
      as.double(rc_sigmas),
      interpolated_calendar_age_start,
      interpolated_rc_age,
      interpolated_rc_sig)
    cluster_identifiers <- DPMM_update$cluster_ids
    phi <- DPMM_update$phi
    tau <- DPMM_update$tau
    observations_per_cluster <- DPMM_update$observations_per_cluster
    mu_phi <- DPMM_update$mu_phi
    calendar_ages <- DPMM_update$calendar_ages
    alpha <- DPMM_update$alpha

    if (iter %% n_thin == 0) {
      output_index <- output_index + 1
      phi_out[[output_index]] <- phi
      tau_out[[output_index]] <- tau
      cluster_identifiers_out[[output_index]] <- cluster_identifiers
      observations_per_cluster_out[[output_index]] <- observations_per_cluster
      calendar_ages_out[output_index, ] <- calendar_ages
      alpha_out[output_index] <- alpha
      mu_phi_out[output_index] <- mu_phi
      n_clust_out[output_index] <- max(cluster_identifiers)
      densities_out[output_index, ] <- FindInstantPredictiveDensityPolyaUrn(
        densities_cal_age_sequence,
        as.integer(observations_per_cluster),
        phi,
        tau,
        alpha,
        mu_phi,
        num_observations,
        lambda,
        nu1,
        nu2)
    }
  }

  density_data <- list(densities = densities_out, calendar_ages = densities_cal_age_sequence)

  return_list <- list(
    cluster_identifiers = cluster_identifiers_out,
    phi = phi_out,
    tau = tau_out,
    observations_per_cluster = observations_per_cluster_out,
    calendar_ages = calendar_ages_out,
    alpha = alpha_out,
    mu_phi = mu_phi_out,
    n_clust = n_clust_out,
    update_type="Polya Urn",
    input_data = input_data,
    input_parameters = input_parameters,
    density_data = density_data)

  if (show_progress) close(progress_bar)
  return(return_list)
}
