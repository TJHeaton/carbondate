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
#' the length of `c14_determinations`.
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
#' # using a saved example data set. Note iterations are kept very small here
#' # for a faster run time.
#' PolyaUrnBivarDirichlet(kerr$c14_ages, kerr$c14_sig, intcal20, n_iter=100, n_thin=10)
PolyaUrnBivarDirichlet <- function(
    c14_determinations,
    c14_sigmas,
    calibration_curve,
    n_iter = 100,
    n_thin = 10,
    slice_width = max(1000, diff(range(c14_determinations)) / 2),
    slice_multiplier = 10,
    n_clust = min(10, length(c14_determinations)),
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
  num_observations <- length(c14_determinations)

  arg_check <- checkmate::makeAssertCollection()

  .CheckInputData(
    arg_check, c14_determinations, c14_sigmas, calibration_curve)
  .CheckDpmmParameters(
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
  .CheckSliceParameters(arg_check, slice_width, slice_multiplier)

  checkmate::reportAssertions(arg_check)

  ##############################################################################
  ## Interpolate cal curve onto single year grid to speed up updating thetas
  integer_cal_year_curve <- InterpolateCalibrationCurve(NA, calibration_curve)
  interpolated_calendar_age_start <- integer_cal_year_curve$calendar_age[1]
  interpolated_c14_age <- integer_cal_year_curve$c14_age
  interpolated_c14_sig <- integer_cal_year_curve$c14_sig

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

  if (sensible_initialisation) {
    initial_probabilities <- mapply(
      .ProbabilitiesForSingleDetermination,
      c14_determinations,
      c14_sigmas,
      MoreArgs = list(calibration_curve=integer_cal_year_curve))
    indices_of_max_probability = apply(initial_probabilities, 2, which.max)

    calendar_ages <- integer_cal_year_curve$calendar_age[indices_of_max_probability]
    maxrange <- max(calendar_ages) - min(calendar_ages)

    mu_phi <- stats::median(calendar_ages)
    A <- stats::median(calendar_ages)
    B <- 1 / (maxrange)^2

    tempspread <- 0.1 * stats::mad(calendar_ages)
    tempprec <- 1/(tempspread)^2

    lambda <- (100 / maxrange)^2
    nu1 <- 0.25
    nu2 <- nu1 / tempprec

    alpha_shape <- 1
    alpha_rate <- 1
  }

  alpha <- 0.0001

  tau <- rep(1 / (diff(range(c14_determinations)) / 4)^2, n_clust)
  phi <- stats::rnorm(
    n_clust, mean = mu_phi, sd = diff(range(c14_determinations)) / 2)

  ##############################################################################
  # Save input data and parameters
  input_data = list(
    c14_determinations = c14_determinations,
    c14_sigmas = c14_sigmas,
    calibration_curve_name = deparse(substitute(calibration_curve)))
  input_parameters = list(
    lambda = lambda,
    nu1 = nu1,
    nu2 = nu2,
    A = A,
    B = B,
    alpha_shape = alpha_shape,
    alpha_rate = alpha_rate,
    mu_phi = mu_phi,
    slice_width = slice_width,
    slice_multiplier = slice_multiplier)

  ##############################################################################
  # Create storage for output
  n_out = floor(n_iter / n_thin) + 1

  phi_out <- list(phi)
  tau_out <- list(tau)
  cluster_identifiers_out <- list(as.integer(cluster_identifiers))
  observations_per_cluster_out <- list(as.integer(observations_per_cluster))
  calendar_ages_out <- matrix(NA, nrow = n_out, ncol = num_observations)
  alpha_out <- rep(NA, length = n_out)
  mu_phi_out <- rep(NA, length = n_out)
  n_clust_out <- rep(NA, length = n_out)

  output_index <- 1
  calendar_ages_out[output_index, ] <- calendar_ages
  alpha_out[output_index] <- alpha
  mu_phi_out[output_index] <- mu_phi
  n_clust_out[output_index] <- length(unique(cluster_identifiers))

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
      as.double(c14_determinations),
      as.double(c14_sigmas),
      interpolated_calendar_age_start,
      interpolated_c14_age,
      interpolated_c14_sig)
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
    }
  }
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
    input_parameters = input_parameters)

  if (show_progress) close(progress_bar)
  return(return_list)
}
