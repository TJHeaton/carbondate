#' Calibration of a set of individual radiocarbon samples by performing the
#' Gibbs sampler using a DPMM
#'
#'
#' @description This function takes as an input a set of radiocarbon determinations and
#' associated 1-sigma uncertainties, as well as the calibration curve which
#' should be used, and returns output data that can be sampled to estimate the
#' joint calendar age density and cluster.
#'
#' @details This method considers both the mean and the variance of the clusters
#' to be unknown. \[TODO Do we want to include more detail about the algorithm
#' here? Or refer to the paper.\]
#'
#' @inheritParams WalkerBivarDirichlet
#'
#' @return A list with 11 items. The first 7 items contain output data, each of
#' which have one dimension of size \eqn{n_{\textrm{out}} =
#' \textrm{floor}( n_{\textrm{iter}}/n_{\textrm{thin}}) + 1}, each row storing
#' the result from every \eqn{n_{\textrm{thin}}}th iteration:
#'
#' \describe{
#'  \item{`cluster_identifiers`}{An \eqn{n_{\textrm{out}}} by
#'     \eqn{n_{\textrm{obs}}} integer matrix. Gives the cluster allocation
#'      - an integer between 1 and n_clust - for each observation.}
#'  \item{`alpha`}{A double vector of length \eqn{n_{\textrm{out}}} giving the
#'      DP concentration parameter.}
#'  \item{`n_clust`}{An integer vector of length \eqn{n_{\textrm{out}}} giving
#'      the number of clusters.}
#'  \item{`phi`}{A list of length \eqn{n_{\textrm{out}}} each entry giving
#'      a vector of length n_clust of the cluster means \eqn{\phi_j}.}
#'  \item{`tau`}{A list of length \eqn{n_{\textrm{out}}} each entry giving
#'      a vector of length n_clust of the cluster uncertainties \eqn{\tau_j}.}
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
#'  \item{`update_type`}{A string that always has the value "neal"}
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
#' # using a saved example data set
#' PolyaUrnBivarDirichlet(kerr$c14_ages, kerr$c14_sig, intcal20)
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

  .check_input_data(
    arg_check, c14_determinations, c14_sigmas, calibration_curve)
  .check_dpmm_parameters(
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
  .check_iteration_parameters(arg_check, n_iter, n_thin)
  .check_slice_parameters(arg_check, slice_width, slice_multiplier)

  checkmate::reportAssertions(arg_check)

  ##############################################################################
  # Initialise parameters
  all_clusters_represented <- FALSE
  while (!all_clusters_represented) {
    cluster_identifiers <- sample(1:n_clust, num_observations, replace = TRUE)
    all_clusters_represented <- (length(unique(cluster_identifiers)) == n_clust)
  }

  if (sensible_initialisation) {
    initial_probabilities <- mapply(
      CalibrateSingleDetermination,
      c14_determinations,
      c14_sigmas,
      MoreArgs = list(calibration_curve=calibration_curve))
    indices_of_max_probability = apply(initial_probabilities, 2, which.max)

    calendar_ages <- calibration_curve$calendar_age[indices_of_max_probability]
    maxrange <- max(calendar_ages) - min(calendar_ages)

    mu_phi <- stats::median(calendar_ages)
    A <- stats::median(calendar_ages)
    B <- 1 / (maxrange)^2

    tempspread <- 0.1 * mad(calendar_ages)
    tempprec <- 1/(tempspread)^2

    lambda <- (100 / maxrange)^2
    nu1 <- 0.25
    nu2 <- nu1 / tempprec

    alpha_shape <- 1
    alpha_rate <- 1
  }

  alpha <- 0.0001

  # ADD COMMENT
  # tau <- rep(1 / (diff(range(c14_determinations)) / 4)^2, n_clust)
  tau <- rep(n_clust, 1 / (diff(range(c14_determinations)) / 4)^2)
  phi <- stats::rnorm(
    n_clust, mean = mu_phi, sd = diff(range(c14_determinations)) / 2)

  ##############################################################################
  # Create storage for output
  n_out = floor(n_iter / n_thin) + 1

  phi_out <- list(phi)
  tau_out <- list(tau)
  cluster_identifiers_out <- matrix(NA, nrow = n_out, ncol = num_observations)
  calendar_ages_out <- matrix(NA, nrow = n_out, ncol = num_observations)
  alpha_out <- rep(NA, length = n_out)
  mu_phi_out <- rep(NA, length = n_out)
  n_clust_out <- rep(NA, length = n_out)

  output_index <- 1
  cluster_identifiers_out[output_index, ] <- cluster_identifiers
  calendar_ages_out[output_index, ] <- calendar_ages
  alpha_out[output_index] <- alpha
  mu_phi_out[output_index] <- mu_phi
  n_clust_out[output_index] <- length(unique(cluster_identifiers))

  ##############################################################################
  ## Interpolate cal curve onto single year grid to speed up updating thetas
  integer_cal_year_curve <- InterpolateCalibrationCurve(
    1:pkg.globals$MAX_YEAR_BP, calibration_curve)
  interpolated_c14_age <- integer_cal_year_curve$c14_age
  interpolated_c14_sig <- integer_cal_year_curve$c14_sig

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
    for (i in 1:num_observations) {
      newclusters <- .BivarUpdateClusterIdentifier(
        i,
        c = cluster_identifiers,
        phi = phi,
        tau = tau,
        theta = calendar_ages[i],
        lambda = lambda,
        nu1 = nu1,
        nu2 = nu2,
        mu_phi = mu_phi,
        alpha = alpha)
      cluster_identifiers <- newclusters$c
      phi <- newclusters$phi
      tau <- newclusters$tau
      if (max(cluster_identifiers) != length(phi)) stop("Lengths do not match")
    }

    for (j in 1:length(phi)) {
      GibbsParams <- .UpdatePhiTau(
        theta = calendar_ages[cluster_identifiers == j],
        mu_phi = mu_phi,
        lambda = lambda,
        nu1 = nu1,
        nu2 = nu2)
      phi[j] <- GibbsParams$phi
      tau[j] <- GibbsParams$tau
    }
    mu_phi <- .UpdateMuPhi(phi = phi, tau = tau, lambda = lambda, A = A, B = B)

    for (k in 1:num_observations) {
      calendar_ages[k] <- .SliceSample(
        TARGET = .ThetaLogLikelihood,
        x0 = calendar_ages[k],
        slice_width = slice_width,
        slice_multiplier = slice_multiplier,
        type = "log",
        prmean = phi[cluster_identifiers[k]],
        prsig = 1 / sqrt(tau[cluster_identifiers[k]]),
        c14obs = c14_determinations[k],
        c14sig = c14_sigmas[k],
        mucalallyr = interpolated_c14_age,
        sigcalallyr = interpolated_c14_sig)
    }

    alpha <- .UpdateAlphaGammaPrior(
        cluster_identifiers, alpha, prshape = alpha_shape, prrate = alpha_rate)

    if (iter %% n_thin == 0) {
      output_index <- output_index + 1
      phi_out[[output_index]] <- phi
      tau_out[[output_index]] <- tau
      cluster_identifiers_out[output_index, ] <- cluster_identifiers
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
    calendar_ages = calendar_ages_out,
    alpha = alpha_out,
    mu_phi = mu_phi_out,
    n_clust = n_clust_out,
    update_type="neal",
    lambda = lambda,
    nu1 = nu1,
    nu2 = nu2)

  if (show_progress) close(progress_bar)
  return(return_list)
}
