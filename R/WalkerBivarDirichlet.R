#' Calibration of a set of individual radiocarbon samples using Walker updating
#' of the DPMM
#'
#' This function takes as an input a set of radiocarbon determinations and
#' associated 1-sigma uncertainties, as well as the calibration curve which
#' should be used, and returns output data that can be sampled to estimate the
#' joint calendar age density and cluster.
#'
#' @inheritParams FindSummedProbabilityDistribution
#' @param n_iter  The number of MCMC iterations (optional). Default is 100.
#' @param n_thin  How much to thin the output (optional). 1 is no thinning,
#' a larger number is more thinning. Default is 10. Must choose an integer more
#' than 1 and not too close to `n_iter`, since after burn-in there are
#' \eqn{(n_{\textrm{iter}}/n_{\textrm{thin}})/2} samples from posterior to
#' potentially use.
#' @param slice_width  Parameter for slice sampling (optional). Default is
#' either 1000 or half the width of the c14 determinations, whichever is larger.
#' @param slice_multiplier  Integer parameter for slice sampling (optional).
#' Default is 10. Limits the slice size to `slice_multiplier * slice_width`.
#' @param n_clust The initial number of clusters (optional). Must
#' be less than the length of `c14_determinations`. Default is 10 or the length
#' of `c14_determinations` if it is less than 10.
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
#' `c14_determinations`.  Required if `sensible_initialisation` is `FALSE`.
#' @param use_cpp DEV ONLY: Use new cpp functions (should be faster).
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
#' the length of `c14_determinations`.
#'
#' The remaining items give information about input data, input parameters (or
#' those calculated using `sensible_initialisation`) and update_type
#'
#' \describe{
#'  \item{`update_type`}{A string that always has the value "Walker".}
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
#' WalkerBivarDirichlet(kerr$c14_ages, kerr$c14_sig, intcal20, n_iter=100, n_thin=10)
WalkerBivarDirichlet <- function(
    c14_determinations,
    c14_sigmas,
    calibration_curve,
    n_iter = 100,
    n_thin = 10,
    slice_width = max(1000, diff(range(c14_determinations)) / 2),
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
    n_clust = min(10, length(c14_determinations)),
    use_cpp = TRUE) {

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
  # Initialise parameters
  if (sensible_initialisation) {
    initial_probabilities <- mapply(
      .ProbabilitiesForSingleDetermination,
      c14_determinations,
      c14_sigmas,
      MoreArgs = list(calibration_curve=calibration_curve))
    indices_of_max_probability = apply(initial_probabilities, 2, which.max)

    calendar_ages <- calibration_curve$calendar_age[indices_of_max_probability]
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

  # do not allow very small values of alpha as this causes crashes
  alpha <- 2

  tau <- stats::rgamma(n_clust, shape = nu1, rate = nu2)
  phi <- stats::rnorm(n_clust, mean = mu_phi, sd = 1 / sqrt(lambda * tau))

  v <- stats::rbeta(n_clust, 1, alpha)
  weight <- v * c(1, cumprod(1 - v)[-n_clust])
  cluster_identifiers <- sample(1:n_clust, num_observations, replace = TRUE)

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
  w_out <- list(weight)
  cluster_identifiers_out <- matrix(NA, nrow = n_out, ncol = num_observations)
  alpha_out <- rep(NA, length = n_out)
  n_clust_out <- rep(NA, length = n_out)
  mu_phi_out <- rep(NA, length = n_out)
  theta_out <- matrix(NA, nrow = n_out, ncol = num_observations)

  output_index <- 1
  cluster_identifiers_out[output_index, ] <- cluster_identifiers
  alpha_out[output_index] <- alpha
  n_clust_out[output_index] <- length(unique(cluster_identifiers))
  mu_phi_out[output_index] <- mu_phi
  theta_out[output_index, ] <- calendar_ages

  ##############################################################################
  ## Interpolate cal curve onto single year grid to speed up updating thetas
  integer_cal_year_curve <- InterpolateCalibrationCurve(
    1:pkg.globals$MAX_YEAR_BP, calibration_curve)
  interpolated_c14_age <- integer_cal_year_curve$c14_age
  interpolated_c14_sig <- integer_cal_year_curve$c14_sig

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
    DPMM_update <- .DPWalkerUpdate(
      theta = calendar_ages,
      w = weight,
      v = v,
      delta = cluster_identifiers,
      phi = phi,
      tau = tau,
      n_clust = n_clust,
      alpha = alpha,
      mu_phi = mu_phi,
      lambda = lambda,
      nu1 = nu1,
      nu2 = nu2)
    weight <- DPMM_update$w
    cluster_identifiers <- DPMM_update$delta
    phi <- DPMM_update$phi
    tau <- DPMM_update$tau
    v <- DPMM_update$v
    n_clust <- DPMM_update$n_clust

    alpha <- .WalkerUpdateAlpha(
      delta = cluster_identifiers,
      alpha = alpha,
      prshape = alpha_shape,
      prrate = alpha_rate)
    mu_phi <- .UpdateMuPhi(phi = phi, tau = tau, lambda = lambda, A = A, B = B)

    if (use_cpp) {
      calendar_ages = Update_calendar_ages_cpp(
        num_observations,
        as.double(calendar_ages),
        slice_width,
        slice_multiplier,
        as.integer(cluster_identifiers),
        as.double(phi),
        as.double(tau),
        as.double(c14_determinations),
        as.double(c14_sigmas),
        as.double(interpolated_c14_age),
        as.double(interpolated_c14_sig))
    } else {
      for (k in 1:num_observations) {
        calendar_ages[k] <- .SliceSample(
          x0 = calendar_ages[k],
          w = slice_width,
          m = slice_multiplier,
          prmean = phi[cluster_identifiers[k]],
          prsig = 1 / sqrt(tau[cluster_identifiers[k]]),
          c14obs = c14_determinations[k],
          c14sig = c14_sigmas[k],
          mucalallyr = interpolated_c14_age,
          sigcalallyr = interpolated_c14_sig)
      }
    }

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
    }
  }
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
    input_parameters = input_parameters)
  if (show_progress) close(progress_bar)
  return(return_list)
}
