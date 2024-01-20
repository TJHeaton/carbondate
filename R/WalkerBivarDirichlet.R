#' Calibrate and Summarise Multiple Radiocarbon Samples via
#' a Bayesian Non-Parametric DPMM (with Walker Updating)
#'
#' @description
#' This function calibrates sets of multiple radiocarbon (\eqn{{}^{14}}C)
#' determinations, and simultaneously summarises the resultant calendar age information.
#' This is achieved using Bayesian non-parametric density estimation,
#' providing a statistically-rigorous alternative to summed probability
#' distributions (SPDs).
#'
#' It takes as an input a set of \eqn{{}^{14}}C determinations and associated \eqn{1\sigma}
#' uncertainties, as well as the radiocarbon age calibration curve to be used. The samples
#' are assumed to arise from an (unknown) shared calendar age distribution \eqn{f(\theta)} that
#' we would like to estimate, alongside performing calibration of each sample.
#'
#' The function models the underlying distribution \eqn{f(\theta)} as a Dirichlet process
#' mixture model (DPMM), whereby the samples are considered to arise from an unknown number of
#' distinct clusters. Fitting is achieved via MCMC.
#'
#' It returns estimates for the calendar age of each individual radiocarbon sample; and broader
#' output (the weights, means and variances of the underpinning calendar age clusters)
#' that can be used by other library functions to provide a predictive estimate of the
#' shared calendar age density \eqn{f(\theta)}.
#'
#' For more information read the vignette: \cr
#' \code{vignette("Non-parametric-summed-density", package = "carbondate")}
#'
#' \strong{Note:} The library provides two slightly-different update schemes for the MCMC. In this particular function, updating of the DPMM is achieved by slice sampling
#' (Walker 2007). We recommend use of the alternative to this, implemented at [carbondate::PolyaUrnBivarDirichlet]
#'
#' \strong{Reference:} \cr
#' Heaton, TJ. 2022. “Non-parametric Calibration of Multiple Related Radiocarbon
#' Determinations and their Calendar Age Summarisation.” \emph{Journal of the Royal Statistical
#' Society Series C: Applied Statistics} \strong{71} (5):1918-56. https://doi.org/10.1111/rssc.12599. \cr
#' Walker, SG. 2007. “Sampling the Dirichlet Mixture Model with Slices.” \emph{Communications in
#' Statistics - Simulation and Computation} \strong{36} (1):45-54. https://doi.org/10.1080/03610910601096262.
#'
#' @inheritParams FindSummedProbabilityDistribution
#' @param n_iter  The number of MCMC iterations (optional). Default is 100,000.
#' @param n_thin  How much to thin the MCMC output (optional). Will store every
#' `n_thin`\eqn{{}^\textrm{th}} iteration. 1 is no thinning, while a larger number will result
#' in more thinning. Default is 10. Must choose an integer greater than 1. Overall
#' number of MCMC realisations stored will be \eqn{n_{\textrm{out}} =
#' \textrm{floor}( n_{\textrm{iter}}/n_{\textrm{thin}}) + 1} so do not choose
#'  `n_thin` too large to ensure there are enough samples from the posterior to
#'  use for later inference.
#' @param use_F14C_space If `TRUE` (default) the calculations within the function are carried
#' out in F\eqn{{}^{14}}C space. If `FALSE` they are carried out in \eqn{{}^{14}}C
#' age space. We recommend selecting `TRUE` as, for very old samples, calibrating in
#' F\eqn{{}^{14}}C space removes affect of potential asymmetry in radiocarbon ages.
#' \emph{Note:} This flag can be set independently of the format/scale on which
#' `rc_determinations` were originally provided.
#' @param slice_width Parameter for slice sampling (optional). If not given a value
#' is chosen intelligently based on the spread of the initial calendar ages.
#' Must be given if `sensible_initialisation` is `FALSE`.
#' @param slice_multiplier  Integer parameter for slice sampling (optional).
#' Default is 10. Limits the slice size to `slice_multiplier * slice_width`.
#' @param n_clust The number of clusters with which to initialise the sampler (optional). Must
#' be less than the length of `rc_determinations`. Default is 10 or the length
#' of `rc_determinations` if that is less than 10.
#' @param show_progress Whether to show a progress bar in the console during
#' execution. Default is `TRUE`.
#' @param sensible_initialisation Whether to use sensible values to initialise the sampler
#' and an automated (adaptive) prior on \eqn{\mu_{\phi}} and (A, B) that is informed
#' by the observed `rc_determinations`. If this is `TRUE` (the recommended default), then
#' all the remaining arguments below are ignored.
#' @param lambda,nu1,nu2  Hyperparameters for the prior on the mean
#' \eqn{\phi_j} and precision \eqn{\tau_j} of each individual calendar age
#' cluster \eqn{j}:
#' \deqn{(\phi_j, \tau_j)|\mu_{\phi} \sim
#' \textrm{NormalGamma}(\mu_{\phi}, \lambda, \nu_1, \nu_2)} where
#' \eqn{\mu_{\phi}} is the overall cluster centering. Required if
#' `sensible_initialisation` is `FALSE`.
#' @param A,B  Prior on \eqn{\mu_{\phi}} giving the mean and precision of the
#' overall centering \eqn{\mu_{\phi} \sim N(A, B^{-1})}. Required if `sensible_initialisation` is `FALSE`.
#' @param mu_phi Initial value of the overall cluster centering \eqn{\mu_{\phi}}.
#' Required if `sensible_initialisation` is `FALSE`.
#' @param alpha_shape,alpha_rate Shape and rate hyperparameters that specify
#' the prior for the Dirichlet Process (DP) concentration, \eqn{\alpha}. This
#' concentration \eqn{\alpha} determines the number of clusters we
#' expect to observe among our \eqn{n} sampled objects. The model
#' places a prior on
#' \eqn{\alpha \sim \Gamma(\eta_1, \eta_2)}, where \eqn{\eta_1, \eta_2} are
#' the `alpha_shape` and `alpha_rate`. A small \eqn{\alpha} means the DPMM is
#' more concentrated (i.e. we expect fewer calendar age clusters) while a large alpha means it is less
#' less concentrated (i.e. many clusters). Required if `sensible_initialisation` is `FALSE`.
#' @param calendar_ages  The initial estimate for the underlying calendar ages
#' (optional). If supplied, it must be a vector with the same length as
#' `rc_determinations`.  Required if `sensible_initialisation` is `FALSE`.
#'
#' @return A list with 11 items. The first 8 items contain output of the model, each of
#' which has one dimension of size \eqn{n_{\textrm{out}} =
#' \textrm{floor}( n_{\textrm{iter}}/n_{\textrm{thin}}) + 1}. The rows in these items store
#' the state of the MCMC from every \eqn{n_{\textrm{thin}}}\eqn{{}^\textrm{th}} iteration:
#'
#' \describe{
#'  \item{`cluster_identifiers`}{An \eqn{n_{\textrm{out}}} by
#'     \eqn{n_{\textrm{obs}}} integer matrix. Provides the cluster allocation
#'      (an integer between 1 and `n_clust`) for each observation on the relevant MCMC iteration.
#'      Information on the state of these calendar age clusters (means, precisions, and weights)
#'      can be found in the other output items.}
#'  \item{`alpha`}{A double vector of length \eqn{n_{\textrm{out}}} giving the Dirichlet Process
#'     concentration parameter \eqn{\alpha}.}
#'  \item{`n_clust`}{An integer vector of length \eqn{n_{\textrm{out}}} giving
#'      the current number of clusters in the model.}
#'  \item{`phi`}{A list of length \eqn{n_{\textrm{out}}} each entry giving
#'      a vector of the means of the current calendar age clusters \eqn{\phi_j}.}
#'  \item{`tau`}{A list of length \eqn{n_{\textrm{out}}} each entry giving
#'      a vector of the precisions of the current calendar age clusters \eqn{\tau_j}.}
#'  \item{`weight`}{A list of length \eqn{n_{\textrm{out}}} each entry giving
#'      the mixing weights of each calendar age cluster.}
#'  \item{`calendar_ages`}{An \eqn{n_{\textrm{out}}} by \eqn{n_{\textrm{obs}}}
#'     integer matrix. Gives the current estimate for the calendar age of each individual
#'     observation.}
#'  \item{`mu_phi`}{A vector of length \eqn{n_{\textrm{out}}} giving the overall
#'      centering \eqn{\mu_{\phi}} of the calendar age clusters.}
#' }
#' where \eqn{n_{\textrm{obs}}} is the number of radiocarbon observations, i.e.,
#' the length of `rc_determinations`.
#'
#' The remaining items give information about the input data, input parameters (or
#' those calculated using `sensible_initialisation`) and the update_type
#'
#' \describe{
#'  \item{`update_type`}{A string that always has the value "Walker".}
#'  \item{`input_data`}{A list containing the \eqn{{}^{14}}C data used, and the name of
#'  the calibration curve used.}
#'  \item{`input_parameters`}{A list containing the values of the fixed
#'  hyperparameters `lambda`, `nu1`, `nu2`, `A`, `B`, `alpha_shape`, and
#'  `alpha_rate`, and the slice parameters `slice_width` and
#'  `slice_multiplier`.}
#' }
#'
#' @export
#'
#' @seealso
#' [carbondate::PolyaUrnBivarDirichlet] for our preferred MCMC method to update the Bayesian DPMM
#' (otherwise an identical model); and [carbondate::PlotCalendarAgeDensityIndividualSample],
#' [carbondate::PlotPredictiveCalendarAgeDensity] and [carbondate::PlotNumberOfClusters]
#' to access the model output and estimate the calendar age information. \cr \cr
#' See also [carbondate::PPcalibrate] for an an alternative (similarly rigorous) approach to
#' calibration and summarisation of related radiocarbon determinations using a variable-rate Poisson process
#'
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
