#' Performs the Gibbs sampler to estimate the density of
#' a set of observations using a Dirichlet Mixture Density
#'
#' THIS CONSIDERS BOTH MEAN AND VARIANCE OF THE CLUSTERS TO BE UNKNOWN
#'
#' @param c14_determinations A vector containing the radiocarbon determinations.
#' @param c14_uncertainties A vector containing the radiocarbon determination
#' uncertainties. Must be the same length as `c14_determinations`.
#' @param calibration_curve A dataframe which should contain one column entitled
#' c14_age and one column entitled c14_sig.
#' This format matches [carbondate::intcal20].
#' @param lambda,nu1,nu2  Hyperparameters for the prior on the means
#' \eqn{\phi_j} and precision \eqn{\tau_j} of each individual calendar age
#' cluster \eqn{j}.
#' \deqn{(\phi_j, \tau_j)|\mu_{\phi} \sim
#' \textrm{NormalGamma}(\mu_{\phi}, \lambda, \nu_1, \nu_2)} where
#' \eqn{\mu_{\phi}} is the overall cluster centering.
#' @param A,B  Prior on \eqn{\mu_{\phi}} giving the mean and precision of the
#' overall centering \eqn{\mu_{\phi} \sim N(A, B^{-1})} i.e.
#' B small is uninformative.
#' @param alpha_type The type of prior on alpha - choose `"lognorm"` or
#' `"gamma"`.  Default is `"gamma"`.
#' @param alpha_mu,alpha_sigma Hyperparameters for the mean and sd on prior
#' for DP concentration, \eqn{\alpha}, determining the number of clusters we
#' expect to observe among our n sampled objects.
#' \eqn{\alpha \sim \textrm{Lognormal}(\mu, \sigma^2)} where \eqn{\mu, \sigma}
#' are the `alpha_mu` and `alpha_sigma`. Note these are only used if
#' `alpha_type` is `"lognorm"`.
#' @param alpha_shape,alpha_rate Hyperparameters for the shape and rate on prior
#' for DP concentration, \eqn{\alpha}, determining the number of clusters we
#' expect to observe among our n sampled objects.
#' \eqn{\alpha \sim \Gamma(\eta_1, \eta_2)} where \eqn{\eta_1, \eta_2} are
#' the `alpha_shape` and `alpha_rate`. Note these are only used if
#' `alpha_type` is `"gamma"`. \[TODO WE DON'T ACTUALLY USE THESE\]
#' @param n_iter  The number of MCMC iterations (optional). Default is 100.
#' @param n_thin  How much to thin the output (optional). 1 is no thinning,
#' a larger number is more thinning. Default is 10. Must choose an integer more
#' than 1 and not too close to `n_iter`, since after burn-in there are
#' \eqn{(n_{\textrm{iter}}/n_{\textrm{thin}})/2} samples from posterior to
#' potentially use.
#' @param calendar_ages The initial estimate for the underlying calendar ages
#' (optional). If supplied it must be a vector with the same length as
#' `c14_determinations`. Will be overridden if `sensible_initialisation` is
#' `TRUE`.
#' @param slice_width  Parameter for slice sampling (optional). Default is 200.
#' @param slice_multiplier  Integer parameter for slice sampling (optional).
#' Default is 50. Limits the slice size to `slice_multiplier * slice_width`.
#' @param n_clust The initial number of clusters (optional). Default is 10.
#' @param sensible_initialisation Whether to use sensible start values and
#' adaptive prior on \eqn{\mu_{\phi}} and  (A, B).
#' If this is `TRUE` (the default), then `calendar_ages`, `A` and `B` will be
#' overridden from any values passed in the arguments.
#' @param show_progress Whether to show a progress bar in the console during
#' execution. Default is `TRUE`.
#'
#' @return A list, containing 8 items, each of which having one dimension of
#' size \eqn{n_{\textrm{out}} = \textrm{floor}( n_{\textrm{iter}}/
#' n_{\textrm{thin}}) + 1}, each row storing the result from every
#' \eqn{n_{\textrm{thin}}}th iteration:
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
#'      a vector of length n_clust of the cluster means \eqn{\phi_j}.}
#'  \item{`tau`}{A list of length \eqn{n_{\textrm{out}}} each entry giving
#'      a vector of length n_clust of the cluster uncertainties \eqn{\tau_j}.}
#'  \item{`calendar_ages`}{An \eqn{n_{\textrm{out}}} by \eqn{n_{\textrm{obs}}}
#'     integer matrix. Gives the calendar age for each observation.}
#'  \item{`mu_phi`}{A vector of length \eqn{n_{\textrm{out}}} giving the overall
#'      centering \eqn{\mu_{\phi}} of the clusters.}
#'  \item{`update_type`}{A string that always has the value "neal". This is
#'      needed if the output is used as an input for post-processing functions.}
#' }
#' where \eqn{n_{\textrm{obs}}} is the number of radiocarbon observations i.e.
#' the length of `c14_determinations`.
#' @export
#'
BivarGibbsDirichletwithSlice <- function(
    c14_determinations,
    c14_uncertainties,
    calibration_curve,
    lambda,
    nu1,
    nu2,
    A,
    B,
    alpha_type = "gamma",
    alpha_mu = NA,
    alpha_sigma = NA,
    alpha_shape = NA,
    alpha_rate = NA,
    n_iter = 100,
    n_thin = 10,
    calendar_ages = NA,
    slice_width = 200,
    slice_multiplier = 50,
    n_clust = 10,
    sensible_initialisation = TRUE,
    show_progress = TRUE) {

  ##############################################################################
  # Check input parameters

  # TODO

  ##############################################################################
  # Initialise parameters
  num_observations <- length(c14_determinations)

  all_clusters_represented <- FALSE
  while (!all_clusters_represented) {
    cluster_identifiers <- sample(1:n_clust, num_observations, replace = TRUE)
    all_clusters_represented <- (length(unique(cluster_identifiers)) == n_clust)
  }

  if (sensible_initialisation) {
    initial_probabilities <- mapply(
      CalibrateSingleDetermination,
      c14_determinations,
      c14_uncertainties,
      MoreArgs = list(calibration_curve=calibration_curve))
    indices_of_max_probability = apply(initial_probabilities, 2, which.max)
    calendar_ages <- calibration_curve$calendar_age[indices_of_max_probability]
    mu_phi <- stats::median(calendar_ages)
    A <- stats::median(calendar_ages)
    B <- 1 / (max(calendar_ages) - min(calendar_ages))^2
  } else {
    scale_val <- 8267 / 8033
    mu_phi <- mean(c14_determinations) * scale_val
    if (is.na(calendar_ages[1])) calendar_ages <- c14_determinations * scale_val
  }

  alpha <- switch(
    alpha_type,
    lognorm = exp(stats::rnorm(1, alpha_mu, sd = alpha_sigma)),
    gamma = 0.0001, # stats::rgamma(1, shape = alpha_shape, rate = alpha_rate),
    stop("Unknown form for prior on gamma"))

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
        c14sig = c14_uncertainties[k],
        mucalallyr = interpolated_c14_age,
        sigcalallyr = interpolated_c14_sig)
    }

    alpha <- switch(
      alpha_type,
      lognorm = .UpdateAlphaLognormPrior(
        cluster_identifiers, alpha, mualpha = alpha_mu, sigalpha = alpha_sigma),
      gamma = .UpdateAlphaGammaPrior(
        cluster_identifiers, alpha, prshape = alpha_shape, prrate = alpha_rate),
      stop("Unknown form for prior on gamma"))

    if (iter %% n_thin == 0) {
      output_index <- output_index + 1
      phi_out[[output_index]] <- phi
      tau_out[[output_index]] <- tau
      cluster_identifiers_out[output_index, ] <- cluster_identifiers
      calendar_ages_out[output_index, ] <- calendar_ages
      alpha_out[output_index] <- alpha
      mu_phi_out[output_index] <- mu_phi
      n_clust_out[output_index] <- length(unique(cluster_identifiers))
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
    update_type="neal")

  if (show_progress) close(progress_bar)
  return(return_list)
}
