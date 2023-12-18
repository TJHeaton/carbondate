#' Title
#'
#' @param rc_determinations A vector of observed radiocarbon determinations
#'
#' @param rc_sigmas A vector of the radiocarbon determinations
#' uncertainties (1-sigma). Must be the same length as `rc_determinations`.
#'
#' @param calibration_curve A dataframe which must contain one column `calendar_age_BP`, and also
#' columns `c14_age` and `c14_sig` or `f14c` and `f14c_sig` (or both sets).
#' This format matches the curves supplied with this package
#'
#' @param F14C_inputs `TRUE` if the provided rc_determinations are F14C concentrations and `FALSE`
#' if they are radiocarbon age BP. Defaults to `FALSE`.
#'
#' @param n_iter  The number of MCMC iterations (optional). Default is 100,000.
#'
#' @param n_thin  How much to thin the output (optional). 1 is no thinning,
#' a larger number is more thinning. Default is 10. Must choose an integer more
#' than 1 and not too close to `n_iter`, to ensure there are enought samples from
#' posterior to potentially use.
#'
#' @param use_F14C_space Whether the calculations are carried out in F14C space (default is TRUE).
#' If FALSE, calculations are carried out in 14C yr BP space.
#'
#' @param show_progress Whether to show a progress bar in the console during
#' execution. Default is `TRUE`.
#'
#' @param calendar_age_range Minimum and maximum calendar ages permitted
#' for the calendar ages of the samples, i.e. range_1 < theta < range_2.
#'
#' @param calendar_grid_resolution The spacing of the calendar age grid on which to consider
#' the ages of the samples, e.g. t, t + resolution, t + 2 * resolution, ..
#'
#' @param prior_h_shape,prior_h_rate Prior for Poisson Process rate height in any interval
#' rate_h ~ Gamma(shape = prior_h_shape, rate = prior_h_rate)
#' If choose sensible_initialisation then prior_h_shape is selected adaptively (internally)
#' to match n_observations (and selecting prior_h_rate = 0.1 to provide diffuse prior)
#'
#' @param prior_n_internal_changepoints_lambda Prior on Poisson parameter specifying
#' n_internal_changepoints ~ Po(prior_n_internal_changepoints_lambda)
#'
#' @param k_max_internal_changepoints Maximum permitted number of internal changepoints
#'
#' @param rescale_factor_rev_jump Factor weighting probability of dimension change
#' in the reversible jump update step for poisson process rate_h and rate_s
#'
#' @param bounding_range_prob_cutoff Probability cut-off for choosing the bounds for the
#' potential calendar ages for the observations
#'
#' @param default_prior_h_rate Default value for rate parameter in Gamma prior for rate
#' if not specifically passed by user alongside prior_h_shape.
#' Important to be able to set independently as determines the diffusivity of the prior
#' Var(rate_h) = (1/default_prior_h_rate) * (n_determinations / calendar_age_interval_length)
#' Smaller values of default_prior_h_rate will mean more diverse prior
#'
#' @param initial_n_internal_changepoints Number of internal changepoints to initialise
#' the sampler with. The default is 10 (so diffuse).
#' Will place them randomly within calendar interval
#'
#' @param grid_extension_factor If you adaptively select the calendar_age_range from
#' rc_determinations, how far you wish to extend the grid beyond this adaptive minimum and maximum
#' The final range will be extended (eqaully on both sides) to cover (1 + grid_extension_factor) * (calendar_age_range)
#'
#' @param use_fast Adding a flag to allow trimming the likelihoods based on prob cutoff
#'
#' @return TODO
#' @export
#'
#' @examples # TODO
PPcalibrate <- function(
    rc_determinations,
    rc_sigmas,
    calibration_curve,
    F14C_inputs = FALSE,
    n_iter = 1e5,
    n_thin = 10,
    use_F14C_space = TRUE,
    show_progress = TRUE,
    calendar_age_range = NA,
    calendar_grid_resolution = 1,
    prior_h_shape = NA,
    prior_h_rate = NA,
    prior_n_internal_changepoints_lambda = 10,
    k_max_internal_changepoints = 30,
    rescale_factor_rev_jump = 0.9,
    bounding_range_prob_cutoff = 0.005,
    default_prior_h_rate = 0.1,
    initial_n_internal_changepoints = 10,
    grid_extension_factor = 0.1,
    use_fast = TRUE) {

  # TODO - Check both prior_h_shape and prior_h_rate specified (or both NA)
  # TODO - Check initial_n_internal_changepoints < k_max
  # TODO - Check parameters are within required bounds e.g. +ve/-ve

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

  # Find initial calendar_age_range for Poisson Process
  # Have converted rc_determinations to use_F14C_space already
  bounds_calendar_range <- .FindBoundingCalendarRange(
    rc_determinations = rc_determinations,
    rc_sigmas = rc_sigmas,
    calibration_curve = calibration_curve,
    F14C_inputs = use_F14C_space,
    prob_cutoff = bounding_range_prob_cutoff)

  if(any(is.na(calendar_age_range))) {
    min_potential_calendar_age <- min(bounds_calendar_range)
    max_potential_calendar_age <- max(bounds_calendar_range)

    # Extend by grid_extension_factor
    mean_potential_calendar_age <- mean(bounds_calendar_range)
    min_potential_calendar_age <- (
      min_potential_calendar_age - grid_extension_factor * (
        mean_potential_calendar_age - min_potential_calendar_age
      )
    )
    max_potential_calendar_age <- (
      max_potential_calendar_age + grid_extension_factor * (
        max_potential_calendar_age - mean_potential_calendar_age
      )
    )

    # Extend to discretised calendar_grid_resolution
    min_potential_calendar_age <- (
      floor(min_potential_calendar_age/calendar_grid_resolution) *
        calendar_grid_resolution
    )
    max_potential_calendar_age <- (
      ceiling(max_potential_calendar_age/calendar_grid_resolution) *
        calendar_grid_resolution
    )

    # Do not extend beyond calibration curve calendar limits
    min_potential_calendar_age <- max(min_potential_calendar_age,
                                      min(calibration_curve$calendar_age_BP))
    max_potential_calendar_age <- min(max_potential_calendar_age,
                                      max(calibration_curve$calendar_age_BP))
  } else {
    ## Create calendar_age_grid covering potential calendar ages
    min_potential_calendar_age <- min(calendar_age_range)
    max_potential_calendar_age <- max(calendar_age_range)
    # Provide warning if these are narrower than the estimated bounds
    if(min_potential_calendar_age > min(bounds_calendar_range)) {
      warning("The minimum of your selected calendar age range may not cover the age range of the samples. Consider reducing the minimum age range." , immediate. = TRUE)
    }
    if(max_potential_calendar_age < max(bounds_calendar_range)) {
      warning("The maximum of your selected calendar age range may not cover the age range of the samples. COnsider increasing the maximum age range." , immediate. = TRUE)
    }
  }

  calendar_age_grid <- seq(
    min_potential_calendar_age,
    max_potential_calendar_age,
    by = calendar_grid_resolution) # May end before max_potential_calendar_age

  # Ensure end of calendar_age_grid extends at least to max_potential_calendar_age
  # If not extend calendar_age_grid and adjust max_potential_calendar_age to match
  if((max(calendar_age_grid) != max_potential_calendar_age)) {
    if( !any(is.na(calendar_age_range)) ) { # Only inform if have user selected values
      warning("Extending calendar age range for Poisson Process as selected grid resolution doesn't divide age range.", immediate. = TRUE)
    }
    max_potential_calendar_age <- (
      calendar_age_grid[length(calendar_age_grid)] + calendar_grid_resolution
    )
    cc <- c(calendar_age_grid,
                           max_potential_calendar_age)
  }

  calendar_age_interval_length <- max_potential_calendar_age - min_potential_calendar_age
  n_determinations <- length(rc_determinations)

  initial_estimate_mean_rate <- n_determinations / calendar_age_interval_length

  if(is.na(prior_h_shape)) {
    # Choose exponential distribution
    # and match mean with initial_estimate_mean_rate above
    prior_h_shape <- 1
    prior_h_rate <- 1 / initial_estimate_mean_rate
  }

  #############
  ## Create initial change points and heights for Poisson process rate
  ## Randomly place initial_n_internal_changepoints
  initial_rate_s <- sort(
    c(
      min_potential_calendar_age,
      stats::runif(initial_n_internal_changepoints,
            min = min_potential_calendar_age,
            max = max_potential_calendar_age),
      max_potential_calendar_age
    )
  )

  # Create diverse initial_rate_h centred on initial_estimate_mean_rate
  # Requires safe setting as as will crash if any sampled rate_h values are zero
  initial_rate_h <- rep(initial_estimate_mean_rate,
                        times = initial_n_internal_changepoints + 1)
  initial_rate_h_log_multiplier <- stats::runif(n = initial_n_internal_changepoints + 1,
                                                min = -1,
                                                max = 1)
  initial_rate_h <- exp(initial_rate_h_log_multiplier) * initial_rate_h

  initial_integrated_rate <- .FindIntegral(
    rate_s = initial_rate_s,
    rate_h = initial_rate_h
  )


  ##############################################################################
  # Save input parameters
  input_parameters <- list(
    pp_cal_age_range = c(min_potential_calendar_age,
                         max_potential_calendar_age),
    prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
    k_max_internal_changepoints = k_max_internal_changepoints,
    prior_h_shape = prior_h_shape,
    prior_h_rate = prior_h_rate,
    rescale_factor_rev_jump = rescale_factor_rev_jump,
    calendar_age_grid = calendar_age_grid,
    calendar_grid_resolution = calendar_grid_resolution,
    n_iter = n_iter,
    n_thin = n_thin)

  ####################################
  ## Create matrix of calendar_likelihoods (stored in main as not updated throughout samples)
  ## Different from .ProbabilitiesForSingleDetermination as not normalised
  ## and can pass theta on different grid to calibration_curve
  likelihood_calendar_ages_from_calibration_curve <- mapply(
    .CalendarAgeLikelihoodGivenCurve,
    rc_determinations,
    rc_sigmas,
    MoreArgs = list(
      theta = calendar_age_grid,
      F14C_inputs = use_F14C_space,
      calibration_curve = calibration_curve)
  )

  if (use_fast) {
    trimmed_likelihood_calendar_ages_from_calibration_curve <- apply(
      likelihood_calendar_ages_from_calibration_curve,
      MARGIN = 2,
      FUN = .FindTrimmedVectorAndIndices,
      prob_cutoff = bounding_range_prob_cutoff,
      simplify = FALSE)
  }

  num_observations <- length(rc_determinations)

  # Set starting values to be initialised ones
  rate_s <- initial_rate_s
  rate_h <- initial_rate_h
  integrated_rate <- initial_integrated_rate

  prob_move <- .FindMoveProbability(
    prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
    k_max_internal_changepoints = k_max_internal_changepoints,
    rescale_factor = rescale_factor_rev_jump)

  ####################################
  # Create storage for output
  n_out <- floor(n_iter / n_thin) + 1

  rate_s_out <- list(rate_s)
  rate_h_out <- list(rate_h)
  n_internal_changes <- rep(NA, length = n_out)
  theta_out <- matrix(NA, nrow = n_out, ncol = num_observations)

  output_index <- 1
  n_internal_changes[output_index] <- length(rate_h) - 1

  ## Store calendar_ages given initial_rate_s and initial_rate_h
  ## (sample from exactly using Gibbs)
  if (use_fast) {
    calendar_ages <- .TrimmedUpdateCalendarAgesGibbs(
      trimmed_likelihood_calendar_ages_from_calibration_curve = trimmed_likelihood_calendar_ages_from_calibration_curve,
      calendar_age_grid = calendar_age_grid,
      rate_s = rate_s,
      rate_h = rate_h
    )
  } else {
    calendar_ages <- UpdateCalendarAgesGibbs(
      likelihood_calendar_ages_from_calibration_curve = likelihood_calendar_ages_from_calibration_curve,
      calendar_age_grid = calendar_age_grid,
      rate_s = rate_s,
      rate_h = rate_h
    )
  }

  theta_out[output_index, ] <- calendar_ages



  #####################################
  # Perform MCMC - RJMCMC within Gibbs
  # Consist of iterating between:
  #    i) Updating calendar_ages given lambda (rate_s, rate_h) and rc_determinations
  #    ii) Updating lambda (rate_s, rate_h) given calendar_ages using RJ MCMC

  if (show_progress) {
    progress_bar <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
  }


  for(iter in 1:n_iter) {

    ## Step 1: Update calendar_ages given rate_s and rate_h (sample from exactly using Gibbs)
    if (use_fast) {
      calendar_ages <- .TrimmedUpdateCalendarAgesGibbs(
        trimmed_likelihood_calendar_ages_from_calibration_curve = trimmed_likelihood_calendar_ages_from_calibration_curve,
        calendar_age_grid = calendar_age_grid,
        rate_s = rate_s,
        rate_h = rate_h
      )
    } else {
      calendar_ages <- UpdateCalendarAgesGibbs(
        likelihood_calendar_ages_from_calibration_curve = likelihood_calendar_ages_from_calibration_curve,
        calendar_age_grid = calendar_age_grid,
        rate_s = rate_s,
        rate_h = rate_h
      )
    }

    ## Step 2: Update rate_s and rate_h given current calendar_ages (using RJMCMC)
    updated_poisson_process <- UpdatePoissonProcessRateRevJump(
      theta = calendar_ages,
      rate_s = rate_s,
      rate_h = rate_h,
      integrated_rate = integrated_rate,
      prior_h_shape = prior_h_shape,
      prior_h_rate = prior_h_rate,
      prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
      prob_move = prob_move
    )

    rate_s <- updated_poisson_process$rate_s
    rate_h <- updated_poisson_process$rate_h
    integrated_rate <- updated_poisson_process$integrated_rate

    # Store output
    if (iter %% n_thin == 0) {
      output_index <- output_index + 1
      n_internal_changes[output_index] <- length(rate_h) - 1
      rate_h_out[[output_index]] <- rate_h
      rate_s_out[[output_index]] <- rate_s
      theta_out[output_index, ] <- calendar_ages
    }

    if (show_progress) {
      if (iter %% 100 == 0) {
        utils::setTxtProgressBar(progress_bar, iter)
      }
    }
  }

  return_list <- list(
    rate_s = rate_s_out,
    rate_h = rate_h_out,
    calendar_ages = theta_out,
    n_internal_changes = n_internal_changes,
    update_type = "RJPP",
    input_data = input_data,
    input_parameters = input_parameters)

  if (show_progress) close(progress_bar)
  return(return_list)
}

## TODO - Think about how to choose min and max cut-off currently ends where there is an observation so will have higher rate
