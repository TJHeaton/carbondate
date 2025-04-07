#' Plot Individual Realisations of Posterior Rate of Sample Occurrence for Poisson Process Model
#'
#' @description
#' Given output from the Poisson process fitting function [carbondate::PPcalibrate] plot
#' individual realisations from the MCMC for the rate of sample occurrence (i.e., realisations
#' of the underlying Poisson process rate \eqn{\lambda(t)}), on a given calendar age grid
#' (provided in cal yr BP). Specify either `n_realisations` if you want to select a random set
#' of realisations, or `realisations` if you want to provide a vector of specific realisations.
#'
#' @inheritParams PlotPosteriorMeanRate
#' @param n_realisations Number of randomly sampled realisations to be drawn from MCMC posterior
#' and plotted. Default is 10.
#' @param plot_realisations_colour The colours to be used to plot the individual realisations.
#' Default is greyscale (otherwise should have same length as number of realisations).
#' @param realisations Specific indices of realisations (in thinned version) to plot if user does not
#' want to sample realisations randomly). If specified will override `n_realisations`.
#' @param interval_width The confidence intervals to show for the
#' calibration curve. Choose from one of `"1sigma"` (68.3%),
#' `"2sigma"` (95.4%) and `"bespoke"`. Default is `"2sigma"`.
#' @param bespoke_probability The probability to use for the confidence interval
#' if `"bespoke"` is chosen above. E.g., if 0.95 is chosen, then the 95% confidence
#' interval is calculated. Ignored if `"bespoke"` is not chosen.
#'
#' @return None
#'
#' @export
#'
#' @examples
#' #' # NOTE: All these examples are shown with a small n_iter and n_posterior_samples
#' # to speed up execution.
#' # Try n_iter and n_posterior_samples as the function defaults.
#'
#' pp_output <- PPcalibrate(
#'     pp_uniform_phase$c14_age,
#'     pp_uniform_phase$c14_sig,
#'     intcal20,
#'     n_iter = 1000,
#'     show_progress = FALSE)
#'
#' # Plot 10 random realisations in greyscale
#' PlotRateIndividualRealisation(
#'     pp_output,
#'     n_realisations = 10)
#'
#' # Plot three random realisations with specific colours
#' PlotRateIndividualRealisation(
#'     pp_output,
#'     n_realisations = 3,
#'     plot_realisations_colour = c("red", "green", "purple"))
#'
#' # Plot some specific realisations
#' PlotRateIndividualRealisation(
#'     pp_output,
#'     realisations = c(60, 73, 92),
#'     plot_realisations_colour = c("red", "green", "purple"))
PlotRateIndividualRealisation <- function(
    output_data,
    n_realisations = 10,
    plot_realisations_colour = NULL,
    realisations = NULL,
    calibration_curve = NULL,
    plot_14C_age = TRUE,
    plot_cal_age_scale = "BP",
    interval_width = "2sigma",
    bespoke_probability = NA,
    denscale = 3,
    resolution = 1,
    n_burn = NA,
    n_end = NA,
    plot_pretty = TRUE,
    plot_lwd = 2) {


  n_iter <- output_data$input_parameters$n_iter
  n_thin <- output_data$input_parameters$n_thin
  n_real <- length(output_data$rate_s)

  arg_check <- .InitializeErrorList()
  .CheckOutputData(arg_check, output_data, "RJPP")
  .CheckInteger(arg_check, n_realisations)
  .CheckRealisations(arg_check, realisations, lower = 1, upper = n_real)
  .CheckCalibrationCurveFromOutput(arg_check, output_data, calibration_curve)
  .CheckFlag(arg_check, plot_14C_age)
  .CheckChoice(arg_check, plot_cal_age_scale, c("BP", "AD", "BC"))
  .CheckIntervalWidth(arg_check, interval_width, bespoke_probability)
  .CheckNumber(arg_check, denscale, lower = 0)
  .CheckNumber(arg_check, resolution, lower = 0.01)
  .ReportErrors(arg_check)

  n_burn <- .SetNBurn(n_burn, n_iter, n_thin)
  n_end <- .SetNEnd(n_end, n_iter, n_thin)

  # Ensure revert to main environment par on exit of function
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  # Set nice plotting parameters
  if(plot_pretty) {
    graphics::par(
      mgp = c(3, 0.7, 0),
      xaxs = "i",
      yaxs = "i",
      mar = c(5, 4.5, 4, 2) + 0.1,
      las = 1)
  }

  if (is.null(calibration_curve)) {
    calibration_curve <- get(output_data$input_data$calibration_curve_name)
  }
  rc_determinations <- output_data$input_data$rc_determinations
  rc_sigmas <- output_data$input_data$rc_sigmas
  F14C_inputs <-output_data$input_data$F14C_inputs

  if (plot_14C_age == TRUE) {
    calibration_curve <- .AddC14ageColumns(calibration_curve)
    if (F14C_inputs == TRUE) {
      converted <- .ConvertF14cTo14Cage(rc_determinations, rc_sigmas)
      rc_determinations <- converted$c14_age
    }
  } else {
    calibration_curve <- .AddF14cColumns(calibration_curve)
    if (F14C_inputs == FALSE) {
      converted <- .Convert14CageToF14c(rc_determinations, rc_sigmas)
      rc_determinations <- converted$f14c
    }
  }

  ##############################################################################
  # Initialise plotting parameters
  calibration_curve_colour <- "blue"
  calibration_curve_bg <- grDevices::rgb(0, 0, 1, .3)
  output_colour <- "purple"

  start_age <- ceiling(min(output_data$rate_s[[1]]) / resolution) * resolution
  end_age <- floor(max(output_data$rate_s[[1]]) / resolution) * resolution
  if (end_age == max(output_data$rate_s[[1]])) {
    # Removes issue of sequence coinciding with end changepoint
    end_age <- end_age - resolution
  }

  calendar_age_sequence <- seq(from = start_age, to = end_age, by = resolution)
  xlim <- rev(range(calendar_age_sequence))

  ##############################################################################

  # Sample realisations if not specific
  if(is.null(realisations)) {
    realisations <- sample((n_burn + 1):n_end,
                           n_realisations,
                           replace = ((n_end - n_burn) < n_realisations))
  }
  n_realisations <- length(realisations)

  if(!is.null(plot_realisations_colour) && (length(plot_realisations_colour) != n_realisations)) {
    warning("Plot_realisations_colour is not the same length as the number of realisations that you have selected, overriding colour choice to be greyscale",
            call. = FALSE)
  }

  if(is.null(plot_realisations_colour) || length(plot_realisations_colour) != n_realisations) {
    plot_realisations_colour <- rep(grey(0.4, alpha = 0.6), n_realisations)
  }

  # Calculate rate for each realisation
  rate <- matrix(NA, nrow = n_realisations, ncol = length(calendar_age_sequence))
  for (i in 1:n_realisations) {
    ind <- realisations[i]
    rate[i,] <- stats::approx(
      x = output_data$rate_s[[ind]],
      y = c(output_data$rate_h[[ind]], 0),
      xout = calendar_age_sequence,
      method = "constant")$y
  }


  ylim_rate <- c(0, denscale * max(rate))

  ##############################################################################
  # Plot curves
  plot_title <- expression(paste("Posterior realisations of rate ", lambda, "(t)"))

  .PlotCalibrationCurveAndInputData(
    plot_cal_age_scale,
    xlim,
    calibration_curve,
    rc_determinations,
    plot_14C_age,
    calibration_curve_colour,
    calibration_curve_bg,
    interval_width,
    bespoke_probability,
    title = plot_title)

  .SetUpDensityPlot(plot_cal_age_scale, xlim, ylim_rate)

  # Plot realisations
  cal_age <- .ConvertCalendarAge(plot_cal_age_scale, calendar_age_sequence)
  for(i in 1:n_realisations) {
    graphics::lines(cal_age, rate[i,], col = plot_realisations_colour[i], lwd = plot_lwd)
  }

  .AddRealisationLegendToRatePlot(
    output_data,
    interval_width,
    bespoke_probability,
    calibration_curve_colour,
    output_colour = plot_realisations_colour[1])

}


.AddRealisationLegendToRatePlot <- function(
    output_data,
    interval_width,
    bespoke_probability,
    calibration_curve_colour,
    output_colour) {

  ci_label <- switch(
    interval_width,
    "1sigma" = expression(paste("1", sigma, " interval")),
    "2sigma"  = expression(paste("2", sigma, " interval")),
    "bespoke" = paste0(round(100 * bespoke_probability), "% interval"))

  legend_labels <- c(
    gsub("intcal", "IntCal",
         gsub("shcal", "SHCal",
              output_data$input_data$calibration_curve_name)), # Both IntCal and SHCal
    ci_label,
    "Realisations")
  lty <- c(1, 2, 1)
  pch <- c(NA, NA, NA)
  col <- c(calibration_curve_colour, calibration_curve_colour, output_colour)

  graphics::legend(
    "topright", legend = legend_labels, lty = lty, pch = pch, col = col)
}


