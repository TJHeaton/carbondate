#' Plot Posterior Mean Rate of Sample Occurrence for Poisson Process Model
#'
#' @description
#' Given output from the Poisson process fitting function [carbondate::PPcalibrate] calculate
#' and plot the posterior mean rate of sample occurrence (i.e., the underlying Poisson process
#' rate \eqn{\lambda(t)}) together with specified probability intervals, on a given calendar age grid
#' (provided in cal yr BP).
#'
#' Will show the original set of radiocarbon determinations (those you are modelling/summarising),
#' the chosen calibration curve, and the estimated posterior rate of occurrence \eqn{\lambda(t)} on the same plot.
#' Can also optionally show the posterior mean of each individual sample's calendar age estimate.
#'
#' An option is also provided to calculate the posterior mean rate \emph{conditional}
#' upon the number of internal changepoints within the period under study (if this is specified as an input
#' to the function).
#'
#' \strong{Note:} If all you are interested in is the value of the posterior mean rate
#' on a grid, without an accompanying plot, you can use
#' [carbondate::FindPosteriorMeanRate] instead.
#'
#' For more information read the vignette: \cr
#' \code{vignette("Poisson-process-modelling", package = "carbondate")}
#'
#' @param output_data The return value from the updating function
#' [carbondate::PPcalibrate]. Optionally, the output data can have an extra list item
#' named `label` which is used to set the label on the plot legend.
#' @param n_posterior_samples Number of samples it will draw, after having removed `n_burn`,
#' from the (thinned) MCMC realisations stored in `output_data` to estimate the
#' rate \eqn{\lambda(t)}. These samples may be repeats if the number of, post burn-in,
#' realisations is less than `n_posterior_samples`. If not given, 5000 is used.
#' @param n_changes (Optional) If wish to condition calculation of the posterior mean on
#' the number of internal changepoints. In this function, if specified, `n_changes` must
#' be a single integer.
#' @param calibration_curve This is usually not required since the name of the
#' calibration curve variable is saved in the output data. However, if the
#' variable with this name is no longer in your environment then you should pass
#' the calibration curve here. If provided, this should be a dataframe which
#' should contain at least 3 columns entitled `calendar_age`, `c14_age` and `c14_sig`.
#' This format matches [carbondate::intcal20].
#' @param plot_14C_age Whether to use the radiocarbon age (\eqn{{}^{14}}C yr BP) as
#' the units of the y-axis in the plot. Defaults to `TRUE`. If `FALSE` uses
#' F\eqn{{}^{14}}C concentration instead.
#' @param plot_cal_age_scale (Optional) The calendar scale to use for the x-axis. Allowed values are
#' "BP", "AD" and "BC". The default is "BP" corresponding to plotting in cal yr BP.
#' @param show_individual_means (Optional) Whether to calculate and show the mean posterior
#' calendar age estimated for each individual \eqn{{}^{14}}C sample on the plot as a rug on
#' the x-axis. Default is `TRUE`.
#' @param show_confidence_intervals Whether to show the pointwise confidence intervals
#' (at chosen probability level) on the plot. Default is `TRUE`.
#' @param interval_width The confidence intervals to show for both the
#' calibration curve and the predictive density. Choose from one of `"1sigma"` (68.3%),
#' `"2sigma"` (95.4%) and `"bespoke"`. Default is `"2sigma"`.
#' @param bespoke_probability The probability to use for the confidence interval
#' if `"bespoke"` is chosen above. E.g., if 0.95 is chosen, then the 95% confidence
#' interval is calculated. Ignored if `"bespoke"` is not chosen.
#' @param denscale (Optional) Whether to scale the vertical range of the Poisson process mean rate plot
#' relative to the calibration curve plot. Default is 3 which means
#' that the maximum of the mean rate will be at 1/3 of the height of the plot.
#' @param resolution The distance between calendar ages at which to calculate the value of the rate
#' \eqn{\lambda(t)}. These ages will be created on a regular grid that automatically covers
#' the calendar period specified in `output_data`. Default is 1.
#' @param n_burn The number of MCMC iterations that should be discarded as burn-in (i.e.,
#' considered to be occurring before the MCMC has converged). This relates to the number
#' of iterations (`n_iter`) when running the original update functions (not the thinned `output_data`).
#' Any MCMC iterations before this are not used in the calculations. If not given, the first half of the
#' MCMC chain is discarded. Note: The maximum value that the function
#' will allow is `n_iter - 100 * n_thin` (where `n_iter` and `n_thin` are the arguments that were given to
#' [carbondate::PPcalibrate]) which would leave only 100 of the (thinned) values in `output_data`.
#' @param n_end The last iteration in the original MCMC chain to use in the calculations. Assumed to be the
#' total number of iterations performed, i.e. `n_iter`, if not given.
#' @param plot_pretty logical, defaulting to `TRUE`. If set `TRUE` then will select pretty plotting
#' margins (that create sufficient space for axis titles and rotates y-axis labels). If `FALSE` will
#' implement current user values.
#'
#'
#' @return A list, each item containing a data frame of the `calendar_age_BP`, the `rate_mean`
#' and the confidence intervals for the rate - `rate_ci_lower` and `rate_ci_upper`.
#'
#' @export
#'
#' @examples
#' # NOTE: All these examples are shown with a small n_iter and n_posterior_samples
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
#' # Default plot with 2 sigma interval
#' PlotPosteriorMeanRate(pp_output, n_posterior_samples = 100)
#'
#' # Specify an 80% confidence interval
#' PlotPosteriorMeanRate(
#'     pp_output,
#'     interval_width = "bespoke",
#'     bespoke_probability = 0.8,
#'     n_posterior_samples = 100)
#'
#' # Plot the posterior rate conditional on 2 internal changes
#' PlotPosteriorMeanRate(
#'     pp_output,
#'     n_changes = 2,
#'     interval_width = "bespoke",
#'     bespoke_probability = 0.8,
#'     n_posterior_samples = 100)
PlotPosteriorMeanRate <- function(
    output_data,
    n_posterior_samples = 5000,
    n_changes = NULL,
    calibration_curve = NULL,
    plot_14C_age = TRUE,
    plot_cal_age_scale = "BP",
    show_individual_means = TRUE,
    show_confidence_intervals = TRUE,
    interval_width = "2sigma",
    bespoke_probability = NA,
    denscale = 3,
    resolution = 1,
    n_burn = NA,
    n_end = NA,
    plot_pretty = TRUE) {

  arg_check <- .InitializeErrorList()
  .CheckOutputData(arg_check, output_data, "RJPP")
  .CheckInteger(arg_check, n_posterior_samples, lower = 10)
  .CheckSingleNChanges(arg_check, n_changes)
  .CheckCalibrationCurveFromOutput(arg_check, output_data, calibration_curve)
  .CheckFlag(arg_check, plot_14C_age)
  .CheckChoice(arg_check, plot_cal_age_scale, c("BP", "AD", "BC"))
  .CheckFlag(arg_check, show_individual_means)
  .CheckFlag(arg_check, show_confidence_intervals)
  .CheckIntervalWidth(arg_check, interval_width, bespoke_probability)
  .CheckNumber(arg_check, denscale, lower = 0)
  .CheckNumber(arg_check, resolution, lower = 0.01)
  .ReportErrors(arg_check)

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
  # Calculate means and rate
  posterior_rate <- FindPosteriorMeanRate(
    output_data,
    calendar_age_sequence,
    n_posterior_samples,
    n_changes,
    interval_width,
    bespoke_probability,
    n_burn,
    n_end)

  if (show_individual_means){
    n_iter <- output_data$input_parameters$n_iter
    n_thin <- output_data$input_parameters$n_thin
    n_burn <- .SetNBurn(n_burn, n_iter, n_thin)
    n_end <- .SetNEnd(n_end, n_iter, n_thin)
    calendar_age_means <- apply(output_data$calendar_ages[(n_burn + 1):n_end, ], 2, mean)
  }

  ylim_rate <- c(0, denscale * max(posterior_rate$rate_mean))

  ##############################################################################
  # Plot curves

  if(!is.null(n_changes)) {
    plot_title <- bquote(paste("Estimate of ",
                               lambda,
                               "(t) conditioned on ",
                               .(n_changes),
                               " internal changes"))
  } else {
    plot_title <- expression(paste("Estimate of Poisson process rate ", lambda, "(t)"))
  }

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

  if (show_individual_means) {
    calendar_age_means <- .ConvertCalendarAge(plot_cal_age_scale, calendar_age_means)
    graphics::rug(calendar_age_means, side = 1, quiet = TRUE)
  }
  .PlotRateEstimateOnCurrentPlot(plot_cal_age_scale, posterior_rate, output_colour, show_confidence_intervals)

  .AddLegendToRatePlot(
    output_data,
    show_confidence_intervals,
    interval_width,
    bespoke_probability,
    calibration_curve_colour,
    output_colour)

  invisible(posterior_rate)
}


.PlotRateEstimateOnCurrentPlot <- function(
    plot_cal_age_scale, posterior_rate, output_colour, show_confidence_intervals) {

  cal_age <- .ConvertCalendarAge(plot_cal_age_scale, posterior_rate$calendar_age_BP)

  graphics::lines(cal_age, posterior_rate$rate_mean, col = output_colour)
  if (show_confidence_intervals) {
    graphics::lines(cal_age, posterior_rate$rate_ci_lower, col = output_colour, lty = 2)
    graphics::lines(cal_age, posterior_rate$rate_ci_upper, col = output_colour, lty = 2)
  }
}


.AddLegendToRatePlot <- function(
    output_data,
    show_confidence_intervals,
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
    "Posterior mean rate")
  lty <- c(1, 2, 1)
  pch <- c(NA, NA, NA)
  col <- c(calibration_curve_colour, calibration_curve_colour, output_colour)

  if (show_confidence_intervals) {
    legend_labels <- c(legend_labels, ci_label)
    lty <- c(lty, 2)
    pch <- c(pch, NA)
    col <- c(col, output_colour)
  }

  graphics::legend(
    "topright", legend = legend_labels, lty = lty, pch = pch, col = col)
}

