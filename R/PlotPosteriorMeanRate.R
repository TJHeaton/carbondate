#' Plots the posterior mean rate from the output data
#'
#' Plots the input radiocarbon determinations and calibration curve, with the
#' output posterior mean rate on the same plot. Can also optionally show the
#' individual mean calendar ages.
#'
#' @param output_data The return value from the updating functions
#' [carbondate::PPcalibrate]. Optionally, the output data can have an extra list item
#' named `label` which is used to set the label on the plot legend.
#' @param n_posterior_samples Current number of samples it will draw from this
#' posterior to estimate the calendar age density (possibly repeats). If not
#' given 5000 is used.
#' @param calibration_curve This is usually not required since the name of the
#' calibration curve variable is saved in the output data. However if the
#' variable with this name is no longer in your environment then you should pass
#' the calibration curve here. If provided this should be a dataframe which
#' should contain at least 3 columns entitled calendar_age, c14_age and c14_sig.
#' This format matches [carbondate::intcal20].
#' @param plot_14C_age Whether to use the 14C yr BP as the units of the y-axis.
#' Defaults to TRUE. If FALSE uses F14C concentration instead.
#' @param show_individual_means Whether to calculate and show the individual mean
#' calendar ages on the plot (optional). Default is `FALSE`.
#' @param show_confidence_intervals Whether to show the confidence intervals
#' for the chosen probability on the plot. Default is `TRUE`.
#' @param interval_width The confidence intervals to show for both the
#' calibration curve and the predictive density. Choose from one of `"1sigma"` (68.3%),
#' `"2sigma"` (95.4%) and `"bespoke"`. Default is `"2sigma"`.
#' @param bespoke_probability The probability to use for the confidence interval
#' if `"bespoke"` is chosen above. E.g. if 0.95 is chosen, then the 95% confidence
#' interval is calculated. Ignored if `"bespoke"` is not chosen.
#' @param denscale Whether to scale the vertical range of the density plot
#' relative to the calibration curve plot (optional). Default is 3 which means
#' that the maximum SPD density will be at 1/3 of the height of the plot.
#' @param n_calc Number of points to use when calculating the predictive
#' density. Default is 1001.
#' @param n_burn The number of samples required for burn-in - any samples before this
#' are not used in the calculation. If not given, the first half of the
#' MCMC chain is discarded. Note that the maximum
#' value that can be chosen is `n_iter - 100 * n_thin` (where `n_iter` and `n_thin` are the
#' arguments given to [carbondate::PPcalibrate]).
#' @param n_end The iteration number of the last sample to use. Assumed to be the number of iterations
#' if not given.
#'
#'
#' @return A list, each item containing a data frame of the `calendar_age`, the `rate_mean`
#' and the confidence intervals for the rate - `rate_ci_lower` and `rate_ci_upper`.
#'
#' @export
#'
#' @examples
#' pp_output <- PPcalibrate(
#'     pp_uniform_phase$c14_age, pp_uniform_phase$c14_sig, intcal20, n_iter = 5000, show_progress = FALSE)
#'
#' # Default plot with 2 sigma interval
#' PlotPosteriorMeanRate(pp_output)
#'
#' # Specify an 80% confidence interval
#' PlotPosteriorMeanRate(
#'     pp_output, interval_width = "bespoke", bespoke_probability = 0.8)
PlotPosteriorMeanRate <- function(
    output_data,
    n_posterior_samples = 5000,
    calibration_curve = NULL,
    plot_14C_age = TRUE,
    show_individual_means = TRUE,
    show_confidence_intervals = TRUE,
    interval_width = "2sigma",
    bespoke_probability = NA,
    denscale = 3,
    n_calc = 1001,
    n_burn = NA,
    n_end = NA) {

  n_iter <- output_data$input_parameters$n_iter
  n_thin <- output_data$input_parameters$n_thin

  arg_check <- .InitializeErrorList()
  .CheckOutputData(arg_check, output_data, "RJPP")
  .CheckInteger(arg_check, n_posterior_samples, lower = 10)
  .CheckCalibrationCurveFromOutput(arg_check, output_data, calibration_curve)
  .CheckFlag(arg_check, plot_14C_age)
  .CheckFlag(arg_check, show_individual_means)
  .CheckFlag(arg_check, show_confidence_intervals)
  .CheckIntervalWidth(arg_check, interval_width, bespoke_probability)
  .CheckNumber(arg_check, denscale, lower = 0)
  .CheckInteger(arg_check, n_calc, lower = 20)
  .CheckNBurnAndNEnd(arg_check, n_burn, n_end, n_iter, n_thin)
  .ReportErrors(arg_check)

  n_burn <- .SetNBurn(n_burn, n_iter, n_thin)
  n_end <- .SetNEnd(n_end, n_iter, n_thin)

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

  calendar_age_sequence <- seq(
    from = min(output_data$rate_s[[1]]), to = max(output_data$rate_s[[1]]), length.out = n_calc)
  xlim <- rev(range(calendar_age_sequence))

  calendar_age_sequence[n_calc] <- calendar_age_sequence[n_calc] - 1e-6 * diff(calendar_age_sequence)[1]
  plot_AD <- any(calendar_age_sequence < 0)
  ##############################################################################
  # Calculate means and rate

  if (show_individual_means){
    calendar_age_means <- apply(output_data$calendar_ages[(n_burn + 1):n_end, ], 2, mean)
  }

  indices <- sample((n_burn + 1):n_end, n_posterior_samples, replace = ((n_end - n_burn) < n_posterior_samples))
  rate <- matrix(NA, nrow = n_posterior_samples, ncol = length(calendar_age_sequence))
  for (i in 1:n_posterior_samples) {
    ind <- indices[i]
    rate[i,] <- stats::approx(
      x = output_data$rate_s[[ind]],
      y = c(output_data$rate_h[[ind]], 0),
      xout = calendar_age_sequence,
      method = "constant")$y
  }
  edge_width <- switch(
    interval_width,
    "1sigma" = 1 - stats::pnorm(1),
    "2sigma"  = 1 - stats::pnorm(2),
    "bespoke" = (1 - bespoke_probability)/2
  )
  posterior_rate <- data.frame(
    calendar_age = calendar_age_sequence,
    rate_mean = apply(rate, 2, mean),
    rate_ci_lower = apply(rate, 2, stats::quantile, probs = edge_width),
    rate_ci_upper = apply(rate, 2, stats::quantile, probs = 1 - edge_width)
  )

  ylim_rate <- c(0, denscale * max(posterior_rate$rate_mean))

  ##############################################################################
  # Plot curves

  .PlotCalibrationCurveAndInputData(
    plot_AD,
    xlim,
    calibration_curve,
    rc_determinations,
    plot_14C_age,
    calibration_curve_colour,
    calibration_curve_bg,
    interval_width,
    bespoke_probability,
    title = expression(paste("Estimate of Poisson process rate ", lambda, "(t)")))

  .SetUpDensityPlot(plot_AD, xlim, ylim_rate)

  if (show_individual_means) {
    if (plot_AD) {
      calendar_age_means <- 1950 - calendar_age_means
    }
    graphics::rug(calendar_age_means, side = 1, quiet = TRUE)
  }
  .PlotRateEstimateOnCurrentPlot(plot_AD, posterior_rate, output_colour, show_confidence_intervals)

  .AddLegendToRatePlot(
    output_data,
    show_confidence_intervals,
    interval_width,
    bespoke_probability,
    calibration_curve_colour,
    output_colour)

  invisible(posterior_rate)
}


.PlotRateEstimateOnCurrentPlot <- function(plot_AD, posterior_rate, output_colour, show_confidence_intervals) {

  if (plot_AD) {
    cal_age <- 1950 - posterior_rate$calendar_age
  } else {
    cal_age <- posterior_rate$calendar_age
  }

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
    "1sigma" = expression(paste(sigma, " interval")),
    "2sigma"  = expression(paste("2", sigma, " interval")),
    "bespoke" = paste0(round(100 * bespoke_probability), "% interval"))

  legend_labels <- c(
    gsub("intcal", "IntCal", output_data$input_data$calibration_curve_name),
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

