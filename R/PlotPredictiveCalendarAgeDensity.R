#' Plots the calendar age density of all objects from the output data
#'
#' Plots the input radiocarbon determinations and calibration curve, with the
#' output predicted density on the same plot. Can also optionally show the
#' SPD estimate. Note that if all you are only interested in is the density
#' data, without an accompanying plot, you can use
#' [carbondate::FindPredictiveCalendarAgeDensity] instead.
#'
#' @param output_data The return value from one of the updating functions e.g.
#' [carbondate::WalkerBivarDirichlet] or
#' [carbondate::PolyaUrnBivarDirichlet] or a list, each item containing
#' one of these values. Optionally, the output data can have an extra list item
#' named `label` which is used to set the label on the plot legend.
#' @param n_posterior_samples Current number of samples it will draw from this
#' posterior to estimate the calendar age density (possibly repeats).
#' @param calibration_curve This is usually not required since the name of the
#' calibration curve variable is saved in the output data. However if the
#' variable with this name is no longer in your environment then you should pass
#' the calibration curve here. If provided this should be a dataframe which
#' should contain at least 3 columns entitled calendar_age, c14_age and c14_sig.
#' This format matches [carbondate::intcal20].
#' @param plot_14C_age Whether to use the 14C yr BP as the units of the y-axis.
#' Defaults to TRUE. If FALSE uses F14C concentration instead.
#' @param show_SPD Whether to calculate and show the summed probability
#' distribution on the plot (optional). Default is `FALSE`.
#' @param show_confidence_intervals Whether to show the confidence intervals
#' for the chosen probability on the plot. Default is `TRUE`.
#' @param interval_width The confidence intervals to show for both the
#' calibration curve and the predictive density. Choose from one of `"1sigma"`,
#' `"2sigma"` and `"bespoke"`. Default is `"2sigma"`.
#' @param bespoke_probability The probability to use for the confidence interval
#' if `"bespoke"` is chosen above. E.g. if 0.95 is chosen, then the 95% confidence
#' interval is calculated. Ignored if `"bespoke"` is not chosen.
#' @param true_density The true calendar age density (optional). Expects a data
#' frame with two columns, the first column being calendar age and the right
#' column being the normalized density.
#' @param denscale Whether to scale the vertical range of the density plot
#' relative to the calibration curve plot (optional). Default is 3 which means
#' that the maximum SPD density will be at 1/3 of the height of the plot.
#' @param n_calc Number of points to use when calculating the predictive
#' density. Default is 1001.
#' @param n_burn The number of samples required for burn-in - any samples before this
#' are not used to calculate the density. If not given, the first half of the
#' MCMC chain is discarded. Looking at the output of [carbondate::PlotConvergenceData] can
#' help determine what an appropriate value is for a given output. Note that the maximum
#' value that can be chosen is `n_iter - 100 * n_thin` (where `n_iter` and `n_thin` are the
#' arguments given to [carbondate::WalkerBivarDirichlet] or c[carbondate::PolyaUrnBivarDirichlet]).
#'
#' @return A list, each item containing a data frame of the `calendar_age`, the
#' `density_mean` and the 95% confidence intervals for the density
#' `density_ci_lower` and `density_ci_upper` for each set of output data.
#'
#' @export
#'
#' @examples
#' # Plot results for a single calibration
#' PlotPredictiveCalendarAgeDensity(walker_example_output, 500)
#'
#' # Plot results from a calibration, and add a label
#' new_output = walker_example_output
#' new_output$label = "My plot"
#' PlotPredictiveCalendarAgeDensity(new_output, 500)
#'
#' # Plot results from two calibrations on the same plot, and show the SPD
#' PlotPredictiveCalendarAgeDensity(
#'   list(walker_example_output, polya_urn_example_output), 500, show_SPD = TRUE)
#'
#' # Plot and show the 1-sigma confidence interval
#' PlotPredictiveCalendarAgeDensity(walker_example_output, 500, interval_width = "1sigma")
#'
#' # Plot and show the 80% confidence interval
#' PlotPredictiveCalendarAgeDensity(
#'   walker_example_output, 500, interval_width = "bespoke", bespoke_probability = 0.8)
PlotPredictiveCalendarAgeDensity <- function(
    output_data,
    n_posterior_samples,
    calibration_curve = NULL,
    plot_14C_age = TRUE,
    show_SPD = FALSE,
    show_confidence_intervals = TRUE,
    interval_width = "2sigma",
    bespoke_probability = NA,
    true_density = NULL,
    denscale = 3,
    n_calc = 1001,
    n_burn = NA) {

  ##############################################################################
  # Check input parameters

  arg_check <- checkmate::makeAssertCollection()

  # Treat single output data as a list of length 1
  if (!is.null(output_data$update_type)) output_data = list(output_data)

  .CheckMultipleOutputDataConsistent(output_data)
  num_data = length(output_data)
  for (i in 1:num_data) {
    .CheckOutputData(arg_check, output_data[[i]])
    .CheckCalibrationCurveFromOutput(
      arg_check, output_data[[i]], calibration_curve)
    if (is.null(output_data[[i]]$label)) {
      output_data[[i]]$label <- output_data[[i]]$update_type
    }
  }

  checkmate::assertInt(n_posterior_samples, lower = 10, add = arg_check)
  .CheckIntervalWidth(arg_check, interval_width, bespoke_probability)
  checkmate::assertNumber(denscale, lower = 0, add = arg_check)
  checkmate::assertDataFrame(
    true_density, types = "numeric", null.ok = TRUE, add = arg_check)
  checkmate::reportAssertions(arg_check)

  if (is.null(calibration_curve)) {
    calibration_curve = get(output_data[[1]]$input_data$calibration_curve_name)
  }
  rc_determinations <- output_data[[1]]$input_data$rc_determinations
  rc_sigmas <- output_data[[1]]$input_data$rc_sigmas
  F14C_inputs <-output_data[[1]]$input_data$F14C_inputs

  if (plot_14C_age == TRUE) {
    calibration_curve = .AddC14ageColumns(calibration_curve)
    if (F14C_inputs == TRUE) {
      converted <- .ConvertF14cTo14Cage(rc_determinations, rc_sigmas)
      rc_determinations <- converted$f14c
      rc_sigmas <- converted$f14c_sig
    }
  } else {
    calibration_curve = .AddF14cColumns(calibration_curve)
    if (F14C_inputs == FALSE) {
      converted <- .Convert14CageToF14c(rc_determinations, rc_sigmas)
      rc_determinations <- converted$c14_age
      rc_sigmas <- converted$c14_sig
    }
  }

  ##############################################################################
  # Initialise plotting parameters

  SPD_colour <- grDevices::grey(0.1, alpha = 0.3)
  calibration_curve_colour <- "blue"
  calibration_curve_bg <- grDevices::rgb(0, 0, 1, .3)
  output_colours <- c("purple", "darkgreen", "darkorange2", "deeppink3")
  if (num_data > 4) {
    output_colours <- c(output_colours, grDevices::hcl.colors(n=num_data-4))
  }
  true_density_colour = "red"

  calendar_age_sequence <- .CreateXRangeToPlotDensity(output_data[[1]], n_calc)
  plot_AD = any(calendar_age_sequence < 0)
  ##############################################################################
  # Calculate density distributions

  if (show_SPD){
    SPD = FindSummedProbabilityDistribution(
      calendar_age_range_BP = floor(range(calendar_age_sequence)),
      rc_determinations = rc_determinations,
      rc_sigmas = rc_sigmas,
      F14C_inputs = !plot_14C_age,
      calibration_curve = calibration_curve)
  }

  predictive_density <- list()
  for (i in 1:num_data) {
    predictive_density[[i]] <- FindPredictiveCalendarAgeDensity(
      output_data[[i]],
      calendar_age_sequence,
      n_posterior_samples,
      interval_width,
      bespoke_probability,
      n_burn)
  }
  ##############################################################################
  # Calculate plot scaling
  xlim <- rev(range(calendar_age_sequence))
  ylim_density = c(0, denscale * max(predictive_density[[1]]$density_mean))

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
    bespoke_probability)

  .SetUpDensityPlot(plot_AD, xlim, ylim_density)

  if (show_SPD){
    .PlotSPDEstimateOnCurrentPlot(plot_AD, SPD, SPD_colour, xlim, ylim_density)
  }

  for (i in 1:num_data) {
    .PlotDensityEstimateOnCurrentPlot(
      plot_AD, predictive_density[[i]], output_colours[[i]], show_confidence_intervals)
  }

  if (is.data.frame(true_density)) {
    .PlotTrueDensityOnCurrentPlot(plot_AD, true_density, true_density_colour)
  }

  .AddLegendToDensityPlot(
    output_data,
    show_SPD,
    is.data.frame(true_density),
    show_confidence_intervals,
    interval_width,
    bespoke_probability,
    calibration_curve_colour,
    output_colours,
    SPD_colour,
    true_density_colour)

  invisible(predictive_density)
}


.CreateXRangeToPlotDensity <- function(output_data, n_calc) {
  start_age = min(output_data$calendar_ages, output_data$density_data$calendar_ages[1])
  end_age = max(output_data$calendar_ages, output_data$density_data$calendar_ages[2])
  calendar_age_sequence <- seq(start_age, end_age, length = n_calc)
  return(calendar_age_sequence)
}

.PlotCalibrationCurveAndInputData <- function(
    plot_AD,
    xlim,
    calibration_curve,
    rc_determinations,
    plot_14C_age,
    calibration_curve_colour,
    calibration_curve_bg,
    interval_width,
    bespoke_probability){
  graphics::par(mar = c(5, 4.5, 4, 2) + 0.1, las = 1)
  .PlotCalibrationCurve(
    plot_AD,
    xlim,
    plot_14C_age,
    calibration_curve,
    calibration_curve_colour,
    calibration_curve_bg,
    interval_width,
    bespoke_probability,
    title = "Calendar age density estimate")
  graphics::rug(rc_determinations, side = 2)
}


.PlotCalibrationCurve = function(
    plot_AD,
    xlim,
    plot_14C_age,
    calibration_curve,
    calibration_curve_colour,
    calibration_curve_bg,
    interval_width,
    bespoke_probability,
    title) {

  cal_age_ind_min <- which.min(abs(calibration_curve$calendar_age_BP - min(xlim)))
  cal_age_ind_max <- which.min(abs(calibration_curve$calendar_age_BP - max(xlim)))
  calendar_age_indices <- cal_age_ind_min:cal_age_ind_max

  if (plot_AD) {
    cal_age = 1950 - calibration_curve$calendar_age
    xlim = 1950 - xlim
    x_label = "Calendar Age (AD)"
  } else {
    cal_age = calibration_curve$calendar_age
    x_label = "Calendar Age (cal yr BP)"
  }

  if (plot_14C_age) {
    rc_age = calibration_curve$c14_age
    rc_sig = calibration_curve$c14_sig
    y_label = expression(paste(""^14, "C", " age (yr BP)"))
  } else {
    rc_age = calibration_curve$f14c
    rc_sig = calibration_curve$f14c_sig
    y_label = expression(paste("F ", ""^14, "C"))
  }

  # multiplier for the confidence interval if you have a standard deviation
  zquant <- switch(
    interval_width,
    "1sigma" = 1,
    "2sigma" = 2,
    "bespoke" = - stats::qnorm((1 - bespoke_probability) / 2))
  calibration_curve$ub <- rc_age + zquant * rc_sig
  calibration_curve$lb <- rc_age - zquant * rc_sig

  # calculate the limits for the y axis
  ylim = c(
    min(calibration_curve$lb[calendar_age_indices]),
    max(calibration_curve$ub[calendar_age_indices])
    )
  ylim = ylim + 0.25 * c(-1, 1) * diff(ylim)

  graphics::plot.default(
    cal_age,
    rc_age,
    col = calibration_curve_colour,
    ylim = ylim,
    xlim = xlim,
    xlab = x_label,
    ylab = y_label,
    type = "l",
    main = title)
  graphics::lines(cal_age, calibration_curve$ub, lty = 2, col = calibration_curve_colour)
  graphics::lines(cal_age, calibration_curve$lb, lty = 2, col = calibration_curve_colour)
  graphics::polygon(
    c(rev(cal_age), cal_age),
    c(rev(calibration_curve$lb), calibration_curve$ub),
    col = calibration_curve_bg,
    border = NA)
}


.SetUpDensityPlot <- function(plot_AD, xlim, ylim) {
  graphics::par(new = TRUE)

  if (plot_AD) xlim = 1950 - xlim
  graphics::plot.default(
    c(),
    c(),
    type = "n",
    ylim = ylim,
    xlim = xlim,
    axes = FALSE,
    xlab = NA,
    ylab = NA)
}


.PlotSPDEstimateOnCurrentPlot <- function(plot_AD, SPD, SPD_colour, xlim, ylim) {
  if (plot_AD) {
    cal_age = 1950 - SPD$calendar_age
  } else {
    cal_age = SPD$calendar_age
  }

  graphics::polygon(
    c(cal_age, rev(cal_age)),
    c(SPD$probability, rep(0, length(SPD$probability))),
    border = NA,
    col = SPD_colour)
}


.PlotDensityEstimateOnCurrentPlot <- function(
    plot_AD, predictive_density, output_colour, show_confidence_intervals) {

  if (plot_AD) {
    cal_age = 1950 - predictive_density$calendar_age
  } else {
    cal_age = predictive_density$calendar_age
  }

  graphics::lines(cal_age, predictive_density$density_mean, col = output_colour)
  if (show_confidence_intervals) {
    graphics::lines(cal_age, predictive_density$density_ci_lower, col = output_colour, lty = 2)
    graphics::lines(cal_age, predictive_density$density_ci_upper, col = output_colour, lty = 2)
  }
}


.PlotTrueDensityOnCurrentPlot <- function(plot_AD, true_density, true_density_colour) {
  if (plot_AD) {
    cal_age = 1950 - true_density[[1]]
  } else {
    cal_age = true_density[[1]]
  }

  graphics::lines(cal_age, true_density[[2]], col = true_density_colour)
}


.AddLegendToDensityPlot <- function(
    output_data,
    show_SPD,
    show_true_density,
    show_confidence_intervals,
    interval_width,
    bespoke_probability,
    calibration_curve_colour,
    output_colours,
    SPD_colour,
    true_density_colour) {

  ci_label = switch(
    interval_width,
    "1sigma" = "1 sigma interval",
    "2sigma"  = "2 sigma interval",
    "bespoke" = paste(round(100*bespoke_probability), "% interval", sep = ""))

  legend_labels = c(output_data[[1]]$input_data$calibration_curve_name, ci_label)
  lty = c(1, 2)
  pch = c(NA, NA)
  col = c(calibration_curve_colour, calibration_curve_colour)

  for (i in 1:length(output_data)) {
    legend_labels <- c(legend_labels, output_data[[i]]$label)
    lty <- c(lty, 1)
    pch <- c(pch, NA)
    col <- c(col, output_colours[[i]])

    if (show_confidence_intervals) {
      legend_labels <- c(legend_labels, ci_label)
      lty <- c(lty, 2)
      pch <- c(pch, NA)
      col <- c(col, output_colours[[i]])
    }
  }

  if (show_true_density) {
    legend_labels <- c(legend_labels, "True density")
    lty <- c(lty, 1)
    pch <- c(pch, NA)
    col <- c(col, true_density_colour)
  }

  if (show_SPD) {
    legend_labels <- c(legend_labels, "SPD Estimate")
    lty <- c(lty, -1)
    pch <- c(pch, 15)
    col <- c(col, SPD_colour)
  }

  graphics::legend(
    "topright", legend = legend_labels, lty = lty, pch = pch, col = col)
}
