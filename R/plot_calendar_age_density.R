#' Plots the predicted calendar age density from the output data
#'
#' Plots the input radiocarbon determinations and calibration curve, with the
#' output predicted density on the same plot. Can also optionally show the
#' SPD estimate
#'
#' @param c14_determinations A vector containing the radiocarbon determinations.
#' @param c14_uncertainties A vector containing the radiocarbon determination
#' uncertainties. Must be the same length as `c14_determinations`.
#' @param calibration_curve A dataframe which should contain one column entitled
#' c14_age and one column entitled c14_sig.
#' This format matches [carbondate::intcal20].
#' @param output_data The return value from one of the updating functions e.g.
#' [carbondate::WalkerBivarDirichlet] or
#' [carbondate::BivarGibbsDirichletwithSlice].
#' @param n_posterior_samples Current number of samples it will draw from this
#' posterior to estimate the calendar age density (possibly repeats).
#' @param show_SPD Whether to calculate and show the summed probability function
#' on the plot (optional). Default is `TRUE`.
#' @param show_confidence_intervals Whether to show the 95% confidence intervals
#' of the posterior density on the plot. Default is `TRUE`.
#' @param true_density The true calendar age density (optional). Expects a data
#' frame with two columns, the first column being calendar age and the right
#' column being the normalized density.
#' @param xlimscal,ylimscal Whether to scale the x or y limits (optional).
#' Default is 1. Values more than 1 will increase range of the limits,
#' values less than 1 will decrease the range of the limits.
#' @param denscale Whether to scale the vertical range of the density plot
#' relative to the calibration curve plot (optional). Default is 3 which means
#' that the maximum SPD density will be at 1/3 of the height of the plot.
#' TODO: base on the output data.
#'
#' @return A list, each item containing a data frame of the `calendar_age`, the
#' `density_mean` and the 95% confidence intervals for the density
#' `density_ci_lower` and `density_ci_upper` for each set of output data.
#'
#' @export
#'
#' @examples
#' # First generate some output data  for plotting using two update methods.
#' c14_determinations = c(602, 805, 954)
#' c14_uncertainties = c(35, 34, 45)
#' neal_temp = BivarGibbsDirichletwithSlice(
#'   c14_determinations = c14_determinations,
#'   c14_uncertainties = c14_uncertainties,
#'   calibration_curve = intcal20,
#'   lambda = 0.1,
#'   nu1 = 0.25,
#'   nu2 = 10,
#'   alpha_shape = 1,
#'   alpha_rate = 1)
#'
#' walker_temp = WalkerBivarDirichlet(
#'   c14_determinations = c14_determinations,
#'   c14_uncertainties = c14_uncertainties,
#'   calibration_curve = intcal20,
#'   lambda = 0.1,
#'   nu1 = 0.25,
#'   nu2 = 10,
#'   alpha_shape = 1,
#'   alpha_rate = 1)
#'
#' # Plot results for both on the same graph
#' PlotCalendarAgeDensity(
#'   c14_determinations = c14_determinations,
#'   c14_uncertainties = c14_uncertainties,
#'   calibration_curve = intcal20,
#'   output_data = list(neal_temp, walker_temp),
#'   n_posterior_samples = 20)
#'
PlotCalendarAgeDensity <- function(
    c14_determinations,
    c14_uncertainties,
    calibration_curve,
    output_data,
    n_posterior_samples,
    show_SPD = TRUE,
    show_confidence_intervals = TRUE,
    true_density = NA,
    xlimscal = 1,
    ylimscal = 1,
    denscale = 3) {

  ##############################################################################
  # Check input parameters

  # Treat single output data as a list of length 1
  if (!is.null(output_data$update_type)) output_data = list(output_data)

  num_data = length(output_data)
  for (i in 1:num_data) {
    if (is.null(output_data[[i]]$label)) {
      output_data[[i]]$label <- stringr::str_to_title(
        output_data[[i]]$update_type)
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

  calendar_age_sequence <- .CreateRangeToPlotDensity(output_data[[1]])

  ##############################################################################
  # Calculate density distributions

  if (show_SPD){
    SPD = FindSPD(
      calendar_age_range = floor(range(output_data[[1]]$calendar_ages)),
      c14_determinations = c14_determinations,
      c14_uncertainties = c14_uncertainties,
      calibration_curve = calibration_curve)
  }

  posterior_density <- list()
  for (i in 1:num_data) {
    posterior_density[[i]] <- .FindPosteriorDensityMeanAndCI(
      output_data[[i]], calendar_age_sequence, n_posterior_samples)
  }
  ##############################################################################
  # Calculate plot scaling

  xlim <- .ScaleLimit(rev(range(calendar_age_sequence)), xlimscal)
  ylim_calibration <- .ScaleLimit(
    range(c14_determinations) +
      c(-2, 2) * stats::quantile(c14_uncertainties, 0.9),
    ylimscal)
  ylim_density = c(0, denscale * max(posterior_density[[1]]$density_mean))

  ##############################################################################
  # Plot curves

  .PlotCalibrationCurveAndInputData(
    xlim,
    ylim_calibration,
    calibration_curve,
    c14_determinations,
    calibration_curve_colour,
    calibration_curve_bg)

  .SetUpDensityPlot(xlim, ylim_density)

  if (show_SPD){
    .PlotSPDEstimateOnCurrentPlot(SPD, SPD_colour, xlim, ylim_density)
  }

  for (i in 1:num_data) {
    .PlotDensityEstimateOnCurrentPlot(
      posterior_density[[i]], output_colours[[i]], show_confidence_intervals)
  }

  if (is.data.frame(true_density)) {
    .PlotTrueDensityOnCurrentPlot(true_density, true_density_colour)
  }

  .AddLegendToDensityPlot(
    output_data,
    show_SPD,
    is.data.frame(true_density),
    show_confidence_intervals,
    calibration_curve_colour,
    output_colours,
    SPD_colour,
    true_density_colour)

  invisible(posterior_density)
}


.CreateRangeToPlotDensity <- function(output_data) {
  calendar_age_sequence <- seq(
    floor(min(output_data$calendar_ages, na.rm = TRUE)),
    ceiling(max(output_data$calendar_ages, na.rm = TRUE)),
    by = 1,
  )
  return(calendar_age_sequence)
}


.ScaleLimit <- function(lim, limscal) {
  lim <- lim + c(1, -1) * diff(lim) * (1 - limscal)
  return(lim)
}


.PlotCalibrationCurveAndInputData <- function(
    xlim,
    ylim,
    calibration_curve,
    c14_determinations,
    calibration_curve_colour,
    calibration_curve_bg){
  graphics::par(mar = c(5, 4.5, 4, 2) + 0.1, las = 1)
  graphics::plot.default(
    calibration_curve$calendar_age,
    calibration_curve$c14_age,
    col = calibration_curve_colour,
    ylim = ylim,
    xlim = xlim,
    xlab = "Calendar Age (cal yr BP)",
    ylab = expression(paste(""^14, "C", " age (", ""^14, "C yr BP)")),
    type = "l",
    main = expression(paste("Posterior calendar age density")),
  )
  calibration_curve$ub <- calibration_curve$c14_age +
    1.96 * calibration_curve$c14_sig
  calibration_curve$lb <- calibration_curve$c14_age -
    1.96 * calibration_curve$c14_sig

  graphics::lines(
    calibration_curve$calendar_age,
    calibration_curve$ub,
    lty = 2,
    col = calibration_curve_colour)
  graphics::lines(
    calibration_curve$calendar_age,
    calibration_curve$lb,
    lty = 2,
    col = calibration_curve_colour)
  graphics::polygon(
    c(rev(calibration_curve$calendar_age), calibration_curve$calendar_age),
    c(rev(calibration_curve$lb), calibration_curve$ub),
    col = calibration_curve_bg,
    border = NA)
  graphics::rug(c14_determinations, side = 2)
}


.SetUpDensityPlot <- function(xlim, ylim) {
  graphics::par(new = TRUE)
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


.PlotSPDEstimateOnCurrentPlot <- function(SPD, SPD_colour, xlim, ylim) {
  graphics::polygon(
    c(SPD$calendar_age, rev(SPD$calendar_age)),
    c(SPD$probability, rep(0, length(SPD$probability))),
    border = NA,
    col = SPD_colour)
}


.PlotDensityEstimateOnCurrentPlot <- function(
    posterior_density, output_colour, show_confidence_intervals) {

  graphics::lines(
    posterior_density$calendar_age,
    posterior_density$density_mean,
    col = output_colour)
  if (show_confidence_intervals) {
    graphics::lines(
      posterior_density$calendar_age,
      posterior_density$density_ci_lower,
      col = output_colour,
      lty = 2)
    graphics::lines(
      posterior_density$calendar_age,
      posterior_density$density_ci_upper,
      col = output_colour,
      lty = 2)
  }
}


.PlotTrueDensityOnCurrentPlot <- function(true_density, true_density_colour) {
  graphics::lines(
    true_density[[1]], true_density[[2]], col = true_density_colour)
}


.AddLegendToDensityPlot <- function(
    output_data,
    show_SPD,
    show_true_density,
    show_confidence_intervals,
    calibration_curve_colour,
    output_colours,
    SPD_colour,
    true_density_colour) {
  legend_labels = "IntCal20"
  lty = 1
  pch = NA
  col = calibration_curve_colour

  for (i in 1:length(output_data)) {
    legend_labels <- c(legend_labels, output_data[[i]]$label)
    lty <- c(lty, 1)
    pch <- c(pch, NA)
    col <- c(col, output_colours[[i]])

    if (show_confidence_intervals) {
      legend_labels <- c(
        legend_labels, paste(output_data[[i]]$label, " 95% prob interval"))
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
