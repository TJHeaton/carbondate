#' Plots the posterior calendar ages for an individual determination
#'
#' Once a function has been run to calibrate a set of radiocarbon
#' determinations, the posterior density histogram for a single determination can be
#' plotted using this function. This also plots a kernel density estimation from the histogram data
#' and shows the highest posterior density range for the interval width specified (default 2\eqn{\sigma}).
#'
#' @inheritParams PlotPredictiveCalendarAgeDensity
#' @param ident the determination you want to show the individual posterior
#' calendar age for.
#' @param resolution The distance between histogram breaks for the calendar age density.
#' Must be an integer greater than one.
#' @param interval_width The confidence intervals to show for the
#' calibration curve and for the highest posterior density ranges.
#' Choose from one of "1sigma" (68.3%), "2sigma" (95.4%) and "bespoke". Default is "2sigma".
#' @param bespoke_probability The probability to use for the confidence interval
#' if "bespoke" is chosen above. E.g. if 0.95 is chosen, then the 95% confidence
#' interval is calculated. Ignored if "bespoke" is not chosen.
#' @param show_hpd_ranges Set to `TRUE` to also show the highest posterior range on the plot.
#' These are calculated using [HDInterval::hdi]. Default is `FALSE`.
#' @param show_unmodelled_density Set to `TRUE` to also show the unmodelled density (i.e. the
#' result of [carbondate::CalibrateSingleDetermination]) on the plot. Default is `FALSE`.
#'
#' @export
#'
#' @examples
#' # Plot results for the 10th determination
#' PlotCalendarAgeDensityIndividualSample(10, polya_urn_example_output)
#'
#' # Plot results for the 10th determination and show the unmodelled density and HPD ranges
#' PlotCalendarAgeDensityIndividualSample(
#'     10, polya_urn_example_output, show_hpd_ranges = TRUE, show_unmodelled_density = TRUE)
#'
#' # Now change to showing 1 sigma interval for HPD range and calibration curve
#' PlotCalendarAgeDensityIndividualSample(
#'     10, polya_urn_example_output, interval_width = "1sigma",
#'     show_hpd_ranges = TRUE, show_unmodelled_density = TRUE)
PlotCalendarAgeDensityIndividualSample <- function(
    ident,
    output_data,
    calibration_curve = NULL,
    plot_14C_age = TRUE,
    resolution = 5,
    interval_width = "2sigma",
    bespoke_probability = NA,
    n_burn = NA,
    n_end = NA,
    show_hpd_ranges = FALSE,
    show_unmodelled_density = FALSE) {

  arg_check <- .InitializeErrorList()
  .CheckInteger(arg_check, ident)
  .CheckOutputData(arg_check, output_data,  c("Polya Urn", "Walker", "RJPP"))
  n_iter <- output_data$input_parameters$n_iter
  n_thin <- output_data$input_parameters$n_thin

  .CheckCalibrationCurveFromOutput(arg_check, output_data, calibration_curve)
  .CheckInteger(arg_check, resolution, lower = 1)
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
      rc_sigmas <- converted$c14_sig
    }
  } else {
    calibration_curve <- .AddF14cColumns(calibration_curve)
    if (F14C_inputs == FALSE) {
      converted <- .Convert14CageToF14c(rc_determinations, rc_sigmas)
      rc_determinations <- converted$f14c
      rc_sigmas <- converted$f14c_sig
    }
  }

  calendar_age <- output_data$calendar_ages[(n_burn+1):n_end, ident]
  rc_age <- rc_determinations[ident]
  rc_sig <- rc_sigmas[ident]

  smoothed_density <- stats::density(calendar_age, bw="SJ")

  # Find the calendar age range to plot
  xrange <- range(calendar_age)
  xrange <- xrange + 0.1 * c(-1, 1) * diff(xrange)
  xrange[1] <- floor(xrange[1])
  if (resolution > 1) while (xrange[1] %% resolution != 0) xrange[1] <- xrange[1] - 1

  xrange[2] <- ceiling(xrange[2])
  if (resolution > 1) while (xrange[2] %% resolution != 0) xrange[2] <- xrange[2] + 1

  if (plot_14C_age == FALSE) {
    title <- substitute(
      paste(
        "Posterior of ",
        i^th,
        " determination ",
        f14c_age,
        "\u00B1",
        f14c_sig,
        " F ",
        ""^14,
        "C"),
      list(i = ident, f14c_age = signif(rc_age, 2), f14c_sig = signif(rc_sig, 2)))
  } else {
    title <- substitute(
      paste(
        "Posterior of ",
        i^th,
        " determination ",
        c14_age,
        "\u00B1",
        c14_sig,
        ""^14,
        "C yr BP"),
      list(i = ident, c14_age = rc_age, c14_sig = rc_sig))
  }

  plot_AD <- any(calendar_age < 0)
  graphics::par(xaxs = "i", yaxs = "i")
  .PlotCalibrationCurve(
    plot_AD,
    xlim = rev(xrange),
    plot_14C_age = plot_14C_age,
    calibration_curve = calibration_curve,
    calibration_curve_colour = "blue",
    calibration_curve_bg = grDevices::rgb(0, 0, 1, .3),
    interval_width = interval_width,
    bespoke_probability = bespoke_probability,
    title = title)

  if (plot_AD) {
    calendar_age <- 1950 - calendar_age
    xrange <- 1950 - xrange
  }

  # Plot the 14C determination on the y-axis
  yfromto <- seq(rc_age - 4 * rc_sig, rc_age + 4 * rc_sig, length.out = 100)
  radpol <- cbind(
    c(0, stats::dnorm(yfromto, mean =rc_age, sd = rc_sig), 0),
    c(min(yfromto), yfromto, max(yfromto))
  )
  relative_height <- 0.1
  radpol[, 1] <- radpol[, 1] * (xrange[2] - xrange[1]) / max(radpol[, 1])
  radpol[, 1] <- radpol[, 1] * relative_height
  radpol[, 1] <- xrange[2] - radpol[, 1]
  graphics::polygon(radpol, col = grDevices::rgb(1, 0, 0, .5))

  # Plot the posterior cal age on the x-axis
  graphics::par(new = TRUE, las = 1)
  # Create hist but do not plot - works out sensible ylim
  if (plot_AD) {
    breaks <-seq(xrange[2], xrange[1], by=resolution)
  } else {
    breaks <-seq(xrange[1], xrange[2], by=resolution)
  }
  temphist <- graphics::hist(calendar_age, breaks = breaks, plot = FALSE)

  graphics::hist(
    calendar_age,
    prob = TRUE,
    breaks = breaks,
    xlim = rev(xrange),
    axes = FALSE,
    xlab = NA,
    ylab = NA,
    main = "",
    xaxs = "i",
    ylim = c(0, 3 * max(temphist$density)),
    col = "lavender",
    border = "purple")
  if (show_unmodelled_density) {
    unmodelled <- CalibrateSingleDetermination(
      rc_age, rc_sig, InterpolateCalibrationCurve(breaks, calibration_curve))
    unmodelled$probability <- unmodelled$probability / sum(unmodelled$probability * resolution)
    graphics::polygon(
      c(unmodelled$calendar_age_BP, rev(unmodelled$calendar_age_BP)),
      c(unmodelled$probability, rep(0, length(unmodelled$probability))),
      border = NA,
      col = grDevices::grey(0.1, alpha = 0.25))
  }
  if (show_hpd_ranges) {
    graphics::lines(smoothed_density, lwd=1, col='black', lty=2)
    hpd_probability <- switch(
      interval_width,
      "1sigma" = 0.682,
      "2sigma" = 0.954,
      "bespoke" = bespoke_probability)
    hpd <- FindHPD(smoothed_density$x, smoothed_density$y, hpd_probability)
    graphics::segments(hpd$start_ages, hpd$height, hpd$end_ages, hpd$height, lwd=4, col='black', lend='butt')
  }

  .AddLegendToDensitySamplePlot(
    interval_width,
    bespoke_probability,
    output_data$input_data$calibration_curve_name,
    show_hpd_ranges,
    show_unmodelled_density,
    hpd)
}

.AddLegendToDensitySamplePlot <- function(
    interval_width,
    bespoke_probability,
    calcurve_name,
    show_hpd_ranges,
    show_unmodelled_density,
    hpd){
  ci_label <- switch(
    interval_width,
    "1sigma" = expression(paste(sigma, " interval")),
    "2sigma"  = expression(paste("2", sigma, " interval")),
    "bespoke" = paste0(round(100 * bespoke_probability), "% interval"))

  legend_labels <- c(
    substitute(paste(""^14, "C determination ")),
    gsub("intcal", "IntCal", calcurve_name),
    ci_label)
  lty <- c(-1, 1, 2)
  lwd <- c(-1, 1, 1)
  pch <- c(15, NA, NA)
  col <- c(grDevices::rgb(1, 0, 0, .5), "blue", "blue")

  if (show_unmodelled_density) {
    legend_labels <- c(legend_labels, "Unmodelled density")
    lty <- c(lty, -1)
    lwd <- c(lwd, -1)
    pch <- c(pch, 15)
    col <- c(col, grDevices::grey(0.1, alpha = 0.25))
  }

  legend_labels <- c(legend_labels, "Posterior density")
  lty <- c(lty, 1)
  lwd <- c(lwd, 1)
  pch <- c(pch, NA)
  col <- c(col, "purple")

  if (show_hpd_ranges) {
    hpd_label <- switch(
      interval_width,
      "1sigma" = "68.3% HPD",
      "2sigma"  = "95.4% HPD",
      "bespoke" = paste0(round(100 * bespoke_probability, 3), "% HPD"))
    legend_labels <- c(legend_labels, hpd_label)
    lty <- c(lty, 1)
    lwd <- c(lwd, 2)
    pch <- c(pch, NA)
    col <- c(col, "black")

    for (i in rev(seq_along(hpd$start_ages))) {
      auc_string <- paste0("(", round(hpd$area_under_curve[i] * 100, digits = 1), "%)")
      legend_labels <- c(
        legend_labels,
        paste("   ", round(hpd$end_ages[i]), "-", round(hpd$start_ages[i]), auc_string))
      lty <- c(lty, NA)
      lwd <- c(lwd, NA)
      pch <- c(pch, NA)
      col <- c(col, NA)
    }
  }

  graphics::legend(
    "topright", legend = legend_labels, lty = lty, lwd=lwd, pch = pch, col = col)

}