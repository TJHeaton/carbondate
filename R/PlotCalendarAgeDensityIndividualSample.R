#' Plots the posterior calendar ages for an individual determination
#'
#' Once a function has been run to calibrate a set of radiocarbon
#' determinations, the posterior density for a single determination can be
#' plotted using this function.
#'
#' @inheritParams PlotPredictiveCalendarAgeDensity
#' @param ident the determination you want to show the individual posterior
#' calendar age for.
#' @param resolution The distance between histogram breaks for the calendar age density.
#' Must be an integer greater than one.
#' @param interval_width The confidence intervals to show for the
#' calibration curve. Choose from one of "1sigma", "2sigma" and "bespoke".
#' Default is "2sigma".
#' @param bespoke_probability The probability to use for the confidence interval
#' if "bespoke" is chosen above. E.g. if 0.95 is chosen, then the 95% confidence
#' interval is calculated. Ignored if "bespoke" is not chosen.
#'
#' @return No return value
#' @export
#'
#' @examples
#' # Plot results for the 10th determinations
#' PlotCalendarAgeDensityIndividualSample(10, walker_example_output)
PlotCalendarAgeDensityIndividualSample <- function(
    ident,
    output_data,
    calibration_curve = NULL,
    plot_14C_age = TRUE,
    resolution = 5,
    interval_width = "2sigma",
    bespoke_probability = NA,
    n_burn = NA,
    n_end = NA) {

  arg_check <- checkmate::makeAssertCollection()
  checkmate::assertInt(ident, add = arg_check)
  .CheckOutputData(arg_check, output_data)
  n_iter = output_data$input_parameters$n_iter
  n_thin = output_data$input_parameters$n_thin

  .CheckCalibrationCurveFromOutput(arg_check, output_data, calibration_curve)
  checkmate::assertInt(resolution, na.ok = FALSE, add = arg_check, lower = 1)
  checkmate::assertInt(n_burn, lower = 0, upper = n_iter - 100 * n_thin, na.ok = TRUE)
  checkmate::reportAssertions(arg_check)

  if (is.null(calibration_curve)) {
    calibration_curve = get(output_data$input_data$calibration_curve_name)
  }
  rc_determinations <- output_data$input_data$rc_determinations
  rc_sigmas <- output_data$input_data$rc_sigmas
  F14C_inputs <-output_data$input_data$F14C_inputs

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

  calendar_age <- output_data$calendar_ages[, ident]
  rc_age <- rc_determinations[ident]
  rc_sig <- rc_sigmas[ident]

  n_out <- length(calendar_age)
  if (is.na(n_burn)) {
    n_burn = floor(n_out / 2)
  } else {
    n_burn = floor(n_burn / n_thin)
  }
  if (is.na(n_end)) {
    n_end = n_out
  } else {
    n_end = floor(n_end / n_thin)
  }
  calendar_age <- calendar_age[(n_burn+1):n_end]

  # Find the calendar age range to plot
  xrange <- range(calendar_age)
  xrange[1] = floor(xrange[1])
  if (resolution > 1) while (xrange[1] %% resolution != 0) xrange[1] = xrange[1] - 1

  xrange[2] = ceiling(xrange[2])
  if (resolution > 1) while (xrange[2] %% resolution != 0) xrange[2] = xrange[2] + 1

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


  plot_AD = any(calendar_age < 0)
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
  relative_height = 0.1
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

  diff = diff(xrange)
  xrange = xrange + c(-1,1) * diff * 0.4
  finalhist <- graphics::hist(
    calendar_age,
    prob = TRUE,
    breaks = breaks,
    xlim = rev(xrange),
    axes = FALSE,
    xlab = NA,
    ylab = NA,
    main = "",
    xaxs = "i",
    ylim = c(0, 3 * max(temphist$density)))
  invisible(finalhist)
}
