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
    resolution = 5,
    interval_width = "2sigma",
    bespoke_probability = NA) {

  arg_check <- checkmate::makeAssertCollection()
  checkmate::assertInt(ident, add = arg_check)
  .CheckOutputData(arg_check, output_data)
  .CheckCalibrationCurveFromOutput(arg_check, output_data, calibration_curve)
  checkmate::assertInt(resolution, na.ok = FALSE, add = arg_check, lower = 1)
  checkmate::reportAssertions(arg_check)

  if (is.null(calibration_curve)) {
    calibration_curve = get(output_data$input_data$calibration_curve_name)
  }
  rc_determinations <- output_data$input_data$rc_determinations
  rc_sigmas <- output_data$input_data$rc_sigmas
  F14C_inputs <-output_data$input_data$F14C_inputs

  calendar_age <- output_data$calendar_ages[, ident]
  rc_age <- rc_determinations[ident]
  rc_sig <- rc_sigmas[ident]

  n_out <- length(calendar_age)
  n_burn <- floor(n_out / 2)
  calendar_age <- calendar_age[(n_burn+1):n_out]

  # Find the calendar age range to plot
  xrange <- range(calendar_age)
  xrange[1] = floor(xrange[1])
  if (resolution > 1) while (xrange[1] %% resolution != 0) xrange[1] = xrange[1] - 1

  xrange[2] = ceiling(xrange[2])
  if (resolution > 1) while (xrange[2] %% resolution != 0) xrange[2] = xrange[2] + 1

  cal_age_ind_min <- which.min(abs(calibration_curve$calendar_age - xrange[1]))
  cal_age_ind_max <- which.min(abs(calibration_curve$calendar_age - xrange[2]))
  calendar_age_indices <- cal_age_ind_min:cal_age_ind_max
  if (F14C_inputs) {
    range_f14c <- range(calibration_curve$f14c[calendar_age_indices])
    yrange <- range_f14c + 0.25 * c(-1, 1) * diff(range_f14c)
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
    range_c14 <- range(calibration_curve$c14_age[calendar_age_indices])
    yrange <- range_c14 + 0.25 * c(-1, 1) * diff(range_c14)
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
    ylim = yrange,
    F14C_inputs = F14C_inputs,
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
