
.PlotCalibrationCurveAndInputData <- function(
    plot_AD,
    xlim,
    calibration_curve,
    rc_determinations,
    plot_14C_age,
    calibration_curve_colour,
    calibration_curve_bg,
    interval_width,
    bespoke_probability,
    title = "Calendar age density estimate"){

  # Set nice plotting parameters
  graphics::par(
    mgp = c(3, 0.7, 0),
    mar = c(5, 4.5, 4, 2) + 0.1,
    las = 1)

  .PlotCalibrationCurve(
    plot_AD,
    xlim,
    plot_14C_age,
    calibration_curve,
    calibration_curve_colour,
    calibration_curve_bg,
    interval_width,
    bespoke_probability,
    title = title)
  graphics::rug(rc_determinations, side = 2)
}


.PlotCalibrationCurve <- function(
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
    cal_age <- 1950 - calibration_curve$calendar_age
    xlim <- 1950 - xlim
    x_label <- "Calendar Age (cal AD)"
  } else {
    cal_age <- calibration_curve$calendar_age
    x_label <- "Calendar Age (cal yr BP)"
  }

  if (plot_14C_age) {
    rc_age <- calibration_curve$c14_age
    rc_sig <- calibration_curve$c14_sig
    y_label <- expression(paste("Radiocarbon age (", ""^14, "C ", "yr BP)"))
  } else {
    rc_age <- calibration_curve$f14c
    rc_sig <- calibration_curve$f14c_sig
    y_label <- expression(paste("F ", ""^14, "C"))
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
  ylim <- c(
    min(calibration_curve$lb[calendar_age_indices]),
    max(calibration_curve$ub[calendar_age_indices])
    )
  ylim <- ylim + 0.05 * c(-3, 1) * diff(ylim)

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

  if (plot_AD) xlim <- 1950 - xlim
  graphics::plot.default(
    NULL,
    NULL,
    type = "n",
    ylim = ylim,
    xlim = xlim,
    axes = FALSE,
    xlab = NA,
    ylab = NA)
}

