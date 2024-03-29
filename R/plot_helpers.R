
.PlotCalibrationCurveAndInputData <- function(
    plot_cal_age_scale,
    xlim,
    calibration_curve,
    rc_determinations,
    plot_14C_age,
    calibration_curve_colour,
    calibration_curve_bg,
    interval_width,
    bespoke_probability,
    title = "Calendar age density estimate"){

  .PlotCalibrationCurve(
    plot_cal_age_scale,
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
    plot_cal_age_scale,
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

  xlim <- .ConvertCalendarAge(plot_cal_age_scale, xlim)
  cal_age <- .ConvertCalendarAge(plot_cal_age_scale, calibration_curve$calendar_age)

  if (plot_cal_age_scale == "AD") {
    x_label <- "Calendar Age (cal AD)"
  } else if (plot_cal_age_scale == "BC") {
    x_label <- "Calendar Age (cal BC)"
  } else {
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


.SetUpDensityPlot <- function(plot_cal_age_scale,  xlim, ylim) {

  xlim <- .ConvertCalendarAge(plot_cal_age_scale, xlim)

  graphics::par(new = TRUE)
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

.PlotSPDEstimateOnCurrentPlot <- function(plot_cal_age_scale, calendar_age_BP, probability) {
  SPD_colour <- grDevices::grey(0.1, alpha = 0.3)

  cal_age <- .ConvertCalendarAge(plot_cal_age_scale, calendar_age_BP)

  graphics::polygon(
    c(cal_age, rev(cal_age)),
    c(probability, rep(0, length(probability))),
    border = NA,
    col = SPD_colour)
}

.ConvertCalendarAge <- function(plot_cal_age_scale, cal_age) {
  if (plot_cal_age_scale == "AD") {
    new_cal_age <- 1950 - cal_age
  } else if (plot_cal_age_scale == "BC") {
    new_cal_age <- cal_age - 1949
  } else {
    new_cal_age <- cal_age
  }
  return (new_cal_age)
}
