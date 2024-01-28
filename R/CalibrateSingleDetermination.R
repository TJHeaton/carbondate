#' Calibrate a Single Radiocarbon Determination
#'
#' Uses the supplied calibration curve to calibrate a single radiocarbon
#' determination and uncertainty (expressed either in terms of radiocarbon age, or
#' as an F\eqn{{}^{14}}C concentration) and obtain its calendar age probability
#' density estimate.
#'
#' @param rc_determination A single observed radiocarbon determination
#' provided either as the radiocarbon age (in \eqn{{}^{14}}C yr BP) or the F\eqn{{}^{14}}C concentration.
#' @param rc_sigma The corresponding measurement uncertainty of the radiocarbon determination
#' (must be in the same units as above, i.e., reported as \eqn{{}^{14}}C age or F\eqn{{}^{14}}C)
#' @param calibration_curve A dataframe which must contain one column `calendar_age_BP`, and also
#' columns `c14_age` and `c14_sig` or `f14c` and `f14c_sig` (or both sets).
#' This format matches the curves supplied with this package, e.g., [carbondate::intcal20],
#' [carbondate::intcal13], which contain all 5 columns.
#' @param F14C_inputs `TRUE` if the provided `rc_determination` is an F\eqn{{}^{14}}C
#' concentration and `FALSE` if it is a radiocarbon age. Defaults to `FALSE`.
#' @param resolution The distance between the calendar ages at which to calculate the calendar age probability.
#' Default is 1.
#' @param plot_output `TRUE` if you wish to plot the determination, the calibration curve,
#' and the posterior calibrated age estimate on the same plot. Defaults to `FALSE`
#' @param plot_cal_age_scale Only for usage when `plot_output = TRUE`.
#' The calendar scale to use for the x-axis. Allowed values are "BP", "AD" and "BC". The default
#' is "BP", corresponding to plotting in cal yr BP.
#' @param interval_width Only for usage when `plot_output = TRUE`. The confidence intervals to show for the
#' calibration curve and for the highest posterior density ranges.
#' Choose from one of "1sigma" (68.3%), "2sigma" (95.4%) and "bespoke". Default is "2sigma".
#' @param bespoke_probability The probability to use for the confidence interval
#' if "bespoke" is chosen above. E.g. if 0.95 is chosen, then the 95% confidence
#' interval is calculated. Ignored if "bespoke" is not chosen.
#' @param denscale Whether to scale the vertical range of the calendar age density plot
#' relative to the calibration curve plot (optional). Default is 3 which means
#' that the maximum calendar age density will be at 1/3 of the height of the plot.
#' @param plot_pretty logical, defaulting to `TRUE`. If set `TRUE` then will select pretty plotting
#' margins (that create sufficient space for axis titles and rotates y-axis labels). If `FALSE` will
#' implement current user values.
#'
#' @export
#'
#' @return A data frame with one column `calendar_age_BP` containing the calendar
#' ages, and the other column `probability` containing the probability at that
#' calendar age.
#'
#' @examples
#' # Calibration of a single determination expressed as 14C age BP
#' calib <- CalibrateSingleDetermination(860, 35, intcal20)
#' plot(calib, type = "l", xlim = c(1000, 600))
#'
#' # Incorporating an automated plot to visualise the calibration
#' CalibrateSingleDetermination(860, 35, intcal20, plot_output = TRUE)
#'
#' # Calibration of a single (old) determination expressed as 14C age BP
#' calib <- CalibrateSingleDetermination(31020, 100, intcal20)
#' plot(calib, type = "l", xlim = c(36500, 34500))
#'
#' # Calibration of a single (old) determination expressed as F14C concentration
#' calib <- CalibrateSingleDetermination(
#'     0.02103493, 0.0002618564, intcal20, F14C_inputs = TRUE)
#' plot(calib, type = "l", xlim = c(36500, 34500))
#'
#' # Calibration of a single determination expressed as 14C age BP
#' # against SHCal20 (and creating an automated plot)
#' CalibrateSingleDetermination(1413, 25, shcal20, plot_output = TRUE)
#'
#' # Implementing a bespoke confidence interval level and plot in AD
#' CalibrateSingleDetermination(
#'     1413,
#'     25,
#'     shcal20,
#'     plot_output = TRUE,
#'     plot_cal_age_scale = "AD",
#'     interval_width = "bespoke",
#'     bespoke_probability = 0.8)
#'
#' # Changing denscale (so the calendar age density takes up less space)
#' CalibrateSingleDetermination(
#'     1413,
#'     25,
#'     shcal20,
#'     plot_output = TRUE,
#'     interval_width = "bespoke",
#'     bespoke_probability = 0.8,
#'     denscale = 5)
CalibrateSingleDetermination <- function(
    rc_determination,
    rc_sigma,
    calibration_curve,
    F14C_inputs = FALSE,
    resolution = 1,
    plot_output = FALSE,
    plot_cal_age_scale = "BP",
    interval_width = "2sigma",
    bespoke_probability = NA,
    denscale = 3,
    plot_pretty = TRUE) {

  arg_check <- .InitializeErrorList()
  .CheckNumber(arg_check, rc_determination)
  .CheckNumber(arg_check, rc_sigma)
  .CheckFlag(arg_check, F14C_inputs)
  .CheckFlag(arg_check, plot_output)
  .CheckChoice(arg_check, plot_cal_age_scale, c("BP", "AD", "BC"))
  .CheckCalibrationCurve(arg_check, calibration_curve, NA)
  .CheckNumber(arg_check, denscale, lower = 0)
  .CheckNumber(arg_check, resolution, lower = 0.01)
  .ReportErrors(arg_check)

  calibration_curve_name <- deparse(substitute(calibration_curve))
  calendar_age_range <- range(calibration_curve$calendar_age_BP)
  calibration_curve <- InterpolateCalibrationCurve(
    seq(from = calendar_age_range[1], to = calendar_age_range[2], by = resolution),
    calibration_curve,
    F14C_outputs = F14C_inputs)

  probabilities <- .ProbabilitiesForSingleDetermination(rc_determination, rc_sigma, F14C_inputs, calibration_curve)

  if(plot_output == TRUE) {
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

    .PlotIndependentCalibration(
      rc_determination = rc_determination,
      rc_sigma = rc_sigma,
      calendar_ages = calibration_curve$calendar_age_BP,
      probabilities = probabilities,
      calibration_curve = calibration_curve,
      calibration_curve_name = calibration_curve_name,
      plot_cal_age_scale = plot_cal_age_scale,
      interval_width = interval_width,
      bespoke_probability = bespoke_probability,
      F14C_inputs = F14C_inputs,
      denscale = denscale)
  }

  return_data <- data.frame(calendar_age_BP=calibration_curve$calendar_age_BP, probability=probabilities)
  if (plot_output == TRUE) {
    invisible(return_data)
  } else {
    return(return_data)
  }
}


.ProbabilitiesForSingleDetermination <- function(rc_determination, rc_sigma, F14C_inputs, calibration_curve) {

  calendar_age_diff <- abs(diff(calibration_curve$calendar_age_BP)[1])

  if (F14C_inputs) {
    calibration_curve <- .AddF14cColumns(calibration_curve)
    calcurve_rc_ages <- calibration_curve$f14c
    calcurve_rc_sigs <- calibration_curve$f14c_sig
  } else {
    calibration_curve <- .AddC14ageColumns(calibration_curve)
    calcurve_rc_ages <- calibration_curve$c14_age
    calcurve_rc_sigs <- calibration_curve$c14_sig
  }

  probabilities <- stats::dnorm(rc_determination, mean=calcurve_rc_ages, sd=sqrt(calcurve_rc_sigs^2 + rc_sigma^2))
  probabilities <- probabilities / (sum(probabilities) * calendar_age_diff)
  return(probabilities)
}


.PlotIndependentCalibration <- function(
    rc_determination,
    rc_sigma,
    calendar_ages,
    probabilities,
    calibration_curve,
    calibration_curve_name,
    F14C_inputs,
    plot_cal_age_scale,
    interval_width = "2sigma",
    bespoke_probability = NA,
    denscale = NA,
    show_hpd_ranges = TRUE,
    prob_cutoff = 0.00001) {

  # What domain to plot in
  plot_14C_age <- !F14C_inputs

  if (F14C_inputs == TRUE) {
    calibration_curve <- .AddF14cColumns(calibration_curve)
  } else {
    calibration_curve <- .AddC14ageColumns(calibration_curve)
  }

  rc_age <- rc_determination
  rc_sig <- rc_sigma

  # Find the calendar age range to plot
  cumulativeprobabilities <- cumsum(probabilities) / sum(probabilities)
  calendar_age_bound_1 <- calibration_curve$calendar_age_BP[
    min(which(cumulativeprobabilities > prob_cutoff))]
  calendar_age_bound_2 <- calibration_curve$calendar_age_BP[
    max(which(cumulativeprobabilities <= (1 - prob_cutoff)))]
  xrange_BP <- sort(c(calendar_age_bound_1, calendar_age_bound_2))
  # Extend by 20% either side
  xrange_BP <- xrange_BP + 0.2 * c(-1, 1) * diff(xrange_BP)

  if (plot_14C_age == FALSE) {
    title <- substitute(
      paste(
        "Individual Calibration of ",
        f14c_age,
        "\u00B1",
        f14c_sig,
        " F"^14,
        "C"),
      list(f14c_age = signif(rc_age, 2), f14c_sig = signif(rc_sig, 2)))
  } else {
    title <- substitute(
      paste(
        "Individual Calibration of ",
        c14_age,
        "\u00B1",
        c14_sig,
        ""^14,
        "C yr BP"),
      list(c14_age = round(rc_age), c14_sig = round(rc_sig, 1)))
  }

  .PlotCalibrationCurve(
    plot_cal_age_scale,
    xlim = rev(xrange_BP),
    plot_14C_age = plot_14C_age,
    calibration_curve = calibration_curve,
    calibration_curve_colour = "blue",
    calibration_curve_bg = grDevices::rgb(0, 0, 1, .3),
    interval_width = interval_width,
    bespoke_probability = bespoke_probability,
    title = title)

  calendar_ages <- .ConvertCalendarAge(plot_cal_age_scale, calendar_ages)
  xrange <- .ConvertCalendarAge(plot_cal_age_scale, xrange_BP)

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

  .SetUpDensityPlot(
    plot_cal_age_scale,
    xlim = rev(xrange_BP),
    ylim = c(0, denscale * max(probabilities)))

  graphics::polygon(
    c(calendar_ages, rev(calendar_ages)),
    c(probabilities, rep(0, length(probabilities))),
    border = NA,
    col = grDevices::grey(0.1, alpha = 0.25))

  if (show_hpd_ranges) {
    hpd_probability <- switch(
      interval_width,
      "1sigma" = 0.682,
      "2sigma" = 0.954,
      "bespoke" = bespoke_probability)
    hpd <- FindHPD(as.numeric(calendar_ages), probabilities, hpd_probability)
    if (length(hpd$start_ages) > 0) {
      graphics::segments(hpd$start_ages, hpd$height, hpd$end_ages, hpd$height, lwd=4, col='black', lend='butt')
    } else {
      warning("HPD Interval could not be found")
    }
  }

  .AddLegendToIndividualCalibrationPlot(
    interval_width,
    bespoke_probability,
    calcurve_name = calibration_curve_name,
    show_hpd_ranges = TRUE,
    hpd)

}


.AddLegendToIndividualCalibrationPlot <- function(
    interval_width,
    bespoke_probability,
    calcurve_name,
    show_hpd_ranges,
    hpd){

  ci_label <- switch(
    interval_width,
    "1sigma" = expression(paste("1", sigma, " interval")),
    "2sigma"  = expression(paste("2", sigma, " interval")),
    "bespoke" = paste0(round(100 * bespoke_probability), "% interval"))

  legend_labels <- c(
    substitute(paste(""^14, "C determination ")),
    gsub("intcal", "IntCal", gsub("shcal", "SHCal", calcurve_name)), # Both IntCal and SHCal
    ci_label)
  lty <- c(-1, 1, 2, -1)
  lwd <- c(-1, 1, 1, -1)
  pch <- c(15, NA, NA, 15)
  col <- c(grDevices::rgb(1, 0, 0, .5),
           "blue",
           "blue",
           grDevices::grey(0.1, alpha = 0.25))

  legend_labels <- c(legend_labels, "Calibrated Age")

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

