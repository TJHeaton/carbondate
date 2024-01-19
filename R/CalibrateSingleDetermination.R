#' Calibrate a single radiocarbon determination
#'
#' Uses the supplied calibration curve to take a single radiocarbon
#' determination and uncertainty (expressed as the radiocarbon age BP or F14C concentration) and
#' calculate the calendar age probability density for it.
#'
#' @param rc_determination A single observed radiocarbon determination
#' (\eqn{{}^{14}}C BP age or F14C concentration)
#' @param rc_sigma The uncertainty of the radiocarbon determination in the same units
#' @param calibration_curve A dataframe which must contain one column `calendar_age_BP`, and also
#' columns `c14_age` and `c14_sig` or `f14c` and `f14c_sig` (or both sets).
#' This format matches the curves supplied with this package e.g. [carbondate::intcal20],
#' which contain all 5 columns.
#' @param F14C_inputs `TRUE` if the provided rc_determinations are F14C concentrations and `FALSE`
#' if they are radiocarbon age BP. Defaults to `FALSE`.
#' @param plot_output `TRUE` if you wish to plot the determination, calibration curve, and the
#' posterior calibrated age estimate. Defaults to `FALSE`
#' @param interval_width Only for usage when `plot_output == TRUE`. The confidence intervals to show for the
#' calibration curve and for the highest posterior density ranges.
#' Choose from one of "1sigma" (68.3%), "2sigma" (95.4%) and "bespoke". Default is "2sigma".
#' @param bespoke_probability The probability to use for the confidence interval
#' if "bespoke" is chosen above. E.g. if 0.95 is chosen, then the 95% confidence
#' interval is calculated. Ignored if "bespoke" is not chosen.
#'
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
#' calib <- CalibrateSingleDetermination(860, 35, intcal20, plot_output = TRUE)
#'
#' # Calibration of a single (old) determination expressed as 14C age BP
#' calib <- CalibrateSingleDetermination(
#'     31020, 100, intcal20)
#' plot(calib, type = "l", xlim = c(36500, 34500))
#'
#' # Calibration of a single determination expressed as F14C concentration
#' calib <- CalibrateSingleDetermination(
#'     0.02103493, 0.0002618564, intcal20, F14C_inputs = TRUE)
#' plot(calib, type = "l", xlim = c(36500, 34500))
#'
#' # Calibration of a single determination expressed as 14C age BP
#' # against SHCal20 (and creating an automated plot)
#' calib <- CalibrateSingleDetermination(1413, 25, shcal20, plot_output = TRUE)
#'
#' # Implementing a bespoke confidence interval level
#' calib <- CalibrateSingleDetermination(1413, 25, shcal20,
#'     plot_output = TRUE, interval_width = "bespoke", bespoke_probability = 0.8)
CalibrateSingleDetermination <- function(
    rc_determination, rc_sigma, calibration_curve, F14C_inputs = FALSE,
    plot_output = FALSE, interval_width = "2sigma", bespoke_probability = NA) {

  arg_check <- .InitializeErrorList()
  .CheckNumber(arg_check, rc_determination)
  .CheckNumber(arg_check, rc_sigma)
  .CheckFlag(arg_check, F14C_inputs)
  .CheckCalibrationCurve(arg_check, calibration_curve, F14C_inputs)
  .ReportErrors(arg_check)

  probabilities <- .ProbabilitiesForSingleDetermination(rc_determination, rc_sigma, F14C_inputs, calibration_curve)

  if(plot_output == TRUE) {
    .PlotIndependentCalibration(
      rc_determination = rc_determination,
      rc_sigma = rc_sigma,
      calendar_ages = calibration_curve$calendar_age,
      probabilities = probabilities,
      calibration_curve = calibration_curve,
      calibration_curve_name = deparse(substitute(calibration_curve)),
      interval_width = interval_width,
      bespoke_probability = bespoke_probability,
      F14C_inputs = F14C_inputs)
  }

  return(
    data.frame(
      calendar_age_BP=calibration_curve$calendar_age, probability=probabilities))
}


.ProbabilitiesForSingleDetermination <- function(
    rc_determination, rc_sigma, F14C_inputs, calibration_curve) {

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
  probabilities <- probabilities / sum(probabilities)
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
    interval_width = "2sigma",
    bespoke_probability = NA,
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
  cumulativeprobabilities <- cumsum(probabilities)
  min_potential_calendar_age <- calibration_curve$calendar_age_BP[
    min(which(cumulativeprobabilities > prob_cutoff))]
  max_potential_calendar_age <- calibration_curve$calendar_age_BP[
    max(which(cumulativeprobabilities <= (1 - prob_cutoff)))]
  xrange <- c(max_potential_calendar_age, min_potential_calendar_age)
  # Extend by 20% either side
  xrange <- xrange + 0.2 * c(-1, 1) * diff(xrange)

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

  plot_AD <- FALSE # Plot in calendar year
  graphics::par(xaxs = "i", yaxs = "i")
  graphics::par(mar = c(5, 4.5, 4, 2) + 0.1, las = 1)

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
  temp <- graphics::plot(calendar_ages,
                         probabilities,
                         axes = FALSE,
                         xlab = NA, ylab = NA,
                         xlim = rev(xrange),
                         ylim = c(0, 3* max(probabilities)),
                         type = "n")
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
    hpd <- FindHPD(as.numeric(rev(calendar_ages)), rev(probabilities), hpd_probability)
    graphics::segments(hpd$start_ages, hpd$height, hpd$end_ages, hpd$height, lwd=4, col='black', lend='butt')
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

