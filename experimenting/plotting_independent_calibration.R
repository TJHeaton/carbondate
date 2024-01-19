CalibrateSingleDetermination(
  0.02103493, 9.164975e-05,
  calibration_curve = intcal20,
  F14C_inputs = TRUE,
  plot_output = TRUE,
  interval_width = "1sigma",
  bespoke_probability = 0.8)



# Testing
rc_determination <- 1413
rc_sigma <- 25
calibration_curve <- intcal20
F14C_inputs <- FALSE
plot_output <- TRUE
prob_cutoff <- 0.00001
interval_width <- "2sigma"
show_hpd_ranges <- TRUE
prob_cutoff <- 0.00001


probabilities <- .ProbabilitiesForSingleDetermination(rc_determination, rc_sigma, F14C_inputs, calibration_curve)
calendar_ages <- calibration_curve$calendar_age


if (plot_14C_age == TRUE) {
  calibration_curve <- .AddC14ageColumns(calibration_curve)
  if (F14C_inputs == TRUE) {
    converted <- .ConvertF14cTo14Cage(rc_determination, rc_sigma)
    rc_determination <- converted$c14_age
    rc_sigma <- converted$c14_sig
  }
} else {
  calibration_curve <- .AddF14cColumns(calibration_curve)
  if (F14C_inputs == FALSE) {
    converted <- .Convert14CageToF14c(rc_determination, rc_sigma)
    rc_determination <- converted$f14c
    rc_sigma <- converted$f14c_sig
  }
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
      "F"^14,
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

plot_AD <- any(xrange < 0)
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
  "IntCal20",
  show_hpd_ranges = TRUE,
  hpd)

