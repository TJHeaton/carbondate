#' Find the summed probability distribution (SPD) for a set of radiocarbon observations
#'
#' Takes a set of radiocarbon determinations and uncertainties, independently
#' calibrates each one, and then averages the resultant calendar age estimates
#' to give the SPD estimate. \cr \cr
#' \strong{Important:} This function should not be used for inference as SPDs are not statistically rigorous.
#' Instead use either of:
#' \itemize{
#' \item the Bayesian non-parametric summarisation approaches [carbondate::PolyaUrnBivarDirichlet]
#'   or [carbondate::WalkerBivarDirichlet];
#' \item or the Poisson process rate approach [carbondate::PPcalibrate]
#' }
#' The SPD function provided here is only intended for comparison. We provide a inbuilt plotting option to show the SPD alongside the
#' determinations and the calibration curve. \cr \cr
#'
#' @inheritParams CalibrateSingleDetermination
#' @param calendar_age_range_BP A vector of length 2 with the start and end
#' calendar age BP to calculate the SPD over. These two values must be in order
#' with start value less than end value, e.g., (600, 1700) cal yr BP.
#' The code will extend this range (by 400 cal yrs each side) for conservativeness
#' @param rc_determinations A vector of observed radiocarbon determinations. Can be provided either as
#' \eqn{{}^{14}}C ages (in \eqn{{}^{14}}C yr BP) or as F\eqn{{}^{14}}C concentrations.
#' @param rc_sigmas A vector of the (1-sigma) measurement uncertainties for the
#' radiocarbon determinations. Must be the same length as `rc_determinations` and
#' given in the same units.
#' @param F14C_inputs `TRUE` if the provided `rc_determinations` are F\eqn{{}^{14}}C
#' concentrations and `FALSE` if they are radiocarbon ages. Defaults to `FALSE`.
#' @param plot_output `TRUE` if you wish to plot the determinations, the calibration curve,
#' and the SPD on the same plot. Defaults to `FALSE`
#' @param interval_width Only for usage when `plot_output = TRUE`. The confidence intervals to show for the
#' calibration curve. Choose from one of "1sigma" (68.3%), "2sigma" (95.4%) and "bespoke". Default is "2sigma".
#' @param bespoke_probability The probability to use for the confidence interval
#' if "bespoke" is chosen above. E.g. if 0.95 is chosen, then the 95% confidence
#' interval is calculated. Ignored if "bespoke" is not chosen.
#' @param denscale Whether to scale the vertical range of the calendar age density plot
#' relative to the calibration curve plot (optional). Default is 3 which means
#' that the maximum SPD will be at 1/3 of the height of the plot.
#'
#' @return A data frame with one column `calendar_age_BP` containing the calendar
#' ages, and the other column `probability` containing the probability at that
#' calendar age
#' @export
#'
#' @seealso [carbondate::PolyaUrnBivarDirichlet],
#' [carbondate::WalkerBivarDirichlet] for rigorous non-parametric Bayesian alternatives; and
#' [carbondate::PPcalibrate] for a rigorous variable-rate Poisson process alternative.
#'
#' @examples
#' # An example using 14C age BP and the IntCal 20 curve
#' SPD <- FindSummedProbabilityDistribution(
#'    calendar_age_range_BP=c(400, 1700),
#'    rc_determinations=c(602, 805, 1554),
#'    rc_sigmas=c(35, 34, 45),
#'    calibration_curve=intcal20)
#' plot(SPD, type = "l",
#'    xlim = rev(range(SPD$calendar_age_BP)),
#'    xlab = "Calendar Age (cal yr BP)")
#'
#' # Using the inbuilt plotting features
#' SPD <- FindSummedProbabilityDistribution(
#'    calendar_age_range_BP=c(400, 1700),
#'    rc_determinations=c(602, 805, 1554),
#'    rc_sigmas=c(35, 34, 45),
#'    calibration_curve=intcal20,
#'    plot_output = TRUE,
#'    interval_width = "bespoke",
#'    bespoke_probability = 0.8)
#'
#' # An different example using F14C concentrations and the IntCal 13 curve
#' SPD <- FindSummedProbabilityDistribution(
#'    calendar_age_range_BP=c(400, 2100),
#'    rc_determinations=c(0.8, 0.85, 0.9),
#'    rc_sigmas=c(0.01, 0.015, 0.012),
#'    F14C_inputs=TRUE,
#'    calibration_curve=intcal13)
#' plot(SPD, type = "l",
#'    xlim = rev(range(SPD$calendar_age_BP)),
#'    xlab = "Calendar Age (cal yr BP)")
FindSummedProbabilityDistribution <- function(
    calendar_age_range_BP,
    rc_determinations,
    rc_sigmas,
    calibration_curve,
    F14C_inputs = FALSE,
    plot_output = FALSE,
    interval_width = "2sigma",
    bespoke_probability = NA,
    denscale = 3) {

  arg_check <- .InitializeErrorList()
  .CheckNumberVector(arg_check, calendar_age_range_BP, len = 2)
  .CheckFlag(arg_check, F14C_inputs)
  .CheckFlag(arg_check, plot_output)
  .CheckNumber(arg_check, denscale, lower = 0)
  .CheckInputData(arg_check, rc_determinations, rc_sigmas, F14C_inputs)
  .CheckCalibrationCurve(arg_check, calibration_curve, NA)
  .ReportErrors(arg_check)


  calibration_data_range <- range(calibration_curve$calendar_age)
  range_margin <- 400

  new_calendar_ages <- seq(
    max(calibration_data_range[1], calendar_age_range_BP[1] - range_margin),
    min(calibration_data_range[2], calendar_age_range_BP[2] + range_margin),
    by = 1)

  interpolated_calibration_data <- InterpolateCalibrationCurve(
    new_calendar_ages, calibration_curve, F14C_inputs)

  # Each column of matrix gives probabilities per calendar age for a given determination
  individual_probabilities <- mapply(
    .ProbabilitiesForSingleDetermination,
    rc_determinations,
    rc_sigmas,
    MoreArgs = list(F14C_inputs = F14C_inputs, calibration_curve = interpolated_calibration_data))

  probabilities_per_calendar_age <- apply(individual_probabilities, 1, sum) /
    dim(individual_probabilities)[2]

  SPD <- data.frame(
    calendar_age_BP=new_calendar_ages,
    probability=probabilities_per_calendar_age)

  if(plot_output == TRUE) {
    .PlotSPD(
      rc_determinations = rc_determinations,
      rc_sigmas = rc_sigmas,
      calendar_ages = new_calendar_ages,
      probabilities = probabilities_per_calendar_age,
      calibration_curve = interpolated_calibration_data,
      calibration_curve_name = deparse(substitute(calibration_curve)),
      F14C_inputs = F14C_inputs,
      interval_width = interval_width,
      bespoke_probability = bespoke_probability,
      denscale = denscale)
  }

  return(SPD)
}


.PlotSPD <- function(
    rc_determinations,
    rc_sigmas,
    calendar_ages,
    probabilities,
    calibration_curve,
    calibration_curve_name,
    F14C_inputs,
    interval_width = "2sigma",
    bespoke_probability = NA,
    denscale = NA) {

  # What domain to plot in
  plot_14C_age <- !F14C_inputs

  if (F14C_inputs == TRUE) {
    calibration_curve <- .AddF14cColumns(calibration_curve)
  } else {
    calibration_curve <- .AddC14ageColumns(calibration_curve)
  }

  rc_ages <- rc_determinations
  rc_sigs <- rc_sigmas

  # Calculate calendar age plotting range
  xrange <- range(calibration_curve$calendar_age)

  title <- "Summed Probability Distribution \n(Do Not Use For Inference)"

  plot_AD <- FALSE # Plot in calendar year

  # Change plotting parameters
  opar <- graphics::par(
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 4, 2) + 0.1,
    las = 1)
  # Revert to main environment pars after complete
  on.exit(graphics::par(opar))

  .PlotCalibrationCurveAndInputData(
    plot_AD,
    xlim = rev(xrange),
    rc_determinations = rc_determinations,
    plot_14C_age = plot_14C_age,
    calibration_curve = calibration_curve,
    calibration_curve_colour = "blue",
    calibration_curve_bg = grDevices::rgb(0, 0, 1, .3),
    interval_width = interval_width,
    bespoke_probability = bespoke_probability,
    title = title)

  # Plot the posterior cal age on the x-axis
  .SetUpDensityPlot(plot_AD,
                    xlim = rev(xrange),
                    ylim = c(0, denscale * max(probabilities)))
  graphics::polygon(
    c(calendar_ages, rev(calendar_ages)),
    c(probabilities, rep(0, length(probabilities))),
    border = NA,
    col = grDevices::grey(0.1, alpha = 0.25))

  .AddLegendToSPDPlot(
    interval_width,
    bespoke_probability,
    calcurve_name = calibration_curve_name)

}

.AddLegendToSPDPlot <- function(
    interval_width,
    bespoke_probability,
    calcurve_name){

  ci_label <- switch(
    interval_width,
    "1sigma" = expression(paste("1", sigma, " interval")),
    "2sigma"  = expression(paste("2", sigma, " interval")),
    "bespoke" = paste0(round(100 * bespoke_probability), "% interval"))

  legend_labels <- c(
    gsub("intcal", "IntCal", gsub("shcal", "SHCal", calcurve_name)), # Both IntCal and SHCal
    ci_label)
  lty <- c(1, 2, -1)
  lwd <- c(1, 1, -1)
  pch <- c(NA, NA, 15)
  col <- c("blue",
           "blue",
           grDevices::grey(0.1, alpha = 0.25))

  legend_labels <- c(legend_labels, "SPD")

  graphics::legend(
    "topright", legend = legend_labels, lty = lty, lwd=lwd, pch = pch, col = col)

}










