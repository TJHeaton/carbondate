#' Plot Predictive Estimate of Shared Calendar Age Density from Bayesian Non-Parametric
#' DPMM Output
#'
#' @description
#' Given output from one of the Bayesian non-parametric summarisation functions (either
#' [carbondate::PolyaUrnBivarDirichlet] or [carbondate::WalkerBivarDirichlet]) calculate
#' and plot the predictive (summarised/shared) calendar age density and probability intervals
#' on a given calendar age grid.
#'
#' Will show the original set of radiocarbon determinations (those you are summarising),
#' the chosen calibration curve, and the summarised predictive calendar age density on the
#' same plot. Can also optionally show the SPD estimate.
#'
#' \strong{Note:} If all you are interested in is the estimated value of the predictive density
#' on a grid, without an accompanying plot, you can use
#' [carbondate::FindPredictiveCalendarAgeDensity] instead.
#'
#' For more information read the vignette: \cr
#' \code{vignette("Non-parametric-summed-density", package = "carbondate")}
#'
#' @param output_data The return value from one of the Bayesian non-parametric DPMM functions, e.g.
#' [carbondate::PolyaUrnBivarDirichlet] or
#' [carbondate::WalkerBivarDirichlet], or a list, each item containing
#' one of these return values. Optionally, the output data can have an extra list item
#' named `label` which is used to set the label on the plot legend.
#' @param n_posterior_samples Number of samples it will draw, after having removed `n_burn`,
#' from the (thinned) realisations stored in the DPMM outputs to estimate the
#' predictive calendar age density. These samples may be repeats if the number of, post burn-in,
#' realisations is less than `n_posterior_samples`. If not given, 5000 is used.
#' @param calibration_curve This is usually not required since the name of the
#' calibration curve variable is saved in the output data. However, if the
#' variable with this name is no longer in your environment then you should pass
#' the calibration curve here. If provided, this should be a dataframe which
#' should contain at least 3 columns entitled `calendar_age`, `c14_age` and `c14_sig`.
#' This format matches [carbondate::intcal20].
#' @param plot_14C_age Whether to use the radiocarbon age (\eqn{{}^{14}}C yr BP) as
#' the units of the y-axis in the plot. Defaults to `TRUE`. If `FALSE` uses
#' F\eqn{{}^{14}}C concentration instead.
#' @param plot_cal_age_scale The scale to use for the x-axis. Allowed values are
#' "BP", "AD" and "BC".
#' @param show_SPD Whether to calculate and show the summed probability
#' distribution on the plot (optional). Default is `FALSE`.
#' @param show_confidence_intervals Whether to show the pointwise confidence intervals
#' (at the chosen probability level) on the plot. Default is `TRUE`.
#' @param interval_width The confidence intervals to show for both the
#' calibration curve and the predictive density. Choose from one of `"1sigma"` (68.3%),
#' `"2sigma"` (95.4%) and `"bespoke"`. Default is `"2sigma"`.
#' @param bespoke_probability The probability to use for the confidence interval
#' if `"bespoke"` is chosen above. E.g., if 0.95 is chosen, then the 95% confidence
#' interval is calculated. Ignored if `"bespoke"` is not chosen.
#' @param denscale Whether to scale the vertical range of the summarised calendar age density plot
#' relative to the calibration curve plot (optional). Default is 3 which means
#' that the maximum predictive density will be at 1/3 of the height of the plot.
#' @param resolution The distance between calendar ages at which to calculate the predictive shared
#' density. These ages will be created on a regular grid that automatically covers the
#' calendar period of the given set of \eqn{{}^{14}}C samples. Default is 1.
#' @param n_burn The number of MCMC iterations that should be discarded as burn-in (i.e.,
#' considered to be occurring before the MCMC has converged). This relates to the number
#' of iterations (`n_iter`) when running the original update functions (not the thinned `output_data`).
#' Any MCMC iterations before this are not used in the calculations. If not given, the first half of
#' the MCMC chain is discarded. Note: The maximum value that the function will allow is
#' `n_iter - 100 * n_thin` (where `n_iter` and `n_thin` are the arguments given to
#' [carbondate::PolyaUrnBivarDirichlet] or [carbondate::WalkerBivarDirichlet])
#' which would leave only 100 of the (thinned) values in `output_data`.
#' @param n_end The last iteration in the original MCMC chain to use in the calculations. Assumed to be the
#' total number of iterations performed, i.e. `n_iter`, if not given.
#'
#' @return A list, each item containing a data frame of the `calendar_age`, the
#' `density_mean` and the confidence intervals for the density
#' `density_ci_lower` and `density_ci_upper` for each set of output data.
#'
#' @export
#'
#' @seealso [carbondate::FindPredictiveCalendarAgeDensity] if only interested in the estimated value of
#' the predictive density on a grid; [carbondate::PlotNumberOfClusters] and
#' [carbondate::PlotCalendarAgeDensityIndividualSample] for more plotting functions using DPMM output.
#'
#' @examples
#' # NOTE: All these examples are shown with a small n_iter and n_posterior_samples
#' # to speed up execution.
#' # Try n_iter and n_posterior_samples as the function defaults.
#'
#' polya_urn_output <- PolyaUrnBivarDirichlet(
#'     two_normals$c14_age,
#'     two_normals$c14_sig,
#'     intcal20,
#'     n_iter = 500,
#'     show_progress = FALSE)
#' walker_output <- WalkerBivarDirichlet(
#'     two_normals$c14_age,
#'     two_normals$c14_sig,
#'     intcal20,
#'     n_iter = 500,
#'     show_progress = FALSE)
#'
#' # Plot results for a single calibration
#' PlotPredictiveCalendarAgeDensity(polya_urn_output, n_posterior_samples = 50)
#'
#' # Plot results from two calibrations on the same plot
#' PlotPredictiveCalendarAgeDensity(
#'     list(walker_output, polya_urn_output), n_posterior_samples = 50)
#'
#' # Plot and show the 80% confidence interval, show the SPD, add a custom label
#' polya_urn_output$label = "My plot"
#' PlotPredictiveCalendarAgeDensity(
#'     polya_urn_output,
#'     n_posterior_samples = 50,
#'     interval_width = "bespoke",
#'     bespoke_probability = 0.8,
#'     show_SPD = TRUE)
PlotPredictiveCalendarAgeDensity <- function(
    output_data,
    n_posterior_samples = 5000,
    calibration_curve = NULL,
    plot_14C_age = TRUE,
    plot_cal_age_scale = "BP",
    show_SPD = FALSE,
    show_confidence_intervals = TRUE,
    interval_width = "2sigma",
    bespoke_probability = NA,
    denscale = 3,
    resolution = 1,
    n_burn = NA,
    n_end = NA) {

  ##############################################################################
  # Check input parameters

  arg_check <- .InitializeErrorList()

  # Treat single output data as a list of length 1
  if (!is.null(output_data$update_type)) output_data <- list(output_data)

  .CheckMultipleOutputDataConsistent(arg_check, output_data)
  num_data <- length(output_data)
  for (i in 1:num_data) {
    .CheckOutputData(arg_check, output_data[[i]], c("Polya Urn", "Walker"))
    .CheckCalibrationCurveFromOutput(
      arg_check, output_data[[i]], calibration_curve)
    if (is.null(output_data[[i]]$label)) {
      output_data[[i]]$label <- output_data[[i]]$update_type
    }
  }
  .CheckFlag(arg_check, plot_14C_age)
  .CheckChoice(arg_check, plot_cal_age_scale, c("BP", "AD", "BC"))
  .CheckInteger(arg_check, n_posterior_samples, lower = 10)
  .CheckIntervalWidth(arg_check, interval_width, bespoke_probability)
  .CheckNumber(arg_check, denscale, lower = 0)
  .CheckNumber(arg_check, resolution, lower = 0.01)
  .ReportErrors(arg_check)

  # Ensure revert to main environment par on exit of function
  opar <- graphics::par()[c("mgp", "xaxs", "yaxs", "mar", "las")]
  on.exit(graphics::par(opar))

  if (is.null(calibration_curve)) {
    calibration_curve <- get(output_data[[1]]$input_data$calibration_curve_name)
  }
  rc_determinations <- output_data[[1]]$input_data$rc_determinations
  rc_sigmas <- output_data[[1]]$input_data$rc_sigmas
  F14C_inputs <-output_data[[1]]$input_data$F14C_inputs

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

  ##############################################################################
  # Initialise plotting parameters

  calibration_curve_colour <- "blue"
  calibration_curve_bg <- grDevices::rgb(0, 0, 1, .3)
  output_colours <- c("purple", "darkgreen", "darkorange2", "deeppink3")
  if (num_data > 4) {
    output_colours <- c(output_colours, grDevices::hcl.colors(n=num_data-4, "Roma"))
  }

  calendar_age_sequence <- .CreateXRangeToPlotDensity(output_data[[1]], resolution)
  ##############################################################################
  # Calculate density distributions

  if (show_SPD){
    SPD <- FindSummedProbabilityDistribution(
      calendar_age_range_BP = floor(range(calendar_age_sequence)),
      rc_determinations = rc_determinations,
      rc_sigmas = rc_sigmas,
      F14C_inputs = !plot_14C_age,
      resolution = resolution,
      calibration_curve = calibration_curve)
  }

  predictive_density <- list()
  for (i in 1:num_data) {
    predictive_density[[i]] <- FindPredictiveCalendarAgeDensity(
      output_data[[i]],
      calendar_age_sequence,
      n_posterior_samples,
      interval_width,
      bespoke_probability,
      n_burn,
      n_end)
  }
  ##############################################################################
  # Calculate plot scaling
  xlim <- rev(range(calendar_age_sequence))
  ylim_density <- c(0, denscale * max(predictive_density[[1]]$density_mean))

  ##############################################################################
  # Plot curves

  .PlotCalibrationCurveAndInputData(
    plot_cal_age_scale,
    xlim,
    calibration_curve,
    rc_determinations,
    plot_14C_age,
    calibration_curve_colour,
    calibration_curve_bg,
    interval_width,
    bespoke_probability)

  .SetUpDensityPlot(plot_cal_age_scale, xlim, ylim_density)

  if (show_SPD) {
    .PlotSPDEstimateOnCurrentPlot(plot_cal_age_scale, SPD$calendar_age_BP, SPD$probability)
  }

  for (i in 1:num_data) {
    .PlotDensityEstimateOnCurrentPlot(
      plot_cal_age_scale, predictive_density[[i]], output_colours[[i]], show_confidence_intervals)
  }

  .AddLegendToDensityPlot(
    output_data,
    show_SPD,
    show_confidence_intervals,
    interval_width,
    bespoke_probability,
    calibration_curve_colour,
    output_colours)

  invisible(predictive_density)
}


.CreateXRangeToPlotDensity <- function(output_data, resolution) {
  start_age <- min(output_data$calendar_ages, output_data$density_data$calendar_ages[1])
  start_age <- floor(start_age / resolution) * resolution

  end_age <- max(output_data$calendar_ages, output_data$density_data$calendar_ages[2])
  end_age <- ceiling(end_age / resolution) * resolution

  calendar_age_sequence <- seq(start_age, end_age, by = resolution)
  return(calendar_age_sequence)
}


.PlotDensityEstimateOnCurrentPlot <- function(
    plot_cal_age_scale, predictive_density, output_colour, show_confidence_intervals) {

  cal_age <- .ConvertCalendarAge(plot_cal_age_scale, predictive_density$calendar_age)

  graphics::lines(cal_age, predictive_density$density_mean, col = output_colour)
  if (show_confidence_intervals) {
    graphics::lines(cal_age, predictive_density$density_ci_lower, col = output_colour, lty = 2)
    graphics::lines(cal_age, predictive_density$density_ci_upper, col = output_colour, lty = 2)
  }
}


.AddLegendToDensityPlot <- function(
    output_data,
    show_SPD,
    show_confidence_intervals,
    interval_width,
    bespoke_probability,
    calibration_curve_colour,
    output_colours) {
  SPD_colour <- grDevices::grey(0.1, alpha = 0.3)

  ci_label <- switch(
    interval_width,
    "1sigma" = expression(paste("1", sigma, " interval")),
    "2sigma"  = expression(paste("2", sigma, " interval")),
    "bespoke" = paste0(round(100 * bespoke_probability), "% interval"))

  legend_labels <- c(
    gsub("intcal", "IntCal",
         gsub("shcal", "SHCal",
              output_data[[1]]$input_data$calibration_curve_name)), # Both IntCal and SHCal
    ci_label)
  lty <- c(1, 2)
  pch <- c(NA, NA)
  col <- c(calibration_curve_colour, calibration_curve_colour)

  for (i in seq_along(output_data)) {
    legend_labels <- c(legend_labels, output_data[[i]]$label)
    lty <- c(lty, 1)
    pch <- c(pch, NA)
    col <- c(col, output_colours[[i]])

    if (show_confidence_intervals) {
      legend_labels <- c(legend_labels, ci_label)
      lty <- c(lty, 2)
      pch <- c(pch, NA)
      col <- c(col, output_colours[[i]])
    }
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
