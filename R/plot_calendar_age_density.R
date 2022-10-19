#' Plots the predicted calendar age density from the output data.
#'
#' Plots the input radiocarbon determinations and calibration curve, with the
#' output predicted density on the same plot. Can also optionally show the
#' SPD estimate.
#'
#' @param c14_determinations A vector containing the radiocarbon determinations.
#' @param c14_uncertainties A vector containing the radiocarbon determination
#' uncertainties. Must be the same length as `c14_determinations`.
#' @param calibration_curve A dataframe which should contain one column entitled
#' c14_age and one column entitled c14_sig.
#' This format matches [carbondate::intcal20].
#' @param output_data Data returned from one of the updating functions e.g.
#' [carbondate::WalkerBivarDirichlet] or
#' [carbondate::BivarGibbsDirichletwithSlice].
#' @param n_posterior_samples Current number of samples it will draw from this
#' posterior to estimate the calendar age density (possibly repeats).
#' @param lambda,nu1,nu2  Hyperparameters used in the updating function for the
#' prior on the means \eqn{\phi_j} and precision \eqn{\tau_j} of each individual
#' calendar age cluster \eqn{j}.
#' @param show_SPD Whether to calculate and show the summed probability function
#' on the plot (optional). Default is `TRUE`.
#' @param xlimscal,ylimscal Whether to scale the x or y limits (optional).
#' Default is 1. Values more than 1 will increase range of the limits,
#' values less than 1 will decrease the range of the limits.
#' @param denscale Whether to scale the vertical range of the density plot
#' relative to the calibration curve plot (optional). Default is 3 which means
#' that the maximum SPD density will be at 1/3 of the height of the plot.
#' TODO: base on the output data.
#'
#' @return A dataframe containing the `calendar_age` and the `density_estimate`
#' from the output data.
#'
#' @export
PlotCalendarAgeDensity <- function(
    c14_determinations,
    c14_uncertainties,
    calibration_curve,
    output_data,
    n_posterior_samples,
    lambda,
    nu1,
    nu2,
    show_SPD = TRUE,
    xlimscal = 1,
    ylimscal = 1,
    denscale = 3) {

  # TODO: Check inputs

  SPD_colour <- grDevices::grey(0.1, alpha = 0.5)
  calibration_curve_colour <- "blue"
  calibration_curve_bg <- grDevices::rgb(0, 0, 1, .3)
  output_colour <- "purple"

  calendar_age_sequence <- .CreateRangeToPlotDensity(output_data)

  xlim <- .ScaleLimit(rev(range(calendar_age_sequence)), xlimscal)
  ylim <- .ScaleLimit(
    range(c14_determinations) +
      c(-2, 2) * stats::quantile(c14_uncertainties, 0.9),
    ylimscal)

  .PlotCalibrationCurveAndInputData(
    xlim,
    ylim,
    calibration_curve,
    c14_determinations,
    calibration_curve_colour,
    calibration_curve_bg)

  if (show_SPD){
    SPD = FindSPD(
      calendar_age_range = floor(range(output_data$calendar_ages)),
      c14_determinations = c14_determinations,
      c14_uncertainties = c14_uncertainties,
      calibration_curve = calibration_curve)
    .PlotSPDEstimateOnCurrentPlot(SPD, SPD_colour, denscale, xlim)
  }

  posterior_density_mean = .PlotDensityEstimateOnCurrentPlot(
    output_data,
    output_colour,
    calendar_age_sequence,
    n_posterior_samples,
    lambda,
    nu1,
    nu2)

  .AddLegendToDensityPlot(
    output_data, show_SPD, calibration_curve_colour, output_colour, SPD_colour)

  invisible(
    data.frame(
      calendar_age=calendar_age_sequence,
      density_estimate=posterior_density_mean))
}


.CreateRangeToPlotDensity <- function(output_data) {
  calendar_age_sequence <- seq(
    floor(min(output_data$calendar_ages, na.rm = TRUE)),
    ceiling(max(output_data$calendar_ages, na.rm = TRUE)),
    by = 1,
  )
  return(calendar_age_sequence)
}


.ScaleLimit <- function(lim, limscal) {
  lim <- lim + c(1, -1) * diff(lim) * (1 - limscal)
  return(lim)
}


.PlotCalibrationCurveAndInputData <- function(
    xlim,
    ylim,
    calibration_curve,
    c14_determinations,
    calibration_curve_colour,
    calibration_curve_bg){
  graphics::par(mar = c(5, 4.5, 4, 2) + 0.1, las = 1)
  graphics::plot.default(
    calibration_curve$calendar_age,
    calibration_curve$c14_age,
    col = calibration_curve_colour,
    ylim = ylim,
    xlim = xlim,
    xlab = "Calendar Age (cal yr BP)",
    ylab = expression(paste(""^14, "C", " age (", ""^14, "C yr BP)")),
    type = "l",
    main = expression(paste(""^14, "C Calibration")),
  )
  calibration_curve$ub <- calibration_curve$c14_age +
    1.96 * calibration_curve$c14_sig
  calibration_curve$lb <- calibration_curve$c14_age -
    1.96 * calibration_curve$c14_sig

  graphics::lines(
    calibration_curve$calendar_age,
    calibration_curve$ub,
    lty = 2,
    col = calibration_curve_colour)
  graphics::lines(
    calibration_curve$calendar_age,
    calibration_curve$lb,
    lty = 2,
    col = calibration_curve_colour)
  graphics::polygon(
    c(rev(calibration_curve$calendar_age), calibration_curve$calendar_age),
    c(rev(calibration_curve$lb), calibration_curve$ub),
    col = calibration_curve_bg,
    border = NA)
  graphics::rug(c14_determinations, side = 2)
}


.PlotSPDEstimateOnCurrentPlot <- function(SPD, SPD_colour, denscale, xlim) {
  graphics::par(new = TRUE)
  graphics::plot.default(
    SPD$calendar_age,
    SPD$probability,
    lty = 1,
    col = SPD_colour,
    type = "n",
    ylim = c(0, denscale * max(SPD$prob)),
    xlim = xlim,
    axes = FALSE,
    xlab = NA,
    ylab = NA)
  graphics::polygon(
    c(SPD$calendar_age, rev(SPD$calendar_age)),
    c(SPD$probability, rep(0, length(SPD$probability))),
    border = NA,
    col = SPD_colour)
}


.PlotDensityEstimateOnCurrentPlot <- function(
    output_data,
    output_colour,
    calendar_age_sequence,
    n_posterior_samples,
    lambda,
    nu1,
    donu2) {

  # each column is the density for a particular sample id
  posterier_density_matrix <- .FindDensityPerSampleID(
    output_data, calendar_age_sequence, n_posterior_samples, lambda, nu1, nu2)

  posterier_density_confidence_intervals <- apply(
    posterier_density_matrix, 1, stats::quantile, probs = c(0.025, 0.975))
  posterier_density_mean <- apply(posterier_density_matrix, 1, mean)

  graphics::lines(
    calendar_age_sequence, posterier_density_mean, col = output_colour)
  graphics::lines(
    calendar_age_sequence,
    posterier_density_confidence_intervals[1, ],
    col = output_colour,
    lty = 2)
  graphics::lines(
    calendar_age_sequence,
    posterier_density_confidence_intervals[2, ],
    col = output_colour,
    lty = 2)

  return(posterier_density_mean)
}


.FindDensityPerSampleID <- function(
    output_data, calendar_age_sequence, n_posterior_samples, lambda, nu1, nu2) {
  n_out <- length(output_data$alpha)
  n_burn <- floor(n_out / 2)

  posterior_sample_ids <- sample(
    x = n_burn:n_out,
    size = n_posterior_samples,
    replace = n_posterior_samples > (n_out - n_burn))

  # Create a matrix where each column is the density for a particular sample id
  posterior_density_matrix <- apply(
    matrix(posterior_sample_ids, 1, n_posterior_samples),
    2,
    function(i, output_data, x, lambda, nu1, nu2) {
      if (output_data$update_type == "walker") {
        .FindPredictedDensityWalker(
          x,
          weight = output_data$weight[[i]],
          phi = output_data$phi[[i]],
          tau = output_data$tau[[i]],
          mu_phi = output_data$mu_phi[i],
          lambda = lambda,
          nu1 = nu1,
          nu2 = nu2)
      } else {
        .FindPredictedDensityNeal(
          x,
          cluster_identifiers = output_data$cluster_identifiers[i, ],
          phi = output_data$phi[[i]],
          tau = output_data$tau[[i]],
          alpha = output_data$alpha[i],
          mu_phi = output_data$mu_phi[i],
          lambda = lambda,
          nu1 = nu1,
          nu2 = nu2)
      }
    },
    output_data = output_data,
    x = calendar_age_sequence,
    lambda = lambda,
    nu1 = nu1,
    nu2 = nu2)
  return(posterior_density_matrix)
}

.AddLegendToDensityPlot <- function(
    output_data,
    show_SPD,
    calibration_curve_colour,
    output_colour,
    SPD_colour) {
  legend_labels = "IntCal20"
  lty = 1
  pch = NA
  col = calibration_curve_colour
  if (output_data$update_type == "walker") {
    legend_labels <- c(legend_labels, "Walker DP", "Walker 95% prob interval")
  } else if (output_data$update_type == "neal") {
    legend_labels <- c(legend_labels, "Neal DP", "Neal 95% prob interval")
  } else {
    stop("Unknown update type")
  }
    lty <- c(lty, 1, 2)
    pch <- c(pch, NA, NA)
    col <- c(col, output_colour, output_colour)

  if (show_SPD) {
    legend_labels <- c(legend_labels, "SPD Estimate")
    lty <- c(lty, -1)
    pch <- c(pch, 15)
    col <- c(col, SPD_colour)
  }

  graphics::legend(
    "topright", legend = legend_labels, lty = lty, pch = pch, col = col)
}
