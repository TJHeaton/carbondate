#' Plots the posterior calendar ages for an individual determination
#'
#' Once a function has been run to calibrate a set of radiocarbon
#' determinations, the posterior density for a single determination can be
#' plotted using this function.
#'
#' @param ident the determination you want to show the individual posterior
#' calendar age for.
#' @param c14_determinations A vector containing the radiocarbon determinations.
#' @param c14_uncertainties A vector containing the radiocarbon determination
#' uncertainties. Must be the same length as `c14_determinations`.
#' @param calibration_curve A dataframe which should contain one column entitled
#' c14_age and one column entitled c14_sig.
#' This format matches [carbondate::intcal20].
#' @param output_data Data returned from one of the updating functions e.g.
#' [carbondate::WalkerBivarDirichlet] or
#' [carbondate::BivarGibbsDirichletwithSlice].
#' @param n_breaks The number of breaks in the histogram (optional). If not
#' given then it is set to the max of 100 and 1/tenth of post burn-in chain.
#'
#' @return No return value
#' @export
#'
#' @examples
#' # First generate some output data
#' c14_determinations = c(602, 805, 954)
#' c14_uncertainties = c(35, 34, 45)
#' walker_temp = WalkerBivarDirichlet(
#'   c14_determinations = c14_determinations,
#'   c14_uncertainties = c14_uncertainties,
#'   calibration_curve = intcal20,
#'   lambda = 0.1,
#'   nu1 = 0.25,
#'   nu2 = 10,
#'   alpha_shape = 1,
#'   alpha_rate = 1)
#'
#' # Plot results for the second determinations
#' PlotIndividualCalendarAgeDensity(
#'   ident = 2,
#'   c14_determinations = c14_determinations,
#'   c14_uncertainties = c14_uncertainties,
#'   calibration_curve = intcal20,
#'   output_data = walker_temp,
#'   n_breaks = 10)
#'
PlotIndividualCalendarAgeDensity <- function(
    ident,
    c14_determinations,
    c14_uncertainties,
    calibration_curve,
    output_data,
    n_breaks = NA) {

  calendar_age <- output_data$calendar_ages[, ident]
  c14_age <- c14_determinations[ident]
  c14_sig <- c14_uncertainties[ident]

  n_out <- length(calendar_age)
  n_burn <- floor(n_out / 2)
  calendar_age <- calendar_age[n_burn:n_out]

  if (is.na(n_breaks)) {
    n_breaks <- min(100, floor(length(calendar_age) / 10))
  }

  # Find the calendar age range to plot
  xrange <- range(calendar_age) + c(-1, 1) * 10
  cal_age_ind_min <- which.min(abs(calibration_curve$calendar_age - xrange[1]))
  cal_age_ind_max <- which.min(abs(calibration_curve$calendar_age - xrange[2]))
  calendar_age_indices <- cal_age_ind_min:cal_age_ind_max

  calibration_curve$ub <- calibration_curve$c14_age +
    1.96 * calibration_curve$c14_sig
  calibration_curve$lb <- calibration_curve$c14_age -
    1.96 * calibration_curve$c14_sig
  yrange <- range(
    calibration_curve$ub[calendar_age_indices],
    calibration_curve$lb[calendar_age_indices] - 10)

  graphics::plot(calibration_curve$calendar_age, calibration_curve$c14_age,
       col = "blue",
       ylim = yrange,
       xlim = rev(xrange),
       xlab = "Calendar Age (cal yr BP)",
       ylab = expression(paste(""^14, "C", " age (", ""^14, "C yr BP)")),
       type = "l",
       main = substitute(
         paste(
           "Posterior of ",
           i^th,
           " determination ",
           c14_age,
           "\u00B1",
           c14_sig,
           ""^14,
           "C yr BP"),
         list(i = ident, c14_age = c14_age, c14_sig = c14_sig)),
       xaxs = "i",
       yaxs = "i")
  graphics::lines(
    calibration_curve$calendar_age, calibration_curve$ub, lty = 2, col = "blue")
  graphics::lines(
    calibration_curve$calendar_age, calibration_curve$lb, lty = 2, col = "blue")

  # Plot the 14C determination on the y-axis
  yfromto <- seq(c14_age - 4 * c14_sig, c14_age + 4 * c14_sig, by = 1)
  radpol <- cbind(
    c(0, stats::dnorm(yfromto, mean = c14_age, sd = c14_sig), 0),
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
  temphist <- graphics::hist(calendar_age, breaks = n_breaks, plot = FALSE)
  graphics::hist(
    calendar_age,
    prob = TRUE,
    breaks = n_breaks,
    xlim = rev(xrange),
    axes = FALSE,
    xlab = NA,
    ylab = NA,
    main = "",
    xaxs = "i",
    ylim = c(0, 2.5 * max(temphist$density)))
}
