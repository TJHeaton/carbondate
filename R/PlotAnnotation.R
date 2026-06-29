#' Add Text Annotation to Various Summary Plots
#'
#' @description
#' Having plotted the output of the various modelling approaches, add text annotation at a given location on the plot. It is based on the base [text] function
#'
#' This function/annotation can be applied after having called any of [carbondate::PlotPosteriorMean]
#'
#' @inheritParams graphics::text
#'
#' @param output_plot The plot onto which you wish to add text
#'
#' @return None
#'
#' @export
#'
#' @examples
#' # NOTE: This example is shown with a small n_iter and n_posterior_samples
#' # to speed up execution.
#' # Try n_iter and n_posterior_samples as the function defaults.
#'
#' pp_output <- PPcalibrate(
#'     pp_uniform_phase$c14_age,
#'     pp_uniform_phase$c14_sig,
#'     intcal20,
#'     n_iter = 1000,
#'     show_progress = FALSE)
#'
#' # Default plot with 2 sigma interval
#' posterior_mean_plot <- PlotPosteriorMeanRate(
#'     pp_output,
#'     n_posterior_samples = 100)
#'
#' # Add text to plot
#' AddTextPlot(posterior_mean_plot,
#'     x = 600, y = 500,
#'     labels = expression(paste("600", ""^14, "C ", "yrs BP")),
#'     cex = 0.7,
#'     col = "red")
#'
AddTextPlot <- function(
    output_plot,
    x,
    y,
    labels,
    adj = NULL,
    pos = NULL,
    offset = 0.5,
    vfont = NULL,
    cex = 1,
    col = NULL,
    font = NULL) {

  # Ensure revert to main environment par on exit of function
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  graphics::par(output_plot$plot_par)
  text(x, y, labels, adj, pos, offset, vfont, cex, col, font)
}


#' Add Shading to Various Summary Plots
#'
#' @description
#' Having plotted the output of the various modelling approaches, add a shaded interval or box to the plot. This function/annotation can be applied after having called any of [carbondate::PlotPosteriorMean]
#'
#' @param output_plot The plot onto which you wish to add text
#' @param x_start,x_end The values on the x-axis (calendar age) at which you want the shading to start and end
#' @param y_start,y_end `NULL` if you want the shading to cover the entire y-axis range. Otherwise (if you want to plot a box) the y-axis values where you want the sahding to start and end
#' @param col The color to be used
#' @param alpha The level of transparency for the shading (`0 < alpha < 1`) to see the plot underneath. `0` means fully transparent and `1` means fully opaque.
#'
#' @return None
#'
#' @export
#'
#' @examples
#' # NOTE: This example is shown with a small n_iter and n_posterior_samples
#' # to speed up execution.
#' # Try n_iter and n_posterior_samples as the function defaults.
#'
#' pp_output <- PPcalibrate(
#'     pp_uniform_phase$c14_age,
#'     pp_uniform_phase$c14_sig,
#'     intcal20,
#'     n_iter = 1000,
#'     show_progress = FALSE)
#'
#' # Default plot with 2 sigma interval
#' posterior_mean_plot <- PlotPosteriorMeanRate(
#'     pp_output,
#'     n_posterior_samples = 100)
#'
#' # Add transparent red shaded period to plot
#' AddShadingPlot(posterior_mean_plot,
#'     x_start = 620, x_end = 600,
#'     col = "red")
#'
#' # Add green (more opaque) shaded box to plot
#' AddShadingPlot(posterior_mean_plot,
#'     x_start = 590, x_end = 550,
#'     y_start = 500, y_end = 410,
#'     col = "green", alpha = 0.8)
#'
AddShadingPlot <- function(
    output_plot,
    x_start,
    x_end,
    y_start = NULL,
    y_end = NULL,
    col, alpha = 0.4) {

  # Ensure revert to main environment par on exit of function
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  add.alpha <- function(cols, alpha) {
    grDevices::rgb(t(grDevices::col2rgb(cols)/255),
                   alpha = alpha)
  }

  shading_col <- add.alpha(col, alpha = alpha)

  graphics::par(output_plot$plot_par)
  if(is.null(y_start)) y_start <- output_plot$plot_par$usr[3]
  if(is.null(y_end)) y_end <- output_plot$plot_par$usr[4]

  graphics::polygon(x = c(rep(x_start, 2),
                          rep(x_end, 2)),
                    y = c(y_start,
                          y_end,
                          y_end,
                          y_start),
                    border = NA,
                    col = shading_col)
}


#' Add Straight Lines to Various Summary Plots
#'
#' @description
#' Having plotted the output of the various modelling approaches, add one or more straight lines to the current plot. It is based on the base [abline] function
#'
#'
#' This function/annotation can be applied after having called any of [carbondate::PlotPosteriorMean]
#'
#' @inheritParams graphics::abline
#'
#' @param output_plot The plot onto which you wish to add text
#' @param reg An object with a [coef] method
#'
#' @return None
#'
#' @export
#'
#' @examples
#' # NOTE: This example is shown with a small n_iter and n_posterior_samples
#' # to speed up execution.
#' # Try n_iter and n_posterior_samples as the function defaults.
#'
#' pp_output <- PPcalibrate(
#'     pp_uniform_phase$c14_age,
#'     pp_uniform_phase$c14_sig,
#'     intcal20,
#'     n_iter = 1000,
#'     show_progress = FALSE)
#'
#' # Default plot with 2 sigma interval
#' posterior_mean_plot <- PlotPosteriorMeanRate(
#'     pp_output,
#'     n_posterior_samples = 100)
#'
#' # Add vertical thick red dashed line to plot
#' AddLinePlot(
#'      posterior_mean_plot,
#'      v = 600,
#'      col = "red",
#'      lwd = 2,
#'      lty = 3)
#'
#' # Add narrow horizontal green solid line
#' AddLinePlot(
#'      posterior_mean_plot,
#'      h = 500,
#'      col = "green",
#'      lwd = 1,
#'      lty = 1)
#'
#' # Add light gray grid lines
#' AddLinePlot(posterior_mean_plot,
#'      #' AddLinePlot(
#'      posterior_mean_plot,
#'      h = seq(250, 700, by = 25),
#'      v = seq(400, 650, by = 25),
#'      col = "lightgray",
#'      lty = 3)
AddLinePlot <- function(
    output_plot,
    a = NULL,
    b = NULL,
    h = NULL,
    v = NULL,
    reg = NULL,
    coef = NULL,
    ...) {

  # Ensure revert to main environment par on exit of function
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  graphics::par(output_plot$plot_par)
  abline(a = a, b = b, h = h, v = v, reg = reg, coef = coef, untf = FALSE, ...)
}


