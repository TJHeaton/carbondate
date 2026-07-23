# Test code to work out best way to write text onto plots

add_text_plot <- function(
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

add_text_plot(Temp,
              x = 550,
              y = 400,
              labels = expression(bar(x) == sum(frac(x[i], n), i==1, n)))




add_shading_plot <- function(
    output_plot,
    x_start,
    x_end,
    y_start = NULL,
    y_end = NULL,
    col, alpha = 0.4) {

  # Ensure revert to main environment par on exit of function
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  add.alpha <- function(cols, alpha) grDevices::rgb(t(grDevices::col2rgb(cols)/255),
                                                    alpha = alpha)
  shading_col <- add.alpha(col, alpha = 0.4)

  graphics::par(output_plot$plot_par)
  if(is.null(y_axis_start)) y_axis_start <- output_plot$plot_par$usr[3]
  if(is.null(y_axis_end)) y_axis_end <- output_plot$plot_par$usr[4]

  cat(y_axis_start, "\n")
  cat(y_axis_end, "\n")

  graphics::polygon(x = c(rep(time_start, 2), rep(time_end, 2)),
          y = c(y_axis_start, y_axis_end, y_axis_end, y_axis_start),
          border = NA,
          col = shading_col)
}

add_shading_plot(Temp, time_start = 600, time_end = 500, col = "red")


polygon(x = c(rep(500, 2), rep(600, 2)),
        y = c(700, 200), 2),
        border = NA,
        col = "red")


  # Add shading wherever you want
  period_start <- 4500
  period_end <- 4300

  # Decide on shading colour (and make somewhat transparent)
  shade_cols <- "red"
  add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)
  cols <- add.alpha(shade_cols, alpha = 0.4)

  polygon(x = c(rep(period_start, 2), rep(period_end, 2)),
          y = c(0, 0.82, 0.82, 0), border = NA,
          col = cols)
}



add_text_plot <- function(
    output_data,
    text,
    cal_age,
    height) {

  # Ensure revert to main environment par on exit of function
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  graphic

  # Set nice plotting parameters
  if(plot_pretty) {
    graphics::par(
      new = TRUE,
      mgp = c(3, 0.7, 0),
      xaxs = "i",
      yaxs = "i",
      mar = c(5, 4.5, 4, 2) + 0.1,
      las = 1)
  }

  xlim <- rev(range(two_normals_DPMM[[1]]$calendar_age_BP))
  ylim <- c(0, 3 * max(two_normals_DPMM[[1]]$density_mean))

}


# Store the old plotting parameters so you can change them back at the end
oldpar <- par(no.readonly = TRUE)

# Run a Polya Urn DPMM model (also works for Poisson process)

polya_urn_output <- PolyaUrnBivarDirichlet(
  rc_determinations = two_normals$c14_age,
  rc_sigmas = two_normals$c14_sig,
  calibration_curve=intcal20,
  n_iter = 1e4,
  show_progress = FALSE)

two_normals_DPMM <- PlotPredictiveCalendarAgeDensity(
  output_data = polya_urn_output,
  show_SPD = TRUE)

# Adjust the plotting parameters to those used internally in the DPMM plotting function
# Note: the ylim corresponds to the density (not the 14C age)
par(new = TRUE,
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 4, 2) + 0.1,
    las = 1)
xlim <- rev(range(two_normals_DPMM[[1]]$calendar_age_BP))
ylim <- c(0, 3 * max(two_normals_DPMM[[1]]$density_mean))

# Create an empty plot
plot(
  NULL,
  NULL,
  type = "n",
  ylim = ylim,
  xlim = xlim,
  axes = FALSE,
  xlab = NA,
  ylab = NA,
  xaxs = "i",
  yaxs = "i")

# Add calendar age ticks wherever you want
xticks_minor <- seq(2000, 6000, by = 100)
axis(1, xticks_minor, labels = FALSE, tcl = -0.2)

# Add shading wherever you want
period_start <- 4500
period_end <- 4300

# Decide on shading colour (and make somewhat transparent)
shade_cols <- "red"
add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)
cols <- add.alpha(shade_cols, alpha = 0.4)

polygon(x = c(rep(period_start, 2), rep(period_end, 2)),
        y = c(0, 0.82, 0.82, 0), border = NA,
        col = cols)

# Place 1/5 up image by using y = 0.2 * max(ylim)
text("Shaded Period", x = period_end, y = 0.2 * max(ylim), pos = 4)

# Reset/Change back to old plotting parameters
par(oldpar)
