# A simple example of how to overlay on the DPMM plot (also works for Poisson process)
library(carbondate)

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
par(opar)

