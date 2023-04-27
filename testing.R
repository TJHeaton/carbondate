library(devtools)
load_all()
measurements <- read.csv('../custom calcurve/8836_modern.csv', header = TRUE, sep=",")
HOBS2022 <- read.csv('../custom calcurve/HOBS2022.14c', header = FALSE, skip=2, sep="")[, c(1, 4, 5)]
names(HOBS2022) <- c("calendar_age", "c14_age", "c14_sig")
HOBS2022$calendar_age = 1950 - HOBS2022$calendar_age

set.seed(7)

###############################################################################
# Perform the MCMC update

walker_output <- WalkerBivarDirichlet(
  c14_determinations = measurements$F14C,
  c14_sigmas = measurements$F14C_sd,
  calibration_curve=HOBS2022,
  n_iter = 1e5,
  n_thin = 10,
  slice_width = 70,
  n_clust = 8)

pu_output <- PolyaUrnBivarDirichlet(
  c14_determinations = measurements$F14C,
  c14_sigmas = measurements$F14C_sd,
  calibration_curve=HOBS2022,
  n_iter = 1e5,
  n_thin = 10,
  slice_width = 70,
  n_clust = 8)

###############################################################################
# Plot results

# Create a layout with 2/3 showing the predictive density 1/3 showing the number
# of clusters
layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
layout(mat = layout.matrix, heights = c(1), widths = c(10, 4.5))

PlotPredictiveCalendarAgeDensity(
  list(walker_output, pu_output),
  n_posterior_samples = 5000,
  interval_width="2sigma",
  ylimscal = 2,
  denscale = 6,
  show_SPD = TRUE)

PlotNumberOfClusters(walker_output)

# New plot for a single determination
par(mfrow = c(1,1))

PlotCalendarAgeDensityIndividualSample(7, walker_output, resolution = 1)
