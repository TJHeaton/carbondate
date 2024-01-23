
hobs22 <- read.csv("experimenting/post-bomb/HOBS2022.14c", header = FALSE, skip=2, sep="")[, c(1, 4, 5)]
names(hobs22) <- c(
  "calendar_age_BP",
  "f14c",
  "f14c_sig")

# convert to calBP
hobs22$calendar_age_BP <- 1950 - hobs22$calendar_age_BP

X8836_modern <- read.csv("experimenting/post-bomb/8836_modern.csv")[, c(2, 3)]

FindSummedProbabilityDistribution(
  c(-1, -69), X8836_modern$F14C, X8836_modern$F14C_sd, hobs22, F14C_inputs = TRUE, plot_output = TRUE)

for (i in seq_along(X8836_modern$F14C)) {
  CalibrateSingleDetermination(
    X8836_modern$F14C[i],
    X8836_modern$F14C_sd[i],
    hobs22,
    resolution = 0.1,
    F14C_inputs = TRUE,
    plot_output = TRUE)
}

polya_urn_output <- PolyaUrnBivarDirichlet(
  X8836_modern$F14C, X8836_modern$F14C_sd, hobs22, F14C_inputs = TRUE)

PlotPredictiveCalendarAgeDensity(polya_urn_output, show_SPD = TRUE, plot_14C_age = FALSE, resolution = 0.1)

for (i in seq_along(X8836_modern$F14C)) {
  PlotCalendarAgeDensityIndividualSample(
    i,
    polya_urn_output,
    plot_14C_age = FALSE,
    hist_resolution = 0.25,
    density_resolution = 0.1,
    show_hpd_ranges = TRUE,
    show_unmodelled_density = TRUE)
}
