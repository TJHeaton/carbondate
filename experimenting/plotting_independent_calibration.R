CalibrateSingleDetermination(
  rc_determination = rc_determination,
  rc_sigma = rc_sigma,
  calibration_curve = shcal20,
  F14C_inputs = TRUE,
  plot_output = TRUE,
  interval_width = "bespoke",
  bespoke_probability = 0.8)

FindSummedProbabilityDistribution(
  calendar_age_range_BP=c(600, 1600),
  rc_determinations=c(602, 805, 1554),
  rc_sigmas=c(35, 34, 45),
  calibration_curve=intcal20,
  plot_output = TRUE,
  interval_width = "bespoke",
  bespoke_probability = 0.8)

# Example for vignette to show SPD
FindSummedProbabilityDistribution(
  calendar_age_range_BP=c(2500, 7000),
  rc_determinations= two_normals$c14_age,
  rc_sigmas = two_normals$c14_sig,
  calibration_curve=intcal20,
  plot_output = TRUE)

# Show true underlying calendar age density
weights_true <- c(0.45, 0.55)
cluster_means_true_calBP <- c(3500, 5000)
cluster_precisions_true <- 1 / c(200, 100)^2

# Find and plot true exact density
truedens <- function(t, w, truemean, trueprec) {
  dens <- 0
  for(i in 1:length(w)) {
    dens <- dens + w[i] * dnorm(t, mean = truemean[i], sd = 1/sqrt(trueprec[i]))
  }
  dens
}
curve(truedens(
  x,
  w = weights_true,
  truemean = cluster_means_true_calBP,
  trueprec = cluster_precisions_true),
      from = 2500, to = 7000, n = 401,
      col = "red", add = TRUE)


FindSummedProbabilityDistribution(
  calendar_age_range_BP=c(2500, 7000),
  rc_determinations= two_normals$c14_age,
  rc_sigmas = two_normals$c14_sig,
  calibration_curve=intcal20,
  plot_output = TRUE)


c14_age <- c(602, 805, 1554)
c14_sig <- c(35, 34, 45)
(f14c <- exp(-c14_age / 8033))
(f14c_sig <- f14c * c14_sig / 8033)

########


# An different example using F14C concentrations and the IntCal 13 curve
SPD <- FindSummedProbabilityDistribution(
  calendar_age_range_BP=c(400, 2100),
  rc_determinations=c(0.8, 0.85, 0.9),
  rc_sigmas=c(0.01, 0.015, 0.012),
  F14C_inputs=TRUE,
  calibration_curve=intcal13,
  plot_output = TRUE)
plot(SPD, type = "l",
     xlim = rev(range(SPD$calendar_age_BP)),
     xlab = "Calendar Age (cal yr BP)")


