library(usethis)

set.seed(9)


# Create some observed data from the clusters according to probabilities
num_observations <- 50


# Create vector of delta_r
offset <- -50
offset_sig <- 30
delta_r <- rep(offset, num_observations)
delta_r_sig <- rep(offset_sig, num_observations)

# Choose base function of calendar age
weights_true <- c(0.45, 0.55)
cluster_means_true_calBP <- c(3500, 5000)
cluster_precisions_true <- 1 / c(200, 100)^2
n_clust_true <- length(weights_true)

true_density <- data.frame(x = marine20$calendar_age_BP, y = 0)
true_density$y <- 0
for(i in seq_along(weights_true)) {
  true_density$y <- true_density$y +
    weights_true[i] * dnorm(
      true_density$x, mean = cluster_means_true_calBP[i], sd = 1/sqrt(cluster_precisions_true[i]))
}

cluster_identifiers_true <- sample(
  1:n_clust_true, num_observations, replace = TRUE, prob = weights_true)
calendar_ages_true <- rnorm(
  num_observations,
  mean = cluster_means_true_calBP[cluster_identifiers_true],
  sd = 1 / sqrt(cluster_precisions_true[cluster_identifiers_true]))
hist(calendar_ages_true, breaks = 20, freq = FALSE)
lines(true_density, col="red")

# Create some radiocarbon determinations
calibration_curve_at_true_calendar_ages <- InterpolateCalibrationCurve(
  new_calendar_ages = calendar_ages_true, calibration_curve = marine20)
c14_age_for_true_calendar_ages <- calibration_curve_at_true_calendar_ages$c14_age
c14_sig_for_true_calendar_ages <- calibration_curve_at_true_calendar_ages$c14_sig

# Sample from the true C14 age and C14 sig values
xcalcurve <- rnorm(num_observations, c14_age_for_true_calendar_ages, c14_sig_for_true_calendar_ages)
c14_sig <- rep(25, num_observations)
c14_age <- rnorm(num_observations, mean = xcalcurve, sd = c14_sig)

# Create delta_r realisations
delta_r_real <- rnorm(n = num_observations,
                      mean = delta_r,
                      sd = delta_r_sig)

c14_age <- c14_age + delta_r_real

two_normals_marine <- data.frame(
  c14_age = c14_age,
  c14_sig = c14_sig,
  sample_source = "marine20",
  delta_r = delta_r,
  delta_r_sig = delta_r_sig)

two_normals_marine$f14c <- exp(-two_normals_marine$c14_age / 8033)
two_normals_marine$f14c_sig <- two_normals_marine$f14c * two_normals_marine$c14_sig / 8033

usethis::use_data(two_normals_marine, overwrite = TRUE)
