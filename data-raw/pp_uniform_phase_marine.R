library(usethis)

set.seed(12)

n_observed <- 40
rc_sigmas <- rep(15, n_observed)
offset <- 100
offset_sig <- 20

# Create samples

delta_r <- rep(offset, n_observed)
delta_r_sig <- rep(offset_sig, n_observed)

# Create artificial rc_determinations
observed_age_range <- c(500, 550)
true_theta <- seq(from = observed_age_range[1],
                  to = observed_age_range[2],
                  length = n_observed)


marine_mean <- approx(x = marine20$calendar_age_BP,
                      y = marine20$c14_age,
                      xout = true_theta)$y
marine_sd <- approx(x = marine20$calendar_age_BP,
                    y = marine20$c14_sig,
                    xout = true_theta)$y

# Create determination and delta_r adjustment
rc_determinations <- rnorm(n = n_observed,
                           mean = marine_mean,
                           sd = sqrt(rc_sigmas^2 + marine_sd^2))
delta_r_real <- rnorm(n = n_observed,
                      mean = delta_r,
                      sd = delta_r_sig)

rc_determinations <- rc_determinations + delta_r_real

pp_uniform_phase_marine <- data.frame(c14_age = rc_determinations,
                                             c14_sig = rc_sigmas,
                                             sample_source = "marine20",
                                             delta_r = delta_r,
                                             delta_r_sig = delta_r_sig)

pp_uniform_phase_marine$f14c <- exp(-pp_uniform_phase_marine$c14_age / 8033)
pp_uniform_phase_marine$f14c_sig <- pp_uniform_phase_marine$f14c * pp_uniform_phase_marine$c14_sig / 8033

usethis::use_data(pp_uniform_phase_marine, overwrite = TRUE)

