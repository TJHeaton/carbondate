library(usethis)

set.seed(15)

n_observed <- 40
rc_sigmas <- rep(15, n_observed)

# Create artificial rc_determinations
observed_age_range <- c(500, 550)
true_theta <- seq(from = observed_age_range[1], to = observed_age_range[2], length = n_observed)
intcal_mean <- approx(x = intcal20$calendar_age_BP, y = intcal20$c14_age, xout = true_theta)$y
intcal_sd <- approx(x = intcal20$calendar_age_BP, y = intcal20$c14_sig, xout = true_theta)$y
rc_determinations <- rnorm(n = n_observed, mean = intcal_mean, sd = sqrt(rc_sigmas^2 + intcal_sd^2))

pp_uniform_phase <- data.frame(c14_age = rc_determinations, c14_sig = rc_sigmas)
pp_uniform_phase$f14c <- exp(-pp_uniform_phase$c14_age / 8033)
pp_uniform_phase$f14c_sig <- pp_uniform_phase$f14c * pp_uniform_phase$c14_sig / 8033

usethis::use_data(pp_uniform_phase, overwrite = TRUE)

