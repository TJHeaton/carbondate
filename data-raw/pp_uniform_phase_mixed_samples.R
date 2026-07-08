library(usethis)

set.seed(15)

n_observed <- 40
rc_sigmas <- rep(15, n_observed)
offset <- 100
offset_sig <- 20

# Create samples
curve_year <- "2020"

sample_source <- sample(c("NH", "SH", "Marine"), 40, replace = TRUE)
sample_source <- factor(sample_source,
                        levels = c("NH",
                                   "SH",
                                   "Marine"))

cal_curves <- list(NH = get(paste("intcal", substr(curve_year, 3, 4), sep = "")),
                   SH = get(paste("shcal", substr(curve_year, 3, 4), sep = "")),
                   Marine = get(paste("marine", substr(curve_year, 3, 4), sep = "")))

delta_r <- rep(0, n_observed)
delta_r_sig <- rep(0, n_observed)
delta_r[sample_source == "Marine"] <- offset
delta_r_sig[sample_source == "Marine"] <- offset_sig

# Create artificial rc_determinations
observed_age_range <- c(500, 550)
true_theta <- seq(from = observed_age_range[1], to = observed_age_range[2], length = n_observed)


# Lazy evaluation
cal_curve_mean <- sapply(1:n_observed, function(i) {
  source <- sample_source[i]
  approx(x = cal_curves[[source]]$calendar_age_BP,
         y = cal_curves[[source]]$c14_age,
         xout = true_theta[i])$y
})


cal_curve_sd <- sapply(1:n_observed, function(i) {
  source <- sample_source[i]
  approx(x = cal_curves[[source]]$calendar_age_BP,
         y = cal_curves[[source]]$c14_sig,
         xout = true_theta[i])$y
})

# Create determination and delta_r adjustment
rc_determinations <- rnorm(n = n_observed,
                           mean = cal_curve_mean,
                           sd = sqrt(rc_sigmas^2 + cal_curve_sd^2))
delta_r_real <- rnorm(n = n_observed,
                      mean = delta_r,
                      sd = delta_r_sig)

rc_determinations <- rc_determinations + delta_r_real

pp_uniform_phase_mixed_samples <- data.frame(c14_age = rc_determinations,
                                             c14_sig = rc_sigmas,
                                             sample_source = as.character(sample_source),
                                             delta_r = delta_r,
                                             delta_r_sig = delta_r_sig)
pp_uniform_phase_mixed_samples$f14c <- exp(-pp_uniform_phase$c14_age / 8033)
pp_uniform_phase_mixed_samples$f14c_sig <- pp_uniform_phase$f14c * pp_uniform_phase$c14_sig / 8033

usethis::use_data(pp_uniform_phase_mixed_samples, overwrite = TRUE)

