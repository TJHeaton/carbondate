library(devtools)
load_all()

get_data_kerr <- function(seeds) {
  for (seed in seeds) {
    set.seed(seed)

    wo = WalkerBivarDirichlet(kerr$c14_age, kerr$c14_sig, FALSE, intcal20, n_iter = 2e5)

    save(wo, file=paste("experimenting/convergence/kerr_seed_", seed, ".rda", sep=""))
  }
}


get_data_postbomb <- function(seeds) {
  HOBS2022 = get_hobs_calcurve()
  measurements = read.csv("experimenting/post-bomb/8836_modern.csv", header = TRUE, sep = ",")

  for (seed in seeds) {
    set.seed(seed)

    wo = WalkerBivarDirichlet(measurements$F14C, measurements$F14C_sd, TRUE, HOBS2022, n_iter = 1e5)

    save(wo, file=paste("experimenting/convergence/postbomb_seed_", seed, ".rda", sep=""))
  }
}


plot_data <- function(seeds, prefix, save_file = TRUE) {
  if (save_file) {
    pdf(
      file = paste("experimenting/convergence/", prefix, ".pdf", sep=""), width = 10, height = 10,
      pointsize = 10)
  }
  for (seed in seeds) {
    load(paste("experimenting/convergence/", prefix ,"_seed_", seed, ".rda", sep=""))

    calendar_ages = wo$density_data$calendar_ages
    iters = wo$density_data$iters/1000
    densities = wo$density_data$densities

    layout.matrix <- matrix(c(1, 4, 2, 3), nrow = 2, ncol = 2)
    layout(mat = layout.matrix, heights = c(1, 1), widths = c(1, 1))

    if (prefix == "postbomb") {
      PlotPredictiveCalendarAgeDensity(
        wo, 5000, plot_14C_age = FALSE, denscale = 5, show_SPD = TRUE,
        calibration_curve = get_hobs_calcurve())
      image_calendar_ages = rev(1950 - calendar_ages)
      densities = densities[length(calendar_ages):1, ]
      zlim = c(0, max(wo$density_data$densities) / 6)
      xlab = "Calendar Age (AD)"
    } else {
      PlotPredictiveCalendarAgeDensity(wo, 5000, show_SPD = TRUE)
      image_calendar_ages = calendar_ages
      zlim = c(0, max(wo$density_data$densities))
      xlab = "Calendar Age (yr BP)"
    }

    image(
      x=image_calendar_ages,
      y=iters,
      z=densities,
      xlab=xlab,
      ylab="k iterations",
      col = hcl.colors(256, "Spectral"),
      main = paste("Seed = ", seed),
      zlim = zlim)

    n_out = length(wo$density_data$iters)
    kld_instant = rep(NA, n_out - 1)
    kld_mean = rep(NA, n_out - 1)
    kld_instant_avg = rep(0, (n_out - 1) / 100)
    kld_mean_avg = rep(0, (n_out - 1) / 100)
    iters_avg = rep(NA, (n_out - 1) / 100)
    mean_density = densities[ ,1]
    avg_index = 1
    for (i in 1:(n_out - 1)) {
      kld_instant[i] = .KLD(densities[, 1], densities[, i + 1])
      kld_mean[i] = .KLD(mean_density, densities[, i + 1])

      kld_instant_avg[avg_index] = kld_instant_avg[avg_index] + kld_instant[i]
      kld_mean_avg[avg_index] = kld_mean_avg[avg_index] + kld_mean[i]
      if (i %% 100 == 0) {
        kld_instant_avg[avg_index] = kld_instant_avg[avg_index] / 100
        kld_mean_avg[avg_index] = kld_mean_avg[avg_index] / 100
        iters_avg[avg_index] = iters[i]
        avg_index = avg_index + 1
      }

      mean_density = mean_density + densities[, i + 1]
    }

    plot(iters[-1], log10(kld_instant), ylab = "log10 of KLD with first 10k iters", xlab = "k iterations", type="p", pch=".")
    lines(iters_avg, log10(kld_instant_avg), col="red")
    plot(iters[-1], log10(kld_mean), ylab = "log10 of KLD with mean density", xlab = "k iterations", type="p", pch=".")
    lines(iters_avg, log10(kld_mean_avg), col="red")
  }
  if (save_file) dev.off()
}



compare_data <- function(seed_1, seed_2, prefix, save_file = TRUE) {
  if (save_file) {
    pdf(
      file = paste("experimenting/convergence/", prefix, "_", seed_1, "_", seed_2, ".pdf", sep=""),
      width = 10, height = 10, pointsize = 10)
  }
  load(paste("experimenting/convergence/", prefix ,"_seed_", seed_1, ".rda", sep=""))
  wo_1 = wo
  wo_1$label = paste("Seed = ", seed_1, sep = "")

  load(paste("experimenting/convergence/", prefix ,"_seed_", seed_2, ".rda", sep=""))
  wo_2 = wo
  wo_2$label = paste("Seed = ", seed_2, sep = "")


  calendar_ages = wo_1$density_data$calendar_ages
  iters = wo_1$density_data$iters/1000
  densities_1 = wo_1$density_data$densities
  densities_2 = wo_2$density_data$densities

  layout.matrix <- matrix(c(1, 3, 2, 4), nrow = 2, ncol = 2)
  layout(mat = layout.matrix, heights = c(1, 1), widths = c(1, 1))

  if (prefix == "postbomb") {
    PlotPredictiveCalendarAgeDensity(
      list(wo_1, wo_2), 5000, plot_14C_age = FALSE, denscale = 5, show_SPD = TRUE,
      calibration_curve = get_hobs_calcurve())
    image_calendar_ages = rev(1950 - calendar_ages)
    densities = densities[length(calendar_ages):1, ]
    zlim = c(0, max(wo$density_data$densities) / 6)
    xlab = "Calendar Age (AD)"
  } else {
    PlotPredictiveCalendarAgeDensity(list(wo_1, wo_2), 5000, show_SPD = TRUE)
    image_calendar_ages = calendar_ages
    zlim = c(0, max(wo$density_data$densities))
    xlab = "Calendar Age (yr BP)"
  }

  n_out = length(wo$density_data$iters)

  kld_1 = rep(NA, n_out - 1)
  kld_2 = rep(NA, n_out - 1)
  kld_diff = rep(NA, n_out)
  kld_1_avg = rep(0, (n_out - 1) / 100)
  kld_2_avg = rep(0, (n_out - 1) / 100)
  kld_diff_avg = rep(0, (n_out - 1) / 100)
  iters_avg = rep(NA, (n_out - 1) / 100)
  avg_index = 1
  kld_diff[1] = .KLD(densities_1[, 1], densities_2[, 1])
  for (i in 1:(n_out - 1)) {
    kld_1[i] = .KLD(densities_1[, 1], densities_1[, i + 1])
    kld_2[i] = .KLD(densities_2[, 1], densities_2[, i + 1])
    kld_diff[i + 1] = .KLD(densities_1[, i + 1], densities_2[, i + 1])

    kld_1_avg[avg_index] = kld_1_avg[avg_index] + kld_1[i]
    kld_2_avg[avg_index] = kld_2_avg[avg_index] + kld_2[i]
    kld_diff_avg[avg_index] = kld_diff_avg[avg_index] + kld_diff[i]
    if (i %% 100 == 0) {
      kld_1_avg[avg_index] = kld_1_avg[avg_index] / 100
      kld_2_avg[avg_index] = kld_2_avg[avg_index] / 100
      kld_diff_avg[avg_index] = kld_diff_avg[avg_index] / 100
      iters_avg[avg_index] = iters[i]
      avg_index = avg_index + 1
    }
  }
  plot(
    iters, log10(kld_diff), ylab = "log10 of KLD between different chains",
    xlab = "k iterations", type="p", pch=".", main = "Both chains")
  lines(iters_avg, log10(kld_diff_avg), col="red")
  plot(
    iters[-1], log10(kld_1), ylab = "log10 of KLD with first 10k iters",
    xlab = "k iterations", type="p", pch=".", main = paste("Seed = ", seed_1))
  lines(iters_avg, log10(kld_1_avg), col="red")
  plot(
    iters[-1], log10(kld_2), ylab = "log10 of KLD with first 10k iters",
    xlab = "k iterations", type="p", pch=".", main = paste("Seed = ", seed_2))
  lines(iters_avg, log10(kld_2_avg), col="red")
  if (save_file) dev.off()
}


gellman_rubin <- function(seeds, prefix) {

  M = length(seeds)
  theta = list()
  theta_mean = rep(0, M)
  theta_variance = rep(0, M)

  for (m in 1:M) {
    load(paste("experimenting/convergence/", prefix ,"_seed_", seeds[m], ".rda", sep=""))

    theta[[m]] = wo$calendar_ages[, 10]
  }

  n_out = length(theta[[1]])
  B = rep(0, n_out - 1)
  W = rep(0, n_out - 1)
  for (N in 2:n_out) {
    for (m in 1:M) {
      theta_mean[m] = mean(theta[[m]][1:N])
      theta_variance[m] = var(theta[[m]][1:N])
    }
    overall_mean = mean(theta_mean)
    B[N - 1] = N / (M - 1) * sum((theta_mean - overall_mean) ^ 2)
    W[N - 1] = 1 / M * sum(theta_variance)
  }
  plot(B, type="l")
  plot(W, type="l")
}

gellman_rubin_single <- function(seed, prefix) {

  M = length(seeds)
  theta = list()
  theta_mean = rep(0, M)
  theta_variance = rep(0, M)

  for (m in 1:M) {
    load(paste("experimenting/convergence/", prefix ,"_seed_", seeds[m], ".rda", sep=""))

    theta[[m]] = wo$calendar_ages[, 10]
  }

  n_out = length(theta[[1]])
  B = rep(0, n_out - 1)
  W = rep(0, n_out - 1)
  for (N in 2:n_out) {
    for (m in 1:M) {
      theta_mean[m] = mean(theta[[m]][1:N])
      theta_variance[m] = var(theta[[m]][1:N])
    }
    overall_mean = mean(theta_mean)
    B[N - 1] = N / (M - 1) * sum((theta_mean - overall_mean) ^ 2)
    W[N - 1] = 1 / M * sum(theta_variance)
  }
  plot(B, type="l")
  plot(W, type="l")
}


get_hobs_calcurve <- function() {
  HOBS2022 = read.csv(
    "experimenting/post-bomb/HOBS2022.14c", header=FALSE, skip=2, sep="")[, c(1, 4, 5)]
  names(HOBS2022) = c("calendar_age_BP", "f14c", "f14c_sig")
  HOBS2022$calendar_age_BP = 1950 - HOBS2022$calendar_age_BP
  return(HOBS2022)
}

