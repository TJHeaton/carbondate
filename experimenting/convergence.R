library(devtools)
library(carbondate)

plot_convergence_walker_kerr <- function(seeds, n_iter, n_thin, n_burn) {
  pdf(
    file = paste0("experimenting/convergence_2/kerr_walker_B_", n_iter / 1000, ".pdf"),
    width = 10, height = 5, pointsize = 10)
  for (seed in seeds) {
    print(paste("seed =", seed))
    set.seed(seed)

    filename <- paste0("experimenting/convergence_2/kerr_walker_B_seed_", seed, ".rda")

    if (file.exists(filename)) {
      load(filename)
    } else {
      wo <- WalkerBivarDirichlet(kerr$c14_age, kerr$c14_sig, intcal20, n_iter = n_iter, n_thin = n_thin)
      save(wo, file=filename)
    }

    n_posterior_samples <- 20000
    layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
    layout(mat = layout.matrix, heights = c(1, 1), widths = c(1, 1))

    PlotPredictiveCalendarAgeDensity(wo, n_posterior_samples, n_burn = n_burn)
    title(sub = paste("n_burn =", n_burn, "n_samples = ", n_posterior_samples))
    PlotConvergenceData(wo)
    title(sub = paste("seed =", seed))
  }
  p <- dev.cur()
  dev.off(p)
}

plot_convergence_kerr_shifted <- function(seeds, n_iter, n_thin, n_burn) {
  pdf(
    file = paste0("experimenting/convergence_2/kerr_walker_shifted_", n_iter / 1000, ".pdf"),
    width = 10, height = 5, pointsize = 10)
  for (seed in seeds) {
    print(paste("seed =", seed))
    set.seed(seed)

    filename <- paste0("experimenting/convergence_2/kerr_walker_shifted_", n_iter / 1000, "_seed_", seed, ".rda")

    if (file.exists(filename)) {
      load(filename)
    } else {
      wo <- WalkerBivarDirichlet(kerr$c14_age + 2000, kerr$c14_sig, intcal20, n_iter = n_iter, n_thin = n_thin)
      save(wo, file=filename)
    }

    n_posterior_samples <- 20000
    layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
    layout(mat = layout.matrix, heights = c(1, 1), widths = c(1, 1))

    PlotPredictiveCalendarAgeDensity(wo, n_posterior_samples, n_burn = n_burn)
    title(sub = paste("n_burn =", n_burn, "n_samples = ", n_posterior_samples))
    PlotConvergenceData(wo)
    title(sub = paste("seed =", seed))
  }
  p <- dev.cur()
  dev.off(p)
}

plot_convergence_walker_pu_kerr <- function(seeds, n_iter, n_thin, n_burn) {
  pdf(
    file = paste0("experimenting/convergence_2/kerr_with_clusters_", n_iter / 1000, ".pdf"),
    width = 12, height = 5, pointsize = 10)
  for (seed in seeds) {
    print(paste("seed =", seed))
    set.seed(seed)

    filename <- paste0("experimenting/convergence_2/data/kerr_walker_A_", n_iter / 1000, "_seed_", seed, ".rda")

    if (file.exists(filename)) {
      load(filename)
    } else {
      wo <- WalkerBivarDirichlet(kerr$c14_age, kerr$c14_sig, intcal20, n_iter = n_iter, n_thin = n_thin)
      save(wo, file=filename)
    }

    set.seed(seed)
    filename <- paste0("experimenting/convergence_2/data/kerr_polya_A_", n_iter / 1000, "_seed_", seed, ".rda")
    if (file.exists(filename)) {
      load(filename)
    } else {
      po <- PolyaUrnBivarDirichlet(kerr$c14_age, kerr$c14_sig, intcal20, n_iter = n_iter, n_thin = n_thin)
      save(po, file=filename)
    }

    n_posterior_samples <- 5000
    layout.matrix <- matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2)
    layout(mat = layout.matrix, heights = c(3, 3), widths = c(1.7, 1))

    PlotPredictiveCalendarAgeDensity(list(wo, po), n_posterior_samples, n_burn = n_burn)
    title(sub = paste("n_burn =", n_burn, "n_samples = ", n_posterior_samples))

    PlotConvergenceData(wo)
    title(sub = paste("seed =", seed))

    PlotConvergenceData(po)
    title(sub = paste("seed =", seed))
  }
  p <- dev.cur()
  dev.off(p)
}


plot_all_kerr <- function(seeds, n_iter, n_thin) {
  pdf(
    file = paste0("experimenting/convergence_2/kerr_all_", n_iter / 1000, "_2_sigma.pdf"),
    width = 12, height = 5, pointsize = 10)
  walker <- list()
  pu <- list()
  for (seed in seeds) {
    print(paste("seed =", seed))
    set.seed(seed)

    filename <- paste0("experimenting/convergence_2/data/kerr_walker_A_", n_iter / 1000, "_seed_", seed, ".rda")

    if (file.exists(filename)) {
      load(filename)
    } else {
      wo <- WalkerBivarDirichlet(kerr$c14_age, kerr$c14_sig, intcal20, n_iter = n_iter, n_thin = n_thin)
      save(wo, file=filename)
    }
    wo$label <- paste(seed)
    walker[[seed]] <- wo

    set.seed(seed)
    filename <- paste0("experimenting/convergence_2/data/kerr_polya_A_", n_iter / 1000, "_seed_", seed, ".rda")
    if (file.exists(filename)) {
      load(filename)
    } else {
      po <- PolyaUrnBivarDirichlet(kerr$c14_age, kerr$c14_sig, intcal20, n_iter = n_iter, n_thin = n_thin)
      save(po, file=filename)
    }
    po$label <- paste(seed)
    pu[[seed]] <- po
  }

  for (n_end in seq(1e5, 1e6, by = 1e5)) {
    print(n_end)
    n_posterior_samples <- 1000
    n_burn <- n_end / 2
    layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
    layout(mat = layout.matrix, heights = c(3), widths = c(1, 1))

    PlotPredictiveCalendarAgeDensity(walker, n_posterior_samples, n_burn = n_burn, n_end = n_end, interval_width = "2sigma")
    title(sub = paste("Walker: n_burn =", n_burn, "n_end = ", n_end))

    PlotPredictiveCalendarAgeDensity(pu, n_posterior_samples, n_burn = n_burn, n_end = n_end, interval_width = "2sigma")
    title(sub = paste("Polya: n_burn =", n_burn, "n_end = ", n_end))
  }

  p <- dev.cur()
  dev.off(p)
}

get_data_postbomb <- function(seeds) {
  HOBS2022 <- get_hobs_calcurve()
  measurements <- read.csv("experimenting/post-bomb/8836_modern.csv", header = TRUE, sep = ",")

  for (seed in seeds) {
    set.seed(seed)

    wo <- WalkerBivarDirichlet(measurements$F14C, measurements$F14C_sd, TRUE, HOBS2022, n_iter = 1e5)

    save(wo, file= paste0("experimenting/convergence/postbomb_seed_", seed, ".rda"))
  }
}

gellman_rubin <- function(seeds, prefix) {

  M <- length(seeds)
  theta <- list()

  for (m in 1:M) {
    load(paste0("experimenting/convergence_2/data/", prefix, "_seed_", seeds[m], ".rda"))

    theta[[m]] <- wo$calendar_ages[, 10]
  }
  R <- CalculatePsrf(theta, 10, length(theta[[1]]))
  plot(R, type="l")
}


single_chain_gellman_rubin <- function(output, n_burn) {
  n_obs <- length(output$input_data$rc_determinations)
  n_thin <- output$input_parameters$n_thin
  R <- rep(0, n_obs)

  for (i in 1:n_obs) {
    R[i] <- single_chain_psrf(output$calendar_ages[, i], n_burn, n_thin)
  }

  return (R)
}


single_chain_psrf <- function(theta, n_burn, n_thin) {
  n_out <- length(theta)
  n_burn <- floor(n_burn / n_thin) + 1
  n_chain <- floor((n_out - n_burn) / 3)

  chain_one <- theta[(n_burn + 1):(n_burn + n_chain)]
  chain_two <- theta[(n_burn + n_chain + 1):(n_burn + 2 * n_chain)]
  chain_three <- theta[(n_burn + 2 * n_chain + 1):(n_burn + 3 * n_chain)]

  overall_mean <- mean(theta[n_burn:(n_burn + 3 * n_chain)])

  within_chain_variance <- 0.5 * (var(chain_one) + var(chain_two) + var(chain_three))
  between_chain_variance <- 0.5 * n_chain * (
    (mean(chain_one) - overall_mean)^2
    + (mean(chain_two) - overall_mean)^2
    + (mean(chain_three) - overall_mean)^2
  )

  pooled_variance <- (n_chain - 1) / n_chain * within_chain_variance
  pooled_variance <- pooled_variance + 4 / (3 * n_chain) * between_chain_variance

  return(pooled_variance / within_chain_variance)
}

calculate_psrf <- function(theta) {

  M <- length(theta)
  theta_mean <- rep(0, M)
  theta_variance <- rep(0, M)

  n_out <- length(theta[[1]])
  B <- rep(0, n_out - 10)
  W <- rep(0, n_out - 10)
  V <- rep(0, n_out - 10)
  for (N in 11:n_out) {
    for (m in 1:M) {
      theta_mean[m] <- mean(theta[[m]][1:N])
      theta_variance[m] <- var(theta[[m]][1:N])
    }
    overall_mean <- mean(theta_mean)
    B[N - 10] <- N / (M - 1) * sum((theta_mean - overall_mean) ^ 2)
    W[N - 10] <- 1 / M * sum(theta_variance)

    V[N - 10] <- (N - 1) / N * W[N - 10] + (M + 1) / (M * N) * B[N - 10]
  }
  return(V / W)
}



get_hobs_calcurve <- function() {
  HOBS2022 <- read.csv(
    "experimenting/post-bomb/HOBS2022.14c", header=FALSE, skip=2, sep="")[, c(1, 4, 5)]
  names(HOBS2022) <- c("calendar_age_BP", "f14c", "f14c_sig")
  HOBS2022$calendar_age_BP <- 1950 - HOBS2022$calendar_age_BP
  return(HOBS2022)
}

