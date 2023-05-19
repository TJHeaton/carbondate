
for (seed in 1:20) {
  set.seed(seed)

  wo = WalkerBivarDirichlet(kerr$c14_age, kerr$c14_sig, FALSE, intcal20, n_iter = 5e5)

  layout.matrix <- matrix(c(1, 4, 2, 3), nrow = 2, ncol = 2)
  layout(mat = layout.matrix, heights = c(1, 1), widths = c(1, 1))

  PlotPredictiveCalendarAgeDensity(wo, 5000)

  calendar_ages = wo$density_data$calendar_ages
  iters = wo$density_data$iters/1000
  densities = wo$density_data$densities

  image(
    x=calendar_ages,
    y=iters,
    z=t(densities),
    xlab="Calendar Age BP",
    ylab="k iterations",
    col = hcl.colors(256, "Greens"),
    main = paste("Seed = ", seed))

  n_out = length(wo$density_data$iters)
  klds = rep(NA, n_out - 1)
  plot(calendar_ages, densities[1,], type="l", col="red", xlab="Calendar Age BP", ylab="Predicitive density for 1000 steps")
  for (i in 2:n_out) {
    frac = (i - 1) / (n_out - 1)
    if (frac < 0.5) {
      color = rgb(1 - 2 * frac, 0, 2 * frac)
    } else {
      color = rgb(0, 2 * (frac - 0.5), 1 - 2 * (frac - 0.5))
    }
    lines(calendar_ages, densities[i,], col=color)
    klds[i-1] = .KLD(densities[i - 1, ], densities[i, ])
  }

  plot(iters[-1], log10(klds), ylab = "log10 of Kullback Liebler Divergence", xlab = "k iterations")

}
