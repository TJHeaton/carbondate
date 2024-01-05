## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(carbondate)

## ----calculate_gr_multiple, fig.width=10, fig.height=8, results=FALSE---------
all_outputs <- list()
for (i in 1:3) {
  set.seed(i + 1)
  all_outputs[[i]] <- PolyaUrnBivarDirichlet(kerr$c14_age, kerr$c14_sig, intcal20, n_iter = 1e4)
}
PlotGelmanRubinDiagnosticMultiChain(all_outputs)

## ----calculate_gr, fig.width=10, fig.height=8, results=FALSE------------------
set.seed(3)
output <- PolyaUrnBivarDirichlet(kerr$c14_age, kerr$c14_sig, intcal20, n_iter = 2e4)

PlotGelmanRubinDiagnosticSingleChain(output, n_burn = 5e3)

## ----calculate_polya_kerr, fig.width=10, fig.height=8, results=FALSE----------
all_outputs <- list()
for (i in 1:3) {
  set.seed(i+1)
  all_outputs[[i]] <- PolyaUrnBivarDirichlet(
  rc_determinations = kerr$c14_age,
  rc_sigmas = kerr$c14_sig,
  calibration_curve=intcal20,
  n_iter = 1e4)
  all_outputs[[i]]$label <- paste("Seed =", i)
}
PlotPredictiveCalendarAgeDensity(
  output_data = all_outputs, n_posterior_samples = 500, denscale = 2.5, interval_width = "1sigma")

## ----calculate_polya_normals, fig.width=10, fig.height=8, results=FALSE-------
all_outputs <- list()
for (i in 1:3) {
  set.seed(i + 1)
  all_outputs[[i]] <- PolyaUrnBivarDirichlet(
  rc_determinations = two_normals$c14_age,
  rc_sigmas = two_normals$c14_sig,
  calibration_curve=intcal20,
  n_iter = 1e4)
  all_outputs[[i]]$label <- paste("Seed =", i)
}
PlotPredictiveCalendarAgeDensity(
  output_data = all_outputs, n_posterior_samples = 500, denscale = 2.5, interval_width = "1sigma")

## ----calculate_kld, fig.width=10, fig.height=8, results=FALSE-----------------
set.seed(50)
output <- WalkerBivarDirichlet(
rc_determinations = kerr$c14_age,
rc_sigmas = kerr$c14_sig,
calibration_curve=intcal20,
n_iter = 1e5)

PlotConvergenceData(output)

