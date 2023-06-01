
<!-- README.md is generated from README.Rmd. Please edit that file -->

# carbondate

<!-- badges: start -->

[![R-CMD-check](https://github.com/TJHeaton/carbondate/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TJHeaton/carbondate/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/carbondate)](https://CRAN.R-project.org/package=carbondate)
<!-- badges: end -->

An R package to analyse multiple radiocarbon determinations. It is based
on the original functions available
[here](https://github.com/TJHeaton/NonparametricCalibration) which were
used for “Non-parametric calibration of multiple related radiocarbon
determinations and their calendar age summarisation”
[arXiv](https://arxiv.org/abs/2109.15024).

## Installation

You can install the development version of carbondate from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("TJHeaton/carbondate", build_vignettes = TRUE)
```

## Example

There are 4 example datasets (`two_normals`, `kerr`, `buchanan`,
`armit`) provided, which can be used to try out the calibration
functions. `two_normals` is a small artificial data set of 50
radiocarbon determinations for which the underlying calendar ages were
drawn from a mixture of two normals. It is included simply to give some
quick-to-run examples. The remaining datasets are from real-life data.
The calibration curves IntCal98 through to IntCal20 are also included in
the package. E.g. see below.

``` r
library(carbondate)

walker_output <- WalkerBivarDirichlet(
  rc_determinations = two_normals$c14_age,
  rc_sigmas = two_normals$c14_sig,
  calibration_curve=intcal20)

polya_urn_output <- PolyaUrnBivarDirichlet(
  rc_determinations = two_normals$c14_age,
  rc_sigmas = two_normals$c14_sig,
  calibration_curve=intcal20)
```

Once the calibration has been run, the calendar age density can be
plotted.

``` r
densities = PlotPredictiveCalendarAgeDensity(
  output_data = list(walker_output, polya_urn_output),
  n_posterior_samples = 5000,
  show_SPD = TRUE)
```

<img src="man/figures/README-plot_density-1.png" width="100%" />

For a full example run-through, load the vignette:

``` r
browseVignettes(package="carbondate")
```
