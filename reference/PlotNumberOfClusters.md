# Plot Number of Calendar Age Clusters Estimated in Bayesian Non-Parametric DPMM Output

Given output from one of the Bayesian non-parametric summarisation
functions (either
[PolyaUrnBivarDirichlet](https://tjheaton.github.io/carbondate/reference/PolyaUrnBivarDirichlet.md)
or
[WalkerBivarDirichlet](https://tjheaton.github.io/carbondate/reference/WalkerBivarDirichlet.md))
plot the estimated number of calendar age clusters represented by the
\\{}^{14}\\C samples.

For more information read the vignette:  
[`vignette("Non-parametric-summed-density", package = "carbondate")`](https://tjheaton.github.io/carbondate/articles/Non-parametric-summed-density.md)

## Usage

``` r
PlotNumberOfClusters(output_data, n_burn = NA, n_end = NA)
```

## Arguments

- output_data:

  The return value from one of the Bayesian non-parametric DPMM
  functions, e.g.
  [PolyaUrnBivarDirichlet](https://tjheaton.github.io/carbondate/reference/PolyaUrnBivarDirichlet.md)
  or
  [WalkerBivarDirichlet](https://tjheaton.github.io/carbondate/reference/WalkerBivarDirichlet.md),
  or a list, each item containing one of these return values.
  Optionally, the output data can have an extra list item named `label`
  which is used to set the label on the plot legend.

- n_burn:

  The number of MCMC iterations that should be discarded as burn-in
  (i.e., considered to be occurring before the MCMC has converged). This
  relates to the number of iterations (`n_iter`) when running the
  original update functions (not the thinned `output_data`). Any MCMC
  iterations before this are not used in the calculations. If not given,
  the first half of the MCMC chain is discarded. Note: The maximum value
  that the function will allow is `n_iter - 100 * n_thin` (where
  `n_iter` and `n_thin` are the arguments given to
  [PolyaUrnBivarDirichlet](https://tjheaton.github.io/carbondate/reference/PolyaUrnBivarDirichlet.md)
  or
  [WalkerBivarDirichlet](https://tjheaton.github.io/carbondate/reference/WalkerBivarDirichlet.md))
  which would leave only 100 of the (thinned) values in `output_data`.

- n_end:

  The last iteration in the original MCMC chain to use in the
  calculations. Assumed to be the total number of iterations performed,
  i.e. `n_iter`, if not given.

## Value

None

## See also

[PlotPredictiveCalendarAgeDensity](https://tjheaton.github.io/carbondate/reference/PlotPredictiveCalendarAgeDensity.md)
and
[PlotCalendarAgeDensityIndividualSample](https://tjheaton.github.io/carbondate/reference/PlotCalendarAgeDensityIndividualSample.md)
for more plotting functions using DPMM output.

## Examples

``` r
# NOTE: these examples are shown with a small n_iter to speed up execution.
# When you run ensure n_iter gives convergence (try function default).

polya_urn_output <- PolyaUrnBivarDirichlet(
    two_normals$c14_age,
    two_normals$c14_sig,
    intcal20,
    n_iter = 500,
    n_thin = 2,
    show_progress = FALSE)

PlotNumberOfClusters(polya_urn_output)
```
