# Find Predictive Estimate of Shared Calendar Age Density from Bayesian Non-Parametric DPMM Output

Given output from one of the Bayesian non-parametric summarisation
functions (either
[PolyaUrnBivarDirichlet](https://tjheaton.github.io/carbondate/reference/PolyaUrnBivarDirichlet.md)
or
[WalkerBivarDirichlet](https://tjheaton.github.io/carbondate/reference/WalkerBivarDirichlet.md))
calculate the predictive (summarised/shared) calendar age density and
probability intervals on a given calendar age grid (provided in cal yr
BP).

**Note:** If you want to calculate and plot the result, use
[PlotPredictiveCalendarAgeDensity](https://tjheaton.github.io/carbondate/reference/PlotPredictiveCalendarAgeDensity.md)
instead.

## Usage

``` r
FindPredictiveCalendarAgeDensity(
  output_data,
  calendar_age_sequence,
  n_posterior_samples = 5000,
  interval_width = "2sigma",
  bespoke_probability = NA,
  n_burn = NA,
  n_end = NA
)
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

- calendar_age_sequence:

  A vector containing the calendar age grid (in cal yr BP) on which to
  calculate the predictive (summarised/shared) density.

- n_posterior_samples:

  Number of samples it will draw, after having removed `n_burn`, from
  the (thinned) realisations stored in the DPMM outputs to estimate the
  predictive calendar age density. These samples may be repeats if the
  number of, post burn-in, realisations is less than
  `n_posterior_samples`. If not given, 5000 is used.

- interval_width:

  The confidence intervals to show for both the calibration curve and
  the predictive density. Choose from one of `"1sigma"` (68.3%),
  `"2sigma"` (95.4%) and `"bespoke"`. Default is `"2sigma"`.

- bespoke_probability:

  The probability to use for the confidence interval if `"bespoke"` is
  chosen above. E.g., if 0.95 is chosen, then the 95% confidence
  interval is calculated. Ignored if `"bespoke"` is not chosen.

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

A data frame of the `calendar_age_BP`, the `density_mean` and the
confidence intervals for the density `density_ci_lower` and
`density_ci_upper`.

## See also

[PlotPredictiveCalendarAgeDensity](https://tjheaton.github.io/carbondate/reference/PlotPredictiveCalendarAgeDensity.md)

## Examples

``` r
# NOTE: All these examples are shown with a small n_iter and n_posterior_samples
# to speed up execution.
# Try n_iter and n_posterior_samples as the function defaults.

# First generate output data
polya_urn_output <- PolyaUrnBivarDirichlet(
    two_normals$c14_age,
    two_normals$c14_sig,
    intcal20,
    n_iter = 100,
    show_progress = FALSE)

# Find results for example output, 2-sigma confidence interval (default)
FindPredictiveCalendarAgeDensity(
    polya_urn_output, seq(3600, 4700, length=12), n_posterior_samples = 500)
#>    calendar_age_BP density_mean density_ci_lower density_ci_upper
#> 1             3600 7.148234e-04     5.121023e-04     9.118296e-04
#> 2             3700 5.365339e-04     3.441233e-04     6.651168e-04
#> 3             3800 3.258515e-04     1.965341e-04     4.089631e-04
#> 4             3900 1.623642e-04     9.553998e-05     2.314653e-04
#> 5             4000 6.795576e-05     3.972726e-05     1.143602e-04
#> 6             4100 2.474368e-05     1.156121e-05     4.938953e-05
#> 7             4200 8.377500e-06     3.731912e-06     1.863494e-05
#> 8             4300 3.111155e-06     1.129557e-06     6.587292e-06
#> 9             4400 1.806583e-06     9.453687e-07     3.310225e-06
#> 10            4500 5.692710e-06     9.634150e-07     1.248723e-05
#> 11            4600 4.855400e-05     1.353136e-06     1.380441e-04
#> 12            4700 1.956020e-04     2.171993e-05     3.628677e-04
```
