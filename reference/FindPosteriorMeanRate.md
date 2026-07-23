# Find Posterior Mean Rate of Sample Occurrence for Poisson Process Model

Given output from the Poisson process fitting function
[PPcalibrate](https://tjheaton.github.io/carbondate/reference/PPcalibrate.md)
calculate the posterior mean rate of sample occurrence (i.e., the
underlying Poisson process rate \\\lambda(t)\\) together with specified
probability intervals, on a given calendar age grid (provided in cal yr
BP).

An option is also provided to calculate the posterior mean rate
*conditional* upon the number of internal changepoints within the period
under study (if this is specified as an input to the function).

**Note:** If you want to calculate and plot the result, use
[PlotPosteriorMeanRate](https://tjheaton.github.io/carbondate/reference/PlotPosteriorMeanRate.md)
instead.

For more information read the vignette:  
[`vignette("Poisson-process-modelling", package = "carbondate")`](https://tjheaton.github.io/carbondate/articles/Poisson-process-modelling.md)

## Usage

``` r
FindPosteriorMeanRate(
  output_data,
  calendar_age_sequence,
  n_posterior_samples = 5000,
  n_changes = NULL,
  interval_width = "2sigma",
  bespoke_probability = NA,
  n_burn = NA,
  n_end = NA
)
```

## Arguments

- output_data:

  The return value from the updating function
  [PPcalibrate](https://tjheaton.github.io/carbondate/reference/PPcalibrate.md).
  Optionally, the output data can have an extra list item named `label`
  which is used to set the label on the plot legend.

- calendar_age_sequence:

  A vector containing the calendar age grid (in cal yr BP) on which to
  calculate the posterior mean rate.

- n_posterior_samples:

  Number of samples it will draw, after having removed `n_burn`, from
  the (thinned) MCMC realisations stored in `output_data` to estimate
  the rate \\\lambda(t)\\. These samples may be repeats if the number
  of, post burn-in, realisations is less than `n_posterior_samples`. If
  not given, 5000 is used.

- n_changes:

  (Optional) If wish to condition calculation of the posterior mean on
  the number of internal changepoints. In this function, if specified,
  `n_changes` must be a single integer.

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
  `n_iter` and `n_thin` are the arguments that were given to
  [PPcalibrate](https://tjheaton.github.io/carbondate/reference/PPcalibrate.md))
  which would leave only 100 of the (thinned) values in `output_data`.

- n_end:

  The last iteration in the original MCMC chain to use in the
  calculations. Assumed to be the total number of iterations performed,
  i.e. `n_iter`, if not given.

## Value

A list, each item containing a data frame of the `calendar_age_BP`, the
`rate_mean` and the confidence intervals for the rate - `rate_ci_lower`
and `rate_ci_upper`.

## See also

[PlotPosteriorMeanRate](https://tjheaton.github.io/carbondate/reference/PlotPosteriorMeanRate.md)

## Examples

``` r
# NOTE: All these examples are shown with a small n_iter and n_posterior_samples
# to speed up execution.
# Try n_iter and n_posterior_samples as the function defaults.

pp_output <- PPcalibrate(
    pp_uniform_phase$c14_age,
    pp_uniform_phase$c14_sig,
    intcal20,
    n_iter = 1000,
    show_progress = FALSE)

# Default plot with 2 sigma interval
FindPosteriorMeanRate(pp_output, seq(450, 640, length=10), n_posterior_samples = 100)
#>    calendar_age_BP  rate_mean rate_ci_lower rate_ci_upper
#> 1         450.0000 0.05386204   0.026957329    0.08402529
#> 2         471.1111 0.05206876   0.041401905    0.08402529
#> 3         492.2222 0.05084057   0.036296024    0.06872355
#> 4         513.3333 0.84174298   0.645966331    1.02825854
#> 5         534.4444 0.84174298   0.645966331    1.02825854
#> 6         555.5556 0.02479378   0.007654306    0.10425931
#> 7         576.6667 0.02112564   0.007654306    0.04869605
#> 8         597.7778 0.02098722   0.007654306    0.04869605
#> 9         618.8889 0.02061186   0.006032444    0.04869605
#> 10        640.0000 0.02145073   0.007654306    0.04869605
```
