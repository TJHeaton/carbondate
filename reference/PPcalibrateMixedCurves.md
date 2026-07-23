# Model Occurrence of Multiple Radiocarbon Samples as a Variable-Rate Poisson Process when the Samples Come from a Variety of Environments and Multiple Calibration Curves are Needed

This function calibrates a set of radiocarbon (\\{}^{14}\\C) samples,
and provides an estimate of how the underlying rate at which the samples
occurred varies over calendar time (including any specific changepoints
in the rate). The method can be used as an alternative approach to
summarise calendar age information contained in a set of related
\\{}^{14}\\C samples, enabling inference on the latent *activity* rate
which led to their creation.

Importantly, this function allows for different samples to be calibrated
against different curves For example, if some of the samples come from
Marine environments (and so should be calibrated against the relevant
MarineCal curve with a \\\Delta R\\ offset) while other samples are
atmospheric and are to be calibrated against the relevant IntCal (NH) or
SHCal (SH) calibration curve.

It takes as an input a set of \\{}^{14}\\C determinations and associated
\\1\sigma\\ uncertainties, as well as the radiocarbon age calibration
curve to be used. The (calendar age) occurrence of these radiocarbon
(\\{}^{14}\\C) samples is modelled as a Poisson process. The underlying
rate of this Poisson process \\\lambda(t)\\, which represents the rate
at which the samples occurred, is considered unknown and assumed to vary
over calendar time.

Specifically, the sample occurrence rate \\\lambda(t)\\ is modelled as
piecewise constant, but with an unknown number of changepoints, which
can occur at unknown times. The value of \\\lambda(t)\\ between any set
of changepoints is also unknown. The function jointly calibrates the
given \\{}^{14}\\C samples under this model, and simultaneously provides
an estimate of \\\lambda(t)\\. Fitting is performed using
reversible-jump MCMC within Gibbs.

It returns estimates for the calendar age of each individual radiocarbon
sample in the input set; and broader output on the estimated value of
\\\lambda(t)\\ which can be used by other library functions. Analysis of
the estimated changepoints in the piecewise \\\lambda(t)\\ permits wider
inference on whether the occurrence rate of samples changed
significantly at any particular calendar time and, if so, when and by
how much.

For more information read the vignette:  
[`vignette("Poisson-process-modelling", package = "carbondate")`](https://tjheaton.github.io/carbondate/articles/Poisson-process-modelling.md)

## Usage

``` r
PPcalibrateMixedCurves(
  rc_determinations,
  rc_sigmas,
  sample_source,
  delta_r = NULL,
  delta_r_sig = NULL,
  curve_year = "2020",
  F14C_inputs = FALSE,
  n_iter = 1e+05,
  n_thin = 10,
  use_F14C_space = TRUE,
  show_progress = TRUE,
  calendar_age_range = NA,
  calendar_grid_resolution = 1,
  prior_h_shape = NA,
  prior_h_rate = NA,
  prior_n_internal_changepoints_lambda = 3,
  k_max_internal_changepoints = 30,
  rescale_factor_rev_jump = 0.9,
  bounding_range_prob_cutoff = 0.001,
  initial_n_internal_changepoints = 10,
  grid_extension_factor = 0.1,
  use_fast = TRUE,
  fast_approx_prob_cutoff = 0.001
)
```

## Arguments

- rc_determinations:

  A vector of observed radiocarbon determinations. Can be provided
  either as \\{}^{14}\\C ages (in \\{}^{14}\\C yr BP) or as
  F\\{}^{14}\\C concentrations.

- rc_sigmas:

  A vector of the (1-sigma) measurement uncertainties for the
  radiocarbon determinations. Must be the same length as
  `rc_determinations` and given in the same units.

- sample_source:

  A vector describing the source of each \\{}^{14}\\C sample and
  determining the calibration curve to use. Choices for the values in
  the vector are `"NH"`, `"SH"` and `"Marine"`. Must be the same length
  as `rc_determinations`

- delta_r, delta_r_sig:

  A vector of the \\\Delta R\\ offset and associated 1\\\sigma\\
  uncertainty for each sample. For atmospheric samples, select `0`
  unless you want to model them as offset from the relevant curve.

- curve_year:

  One of `"2020"`, `"2013"`, `"2009"` or `"2004"` specifying the
  calibration curve year to use. Make sure you use the \\\Delta R\\
  appropriate to the calibration curve.

- F14C_inputs:

  `TRUE` if the provided `rc_determinations` are F\\{}^{14}\\C
  concentrations and `FALSE` if they are radiocarbon ages. Defaults to
  `FALSE`.

- n_iter:

  The number of MCMC iterations (optional). Default is 100,000.

- n_thin:

  How much to thin the MCMC output (optional). Will store every
  `n_thin`\\{}^\textrm{th}\\ iteration. 1 is no thinning, while a larger
  number will result in more thinning. Default is 10. Must choose an
  integer greater than 1. Overall number of MCMC realisations stored
  will be \\n\_{\textrm{out}} = \textrm{floor}(
  n\_{\textrm{iter}}/n\_{\textrm{thin}}) + 1\\ so do not choose `n_thin`
  too large to ensure there are enough samples from the posterior to use
  for later inference.

- use_F14C_space:

  If `TRUE` (default) the calculations within the function are carried
  out in F\\{}^{14}\\C space. If `FALSE` they are carried out in
  \\{}^{14}\\C age space. We recommend selecting `TRUE` as, for very old
  samples, calibrating in F\\{}^{14}\\C space removes the potential
  affect of asymmetry in the radiocarbon age uncertainty. *Note:* This
  flag can be set independently of the format/scale on which
  `rc_determinations` were originally provided.

- show_progress:

  Whether to show a progress bar in the console during execution.
  Default is `TRUE`.

- calendar_age_range:

  (Optional) Overall minimum and maximum calendar ages (in cal yr BP)
  permitted for the set of radiocarbon samples, i.e.,
  `calendar_age_range[1]` \< \\\theta_1, \ldots, \theta_n\\ \<
  `calendar_age_range[1]`. This is used to bound the start and end of
  the Poisson process (so no events will be permitted to occur outside
  this interval). If not selected then automated selection will be made
  based on given `rc_determinations` and value of
  `bounding_range_prob_cutoff`

- calendar_grid_resolution:

  The spacing of the calendar age grid on which to restrict the
  potential calendar ages of the samples, e.g., calendar ages of samples
  are limited to being one of
  `t, t + resolution, t + 2 * resolution, ...` Default is 1 (i.e., all
  calendar years in the overall calendar range are considered).
  Primarily used to speed-up code if have large range, when may select
  larger resolution.

- prior_h_shape, prior_h_rate:

  (Optional) Prior for the value of the Poisson Process rate (the height
  `rate_h`) in any specific interval: \$\$\textrm{rate}\\\textrm{h} \sim
  \textrm{Gamma}( \textrm{shape} =
  \textrm{prior}\\\textrm{h}\\\textrm{shape}, \textrm{rate} =
  \textrm{prior}\\\textrm{h}\\\textrm{rate}).\$\$ If they are both `NA`
  then `prior_h_shape` is selected to be 1 (so `rate_h` follows an
  Exponential distribution) and `prior_h_rate` chosen adaptively
  (internally) to match `n_observations`.

- prior_n_internal_changepoints_lambda:

  Prior mean for the number of internal changepoints in the rate
  \\\lambda(t)\\.
  \$\$\textrm{n}\\\textrm{internal}\\\textrm{changepoints} \sim
  \textrm{Po}(\textrm{prior}\\\textrm{n}\\\textrm{internal}\\\textrm{changepoints}\\\textrm{lambda})\$\$

- k_max_internal_changepoints:

  Maximum permitted number of internal changepoints

- rescale_factor_rev_jump:

  Factor weighting probability of dimension change in the reversible
  jump update step for Poisson process `rate_h` and `rate_s`

- bounding_range_prob_cutoff:

  Probability cut-off for choosing the bounds for the potential calendar
  ages for the observations

- initial_n_internal_changepoints:

  Number of internal changepoints to initialise MCMC sampler with. The
  default is 10 (so initialise with diffuse state). Will place these
  initial changepoints uniformly at random within overall calendar age
  range.

- grid_extension_factor:

  If you adaptively select the `calendar_age_range` from
  `rc_determinations`, how far you wish to extend the grid beyond this
  adaptive minimum and maximum. The final range will be extended
  (equally on both sides) to cover
  `(1 + grid_extension_factor) * (calendar_age_range)`

- use_fast, fast_approx_prob_cutoff:

  A flag to allow trimming the calendar age likelihood (i.e., reducing
  the range of calendar ages) for each individual sample to speed up the
  sampler. If `TRUE` (default), for each individual sample, those tail
  calendar ages (in the overall `calendar_age_grid`) with very small
  likelihoods will be trimmed (speeding up the updating of the calendar
  ages). If `TRUE` the probability cut-off to remove the tails is
  `fast_approx_prob_cutoff`.

## Value

A list with 7 items. The first 4 items contain output of the model, each
of which has one dimension of size \\n\_{\textrm{out}} = \textrm{floor}(
n\_{\textrm{iter}}/n\_{\textrm{thin}}) + 1\\. The rows in these items
store the state of the MCMC from every
\\n\_{\textrm{thin}}\\\\{}^\textrm{th}\\ iteration:

- `rate_s`:

  A list of length \\n\_{\textrm{out}}\\ each entry giving the current
  set of (calendar age) changepoint locations in the piecewise-constant
  rate \\\lambda(t)\\.

- `rate_h`:

  A list of length \\n\_{\textrm{out}}\\ each entry giving the current
  set of heights (values for the rate) in each piecewise-constant
  section of \\\lambda(t)\\.

- `calendar_ages`:

  An \\n\_{\textrm{out}}\\ by \\n\_{\textrm{obs}}\\ matrix. Gives the
  current estimate for the calendar age of each individual observation.

- `n_internal_changes`:

  A vector of length \\n\_{\textrm{out}}\\ giving the current number of
  internal changes in the value of \\\lambda(t)\\.

where \\n\_{\textrm{obs}}\\ is the number of radiocarbon observations,
i.e., the length of `rc_determinations`.

The remaining items give information about input data, input parameters
(or those calculated) and update_type

- `update_type`:

  A string that always has the value "RJPP".

- `input_data`:

  A list containing the \\{}^{14}\\C data used, and the name of the
  calibration curve used.

- `input_parameters`:

  A list containing the values of the fixed parameters
  `pp_cal_age_range`, `prior_n_internal_changepoints_lambda`,
  `k_max_internal_changepoints`, `prior_h_shape`, `prior_h_rate`,
  `rescale_factor_rev_jump`, `calendar_age_grid`,
  `calendar_grid_resolution`, `n_iter` and `n_thin`.

## See also

See
[PPcalibrate](https://tjheaton.github.io/carbondate/reference/PPcalibrate.md)
for summarisation of samples against a single calibration curve.

Also
[PlotPosteriorMeanRate](https://tjheaton.github.io/carbondate/reference/PlotPosteriorMeanRate.md),
[PlotNumberOfInternalChanges](https://tjheaton.github.io/carbondate/reference/PlotNumberOfInternalChanges.md),
[PlotPosteriorChangePoints](https://tjheaton.github.io/carbondate/reference/PlotPosteriorChangePoints.md)
and
[PlotPosteriorHeights](https://tjheaton.github.io/carbondate/reference/PlotPosteriorHeights.md)
for plotting of the results; and
[FindPosteriorMeanRate](https://tjheaton.github.io/carbondate/reference/FindPosteriorMeanRate.md)
if you just want the posterior mean rate without a plot.

## Examples

``` r
# NOTE: This example is shown with a small n_iter to speed up execution.
# When you run ensure n_iter gives convergence (try function default).

pp_output_mixed <- PPcalibrateMixedCurves(
    pp_uniform_phase_mixed$c14_age,
    pp_uniform_phase_mixed$c14_sig,
    pp_uniform_phase_mixed$sample_source,
    pp_uniform_phase_mixed$delta_r,
    pp_uniform_phase_mixed$delta_r_sig,
    n_iter = 100,
    show_progress = FALSE)
```
