# Calibrate and Summarise Multiple Radiocarbon Samples via a Bayesian Non-Parametric DPMM (with Polya Urn Updating)

This function calibrates sets of multiple radiocarbon (\\{}^{14}\\C)
determinations, and simultaneously summarises the resultant calendar age
information. This is achieved using Bayesian non-parametric density
estimation, providing a statistically-rigorous alternative to summed
probability distributions (SPDs).

It takes as an input a set of \\{}^{14}\\C determinations and associated
\\1\sigma\\ uncertainties, as well as the radiocarbon age calibration
curve to be used. The samples are assumed to arise from an (unknown)
shared calendar age distribution \\f(\theta)\\ that we would like to
estimate, alongside performing calibration of each sample.

The function models the underlying distribution \\f(\theta)\\ as a
Dirichlet process mixture model (DPMM), whereby the samples are
considered to arise from an unknown number of distinct clusters. Fitting
is achieved via MCMC.

It returns estimates for the calendar age of each individual radiocarbon
sample; and broader output (including the means and variances of the
underpinning calendar age clusters) that can be used by other library
functions to provide a predictive estimate of the shared calendar age
density \\f(\theta)\\.

For more information read the vignette:  
[`vignette("Non-parametric-summed-density", package = "carbondate")`](https://tjheaton.github.io/carbondate/articles/Non-parametric-summed-density.md)

**Note:** The library provides two slightly-different update schemes for
the MCMC. In this particular function, updating of the DPMM is achieved
by a Polya Urn approach (Neal 2000) This is our recommended updating
approach based on testing. The alternative, slice-sampled, approach can
be found at
[WalkerBivarDirichlet](https://tjheaton.github.io/carbondate/reference/WalkerBivarDirichlet.md)

**References:**  
Heaton, TJ. 2022. “Non-parametric Calibration of Multiple Related
Radiocarbon Determinations and their Calendar Age Summarisation.”
*Journal of the Royal Statistical Society Series C: Applied Statistics*
**71** (5):1918-56. https://doi.org/10.1111/rssc.12599.  
Neal, RM. 2000. “Markov Chain Sampling Methods for Dirichlet Process
Mixture Models.” *Journal of Computational and Graphical Statistics*
**9** (2):249 https://doi.org/10.2307/1390653.

## Usage

``` r
PolyaUrnBivarDirichlet(
  rc_determinations,
  rc_sigmas,
  calibration_curve,
  F14C_inputs = FALSE,
  delta_r = NULL,
  delta_r_sig = NULL,
  n_iter = 1e+05,
  n_thin = 10,
  use_F14C_space = TRUE,
  slice_width = NA,
  slice_multiplier = 10,
  n_clust = min(10, length(rc_determinations)),
  show_progress = TRUE,
  sensible_initialisation = TRUE,
  lambda = NA,
  nu1 = NA,
  nu2 = NA,
  A = NA,
  B = NA,
  alpha_shape = NA,
  alpha_rate = NA,
  mu_phi = NA,
  calendar_ages = NA
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

- calibration_curve:

  A dataframe which must contain one column `calendar_age_BP`, and also
  columns `c14_age` and `c14_sig` or `f14c` and `f14c_sig` (or both
  sets). This format matches the curves supplied with this package,
  e.g.,
  [intcal20](https://tjheaton.github.io/carbondate/reference/intcal20.md),
  [intcal13](https://tjheaton.github.io/carbondate/reference/intcal13.md),
  which contain all 5 columns.

- F14C_inputs:

  `TRUE` if the provided `rc_determinations` are F\\{}^{14}\\C
  concentrations and `FALSE` if they are radiocarbon ages. Defaults to
  `FALSE`.

- delta_r, delta_r_sig:

  (Optional) The \\\Delta R\\ offset and associated 1\\\sigma\\
  uncertainty if calibrating a set of marine samples. This offset must
  be a single value that is shared by all the samples.

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

- slice_width:

  Parameter for slice sampling (optional). If not given a value is
  chosen intelligently based on the spread of the initial calendar ages.
  Must be given if `sensible_initialisation` is `FALSE`.

- slice_multiplier:

  Integer parameter for slice sampling (optional). Default is 10. Limits
  the slice size to `slice_multiplier * slice_width`.

- n_clust:

  The number of clusters with which to initialise the sampler
  (optional). Must be less than the length of `rc_determinations`.
  Default is 10 or the length of `rc_determinations` if that is less
  than 10.

- show_progress:

  Whether to show a progress bar in the console during execution.
  Default is `TRUE`.

- sensible_initialisation:

  Whether to use sensible values to initialise the sampler and an
  automated (adaptive) prior on \\\mu\_{\phi}\\ and (A, B) that is
  informed by the observed `rc_determinations`. If this is `TRUE` (the
  recommended default), then all the remaining arguments below are
  ignored.

- lambda, nu1, nu2:

  Hyperparameters for the prior on the mean \\\phi_j\\ and precision
  \\\tau_j\\ of each individual calendar age cluster \\j\\: \$\$(\phi_j,
  \tau_j)\|\mu\_{\phi} \sim \textrm{NormalGamma}(\mu\_{\phi}, \lambda,
  \nu_1, \nu_2)\$\$ where \\\mu\_{\phi}\\ is the overall cluster
  centering. Required if `sensible_initialisation` is `FALSE`.

- A, B:

  Prior on \\\mu\_{\phi}\\ giving the mean and precision of the overall
  centering \\\mu\_{\phi} \sim N(A, B^{-1})\\. Required if
  `sensible_initialisation` is `FALSE`.

- alpha_shape, alpha_rate:

  Shape and rate hyperparameters that specify the prior for the
  Dirichlet Process (DP) concentration, \\\alpha\\. This concentration
  \\\alpha\\ determines the number of clusters we expect to observe
  among our \\n\\ sampled objects. The model places a prior on \\\alpha
  \sim \Gamma(\eta_1, \eta_2)\\, where \\\eta_1, \eta_2\\ are the
  `alpha_shape` and `alpha_rate`. A small \\\alpha\\ means the DPMM is
  more concentrated (i.e. we expect fewer calendar age clusters) while a
  large alpha means it is less less concentrated (i.e. many clusters).
  Required if `sensible_initialisation` is `FALSE`.

- mu_phi:

  Initial value of the overall cluster centering \\\mu\_{\phi}\\.
  Required if `sensible_initialisation` is `FALSE`.

- calendar_ages:

  The initial estimate for the underlying calendar ages (optional). If
  supplied, it must be a vector with the same length as
  `rc_determinations`. Required if `sensible_initialisation` is `FALSE`.

## Value

A list with 10 items. The first 8 items contain output of the model,
each of which has one dimension of size \\n\_{\textrm{out}} =
\textrm{floor}( n\_{\textrm{iter}}/n\_{\textrm{thin}}) + 1\\. The rows
in these items store the state of the MCMC from every
\\n\_{\textrm{thin}}\\\\{}^\textrm{th}\\ iteration:

- `cluster_identifiers`:

  A list of length \\n\_{\textrm{out}}\\ each entry gives the cluster
  allocation (an integer between 1 and `n_clust`) for each observation
  on the relevant MCMC iteration. Information on the state of these
  calendar age clusters (means and precisions) can be found in the other
  output items.

- `alpha`:

  A double vector of length \\n\_{\textrm{out}}\\ giving the Dirichlet
  Process concentration parameter \\\alpha\\.

- `n_clust`:

  An integer vector of length \\n\_{\textrm{out}}\\ giving the current
  number of clusters in the model.

- `phi`:

  A list of length \\n\_{\textrm{out}}\\ each entry giving a vector of
  length `n_clust` of the means of the current calendar age clusters
  \\\phi_j\\.

- `tau`:

  A list of length \\n\_{\textrm{out}}\\ each entry giving a vector of
  length `n_clust` of the precisions of the current calenadar age
  cluster \\\tau_j\\.

- `observations_per_cluster`:

  A list of length \\n\_{\textrm{out}}\\ each entry giving a vector of
  length `n_clust` of the number of observations for that cluster.

- `calendar_ages`:

  An \\n\_{\textrm{out}}\\ by \\n\_{\textrm{obs}}\\ integer matrix.
  Gives the current estimate for the calendar age of each individual
  observation.

- `mu_phi`:

  A vector of length \\n\_{\textrm{out}}\\ giving the overall centering
  \\\mu\_{\phi}\\ of the calendar age clusters.

where \\n\_{\textrm{obs}}\\ is the number of radiocarbon observations
i.e. the length of `rc_determinations`.

The remaining items give information about the input data, input
parameters (or those calculated using `sensible_initialisation`) and the
update_type

- `update_type`:

  A string that always has the value "Polya Urn".

- `input_data`:

  A list containing the \\{}^{14}\\C data used, the name of the
  calibration curve used, and `delta_r` information (if specified).

- `input_parameters`:

  A list containing the values of the fixed hyperparameters `lambda`,
  `nu1`, `nu2`, `A`, `B`, `alpha_shape`, `alpha_rate` and `mu_phi`, and
  the slice parameters `slice_width` and `slice_multiplier`.

## See also

[WalkerBivarDirichlet](https://tjheaton.github.io/carbondate/reference/WalkerBivarDirichlet.md)
for our less-preferred MCMC method to update the Bayesian DPMM
(otherwise an identical model); and
[PlotCalendarAgeDensityIndividualSample](https://tjheaton.github.io/carbondate/reference/PlotCalendarAgeDensityIndividualSample.md),
[PlotPredictiveCalendarAgeDensity](https://tjheaton.github.io/carbondate/reference/PlotPredictiveCalendarAgeDensity.md)
and
[PlotNumberOfClusters](https://tjheaton.github.io/carbondate/reference/PlotNumberOfClusters.md)
to access the model output and estimate the calendar age information.  
  
See also
[PPcalibrate](https://tjheaton.github.io/carbondate/reference/PPcalibrate.md)
for an an alternative (similarly rigorous) approach to calibration and
summarisation of related radiocarbon determinations using a
variable-rate Poisson process

## Examples

``` r
# Note these examples are shown with a small n_iter to speed up execution.
# When you run ensure n_iter gives convergence (try function default).

# Basic usage making use of sensible initialisation to set most values and
# using a saved example data set and the IntCal20 curve.
polya_urn_output <- PolyaUrnBivarDirichlet(
    two_normals$c14_age,
    two_normals$c14_sig,
    intcal20,
    n_iter = 100,
    show_progress = FALSE)

# The radiocarbon determinations can be given as F14C concentrations
polya_urn_output <- PolyaUrnBivarDirichlet(
    two_normals$f14c,
    two_normals$f14c_sig,
    intcal20,
    F14C_inputs = TRUE,
    n_iter = 100,
    show_progress = FALSE)
```
