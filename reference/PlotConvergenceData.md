# Plot KL Divergence of Predictive Density to Assess Convergence of Bayesian Non-Parametric DPMM Sampler

This plots the Kullback-Leibler (KL) divergence between a fixed
(initial/baseline) predictive density and the predictive density
calculated from later individual realisations in the MCMC run of one of
the Bayesian non-parametric summarisation approach. The divergence from
the initial predictive density is plotted as a function of the
realisation/iteration number.

This aims to identify when the divergence, from the initial estimate of
the shared \\f(\theta)\\ to the current estimate, has begun to
stabilise. Hence, to (informally) assess when the MCMC chain has
converged to equilibrium for the shared, underlying, predictive
\\f(\theta)\\.

For more information read the vignette:  
[`vignette("determining-convergence", package = "carbondate")`](https://tjheaton.github.io/carbondate/articles/determining-convergence.md)

## Usage

``` r
PlotConvergenceData(output_data, n_initial = NA)
```

## Arguments

- output_data:

  The return value from one of the Bayesian non-parametric DPMM
  summarisation functions, i.e.,
  [PolyaUrnBivarDirichlet](https://tjheaton.github.io/carbondate/reference/PolyaUrnBivarDirichlet.md)
  or
  [WalkerBivarDirichlet](https://tjheaton.github.io/carbondate/reference/WalkerBivarDirichlet.md).

- n_initial:

  The number of (thinned) realisations to use for the 'initial'
  predictive shared density. This predictive density is then compared
  with the predictive obtained at each subsequent realisation in the
  (thinned) DPMM output. If not specified, then the minimum of 1000
  realisations, or 1 / 10 of the total number of realisations, will be
  used.

## Value

None

## Examples

``` r
# Plot results for the example two_normal data
# NOTE: This does not show meaningful results as n_iter
# is too small. Try increasing n_iter to 1e5.
polya_urn_output <- PolyaUrnBivarDirichlet(
    two_normals$c14_age,
    two_normals$c14_sig,
    intcal20,
    n_iter = 500,
    show_progress = FALSE)
PlotConvergenceData(polya_urn_output)
```
