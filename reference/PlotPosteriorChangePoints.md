# Plot Calendar Ages of Changes in Rate of Sample Occurrence for Poisson Process Model

Given output from the Poisson process fitting function
[PPcalibrate](https://tjheaton.github.io/carbondate/reference/PPcalibrate.md),
plot the posterior density estimates for the calendar ages at which
there are internal changepoints in the rate of sample occurrence
\\\lambda(t)\\. These density estimates are calculated **conditional**
upon the number of internal changepoints within the period under study
(which is specified as an input to the function).

Having conditioned on the number of changes, `n_change`, the code will
extract all realisations from the the posterior of the MCMC sampler
which have that number of internal changepoints in the estimate of
\\\lambda(t)\\. It will then provide density estimates for the (ordered)
calendar ages of those internal changepoints. These density estimates
are obtained using a Gaussian kernel.

**Note: These graphs will become harder to interpret as the specified
number of changepoints increases**

For more information read the vignette:  
[`vignette("Poisson-process-modelling", package = "carbondate")`](https://tjheaton.github.io/carbondate/articles/Poisson-process-modelling.md)

## Usage

``` r
PlotPosteriorChangePoints(
  output_data,
  n_changes = c(1, 2, 3),
  plot_cal_age_scale = "BP",
  n_burn = NA,
  n_end = NA,
  kernel_bandwidth = NA,
  plot_lwd = 2
)
```

## Arguments

- output_data:

  The return value from the updating function
  [PPcalibrate](https://tjheaton.github.io/carbondate/reference/PPcalibrate.md).
  Optionally, the output data can have an extra list item named `label`
  which is used to set the label on the plot legend.

- n_changes:

  Number of internal changepoints to condition on, and plot for. A
  vector which can contain at most 4 elements, with values in the range
  1 to 6. If not given, then `c(1, 2, 3)` will be used.

- plot_cal_age_scale:

  (Optional) The calendar scale to use for the x-axis. Allowed values
  are "BP", "AD" and "BC". The default is "BP" corresponding to plotting
  in cal yr BP.

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

- kernel_bandwidth:

  (Optional) The bandwidth used for the (Gaussian) kernel smoothing of
  the calendar age densities. If not given, then 1/50th of the overall
  calendar age range will be used.

- plot_lwd:

  The line width to use when plotting the posterior mean (and confidence
  intervals). Default is 2 (to add emphasis).

## Value

A list with single element `plot_par` that contains the
plotting/graphical parameters of the plot to allow for
editing/annotation.

## See also

For annotating the plot, see
[AddTextPlot](https://tjheaton.github.io/carbondate/reference/AddTextPlot.md),
[AddLinePlot](https://tjheaton.github.io/carbondate/reference/AddLinePlot.md)
and
[AddShadingPlot](https://tjheaton.github.io/carbondate/reference/AddShadingPlot.md)

## Examples

``` r
# NOTE: This example is shown with a small n_iter to speed up execution.
# Try n_iter and n_posterior_samples as the function defaults.

pp_output <- PPcalibrate(
    pp_uniform_phase$c14_age,
    pp_uniform_phase$c14_sig,
    intcal20,
    n_iter = 1000,
    show_progress = FALSE)

# Plot the posterior change points for only 2 or 3 internal changes
PlotPosteriorChangePoints(pp_output, n_changes = c(2, 3))


# Changing the calendar age plotting scale to cal AD
PlotPosteriorChangePoints(pp_output, n_changes = c(2, 3),
    plot_cal_age_scale = "AD")



# Creating a changepoint plot and then adding annotations
posterior_changepoint_plot <- PlotPosteriorChangePoints(pp_output, n_changes = c(2, 3))
# Note: Assigning plot to a variable (with <-) is only needed for annotation.

# Add annotations to plot
AddShadingPlot(posterior_changepoint_plot,
    x_start = 640, x_end = 620,
    col = "red")

AddLinePlot(
     posterior_changepoint_plot,
     v = 600,
     col = "purple",
     lwd = 1,
     lty = 2)

AddTextPlot(posterior_changepoint_plot,
    x = 600, y = 0.04,
    labels = expression(paste("600 cal yrs BP")),
    cex = 0.7,
    pos = 4,
    offset = 0.2,
    col = "purple")

```
