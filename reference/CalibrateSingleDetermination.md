# Calibrate a Single Radiocarbon Determination

Uses the supplied calibration curve to calibrate a single radiocarbon
determination and uncertainty (expressed either in terms of radiocarbon
age, or as an F\\{}^{14}\\C concentration) and obtain its calendar age
probability density estimate.

## Usage

``` r
CalibrateSingleDetermination(
  rc_determination,
  rc_sigma,
  calibration_curve,
  delta_r = NULL,
  delta_r_sig = NULL,
  F14C_inputs = FALSE,
  resolution = 1,
  plot_output = FALSE,
  plot_cal_age_scale = "BP",
  interval_width = "2sigma",
  bespoke_probability = NA,
  denscale = 3,
  plot_pretty = TRUE
)
```

## Arguments

- rc_determination:

  A single observed radiocarbon determination provided either as the
  radiocarbon age (in \\{}^{14}\\C yr BP) or the F\\{}^{14}\\C
  concentration.

- rc_sigma:

  The corresponding measurement uncertainty of the radiocarbon
  determination (must be in the same units as above, i.e., reported as
  \\{}^{14}\\C age or F\\{}^{14}\\C)

- calibration_curve:

  A dataframe which must contain one column `calendar_age_BP`, and also
  columns `c14_age` and `c14_sig` or `f14c` and `f14c_sig` (or both
  sets). This format matches the curves supplied with this package,
  e.g.,
  [intcal20](https://tjheaton.github.io/carbondate/reference/intcal20.md),
  [intcal13](https://tjheaton.github.io/carbondate/reference/intcal13.md),
  which contain all 5 columns.

- delta_r, delta_r_sig:

  For marine curve calibration only. The \\\Delta R\\ offset and
  associated 1\\\sigma\\ uncertainty on the regional offset to the
  marine calibration curve (must be given in \\{}^{14}\\C yrs)

- F14C_inputs:

  `TRUE` if the provided `rc_determination` is an F\\{}^{14}\\C
  concentration and `FALSE` if it is a radiocarbon age. Defaults to
  `FALSE`.

- resolution:

  The distance between the calendar ages at which to calculate the
  calendar age probability. Default is 1.

- plot_output:

  `TRUE` if you wish to plot the determination, the calibration curve,
  and the posterior calibrated age estimate on the same plot. Defaults
  to `FALSE`

- plot_cal_age_scale:

  Only for usage when `plot_output = TRUE`. The calendar scale to use
  for the x-axis. Allowed values are "BP", "AD" and "BC". The default is
  "BP", corresponding to plotting in cal yr BP.

- interval_width:

  Only for usage when `plot_output = TRUE`. The confidence intervals to
  show for the calibration curve and for the highest posterior density
  ranges. Choose from one of "1sigma" (68.3%), "2sigma" (95.4%) and
  "bespoke". Default is "2sigma".

- bespoke_probability:

  The probability to use for the confidence interval if "bespoke" is
  chosen above. E.g. if 0.95 is chosen, then the 95% confidence interval
  is calculated. Ignored if "bespoke" is not chosen.

- denscale:

  Whether to scale the vertical range of the calendar age density plot
  relative to the calibration curve plot (optional). Default is 3 which
  means that the maximum calendar age density will be at 1/3 of the
  height of the plot.

- plot_pretty:

  logical, defaulting to `TRUE`. If set `TRUE` then will select pretty
  plotting margins (that create sufficient space for axis titles and
  rotates y-axis labels). If `FALSE` will implement current user values.

## Value

If `plot = FALSE` then a dataframe with with one column
`calendar_age_BP` containing the calendar ages, and the other column
`probability` containing the probability at that calendar age.

If `plot = TRUE` then returns a list, the first element
`posterior_cal_age` is as above. The second list element, `plot_par`,
contains the plotting/graphical parameters of the plot to allow for
editing/annotation.

## See also

For annotating the plot, see
[AddTextPlot](https://tjheaton.github.io/carbondate/reference/AddTextPlot.md),
[AddLinePlot](https://tjheaton.github.io/carbondate/reference/AddLinePlot.md)
and
[AddShadingPlot](https://tjheaton.github.io/carbondate/reference/AddShadingPlot.md)

## Examples

``` r
# Calibration of a single determination expressed as 14C age BP
calib <- CalibrateSingleDetermination(860, 35, intcal20)
plot(calib, type = "l", xlim = c(1000, 600))


# Incorporating an automated plot to visualise the calibration
CalibrateSingleDetermination(860, 35, intcal20, plot_output = TRUE)


# Calibration of a single (old) determination expressed as 14C age BP
calib <- CalibrateSingleDetermination(31020, 100, intcal20)
plot(calib, type = "l", xlim = c(36500, 34500))


# Calibration of a single (old) determination expressed as F14C concentration
calib <- CalibrateSingleDetermination(
    0.02103493, 0.0002618564, intcal20, F14C_inputs = TRUE)
plot(calib, type = "l", xlim = c(36500, 34500))


# Calibration of a single determination expressed as 14C age BP
# against SHCal20 (and creating an automated plot)
CalibrateSingleDetermination(1413, 25, shcal20, plot_output = TRUE)


# Implementing a bespoke confidence interval level and plot in AD
CalibrateSingleDetermination(
    1413,
    25,
    shcal20,
    plot_output = TRUE,
    plot_cal_age_scale = "AD",
    interval_width = "bespoke",
    bespoke_probability = 0.8)


# Changing denscale (so the calendar age density takes up less space)
CalibrateSingleDetermination(
    1413,
    25,
    shcal20,
    plot_output = TRUE,
    interval_width = "bespoke",
    bespoke_probability = 0.8,
    denscale = 5)


# Annotating a plot
# Assign plot to a variable (with <-):
individual_calibration_plot <- CalibrateSingleDetermination(860, 35, intcal20, plot_output = TRUE)

AddLinePlot(individual_calibration_plot,
    v = 850,
    col = "purple",
    lwd = 1,
    lty = 2)

AddTextPlot(individual_calibration_plot,
    x = 850, y = 750,
    labels = expression(paste("850 cal yrs BP")),
    cex = 0.7,
    pos = 2,
    offset = 0.2,
    col = "purple")
```
