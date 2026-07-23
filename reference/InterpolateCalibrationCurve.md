# Interpolate a calibration curve at a set of calendar ages

Interpolate a calibration curve at a set of calendar ages

## Usage

``` r
InterpolateCalibrationCurve(
  new_calendar_ages_BP,
  calibration_curve,
  F14C_outputs = NA
)
```

## Arguments

- new_calendar_ages_BP:

  A scalar or vector containing calendar ages (in cal yr BP) at which to
  interpolate the values (both the means and uncertainties) of the given
  calibration curve. If not provided (and `NA` is given), will use the
  range from the minimum calendar age to the maximum calendar age of the
  original calibration curve spaced by 1.

- calibration_curve:

  A dataframe which must contain one column `calendar_age_BP`, and also
  columns `c14_age` and `c14_sig` or `f14c` and `f14c_sig` (or both
  sets). This format matches the curves supplied with this package,
  e.g.,
  [intcal20](https://tjheaton.github.io/carbondate/reference/intcal20.md),
  [intcal13](https://tjheaton.github.io/carbondate/reference/intcal13.md),
  which contain all 5 columns.

- F14C_outputs:

  `TRUE` if only F\\{}^{14}\\C concentrations are required, `FALSE` if
  only the radiocarbon ages (in \\{}^{14}\\C yrs BP) are required and
  `NA` if both are required for the new curve.

## Value

A new dataframe with entries for the interpolated `c14_age`, and
`c14_sig`, `f14c` and `f14c_sig` values at the `calendar_age_BP` values
that were given in `new_calendar_ages_BP`.

## Examples

``` r

# Interpolate intcal20 at a single calendar age. Generates both 14C ages and F14C scales.
InterpolateCalibrationCurve(51020, intcal20)
#>   calendar_age_BP c14_age c14_sig        f14c     f14c_sig
#> 1           51020   48225     413 0.002470435 0.0001270123

# Interpolate intcal20 at two calendar ages. Generates F14C estimates only.
InterpolateCalibrationCurve(c(51017, 51021), intcal20, TRUE)
#>   calendar_age_BP        f14c     f14c_sig
#> 1           51017 0.002470619 0.0001267448
#> 2           51021 0.002470389 0.0001271021

# Interpolate intcal20 at two calendar ages. Generate 14C age estimates (cal yr BP) only.
InterpolateCalibrationCurve(c(51017, 51021), intcal20, FALSE)
#>   calendar_age_BP  c14_age c14_sig
#> 1           51017 48224.40   412.1
#> 2           51021 48225.15   413.3

# Interpolate intcal20 at every integer calendar age within the range of dates
# (for intcal20 this is 0 to 55000 cal yr BP), and create estimates for both radiocarbon scales.
cal_curve <- InterpolateCalibrationCurve(NA, intcal20)
```
