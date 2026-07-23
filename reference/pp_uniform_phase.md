# Example artificial data - Uniform Phase

40 simulated radiocarbon determinations for which the underlying
calendar ages are drawn (uniformly at random) from the period 550–500
cal yr BP. \$\$f(\theta) = U\[550, 500\]\$\$ The observational
uncertainty of each determination is set to be 15 \\{}^{14}\\C yrs.  
  
The corresponding \\{}^{14}\\C ages are then simulated based upon the
IntCal20 calibration curve (convolved with the 15 \\{}^{14}\\C yr
measurement uncertainty): \$\$X_i \| \theta_i \sim N(m(\theta_i),
\rho(\theta_i)^2 + 15^2),\$\$ where \\m(\theta_i)\\ and
\\\rho(\theta_i)\\ are the IntCal20 pointwise means and uncertainties.  
  
This dataset matches that used in the package vignette to illustrate the
Poisson process modelling.

## Usage

``` r
pp_uniform_phase
```

## Format

### `pp_uniform_phase`

A data frame with 40 rows and 4 columns:

- c14_age:

  The simulated \\{}^{14}\\C age (in \\{}^{14}\\C yr BP)

- c14_sig:

  The (fixed) \\{}^{14}\\C age measurement uncertainty used in the
  simulation (set at 15 \\{}^{14}\\C yrs)

- f14c:

  The corresponding simulated values of F\\{}^{14}\\C concentration

- f14c_sig:

  The (fixed) corresponding F\\{}^{14}\\C measurement uncertainty used
  in the simulation
