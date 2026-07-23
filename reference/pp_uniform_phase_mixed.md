# Example artificial data requiring multiple calibration curves - Uniform Phase

40 simulated radiocarbon determinations that are set up to need
calibration against multiple curves (and with \\\Delta R\\ offsets for
the marine samples). As for
[pp_uniform_phase](https://tjheaton.github.io/carbondate/reference/pp_uniform_phase.md),
the samples are simulated with underlying calendar ages drawn (uniformly
at random) from the period 550–500 cal yr BP \$\$f(\theta) = U\[550,
500\]\$\$ however in this case they are assumed to come from different
environments - covering the Northern Hemisphere (NH - IntCal20),
Southern Hemisphere (HS - SHCal20), and marine/surface ocean (Marine20).
The observational uncertainty of each determination is set to be 15
\\{}^{14}\\C yrs.  
  
The corresponding \\{}^{14}\\C ages are then simulated based upon the
relevant (environment based) calibration curve (convolved with the 15
\\{}^{14}\\C yr measurement uncertainty): \$\$X_i \| \theta_i \sim
N(m(\theta_i), \rho(\theta_i)^2 + 15^2),\$\$ where \\m(\theta_i)\\ and
\\\rho(\theta_i)\\ are the calibration curve pointwise means and
uncertainties.  
  
For the marine samples, we have incorporated into the simulation a
\\\Delta R\\ of 100 \\\pm\\ 20 \\{}^{14}\\C yr.

## Usage

``` r
pp_uniform_phase_mixed
```

## Format

### `pp_uniform_phase_mixed`

A data frame with 40 rows and 7 columns:

- c14_age:

  The simulated \\{}^{14}\\C age (in \\{}^{14}\\C yr BP)

- c14_sig:

  The (fixed) \\{}^{14}\\C age measurement uncertainty used in the
  simulation (set at 15 \\{}^{14}\\C yrs)

- sample_source:

  The environmental source of the sample (i.e., NH, SH, or Marine). This
  determines the calibration curve that needs to be used

- delta_r:

  The \\\Delta R\\ offset to the relevant calibration curve, set at 0
  \\{}^{14}\\C yrs for all atmospheric (NH or SH) samples and 100
  \\{}^{14}\\C yrs for marine samples

- delta_r_sig:

  The (1\\\sigma\\) uncertainty in the \\\Delta R\\ offset. Set at 20
  \\{}^{14}\\C yrs for marine samples

- f14c:

  The corresponding simulated values of F\\{}^{14}\\C concentration

- f14c_sig:

  The (fixed) corresponding F\\{}^{14}\\C measurement uncertainty used
  in the simulation
