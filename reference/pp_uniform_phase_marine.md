# Example artificial marine data - Uniform Phase

40 simulated marine radiocarbon determinations that are created as
though they were marine samples (with a \\\Delta R\\ offset of 100
\\\pm\\ 20 \\{}^{14}\\C yr). As for
[pp_uniform_phase](https://tjheaton.github.io/carbondate/reference/pp_uniform_phase.md),
the samples are simulated with underlying calendar ages drawn (uniformly
at random) from the period 550–500 cal yr BP \$\$f(\theta) = U\[550,
500\]\$\$ The corresponding \\{}^{14}\\C ages are then simulated based
upon the Marine20 calibration curve with a \\\Delta R\\ of 100 \\\pm\\
20 \\{}^{14}\\C yr. The observational uncertainty of each determination
is set to be 15 \\{}^{14}\\C yrs.  
  
The corresponding \\{}^{14}\\C ages are then simulated based upon the
Marine20 calibration curve (convolved with the 15 \\{}^{14}\\C yr
measurement uncertainty): \$\$X_i \| \theta_i \sim N(m(\theta_i),
\rho(\theta_i)^2 + 15^2),\$\$ where \\m(\theta_i)\\ and
\\\rho(\theta_i)\\ are the calibration curve pointwise means and
uncertainties. The simulated \\\Delta R\\ of 100 \\\pm\\ 20 \\{}^{14}\\C
yr is then added.  
  

## Usage

``` r
pp_uniform_phase_marine
```

## Format

### `pp_uniform_phase_marine`

A data frame with 40 rows and 7 columns:

- c14_age:

  The simulated \\{}^{14}\\C age (in \\{}^{14}\\C yr BP)

- c14_sig:

  The (fixed) \\{}^{14}\\C age measurement uncertainty used in the
  simulation (set at 15 \\{}^{14}\\C yrs)

- sample_source:

  A character vector to clarify that the data were simultaed using the
  Marine20 calibration curve

- delta_r:

  The \\\Delta R\\ offset, set at 100 \\{}^{14}\\C yrs for marine
  samples

- delta_r_sig:

  The (1\\\sigma\\) uncertainty in the \\\Delta R\\ offset. Set at 20
  \\{}^{14}\\C yrs

- f14c:

  The corresponding simulated values of F\\{}^{14}\\C concentration

- f14c_sig:

  The (fixed) corresponding F\\{}^{14}\\C measurement uncertainty used
  in the simulation
