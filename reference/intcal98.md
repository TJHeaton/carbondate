# IntCal98 calibration curve

The IntCal98 Northern Hemisphere radiocarbon age calibration curve on a
calendar grid spanning from 24,000–0 cal yr BP (Before Present, 0 cal yr
BP corresponds to 1950 CE).  
  
*Note:* This dataset provides \\{}^{14}\\C ages and F\\{}^{14}\\C values
on a calendar age grid. This is a different format from oxcal/calib .14c
files which give the \\{}^{14}\\C ages and \\{\Delta}^{14}\\C values.  
  
**Reference:**  
M. Stuiver, P. J. Reimer, E. Bard, J. W. Beck, G. S. Burr, K. A. Hughen,
B. Kromer, F. G. McCormac, J. v. d. Plicht and M. Spurk. 1998. INTCAL98
Radiocarbon Age Calibration, 24,000–0 cal BP. *Radiocarbon*
**40**(3):1041-1083 https://doi.org/10.1017/S0033822200019123.  
  

## Usage

``` r
intcal98
```

## Format

### `intcal98`

A data frame with 1,538 rows and 5 columns providing the IntCal98
radiocarbon age calibration curve on a calendar grid spanning from
24,000–0 cal yr BP:

- calendar_age:

  The calendar age (in cal yr BP)

- c14_age:

  The \\{}^{14}\\C age (in \\{}^{14}\\C yr BP)

- c14_sig:

  The (1\\\sigma\\) uncertainty in the \\{}^{14}\\C age

- f14c:

  The \\{}^{14}\\C age expressed as F\\{}^{14}\\C concentration

- f14c_sig:

  The (1\\\sigma\\) uncertainty in the F\\{}^{14}\\C concentration

## Source

https://doi.org/10.1017/S0033822200019123
