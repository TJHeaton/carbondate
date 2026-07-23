# SHCal04 calibration curve

The SHCal04 Southern Hemisphere radiocarbon age calibration curve on a
calendar grid spanning from 11,000–0 cal yr BP (Before Present, 0 cal yr
BP corresponds to 1950 CE).  
  
*Note:* This dataset provides \\{}^{14}\\C ages and F\\{}^{14}\\C values
on a calendar age grid. This is a different format from oxcal/calib .14c
files which give the \\{}^{14}\\C ages and \\{\Delta}^{14}\\C values.  
  
**Reference:**  
FG McCormac, AG Hogg, PG Blackwell, CE Buck, TFG Higham, and PJ Reimer
2004. SHCal04 Southern Hemisphere Calibration 0–11.0 cal kyr BP.
*Radiocarbon* **46**(3):1087–1092
https://doi.org/10.1017/S0033822200033014.  
  

## Usage

``` r
shcal04
```

## Format

### `shcal04`

A data frame with 2,202 rows and 5 columns providing the SHCal04
radiocarbon age calibration curve on a calendar grid spanning from
11,000–0 cal yr BP:

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

http://doi.org/10.1017/S0033822200033014
