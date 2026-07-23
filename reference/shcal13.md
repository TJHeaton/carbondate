# SHCal13 calibration curve

The SHCal13 Southern Hemisphere radiocarbon age calibration curve on a
calendar grid spanning from 50,000–0 cal yr BP (Before Present, 0 cal yr
BP corresponds to 1950 CE).  
  
*Note:* This dataset provides \\{}^{14}\\C ages and F\\{}^{14}\\C values
on a calendar age grid. This is a different format from oxcal/calib .14c
files which give the \\{}^{14}\\C ages and \\{\Delta}^{14}\\C values.  
  
**Reference:**  
Alan G Hogg, Quan Hua, Paul G Blackwell, Caitlin E Buck, Thomas P
Guilderson, Timothy J Heaton, Mu Niu, Jonathan G Palmer, Paula J Reimer,
Ron W Reimer, Christian S M Turney, Susan R H Zimmerman. 2013. SHCal13
Southern Hemisphere Calibration, 0-50,000 Years cal BP. *Radiocarbon*
**55**(4):1889–1903 https://doi.org/10.2458/azu_js_rc.55.16783.  
  

## Usage

``` r
shcal13
```

## Format

### `shcal13`

A data frame with 5,141 rows and 5 columns providing the SHCal13
radiocarbon age calibration curve on a calendar grid spanning from
50,000–0 cal yr BP:

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

http://doi.org/10.2458/azu_js_rc.55.16783
