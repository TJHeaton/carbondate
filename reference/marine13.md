# Marine13 calibration curve

The Marine13 marine radiocarbon age calibration curve on a calendar grid
spanning from 55,000–0 cal yr BP (Before Present, 0 cal yr BP
corresponds to 1950 CE).  
  
*Note:* This dataset provides \\{}^{14}\\C ages and F\\{}^{14}\\C values
on a calendar age grid. This is a different format from oxcal/calib .14c
files which give the \\{}^{14}\\C ages and \\{\Delta}^{14}\\C values.  
  
**Reference:**  
Reimer PJ, Bard E, Bayliss A, Beck JW, Blackwell PG, Bronk Ramsey C,
Buck CE, Cheng H, Edwards RL, Friedrich M, Grootes PM, Guilderson TP,
Haflidason H, Hajdas I, Hatt? C, Heaton TJ, Hogg AG, Hughen KA, Kaiser
KF, Kromer B, Manning SW, Niu M, Reimer RW, Richards DA, Scott EM,
Southon JR, Turney CSM, van der Plicht J. 2013. IntCal13 and Marine13
radiocarbon age calibration curves 0–50000 years calBP. *Radiocarbon*
**55**(4) https://doi.org/10.2458/azu_js_rc.55.16947.  
  

## Usage

``` r
marine13
```

## Format

### `marine13`

A data frame with 4,801 rows and 5 columns providing the Marine13
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

http://doi.org/10.2458/azu_js_rc.55.16947
