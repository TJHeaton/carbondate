# SHCal20 calibration curve

The SHCal20 Southern Hemisphere radiocarbon age calibration curve on a
calendar grid spanning from 55,000–0 cal yr BP (Before Present, 0 cal yr
BP corresponds to 1950 CE).  
  
*Note:* This dataset provides \\{}^{14}\\C ages and F\\{}^{14}\\C values
on a calendar age grid. This is a different format from oxcal/calib .14c
files which give the \\{}^{14}\\C ages and \\{\Delta}^{14}\\C values.  
  
**Reference:**  
Hogg AG, Heaton TJ, Hua Q, Palmer JG, Turney CSM, Southon J, Bayliss A,
Blackwell PG, Boswijk G, Bronk Ramsey C, Pearson C, Petchey F, Reimer P,
Reimer R, Wacker L. 2020. SHCal20 Southern Hemisphere calibration,
0–55,000 years cal BP. *Radiocarbon* **62**
https://doi.org/10.1017/RDC.2020.59.  
  

## Usage

``` r
shcal20
```

## Format

### `shcal20`

A data frame with 9,501 rows and 5 columns providing the SHCal20
radiocarbon age calibration curve on a calendar grid spanning from
55,000–0 cal yr BP:

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

http://doi.org/10.1017/RDC.2020.59
