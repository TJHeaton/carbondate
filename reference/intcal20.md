# IntCal20 calibration curve

The IntCal20 Northern Hemisphere radiocarbon age calibration curve on a
calendar grid spanning from 55,000–0 cal yr BP (Before Present, 0 cal yr
BP corresponds to 1950 CE).  
  
*Note:* This dataset provides \\{}^{14}\\C ages and F\\{}^{14}\\C values
on a calendar age grid. This is a different format from oxcal/calib .14c
files which give the \\{}^{14}\\C ages and \\{\Delta}^{14}\\C values.  
  
**Reference:**  
Reimer PJ, Austin WEN, Bard E, Bayliss A, Blackwell PG, Bronk Ramsey C,
Butzin M, Cheng H, Edwards RL, Friedrich M, Grootes PM, Guilderson TP,
Hajdas I, Heaton TJ, Hogg AG, Hughen KA, Kromer B, Manning SW, Muscheler
R, Palmer JG, Pearson C, van der Plicht J, Reimer RW, Richards DA, Scott
EM, Southon JR, Turney CSM, Wacker L, Adolphi F, Büntgen U, Capano M,
Fahrni S, Fogtmann-Schulz A, Friedrich R, Köhler P, Kudsk S, Miyake F,
Olsen J, Reinig F, Sakamoto M, Sookdeo A, Talamo S. 2020. The IntCal20
Northern Hemisphere radiocarbon age calibration curve (0-55 cal kBP).
*Radiocarbon* **62** https://doi.org/10.1017/RDC.2020.41.

## Usage

``` r
intcal20
```

## Format

### `intcal20`

A data frame with 9,501 rows and 5 columns providing the IntCal20
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

http://doi.org/10.1017/RDC.2020.41
