#' IntCal20 calibration curve
#'
#' The IntCal20 Northern Hemisphere radiocarbon age calibration curve
#' on a calendar grid spanning from 55,000--0 cal yr BP
#' (Before Present, 0 cal yr BP corresponds to 1950 CE). \cr \cr
#' \emph{Note:} This dataset provides \eqn{{}^{14}}C ages and F\eqn{{}^{14}}C values
#' on a calendar age grid. This is a different format from oxcal/calib .14c files
#' which give the \eqn{{}^{14}}C ages and \eqn{{\Delta}^{14}}C values.\cr \cr
#' \strong{Reference:} \cr Reimer PJ, Austin WEN, Bard E, Bayliss A, Blackwell PG, Bronk Ramsey C,
#' Butzin M, Cheng H, Edwards RL, Friedrich M, Grootes PM, Guilderson TP,
#' Hajdas I, Heaton TJ, Hogg AG, Hughen KA, Kromer B, Manning SW, Muscheler R,
#' Palmer JG, Pearson C, van der Plicht J, Reimer RW, Richards DA, Scott EM,
#' Southon JR, Turney CSM, Wacker L, Adolphi F, Büntgen U, Capano M, Fahrni S,
#' Fogtmann-Schulz A, Friedrich R, Köhler P, Kudsk S, Miyake F, Olsen J,
#' Reinig F, Sakamoto M, Sookdeo A, Talamo S. 2020.
#' The IntCal20 Northern Hemisphere radiocarbon age calibration curve
#' (0-55 cal kBP). \emph{Radiocarbon} \strong{62} https://doi.org/10.1017/RDC.2020.41.
#'
#' @format ## `intcal20`
#' A data frame with 9,501 rows and 5 columns providing the IntCal20 radiocarbon age
#' calibration curve on a calendar grid spanning from 55,000--0 cal yr BP:
#' \describe{
#'   \item{calendar_age}{The calendar age (in cal yr BP)}
#'   \item{c14_age}{The \eqn{{}^{14}}C age (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The (1\eqn{\sigma}) uncertainty in the \eqn{{}^{14}}C age}
#'   \item{f14c}{The \eqn{{}^{14}}C age expressed as F\eqn{{}^{14}}C concentration}
#'   \item{f14c_sig}{The (1\eqn{\sigma}) uncertainty in the F\eqn{{}^{14}}C concentration}
#' }
#' @source http://doi.org/10.1017/RDC.2020.41
"intcal20"

#' IntCal13 calibration curve
#'
#' The IntCal13 Northern Hemisphere radiocarbon age calibration curve
#' on a calendar grid spanning from 50,000--0 cal yr BP
#' (Before Present, 0 cal yr BP corresponds to 1950 CE). \cr \cr
#' \emph{Note:} This dataset provides \eqn{{}^{14}}C ages and F\eqn{{}^{14}}C values
#' on a calendar age grid. This is a different format from oxcal/calib .14c files
#' which give the \eqn{{}^{14}}C ages and \eqn{{\Delta}^{14}}C values.\cr \cr
#' \strong{Reference:} \cr Reimer PJ, Bard E, Bayliss A, Beck JW, Blackwell PG, Bronk Ramsey C, Buck CE,
#' Cheng H, Edwards RL, Friedrich M, Grootes PM, Guilderson TP, Haflidason H,
#' Hajdas I, Hatt? C, Heaton TJ, Hogg AG, Hughen KA, Kaiser KF, Kromer B,
#' Manning SW, Niu M, Reimer RW, Richards DA, Scott EM, Southon JR, Turney CSM,
#' van der Plicht J. 2013.
#' IntCal13 and Marine13 radiocarbon age calibration curves 0--50000 years calBP.
#' \emph{Radiocarbon} \strong{55}(4) https://doi.org/10.2458/azu_js_rc.55.16947. \cr \cr
#'
#' @format ## `intcal13`
#' A data frame with 5,141 rows and 5 columns providing the IntCal13 radiocarbon age
#' calibration curve on a calendar grid spanning from 50,000--0 cal yr BP:
#' \describe{
#'   \item{calendar_age}{The calendar age (in cal yr BP)}
#'   \item{c14_age}{The \eqn{{}^{14}}C age (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The (1\eqn{\sigma}) uncertainty in the \eqn{{}^{14}}C age}
#'   \item{f14c}{The \eqn{{}^{14}}C age expressed as F\eqn{{}^{14}}C concentration}
#'   \item{f14c_sig}{The (1\eqn{\sigma}) uncertainty in the F\eqn{{}^{14}}C concentration}
#' }
#' @source http://doi.org/10.2458/azu_js_rc.55.16947
"intcal13"

#' IntCal09 calibration curve
#'
#' The IntCal09 Northern Hemisphere radiocarbon age calibration curve
#' on a calendar grid spanning from 50,000--0 cal yr BP
#' (Before Present, 0 cal yr BP corresponds to 1950 CE). \cr \cr
#' \emph{Note:} This dataset provides \eqn{{}^{14}}C ages and F\eqn{{}^{14}}C values
#' on a calendar age grid. This is a different format from oxcal/calib .14c files
#' which give the \eqn{{}^{14}}C ages and \eqn{{\Delta}^{14}}C values.\cr \cr
#' \strong{Reference:} \cr PJ Reimer, MGL Baillie, E Bard, A Bayliss, JW Beck, PG Blackwell,
#' C Bronk Ramsey, CE Buck, GS Burr, RL Edwards, M Friedrich, PM Grootes,
#' TP Guilderson, I Hajdas, TJ Heaton, AG Hogg, KA Hughen, KF Kaiser, B Kromer,
#' FG McCormac, SW Manning, RW Reimer, DA Richards, JR Southon, S Talamo,
#' CSM Turney, J van der Plicht, CE Weyhenmeyer. 2009.
#' IntCal09 and Marine09 Radiocarbon Age Calibration Curves, 0--50,000 Years cal BP
#' \emph{Radiocarbon} \strong{51}(4):1111-1150 https://doi.org/10.1017/S0033822200034202. \cr \cr
#'
#' @format ## `intcal09`
#' A data frame with 3,521 rows and 5 columns providing the IntCal09 radiocarbon age
#' calibration curve on a calendar grid spanning from 50,000--0 cal yr BP:
#' \describe{
#'   \item{calendar_age}{The calendar age (in cal yr BP)}
#'   \item{c14_age}{The \eqn{{}^{14}}C age (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The (1\eqn{\sigma}) uncertainty in the \eqn{{}^{14}}C age}
#'   \item{f14c}{The \eqn{{}^{14}}C age expressed as F\eqn{{}^{14}}C concentration}
#'   \item{f14c_sig}{The (1\eqn{\sigma}) uncertainty in the F\eqn{{}^{14}}C concentration}
#' }
#' @source http://doi.org/10.1017/S0033822200034202
"intcal09"

#' IntCal04 calibration curve
#'
#' The IntCal04 Northern Hemisphere radiocarbon age calibration curve
#' on a calendar grid spanning from 26,000--0 cal yr BP
#' (Before Present, 0 cal yr BP corresponds to 1950 CE). \cr \cr
#' \emph{Note:} This dataset provides \eqn{{}^{14}}C ages and F\eqn{{}^{14}}C values
#' on a calendar age grid. This is a different format from oxcal/calib .14c files
#' which give the \eqn{{}^{14}}C ages and \eqn{{\Delta}^{14}}C values.\cr \cr
#' \strong{Reference:} \cr PJ Reimer, MGL Baillie, E Bard, A Bayliss, JW Beck, C Bertrand, PG Blackwell,
#' CE Buck, G Burr, KB Cutler, PE Damon, RL Edwards, RG Fairbanks, M Friedrich,
#' TP Guilderson, KA Hughen, B Kromer, FG McCormac, S Manning, C Bronk Ramsey,
#' RW Reimer, S Remmele, JR Southon, M Stuiver, S Talamo, FW Taylor,
#' J van der Plicht, and CE Weyhenmeyer. 2004.
#' Intcal04 Terrestrial Radiocarbon Age Calibration, 0--26 Cal Kyr BP.
#' \emph{Radiocarbon} \strong{46}(3):1029-1058 https://doi.org/10.1017/S0033822200032999. \cr \cr
#'
#' @format ## `intcal04`
#' A data frame with 3,301 rows and 5 columns providing the IntCal04 radiocarbon age
#' calibration curve on a calendar grid spanning from 26,000--0 cal yr BP:
#' \describe{
#'   \item{calendar_age}{The calendar age (in cal yr BP)}
#'   \item{c14_age}{The \eqn{{}^{14}}C age (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The (1\eqn{\sigma}) uncertainty in the \eqn{{}^{14}}C age}
#'   \item{f14c}{The \eqn{{}^{14}}C age expressed as F\eqn{{}^{14}}C concentration}
#'   \item{f14c_sig}{The (1\eqn{\sigma}) uncertainty in the F\eqn{{}^{14}}C concentration}
#' }
#' @source https://doi.org/10.1017/S0033822200032999
"intcal04"

#' IntCal98 calibration curve
#'
#' The IntCal98 Northern Hemisphere radiocarbon age calibration curve
#' on a calendar grid spanning from 24,000--0 cal yr BP
#' (Before Present, 0 cal yr BP corresponds to 1950 CE). \cr \cr
#' \emph{Note:} This dataset provides \eqn{{}^{14}}C ages and F\eqn{{}^{14}}C values
#' on a calendar age grid. This is a different format from oxcal/calib .14c files
#' which give the \eqn{{}^{14}}C ages and \eqn{{\Delta}^{14}}C values.\cr \cr
#' \strong{Reference:} \cr M. Stuiver, P. J. Reimer, E. Bard, J. W. Beck, G. S. Burr, K. A. Hughen,
#' B. Kromer, F. G. McCormac, J. v. d. Plicht and M. Spurk. 1998.
#' INTCAL98 Radiocarbon Age Calibration, 24,000--0 cal BP.
#' \emph{Radiocarbon} \strong{40}(3):1041-1083 https://doi.org/10.1017/S0033822200019123. \cr \cr
#'
#' @format ## `intcal98`
#' A data frame with 1,538 rows and 5 columns providing the IntCal98 radiocarbon age
#' calibration curve on a calendar grid spanning from 24,000--0 cal yr BP:
#' \describe{
#'   \item{calendar_age}{The calendar age (in cal yr BP)}
#'   \item{c14_age}{The \eqn{{}^{14}}C age (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The (1\eqn{\sigma}) uncertainty in the \eqn{{}^{14}}C age}
#'   \item{f14c}{The \eqn{{}^{14}}C age expressed as F\eqn{{}^{14}}C concentration}
#'   \item{f14c_sig}{The (1\eqn{\sigma}) uncertainty in the F\eqn{{}^{14}}C concentration}
#' }
#' @source https://doi.org/10.1017/S0033822200019123
"intcal98"

#' SHCal20 calibration curve
#'
#' The SHCal20 Southern Hemisphere radiocarbon age calibration curve
#' on a calendar grid spanning from 55,000--0 cal yr BP
#' (Before Present, 0 cal yr BP corresponds to 1950 CE). \cr \cr
#' \emph{Note:} This dataset provides \eqn{{}^{14}}C ages and F\eqn{{}^{14}}C values
#' on a calendar age grid. This is a different format from oxcal/calib .14c files
#' which give the \eqn{{}^{14}}C ages and \eqn{{\Delta}^{14}}C values.\cr \cr
#' \strong{Reference:} \cr Hogg AG, Heaton TJ, Hua Q, Palmer JG, Turney CSM, Southon J, Bayliss A, Blackwell PG,
#' Boswijk G, Bronk Ramsey C, Pearson C, Petchey F, Reimer P, Reimer R, Wacker L.
#' 2020.
#' SHCal20 Southern Hemisphere calibration, 0--55,000 years cal BP.
#' \emph{Radiocarbon} \strong{62} https://doi.org/10.1017/RDC.2020.59. \cr \cr
#'
#' @format ## `shcal20`
#' A data frame with 9,501 rows and 5 columns providing the SHCal20 radiocarbon age
#' calibration curve on a calendar grid spanning from 55,000--0 cal yr BP:
#' \describe{
#'   \item{calendar_age}{The calendar age (in cal yr BP)}
#'   \item{c14_age}{The \eqn{{}^{14}}C age (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The (1\eqn{\sigma}) uncertainty in the \eqn{{}^{14}}C age}
#'   \item{f14c}{The \eqn{{}^{14}}C age expressed as F\eqn{{}^{14}}C concentration}
#'   \item{f14c_sig}{The (1\eqn{\sigma}) uncertainty in the F\eqn{{}^{14}}C concentration}
#' }
#' @source http://doi.org/10.1017/RDC.2020.59
"shcal20"


#' SHCal13 calibration curve
#'
#' The SHCal13 Southern Hemisphere radiocarbon age calibration curve
#' on a calendar grid spanning from 50,000--0 cal yr BP
#' (Before Present, 0 cal yr BP corresponds to 1950 CE). \cr \cr
#' \emph{Note:} This dataset provides \eqn{{}^{14}}C ages and F\eqn{{}^{14}}C values
#' on a calendar age grid. This is a different format from oxcal/calib .14c files
#' which give the \eqn{{}^{14}}C ages and \eqn{{\Delta}^{14}}C values.\cr \cr
#' \strong{Reference:} \cr Alan G Hogg, Quan Hua, Paul G Blackwell, Caitlin E Buck, Thomas P Guilderson,
#' Timothy J  Heaton, Mu Niu, Jonathan G Palmer, Paula J Reimer, Ron W Reimer,
#' Christian S M Turney, Susan R H Zimmerman. 2013.
#' SHCal13 Southern Hemisphere Calibration, 0-50,000 Years cal BP.
#' \emph{Radiocarbon} \strong{55}(4):1889--1903 https://doi.org/10.2458/azu_js_rc.55.16783. \cr \cr
#'
#' @format ## `shcal13`
#' A data frame with 5,141 rows and 5 columns providing the SHCal13 radiocarbon age
#' calibration curve on a calendar grid spanning from 50,000--0 cal yr BP:
#' \describe{
#'   \item{calendar_age}{The calendar age (in cal yr BP)}
#'   \item{c14_age}{The \eqn{{}^{14}}C age (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The (1\eqn{\sigma}) uncertainty in the \eqn{{}^{14}}C age}
#'   \item{f14c}{The \eqn{{}^{14}}C age expressed as F\eqn{{}^{14}}C concentration}
#'   \item{f14c_sig}{The (1\eqn{\sigma}) uncertainty in the F\eqn{{}^{14}}C concentration}
#' }
#' @source http://doi.org/10.2458/azu_js_rc.55.16783
"shcal13"


#' SHCal04 calibration curve
#'
#' The SHCal04 Southern Hemisphere radiocarbon age calibration curve
#' on a calendar grid spanning from 11,000--0 cal yr BP
#' (Before Present, 0 cal yr BP corresponds to 1950 CE). \cr \cr
#' \emph{Note:} This dataset provides \eqn{{}^{14}}C ages and F\eqn{{}^{14}}C values
#' on a calendar age grid. This is a different format from oxcal/calib .14c files
#' which give the \eqn{{}^{14}}C ages and \eqn{{\Delta}^{14}}C values.\cr \cr
#' \strong{Reference:} \cr FG McCormac, AG Hogg, PG Blackwell, CE Buck, TFG Higham, and PJ Reimer
#' 2004.
#' SHCal04 Southern Hemisphere Calibration 0--11.0 cal kyr BP.
#' \emph{Radiocarbon} \strong{46}(3):1087--1092 https://doi.org/10.1017/S0033822200033014. \cr \cr
#'
#' @format ## `shcal04`
#' A data frame with 2,202 rows and 5 columns providing the SHCal04 radiocarbon age
#' calibration curve on a calendar grid spanning from 11,000--0 cal yr BP:
#' \describe{
#'   \item{calendar_age}{The calendar age (in cal yr BP)}
#'   \item{c14_age}{The \eqn{{}^{14}}C age (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The (1\eqn{\sigma}) uncertainty in the \eqn{{}^{14}}C age}
#'   \item{f14c}{The \eqn{{}^{14}}C age expressed as F\eqn{{}^{14}}C concentration}
#'   \item{f14c_sig}{The (1\eqn{\sigma}) uncertainty in the F\eqn{{}^{14}}C concentration}
#' }
#' @source http://doi.org/10.1017/S0033822200033014
"shcal04"


#' Marine20 calibration curve
#'
#' The Marine20 marine radiocarbon age calibration curve
#' on a calendar grid spanning from 55,000--0 cal yr BP
#' (Before Present, 0 cal yr BP corresponds to 1950 CE). \cr \cr
#' \emph{Note:} This dataset provides \eqn{{}^{14}}C ages and F\eqn{{}^{14}}C values
#' on a calendar age grid. This is a different format from oxcal/calib .14c files
#' which give the \eqn{{}^{14}}C ages and \eqn{{\Delta}^{14}}C values.\cr \cr
#' \strong{Reference:} \cr Heaton TJ, Köhler P, Butzin M, Bard E, Reimer RW,
#'  Austin WEN, Bronk Ramsey C, Grootes PM, Hughen KA, Kromer B, Reimer PJ,
#'  Adkins J, Burke A, Cook MA, Olsen J, Skinner L. 2020 Marine20 - the marine radiocarbon age calibration curve (0–55,000 cal BP).
#' \emph{Radiocarbon} \strong{62} https://doi.org/10.1017/RDC.2020.68. \cr \cr
#'
#' @format ## `marine20`
#' A data frame with 5,501 rows and 5 columns providing the Marine20 radiocarbon age
#' calibration curve on a calendar grid spanning from 55,000--0 cal yr BP:
#' \describe{
#'   \item{calendar_age}{The calendar age (in cal yr BP)}
#'   \item{c14_age}{The \eqn{{}^{14}}C age (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The (1\eqn{\sigma}) uncertainty in the \eqn{{}^{14}}C age}
#'   \item{f14c}{The \eqn{{}^{14}}C age expressed as F\eqn{{}^{14}}C concentration}
#'   \item{f14c_sig}{The (1\eqn{\sigma}) uncertainty in the F\eqn{{}^{14}}C concentration}
#' }
#' @source http://doi.org/10.1017/RDC.2020.68
"marine20"


#' Marine13 calibration curve
#'
#' The Marine13 marine radiocarbon age calibration curve
#' on a calendar grid spanning from 55,000--0 cal yr BP
#' (Before Present, 0 cal yr BP corresponds to 1950 CE). \cr \cr
#' \emph{Note:} This dataset provides \eqn{{}^{14}}C ages and F\eqn{{}^{14}}C values
#' on a calendar age grid. This is a different format from oxcal/calib .14c files
#' which give the \eqn{{}^{14}}C ages and \eqn{{\Delta}^{14}}C values.\cr \cr
#' \strong{Reference:} \cr Reimer PJ, Bard E, Bayliss A, Beck JW, Blackwell PG, Bronk Ramsey C, Buck CE,
#' Cheng H, Edwards RL, Friedrich M, Grootes PM, Guilderson TP, Haflidason H,
#' Hajdas I, Hatt? C, Heaton TJ, Hogg AG, Hughen KA, Kaiser KF, Kromer B,
#' Manning SW, Niu M, Reimer RW, Richards DA, Scott EM, Southon JR, Turney CSM,
#' van der Plicht J. 2013.
#' IntCal13 and Marine13 radiocarbon age calibration curves 0--50000 years calBP.
#' \emph{Radiocarbon} \strong{55}(4) https://doi.org/10.2458/azu_js_rc.55.16947. \cr \cr
#'
#' @format ## `marine13`
#' A data frame with 4,801 rows and 5 columns providing the Marine13 radiocarbon age
#' calibration curve on a calendar grid spanning from 50,000--0 cal yr BP:
#' \describe{
#'   \item{calendar_age}{The calendar age (in cal yr BP)}
#'   \item{c14_age}{The \eqn{{}^{14}}C age (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The (1\eqn{\sigma}) uncertainty in the \eqn{{}^{14}}C age}
#'   \item{f14c}{The \eqn{{}^{14}}C age expressed as F\eqn{{}^{14}}C concentration}
#'   \item{f14c_sig}{The (1\eqn{\sigma}) uncertainty in the F\eqn{{}^{14}}C concentration}
#' }
#' @source http://doi.org/10.2458/azu_js_rc.55.16947
"marine13"


#' Marine09 calibration curve
#'
#' The Marine09 marine radiocarbon age calibration curve
#' on a calendar grid spanning from 50,000--0 cal yr BP
#' (Before Present, 0 cal yr BP corresponds to 1950 CE). \cr \cr
#' \emph{Note:} This dataset provides \eqn{{}^{14}}C ages and F\eqn{{}^{14}}C values
#' on a calendar age grid. This is a different format from oxcal/calib .14c files
#' which give the \eqn{{}^{14}}C ages and \eqn{{\Delta}^{14}}C values.\cr \cr
#' \strong{Reference:} \cr PJ Reimer, MGL Baillie, E Bard, A Bayliss, JW Beck, PG Blackwell,
#' C Bronk Ramsey, CE Buck, GS Burr, RL Edwards, M Friedrich, PM Grootes,
#' TP Guilderson, I Hajdas, TJ Heaton, AG Hogg, KA Hughen, KF Kaiser, B Kromer,
#' FG McCormac, SW Manning, RW Reimer, DA Richards, JR Southon, S Talamo,
#' CSM Turney, J van der Plicht, CE Weyhenmeyer. 2009.
#' IntCal09 and Marine09 Radiocarbon Age Calibration Curves, 0--50,000 Years cal BP
#' \emph{Radiocarbon} \strong{51}(4):1111-1150 https://doi.org/10.1017/S0033822200034202. \cr \cr
#'
#' @format ## `marine09`
#' A data frame with 3,651 rows and 5 columns providing the Marine09 radiocarbon age
#' calibration curve on a calendar grid spanning from 50,000--0 cal yr BP:
#' \describe{
#'   \item{calendar_age}{The calendar age (in cal yr BP)}
#'   \item{c14_age}{The \eqn{{}^{14}}C age (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The (1\eqn{\sigma}) uncertainty in the \eqn{{}^{14}}C age}
#'   \item{f14c}{The \eqn{{}^{14}}C age expressed as F\eqn{{}^{14}}C concentration}
#'   \item{f14c_sig}{The (1\eqn{\sigma}) uncertainty in the F\eqn{{}^{14}}C concentration}
#' }
#' @source http://doi.org/10.1017/S0033822200034202
"marine09"


#' Marine04 calibration curve
#'
#' The Marine04 marine radiocarbon age calibration curve
#' on a calendar grid spanning from 50,000--0 cal yr BP
#' (Before Present, 0 cal yr BP corresponds to 1950 CE). \cr \cr
#' \emph{Note:} This dataset provides \eqn{{}^{14}}C ages and F\eqn{{}^{14}}C values
#' on a calendar age grid. This is a different format from oxcal/calib .14c files
#' which give the \eqn{{}^{14}}C ages and \eqn{{\Delta}^{14}}C values.\cr \cr
#' \strong{Reference:} \cr Hughen KA, Baillie MGL, Bard E, Beck JW,
#'  Bertrand CJH, Blackwell PG, Buck CE, Burr GS, Cutler KB, Damon PE,
#'  Edwards RL, Fairbanks RG, Friedrich M, Guilderson TP, Kromer B,
#'  McCormac G, Manning S, Bronk Ramsey C, Reimer PJ, Reimer RW,
#'  Remmele S, Southon JR, Stuiver M, Talamo S, Taylor FW,
#'  van der Plicht J, Weyhenmeyer CE. 2004. Marine04 marine radiocarbon age calibration, 0-26 cal kyr BP.
#' \emph{Radiocarbon} \strong{46}(3):1059-1086 https://doi.org/10.1017/s0033822200033002. \cr \cr
#'
#' @format ## `marine04`
#' A data frame with 3,301 rows and 5 columns providing the Marine04 radiocarbon age
#' calibration curve on a calendar grid spanning from 50,000--0 cal yr BP:
#' \describe{
#'   \item{calendar_age}{The calendar age (in cal yr BP)}
#'   \item{c14_age}{The \eqn{{}^{14}}C age (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The (1\eqn{\sigma}) uncertainty in the \eqn{{}^{14}}C age}
#'   \item{f14c}{The \eqn{{}^{14}}C age expressed as F\eqn{{}^{14}}C concentration}
#'   \item{f14c_sig}{The (1\eqn{\sigma}) uncertainty in the F\eqn{{}^{14}}C concentration}
#' }
#' @source http://doi.org/10.1017/s0033822200033002
"marine04"


#' Example artificial data - Mixture of Normal Phases
#'
#' 50 simulated radiocarbon determinations for which the underlying calendar ages are
#' drawn from a mixture of two normals:
#' \deqn{f(\theta) = 0.45 N(3500, 200^2) + 0.55 N(5000, 100^2) }
#' i.e., a mixture of a normal centred around 3500 cal yr BP; and another
#' (slightly more concentrated/narrower) normal centred around 5000 cal yr BP. \cr \cr
#' The corresponding 50 radiocarbon ages were then simulated using the IntCal20 calibration curve
#' incorporating both the uncertainty in the calibration curve and a hypothetical measurement
#' uncertainty:
#' \deqn{X_i | \theta_i \sim N(m(\theta_i), \rho(\theta_i)^2 + \sigma_{i,\textrm{lab}}^2),}
#' where \eqn{m(\theta_i)} and \eqn{\rho(\theta_i)} are the IntCal20 pointwise
#' means and uncertainties; and \eqn{\sigma_{i,\textrm{lab}}}, the simulated
#' laboratory measurement uncertainty, was fixed at a common value of 25 \eqn{{}^{14}}C yrs. \cr \cr
#' This dataset is included simply to give some quick-to-run examples.
#'
#' @examples
#' # Plotting calendar age density underlying two_normals
#' # Useful for comparisons against estimation techniques
#' weights_true <- c(0.45, 0.55)
#' cluster_means_true_calBP <- c(3500, 5000)
#' cluster_precisions_true <- 1 / c(200, 100)^2
#'
#' # Create mixture density
#' truedens <- function(t, w, truemean, trueprec) {
#'   dens <- 0
#'   for(i in 1:length(w)) {
#'     dens <- dens + w[i] * dnorm(t, mean = truemean[i], sd = 1/sqrt(trueprec[i]))
#'   }
#'   dens
#' }
#'
#' # Visualise mixture
#' curve(truedens(
#'   x,
#'   w = weights_true,
#'   truemean = cluster_means_true_calBP,
#'   trueprec = cluster_precisions_true),
#'   from = 2500, to = 7000, n = 401,
#'   xlim = c(7000, 2500),
#'   xlab = "Calendar Age (cal yr BP)",
#'   ylab = "Density",
#'   col = "red"
#' )
#'
#' @format ## `two_normals`
#' A data frame with 50 rows and 4 columns:
#' \describe{
#'   \item{c14_age}{The simulated \eqn{{}^{14}}C age (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The (fixed) \eqn{{}^{14}}C age measurement uncertainty used in the simulation (set at 25 \eqn{{}^{14}}C yrs)}
#'   \item{f14c}{The corresponding simulated values of F\eqn{{}^{14}}C concentration}
#'   \item{f14c_sig}{The (fixed) corresponding F\eqn{{}^{14}}C measurement uncertainty used in the simulation}
#' }
"two_normals"


#' Example artificial data - Mixture of Normal Phases Using Marine20 calibration curve
#'
#' 50 simulated marine radiocarbon determinations for which the underlying calendar ages are
#' drawn from a mixture of two normals:
#' \deqn{f(\theta) = 0.45 N(3500, 200^2) + 0.55 N(5000, 100^2) }
#' i.e., a mixture of a normal centred around 3500 cal yr BP; and another
#' (slightly more concentrated/narrower) normal centred around 5000 cal yr BP. \cr \cr
#' The corresponding 50 radiocarbon ages were then simulated using the Marine20 calibration curve with a \eqn{\Delta R} of \eqn{-50 \pm 30} \eqn{{}^{14}}C yr.
#' The simulated values incorporate both the uncertainty in the calibration curve and a hypothetical measurement
#' uncertainty:
#' \deqn{X_i | \theta_i \sim N(m(\theta_i), \rho(\theta_i)^2 + \sigma_{i,\textrm{lab}}^2),}
#' where \eqn{m(\theta_i)} and \eqn{\rho(\theta_i)} are the Marine20 pointwise
#' means and uncertainties; and \eqn{\sigma_{i,\textrm{lab}}}, the simulated
#' laboratory measurement uncertainty, was fixed at a common value of 25 \eqn{{}^{14}}C yrs. The simulated \eqn{\Delta R} of \eqn{-50 \pm 30} \eqn{{}^{14}}C yrs is then added. \cr \cr
#' This dataset is included simply to give some quick-to-run examples.
#'
#' @examples
#' # Plotting calendar age density underlying two_normals
#' # Useful for comparisons against estimation techniques
#' weights_true <- c(0.45, 0.55)
#' cluster_means_true_calBP <- c(3500, 5000)
#' cluster_precisions_true <- 1 / c(200, 100)^2
#'
#' # Create mixture density
#' truedens <- function(t, w, truemean, trueprec) {
#'   dens <- 0
#'   for(i in 1:length(w)) {
#'     dens <- dens + w[i] * dnorm(t, mean = truemean[i], sd = 1/sqrt(trueprec[i]))
#'   }
#'   dens
#' }
#'
#' # Visualise mixture
#' curve(truedens(
#'   x,
#'   w = weights_true,
#'   truemean = cluster_means_true_calBP,
#'   trueprec = cluster_precisions_true),
#'   from = 2500, to = 7000, n = 401,
#'   xlim = c(7000, 2500),
#'   xlab = "Calendar Age (cal yr BP)",
#'   ylab = "Density",
#'   col = "red"
#' )
#'
#' @format ## `two_normals_marine`
#' A data frame with 50 rows and 7 columns:
#' \describe{
#'   \item{c14_age}{The simulated \eqn{{}^{14}}C age (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The (fixed) \eqn{{}^{14}}C age measurement uncertainty used in the simulation (set at 25 \eqn{{}^{14}}C yrs)}
#'   \item{sample_source}{A character vector to clarify that the data were simultaed using the Marine20 calibration curve}
#'   \item{delta_r}{The \eqn{\Delta R} offset, set at \eqn{-50} \eqn{{}^{14}}C yrs}
#'   \item{delta_r_sig}{The (1\eqn{\sigma}) uncertainty in the \eqn{\Delta R} offset. Set at 30 \eqn{{}^{14}}C yrs}
#'   \item{f14c}{The corresponding simulated values of F\eqn{{}^{14}}C concentration}
#'   \item{f14c_sig}{The (fixed) corresponding F\eqn{{}^{14}}C measurement uncertainty used in the simulation}
#' }
#'
#'
"two_normals_marine"



#' Example artificial data - Uniform Phase
#'
#' 40 simulated radiocarbon determinations for which the underlying calendar ages are
#' drawn (uniformly at random) from the period  550--500 cal yr BP.
#' \deqn{f(\theta) = U[550, 500]}
#' The observational uncertainty of each determination is set to be 15 \eqn{{}^{14}}C yrs. \cr \cr
#' The corresponding \eqn{{}^{14}}C ages are then simulated based upon the IntCal20 calibration curve
#' (convolved with the 15 \eqn{{}^{14}}C yr measurement uncertainty):
#' \deqn{X_i | \theta_i \sim N(m(\theta_i), \rho(\theta_i)^2 + 15^2),}
#' where \eqn{m(\theta_i)} and \eqn{\rho(\theta_i)} are the IntCal20 pointwise
#' means and uncertainties. \cr \cr
#' This dataset matches that used in the package vignette to illustrate the Poisson process modelling.
#'
#' @format ## `pp_uniform_phase`
#' A data frame with 40 rows and 4 columns:
#' \describe{
#'   \item{c14_age}{The simulated \eqn{{}^{14}}C age (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The (fixed) \eqn{{}^{14}}C age measurement uncertainty used in the simulation (set at 15 \eqn{{}^{14}}C yrs)}
#'   \item{f14c}{The corresponding simulated values of F\eqn{{}^{14}}C concentration}
#'   \item{f14c_sig}{The (fixed) corresponding F\eqn{{}^{14}}C measurement uncertainty used in the simulation}
#' }
"pp_uniform_phase"


#' Example artificial marine data - Uniform Phase
#'
#' 40 simulated marine radiocarbon determinations that are created as though they were marine samples (with a \eqn{\Delta R} offset of 100 \eqn{\pm} 20 \eqn{{}^{14}}C yr). As for [carbondate::pp_uniform_phase], the samples are simulated with underlying calendar ages drawn (uniformly at random) from the period  550--500 cal yr BP
#' \deqn{f(\theta) = U[550, 500]}
#' The corresponding \eqn{{}^{14}}C ages are then simulated based upon the Marine20 calibration curve with a \eqn{\Delta R} of 100 \eqn{\pm} 20 \eqn{{}^{14}}C yr. The observational uncertainty of each determination is set to be 15 \eqn{{}^{14}}C yrs. \cr \cr
#' The corresponding \eqn{{}^{14}}C ages are then simulated based upon the Marine20 calibration curve
#' (convolved with the 15 \eqn{{}^{14}}C yr measurement uncertainty):
#' \deqn{X_i | \theta_i \sim N(m(\theta_i), \rho(\theta_i)^2 + 15^2),}
#' where \eqn{m(\theta_i)} and \eqn{\rho(\theta_i)} are the calibration curve pointwise
#' means and uncertainties. The simulated \eqn{\Delta R} of 100 \eqn{\pm} 20 \eqn{{}^{14}}C yr is then added.  \cr \cr
#'
#' @format ## `pp_uniform_phase_marine`
#' A data frame with 40 rows and 7 columns:
#' \describe{
#'   \item{c14_age}{The simulated \eqn{{}^{14}}C age (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The (fixed) \eqn{{}^{14}}C age measurement uncertainty used in the simulation (set at 15 \eqn{{}^{14}}C yrs)}
#'   \item{sample_source}{A character vector to clarify that the data were simultaed using the Marine20 calibration curve}
#'   \item{delta_r}{The \eqn{\Delta R} offset, set at 100 \eqn{{}^{14}}C yrs for marine samples}
#'   \item{delta_r_sig}{The (1\eqn{\sigma}) uncertainty in the \eqn{\Delta R} offset. Set at 20 \eqn{{}^{14}}C yrs}
#'   \item{f14c}{The corresponding simulated values of F\eqn{{}^{14}}C concentration}
#'   \item{f14c_sig}{The (fixed) corresponding F\eqn{{}^{14}}C measurement uncertainty used in the simulation}
#' }
"pp_uniform_phase_marine"


#' Example artificial data requiring multiple calibration curves - Uniform Phase
#'
#' 40 simulated radiocarbon determinations that are set up to need calibration against multiple curves (and with \eqn{\Delta R} offsets for the marine samples). As for [carbondate::pp_uniform_phase], the samples are simulated with underlying calendar ages
#' drawn (uniformly at random) from the period  550--500 cal yr BP
#' \deqn{f(\theta) = U[550, 500]}
#' however in this case they are assumed to come from different environments - covering the Northern Hemisphere (NH - IntCal20), Southern Hemisphere (HS - SHCal20), and marine/surface ocean (Marine20). The observational uncertainty of each determination is set to be 15 \eqn{{}^{14}}C yrs. \cr \cr
#' The corresponding \eqn{{}^{14}}C ages are then simulated based upon the relevant (environment based) calibration curve
#' (convolved with the 15 \eqn{{}^{14}}C yr measurement uncertainty):
#' \deqn{X_i | \theta_i \sim N(m(\theta_i), \rho(\theta_i)^2 + 15^2),}
#' where \eqn{m(\theta_i)} and \eqn{\rho(\theta_i)} are the calibration curve pointwise
#' means and uncertainties. \cr \cr
#' For the marine samples, we have incorporated into the simulation a \eqn{\Delta R} of 100 \eqn{\pm} 20 \eqn{{}^{14}}C yr.
#'
#' @format ## `pp_uniform_phase_mixed`
#' A data frame with 40 rows and 7 columns:
#' \describe{
#'   \item{c14_age}{The simulated \eqn{{}^{14}}C age (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The (fixed) \eqn{{}^{14}}C age measurement uncertainty used in the simulation (set at 15 \eqn{{}^{14}}C yrs)}
#'   \item{sample_source}{The environmental source of the sample (i.e., NH, SH, or Marine). This determines the calibration curve that needs to be used}
#'   \item{delta_r}{The \eqn{\Delta R} offset to the relevant calibration curve, set at 0 \eqn{{}^{14}}C yrs for all atmospheric (NH or SH) samples and 100 \eqn{{}^{14}}C yrs for marine samples}
#'   \item{delta_r_sig}{The (1\eqn{\sigma}) uncertainty in the \eqn{\Delta R} offset. Set at 20 \eqn{{}^{14}}C yrs for marine samples}
#'   \item{f14c}{The corresponding simulated values of F\eqn{{}^{14}}C concentration}
#'   \item{f14c_sig}{The (fixed) corresponding F\eqn{{}^{14}}C measurement uncertainty used in the simulation}
#' }
"pp_uniform_phase_mixed"


#' Example real-life data - Irish Rath
#'
#' 255 radiocarbon determinations collated by Kerr and McCormick related to the
#' building and use of raths in Ireland in the early-medieval period. \cr \cr
#' \strong{Reference:} \cr
#' Kerr, T., McCormick, F., 2014. Statistics, sunspots and settlement:
#' influences on sum of probability curves.
#' \emph{Journal of Archaeological Science} \strong{41}, 493--501.
#'
#' @format ## `kerr`
#' A data frame with 255 rows and 4 columns:
#' \describe{
#'   \item{c14_age}{The observed \eqn{{}^{14}}C ages of the samples (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The uncertainty in the observed \eqn{{}^{14}}C ages reported by the radiocarbon laboratory}
#'   \item{f14c}{The observed F\eqn{{}^{14}}C concentrations}
#'   \item{f14c_sig}{The uncertainty in the observed F\eqn{{}^{14}}C concentrations reported by the radiocarbon laboratory}
#' }
#' @source http://doi.org/10.1016/j.jas.2013.09.002
"kerr"


#' Example real-life data - Palaeo-Indian demography
#'
#' 628 radiocarbon determinations collated by Buchanan et al. representing the
#' ages of distinct archaeological  sites found across Canada and North America
#' during the time of the palaeoindians. \cr \cr
#' \strong{Reference:} \cr
#' Buchanan, B., Collard, M., Edinborough, K., 2008. Paleoindian demography and
#' the extraterrestrial impact hypothesis. \emph{Proceedings of the National Academy
#' of Sciences} \strong{105}, 11651--11654.
#'
#' @format ## `buchanan`
#' A data frame with 628 rows and 4 columns:
#' \describe{
#'   \item{c14_age}{The observed \eqn{{}^{14}}C ages of the samples (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The uncertainty in the observed \eqn{{}^{14}}C ages reported by the radiocarbon laboratory}
#'   \item{f14c}{The observed F\eqn{{}^{14}}C concentrations}
#'   \item{f14c_sig}{The uncertainty in the observed F\eqn{{}^{14}}C concentrations reported by the radiocarbon laboratory}
#' }
#' @source http://doi.org/10.1073/pnas.0803762105
"buchanan"


#' Example real-life data - Population Decline in Iron Age Ireland
#'
#' 2021 radiocarbon determinations collated by Armit et al. from archaeological
#' groups operating in Ireland, to investigate whether a wetter environment
#' around 2700 cal yr BP led to a population collapse. \cr \cr
#' \strong{Reference:} \cr
#' Armit, I., Swindles, G.T., Becker, K., Plunkett, G., Blaauw, M., 2014. Rapid
#' climate change did not cause population collapse at the end of the European
#' Bronze Age. \emph{Proceedings of the National Academy
#' of Sciences} \strong{111}, 17045--17049.
#'
#' @format ## `armit`
#' A data frame with 2021 rows and 4 columns:
#' \describe{
#'   \item{c14_age}{The observed \eqn{{}^{14}}C ages of the samples (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The uncertainty in the observed \eqn{{}^{14}}C ages reported by the radiocarbon laboratory}
#'   \item{f14c}{The observed F\eqn{{}^{14}}C concentrations}
#'   \item{f14c_sig}{The uncertainty in the observed F\eqn{{}^{14}}C concentrations reported by the radiocarbon laboratory}
#' }
#' @source http://doi.org/10.1073/pnas.1408028111
"armit"


#' Example real-life data - Alces in Yukon and Alaska
#'
#' 58 radiocarbon determinations collated by Dale Guthrie, R. related to
#' Alces (moose) in Yukon and Alaska. Samples
#' are restricted to those between 25,000--6000 \eqn{{}^{14}}C yrs BP. \cr \cr
#' \strong{Reference:} \cr
#' Dale Guthrie, R. New carbon dates link climatic change with human colonization
#' and Pleistocene extinctions. \emph{Nature} \strong{441}, 207--209 (2006).
#' https://doi.org/10.1038/nature04604
#'
#' @format ## `alces`
#' A data frame with 58 rows and 7 columns:
#' \describe{
#'   \item{lab_code}{The sample code for the \eqn{{}^{14}}C laboratory}
#'   \item{site_code}{The site/museum code}
#'   \item{location}{The location of the sample}
#'   \item{c14_age}{The observed \eqn{{}^{14}}C ages of the samples (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The uncertainty in the observed \eqn{{}^{14}}C ages reported by the radiocarbon laboratory}
#'   \item{f14c}{The observed F\eqn{{}^{14}}C concentrations}
#'   \item{f14c_sig}{The uncertainty in the observed F\eqn{{}^{14}}C concentrations reported by the radiocarbon laboratory}
#' }
#' @source https://doi.org/10.1038/nature04604
"alces"


#' Example real-life data - Bison in Yukon and Alaska
#'
#' 64 radiocarbon determinations collated by Dale Guthrie, R. related to
#' Bison in Yukon and Alaska. Samples
#' are restricted to those between 25,000--6000 \eqn{{}^{14}}C yrs BP.\cr \cr
#' \strong{Reference:} \cr
#' Dale Guthrie, R. New carbon dates link climatic change with human colonization
#' and Pleistocene extinctions. \emph{Nature} \strong{441}, 207--209 (2006).
#' https://doi.org/10.1038/nature04604
#'
#' @format ## `bison`
#' A data frame with 64 rows and 7 columns:
#' \describe{
#'   \item{lab_code}{The sample code for the \eqn{{}^{14}}C laboratory}
#'   \item{site_code}{The site/museum code}
#'   \item{location}{The location of the sample}
#'   \item{c14_age}{The observed \eqn{{}^{14}}C ages of the samples (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The uncertainty in the observed \eqn{{}^{14}}C ages reported by the radiocarbon laboratory}
#'   \item{f14c}{The observed F\eqn{{}^{14}}C concentrations}
#'   \item{f14c_sig}{The uncertainty in the observed F\eqn{{}^{14}}C concentrations reported by the radiocarbon laboratory}
#' }
#' @source https://doi.org/10.1038/nature04604
"bison"


#' Example real-life data - Cervus in Yukon and Alaska
#'
#' 63 radiocarbon determinations collated by Dale Guthrie, R. related to
#' Cervus (wapiti) in Yukon and Alaska. Samples
#' are restricted to those between 25,000--6000 \eqn{{}^{14}}C yrs BP. \cr \cr
#' \strong{Reference:} \cr
#' Dale Guthrie, R. New carbon dates link climatic change with human colonization
#' and Pleistocene extinctions. \emph{Nature} \strong{441}, 207--209 (2006).
#' https://doi.org/10.1038/nature04604
#'
#' @format ## `cervus`
#' A data frame with 63 rows and 7 columns:
#' \describe{
#'   \item{lab_code}{The sample code for the \eqn{{}^{14}}C laboratory}
#'   \item{site_code}{The site/museum code}
#'   \item{location}{The location of the sample}
#'   \item{c14_age}{The observed \eqn{{}^{14}}C ages of the samples (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The uncertainty in the observed \eqn{{}^{14}}C ages reported by the radiocarbon laboratory}
#'   \item{f14c}{The observed F\eqn{{}^{14}}C concentrations}
#'   \item{f14c_sig}{The uncertainty in the observed F\eqn{{}^{14}}C concentrations reported by the radiocarbon laboratory}
#' }
#' @source https://doi.org/10.1038/nature04604
"cervus"


#' Example real-life data - Equus in Yukon and Alaska
#'
#' 84 radiocarbon determinations collated by Dale Guthrie, R. related to
#' Equus (horse) in Yukon and Alaska. Samples
#' are restricted to those between 25,000--6000 \eqn{{}^{14}}C yrs BP. \cr \cr
#' \strong{Reference:} \cr
#' Dale Guthrie, R. New carbon dates link climatic change with human colonization
#' and Pleistocene extinctions. \emph{Nature} \strong{441}, 207--209 (2006).
#' https://doi.org/10.1038/nature04604
#'
#' @format ## `equus`
#' A data frame with 84 rows and 7 columns:
#' \describe{
#'   \item{lab_code}{The sample code for the \eqn{{}^{14}}C laboratory}
#'   \item{site_code}{The site/museum code}
#'   \item{location}{The location of the sample}
#'   \item{c14_age}{The observed \eqn{{}^{14}}C ages of the samples (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The uncertainty in the observed \eqn{{}^{14}}C ages reported by the radiocarbon laboratory}
#'   \item{f14c}{The observed F\eqn{{}^{14}}C concentrations}
#'   \item{f14c_sig}{The uncertainty in the observed F\eqn{{}^{14}}C concentrations reported by the radiocarbon laboratory}
#' }
#' @source https://doi.org/10.1038/nature04604
"equus"


#' Example real-life data - Humans in Yukon and Alaska
#'
#' 46 radiocarbon determinations collated by Dale Guthrie, R. related to
#' archaeological sites (i.e., evidence of human existence) in Alaska. Samples
#' are restricted to those between 25,000--6000 \eqn{{}^{14}}C yrs BP. \cr \cr
#' \strong{Reference:} \cr
#' Dale Guthrie, R. New carbon dates link climatic change with human colonization
#' and Pleistocene extinctions. \emph{Nature} \strong{441}, 207--209 (2006).
#' https://doi.org/10.1038/nature04604
#'
#' @format ## `human`
#' A data frame with 46 rows and 7 columns:
#' \describe{
#'   \item{lab_code}{The sample code for the \eqn{{}^{14}}C laboratory}
#'   \item{site_code}{The site/museum code}
#'   \item{location}{The location of the sample}
#'   \item{c14_age}{The observed \eqn{{}^{14}}C ages of the samples (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The uncertainty in the observed \eqn{{}^{14}}C ages reported by the radiocarbon laboratory}
#'   \item{f14c}{The observed F\eqn{{}^{14}}C concentrations}
#'   \item{f14c_sig}{The uncertainty in the observed F\eqn{{}^{14}}C concentrations reported by the radiocarbon laboratory}
#' }
#' @source https://doi.org/10.1038/nature04604
"human"


#' Example real-life data - Mammuthus in Yukon and Alaska
#'
#' 117 radiocarbon determinations collated by Dale Guthrie, R. related to
#' Mammuthus (mammoth) in Yukon and Alaska. Samples
#' are restricted to those between 25,000--6000 \eqn{{}^{14}}C yrs BP. \cr \cr
#' \strong{Reference:} \cr
#' Dale Guthrie, R. New carbon dates link climatic change with human colonization
#' and Pleistocene extinctions. \emph{Nature} \strong{441}, 207--209 (2006).
#' https://doi.org/10.1038/nature04604
#'
#' @format ## `mammuthus`
#' A data frame with 117 rows and 7 columns:
#' \describe{
#'   \item{lab_code}{The sample code for the \eqn{{}^{14}}C laboratory}
#'   \item{site_code}{The site/museum code}
#'   \item{location}{The location of the sample}
#'   \item{c14_age}{The observed \eqn{{}^{14}}C ages of the samples (in \eqn{{}^{14}}C yr BP)}
#'   \item{c14_sig}{The uncertainty in the observed \eqn{{}^{14}}C ages reported by the radiocarbon laboratory}
#'   \item{f14c}{The observed F\eqn{{}^{14}}C concentrations}
#'   \item{f14c_sig}{The uncertainty in the observed F\eqn{{}^{14}}C concentrations reported by the radiocarbon laboratory}
#' }
#' @source https://doi.org/10.1038/nature04604
"mammuthus"
