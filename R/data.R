#' IntCal20 calibration curve
#'
#' Atmospheric data from Reimer et al (2020)
#' Reimer et al. 2020
#' Reimer P, Austin WEN, Bard E, Bayliss A, Blackwell PG, Bronk Ramsey C,
#' Butzin M, Cheng H, Edwards RL, Friedrich M, Grootes PM, Guilderson TP,
#' Hajdas I, Heaton TJ, Hogg AG, Hughen KA, Kromer B, Manning SW, Muscheler R,
#' Palmer JG, Pearson C, van der Plicht J, Reimer RW, Richards DA, Scott EM,
#' Southon JR, Turney CSM, Wacker L, Adolphi F, Büntgen U, Capano M, Fahrni S,
#' Fogtmann-Schulz A, Friedrich R, Köhler P, Kudsk S, Miyake F, Olsen J,
#' Reinig F, Sakamoto M, Sookdeo A, Talamo S. 2020.
#' The IntCal20 Northern Hemisphere radiocarbon age calibration curve
#' (0-55 cal kBP). Radiocarbon 62. doi: 10.1017/RDC.2020.41.
#'
#' @format ## `intcal20`
#' A data frame with 9,501 rows and 5 columns:
#' \describe{
#'   \item{calendar_age}{The calendar age (cal yr BP)}
#'   \item{c14_age}{The 14C age (14C yr BP)}
#'   \item{c14_sig}{The uncertainty in the 14C age reported by the radiocarbon laboratory}
#'   \item{Delta14C}{TODO: What does this describe? Do we need it in the package?}
#'   \item{DeltaSigma}{TODO: What does this describe? Do we need it in the package?}
#' }
#' @source <http://doi.org/10.1017/RDC.2020.41>
"intcal20"