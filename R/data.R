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
#'   \item{calendar_age}{The calendar age (yr BP)}
#'   \item{c14_age}{The 14C age (yr BP)}
#'   \item{c14_sig}{The uncertainty in the 14C age reported by the radiocarbon laboratory}
#'   \item{Delta14C}{TODO: What does this describe? Do we need it in the package?}
#'   \item{DeltaSigma}{TODO: What does this describe? Do we need it in the package?}
#' }
#' @source <http://doi.org/10.1017/RDC.2020.41>
"intcal20"


#' Example real-life data - Irish Rath
#'
#' 255 radiocarbon determinations collated by Kerr and McCormick related to the
#' building and use of raths in Ireland in the early-medieval period. Ref:
#' Kerr, T., McCormick, F., 2014. Statistics, sunspots and settlement:
#' influences on sum of probability curves.
#' Journal of Archaeological Science 41, 493 – 501.
#'
#' @format ## `kerr`
#' A data frame with 255 rows and 2 columns:
#' \describe{
#'   \item{c14_age}{The 14C age (yr BP)}
#'   \item{c14_sig}{The uncertainty in the 14C age reported by the radiocarbon laboratory}
#' }
#' @source <http://doi.org/10.1016/j.jas.2013.09.002>
"kerr"


#' Example real-life data - Palaeo-Indian demography
#'
#' 628 radiocarbon determinations collated by Buchanan et al. representing the
#' ages of distinct archaeological  sites found across Canada and North America
#' during the time of the palaeoindians. Ref:
#' Buchanan, B., Collard, M., Edinborough, K., 2008. Paleoindian demography and
#' the extraterrestrial impact hypothesis. Proceedings of the National Academy
#' of Sciences 105, 11651–11654.
#'
#' @format ## `buchanan`
#' A data frame with 628 rows and 2 columns:
#' \describe{
#'   \item{c14_age}{The 14C age (yr BP)}
#'   \item{c14_sig}{The uncertainty in the 14C age reported by the radiocarbon laboratory}
#' }
#' @source <http://doi.org/10.1073/pnas.0803762105>
"buchanan"


#' Example real-life data - Population Decline in Iron Age Ireland
#'
#' 2021 radiocarbon determinations collated by Armit et al. from archaeological
#' groups operating in Ireland, to investigate whether a wetter environment
#' around 2700 cal yr BP led to a population collapse. Ref:
#' Armit, I., Swindles, G.T., Becker, K., Plunkett, G., Blaauw, M., 2014. Rapid
#' climate change did not cause population collapse at the end of the European
#' Bronze Age. Proceedings of the National Academy of Sciences 111, 17045–17049.
#'
#' @format ## `armit`
#' A data frame with 2021 rows and 2 columns:
#' \describe{
#'   \item{c14_age}{The 14C age (yr BP)}
#'   \item{c14_sig}{The uncertainty in the 14C age reported by the radiocarbon laboratory}
#' }
#' @source <http://doi.org/10.1073/pnas.1408028111>
"armit"


#' Example output from Walker calibration
#'
#' This output has been provided for trying out the plotting and summarisation
#' functions without having to first run a calibration. It shows the output from
#' running [carbondate::WalkerBivarDirichlet] with the [carbondate::kerr]
#' radiocarbon input data and the [carbondate::intcal20] calibration curve for
#' 100,000 iterations. To reduce file size a large \eqn{n_{\textrm{thin}} = 100}
#' was used and the first half of the observations where discarded.
#'
#' @format ## `walker_example_output`
#' A list with 11 items, as described in [carbondate::WalkerBivarDirichlet]. In
#' this case \eqn{n_{\textrm{out}}} is 501 and \eqn{n_{\textrm{obs}}} is 255.
"walker_example_output"


#' Example output from Polya Urn calibration
#'
#' This output has been provided for trying out the plotting and summarisation
#' functions without having to first run a calibration. It shows the output from
#' running [carbondate::PolyaUrnBivarDirichlet] with the [carbondate::kerr]
#' radiocarbon input data and the [carbondate::intcal20] calibration curve for
#' 100,000 iterations. To reduce file size a large \eqn{n_{\textrm{thin}} = 100}
#' was used and the first half of the observations where discarded.
#'
#' @format ## `polya_urn_example_output`
#' A list with 10 items, as described in [carbondate::PolyaUrnBivarDirichlet]. In
#' this case \eqn{n_{\textrm{out}}} is 501 and \eqn{n_{\textrm{obs}}} is 255.
"polya_urn_example_output"

