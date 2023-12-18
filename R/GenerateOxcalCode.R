#' Outputs code suitable for running in OxCal from a series of radiocarbon determinations given as 14C age.
#'
#' @inheritParams FindSummedProbabilityDistribution
#' @param model_name The name given to the model in the OxCal code.
#' @param rc_names Optional. The name of each data point - if given it must be the same length of rc_determinations.
#' @param outfile_path Optional. If given the OxCal code will be output to the file at the path given, otherwise it
#' will be output to the terminal.
#' @export
#'
#' @examples
#' # Generate names automatically and outputs to the screen for 14C ages
#' GenerateOxcalCode("My_data", c(1123, 1128, 1135), c(32, 24, 25))
#'
#' # Provide name automatically and outputs to the screen for F14C concentrations
#' GenerateOxcalCode(
#'   "My_data",
#'   c(0.832, 0.850, 0.846),
#'   c(0.004, 0.003, 0.009),
#'   c("P-1", "P-2", "P-3"),
#'   F14C_inputs=TRUE)
#'
GenerateOxcalCode <- function(
  model_name, rc_determinations, rc_sigmas, rc_names = NULL, F14C_inputs = FALSE, outfile_path = NULL) {

  arg_check <- .InitializeErrorList()
  .CheckInputData(arg_check, rc_determinations, rc_sigmas, F14C_inputs)
  if (is.null(rc_names)) {
    rc_names <- seq_along(rc_determinations)
  } else {
    .CheckVector(arg_check, rc_names, len = length(rc_determinations))
  }
  .ReportErrors(arg_check)

  if (!is.null(outfile_path)) sink(outfile_path)
  if (F14C_inputs) {
    oxcal_command <- "R_F14C"
  } else {
    oxcal_command <- "R_Date"
  }

  cat(" Plot()\n")
  cat(" {\n")
  cat(paste0('  NP_Model(\"', model_name, '")\n'))
  cat("  {\n")
  for (i in seq_along(rc_determinations)) {
    cat(paste0('  ', oxcal_command , '(\"', rc_names[i], '",', rc_determinations[i], ',', rc_sigmas[i], ');\n'))
  }
  cat("  };\n")
  cat(" };\n")

  if (!is.null(outfile_path)) sink()
}
