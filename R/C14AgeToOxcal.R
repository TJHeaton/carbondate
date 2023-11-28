#' Outputs code suitable for running in OxCal from a series of radiocarbon determinations given as 14C age.
#'
#' @param model_name The name given to the model in the OxCal code.
#' @param c14_age An array of 14C determinations in yr BP.
#' @param c14_sig An array of the errors for the 14C determinations, must be the same length as c14_age.
#' @param c14_name Optional. The name of each data point - if given it must be the same length of c14_age.
#' @param outfile_path Optional. If given the OxCal code will be output to the file at the path given, otherwise it
#' will be output to the terminal.
#' @export
#'
#' @examples
#' # Generate names automatically and outputs to the screen
#' C14AgeToOxcal("My_data", c(1123, 1128, 1135), c(32, 24, 25))
#'
C14AgeToOxcal <- function(model_name, c14_age, c14_sig, c14_name = NULL, outfile_path = NULL) {

  if (length(c14_age) != length(c14_sig)) cat("error")
  if (is.null(c14_name)) c14_name <- seq_len(length(c14_age))

  if (!is.null(outfile_path)) {
    sink(outfile_path)
  }

  cat(" Plot()\n")
  cat(" {\n")
  cat(paste0('  NP_Model(\"', model_name, '")\n'))
  cat("  {\n")
  for (i in seq_along(c14_age)) {
    cat(paste0('  R_Date(\"', c14_name[i], '",', c14_age[i], ',', c14_sig[i], ');\n'))
  }
  cat("  };\n")
  cat(" };\n")

  if (!is.null(outfile_path)) {
    sink()
  }
}
