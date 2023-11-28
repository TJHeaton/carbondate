#' Outputs text suitable for running in OxCal
#'
#' @param ident the determination you want to show the individual posterior
#' calendar age for.
#' @param resolution The distance between histogram breaks for the calendar age density.
#' Must be an integer greater than one.
#' @param interval_width The confidence intervals to show for the
#' calibration curve and for the highest posterior density ranges.
#' Choose from one of "1sigma" (68.3%), "2sigma" (95.4%) and "bespoke". Default is "2sigma".
#' @param bespoke_probability The probability to use for the confidence interval
#' if "bespoke" is chosen above. E.g. if 0.95 is chosen, then the 95% confidence
#' interval is calculated. Ignored if "bespoke" is not chosen.
#' @param show_hpd_ranges Set to `TRUE` to also show the highest posterior range on the plot.
#' These are calculated using [HDInterval::hdi]. Default is `FALSE`.
#' @param show_unmodelled_density Set to `TRUE` to also show the unmodelled density (i.e. the
#' result of [carbondate::CalibrateSingleDetermination]) on the plot. Default is `FALSE`.
#'
#' @export
#'
#' @examples
#' # Generate names automatically
#' C14AgeToOxcal("My_data", c(1123, 1128, 1135), c(32, 24, 25))
#'
C14AgeToOxcal <- function(model_name, c14_age, c14_sig, c14_name = NULL) {

  if (length(c14_age) != length(c14_sig)) cat("error")
  if (is.null(c14_name)) c14_name = seq_len(length(c14_age))

  cat(" Plot()\n")
  cat(" {\n")
  cat(paste('  NP_Model(\"', model_name ,'")\n', sep=""))
  cat("  {\n")
  for (i in 1:length(c14_age)) {
    cat(paste('  R_Date(\"', c14_name[i] ,'",', c14_age[i],',', c14_sig[i],');\n', sep=""))
  }
  cat("  };\n")
  cat(" };\n")
}
