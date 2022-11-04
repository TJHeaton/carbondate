.check_calibration_curve <- function(arg_check, calibration_curve){
  checkmate::assertDataFrame(
    calibration_curve,
    types = "numeric",
    any.missing = FALSE,
    min.cols = 3,
    col.names = "named",
    add = arg_check)
  checkmate::assertSubset(
    c("calendar_age", "c14_age", "c14_sig"),
    names(calibration_curve),
    .var.name ="calibration_curve required column names",
    add = arg_check)
}


.check_input_data <- function(
    arg_check, c14_determinations, c14_sigmas, calibration_curve){

  checkmate::assertNumeric(
    c14_determinations,
    any.missing = FALSE,
    lower = 0,
    null.ok = FALSE,
    typed.missing = FALSE,
    add = arg_check)
  checkmate::assertNumeric(
    c14_sigmas,
    any.missing = FALSE,
    lower = 0,
    len = length(c14_determinations),
    null.ok = FALSE,
    typed.missing = FALSE,
    add = arg_check)
  .check_calibration_curve(arg_check, calibration_curve)

}


.check_dpmm_parameters <- function(
    arg_check,
    sensible_initialisation,
    num_observations,
    lambda,
    nu1,
    nu2,
    A,
    B,
    alpha_shape,
    alpha_rate,
    mu_phi,
    calendar_ages,
    n_clust){

  if (!sensible_initialisation) {
    checkmate::assertNumber(lambda, add = arg_check)
    checkmate::assertNumber(nu1, add = arg_check)
    checkmate::assertNumber(nu2, add = arg_check)
    checkmate::assertNumber(A, add = arg_check)
    checkmate::assertNumber(B, add = arg_check)
    checkmate::assertNumber(alpha_shape, add = arg_check)
    checkmate::assertNumber(alpha_rate, add = arg_check)
    checkmate::assertNumber(mu_phi, add = arg_check)
    checkmate::assertNumeric(
      calendar_ages,
      len = num_observations,
      any.missing = FALSE,
      add = arg_check)
  } else {
    reason = "sensible_initialisation is TRUE"
    .WarnIfValueOverwritten(lambda, reason)
    .WarnIfValueOverwritten(nu1, reason)
    .WarnIfValueOverwritten(nu2, reason)
    .WarnIfValueOverwritten(A, reason)
    .WarnIfValueOverwritten(B, reason)
    .WarnIfValueOverwritten(alpha_shape, reason)
    .WarnIfValueOverwritten(alpha_rate, reason)
    .WarnIfValueOverwritten(mu_phi, reason)
    .WarnIfValueOverwritten(calendar_ages, reason)
  }
  checkmate::assertInt(
    n_clust, upper = num_observations, add = arg_check)
}


.WarnIfValueOverwritten <- function(var, reason = NULL) {
  varname = deparse(substitute(var))
  if (!is.na(var)) {
    warning_msg = paste("Provided value of", varname, "was overwritten")
    if (!is.null(reason)) {
      warning_msg = paste(warning_msg, "since", reason)
    }
    warning(warning_msg)
  }
}


.check_slice_parameters <- function(
    arg_check, slice_width, slice_multiplier){
  checkmate::assertNumber(slice_width, lower = 1, add = arg_check)
  checkmate::assertNumber(slice_multiplier, lower = 1, add = arg_check)
}


.check_iteration_parameters <- function(arg_check, n_iter, n_thin){
  checkmate::assertInt(n_iter, lower = 1, add = arg_check)
  checkmate::assertInt(n_thin, lower = 1, add = arg_check)
}


.check_output_data <- function(arg_check, output_data) {
  checkmate::assertList(
    output_data, names = "named", min.len = 10, add=arg_check)
  checkmate::assertSubset(
    c(
      "cluster_identifiers",
      "alpha",
      "n_clust",
      "phi",
      "tau",
      "calendar_ages",
      "mu_phi",
      "update_type",
      "input_data",
      "input_parameters"),
    names(output_data),
    .var.name ="output_data required list item names")
  checkmate::assertChoice(output_data$update_type, c("Polya Urn", "Walker"))
  if (output_data$update_type == "walker") {
    checkmate::assertChoice(
      "weight",
      names(output_data),
      .var.name ="output_data required list item names")
  }
}


.CheckCalibrationCurveFromOutput <- function(
    arg_check, output_data, calibration_curve) {
  calibration_curve_name = output_data$input_data$calibration_curve_name
  if (!exists(calibration_curve_name) && is.null(calibration_curve)){
    checkmate::reportAssertions(arg_check)
    cli::cli_abort(
      c(
        paste("Calibration curve", calibration_curve_name, "does not exist\n"),
        "Either ensure a variable with this name exists, or pass the variable
        in the arguments"
      )
    )
  }
  if (!is.null(calibration_curve)) {
    .check_calibration_curve(arg_check, calibration_curve)
  }

}


.CheckMultipleOutputDataConsistent <- function(output_data_list) {
  if (length(output_data_list) == 1) {
    return()
  }
  for (i in 2:length(output_data_list)) {
    if (output_data_list[[1]]$input_data != output_data_list[[i]]$input_data) {
      cli::cli_abort(
        c(
          "Output data is not consistent.",
          "Ensure all output data given in the list comes from the same
          calibration curve and c14 values."
        )
      )
    }
  }
}
