.CheckCalibrationCurve <- function(arg_check, calibration_curve, F14C_inputs){
  checkmate::assertDataFrame(
    calibration_curve,
    types = "numeric",
    any.missing = FALSE,
    min.cols = 3,
    col.names = "named",
    add = arg_check)
  if (is.na(F14C_inputs)) {
    required_column_names = c("calendar_age_BP")
    # TODO: need to check they have one type of the other
  } else if (F14C_inputs == TRUE) {
    required_column_names = c("calendar_age_BP", "f14c", "f14c_sig")
  } else if (F14C_inputs == FALSE) {
    required_column_names = c("calendar_age_BP", "c14_age", "c14_sig")
  }
  checkmate::assertSubset(
    required_column_names,
    names(calibration_curve),
    .var.name ="calibration_curve required column names",
    add = arg_check)
}


.CheckInputData <- function(
    arg_check, rc_determinations, rc_sigmas, calibration_curve, F14C_inputs){

  checkmate::assertNumeric(
    rc_determinations,
    any.missing = FALSE,
    null.ok = FALSE,
    typed.missing = FALSE,
    add = arg_check)
  checkmate::assertNumeric(
    rc_sigmas,
    any.missing = FALSE,
    lower = 0,
    len = length(rc_determinations),
    null.ok = FALSE,
    typed.missing = FALSE,
    add = arg_check)
  .CheckCalibrationCurve(arg_check, calibration_curve, F14C_inputs)

}


.CheckDpmmParameters <- function(
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


.CheckSliceParameters <- function(
    arg_check, slice_width, slice_multiplier){
  checkmate::assertNumber(slice_width, lower = 1, add = arg_check)
  checkmate::assertNumber(slice_multiplier, lower = 1, add = arg_check)
}


.CheckIterationParameters <- function(arg_check, n_iter, n_thin){
  checkmate::assertInt(n_iter, lower = 10, add = arg_check)
  checkmate::assertInt(n_thin, lower = 1, upper = n_iter/10, add = arg_check)
}


.CheckOutputData <- function(arg_check, output_data) {
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
    .CheckCalibrationCurve(arg_check, calibration_curve)
  }

}


.CheckMultipleOutputDataConsistent <- function(output_data_list) {
  if (length(output_data_list) == 1) {
    return()
  }
  for (i in 2:length(output_data_list)) {
    first = output_data_list[[1]]$input_data
    other = output_data_list[[i]]$input_data
    if (
      !identical(first$c14_determinations, other$c14_determinations)
      || !identical(first$sigma, other$sigma)
      || first$F14C_inputs != other$F14C_inputs
      || first$calibration_curve_name != other$calibration_curve_name) {
      cli::cli_abort(
        c(
          "Output data is not consistent.",
          "Ensure all output data given in the list comes from the same
          calibration curve and radiocarbon values, in the same scale."
        )
      )
    }
  }
}


.CheckIntervalWidth = function(arg_check, interval_width, bespoke_probability) {
  checkmate::assertChoice(
    interval_width, c("1sigma", "2sigma", "bespoke"), add = arg_check)
  if (interval_width == "bespoke") {
    checkmate::assertNumber(
      bespoke_probability, lower = 0, upper = 1, add = arg_check)
  } else {
    if (!is.na(bespoke_probability)) {
      cli::cli_warn(
        c(paste("You have chosed an interval width of", interval_width)),
        c(
          paste("The value you have chosen for `bespoke_probability` will
                therefore be ignored")))
    }
  }
}


.CheckCalendarAgeSequence = function(arg_check, calendar_age_sequence) {
  checkmate::assertNumeric(
    calendar_age_sequence,
    unique = TRUE,
    sorted = TRUE,
    any.missing = FALSE,
    add = arg_check)
}
