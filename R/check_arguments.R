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
    alpha_type,
    alpha_shape,
    alpha_rate,
    alpha_mu,
    alpha_sigma,
    calendar_ages,
    n_clust){

  if (!sensible_initialisation) {
    checkmate::assertNumber(A, add = arg_check)
    checkmate::assertNumber(B, add = arg_check)
    checkmate::assertDouble(
      calendar_ages,
      len = num_observations,
      any.missing = FALSE,
      add = arg_check)
  }
  checkmate::assertNumber(lambda, add = arg_check)
  checkmate::assertNumber(nu1, add = arg_check)
  checkmate::assertNumber(nu2, add = arg_check)
  checkmate::assertChoice(alpha_type, c("gamma", "lognorm"))
  if (alpha_type == "gamma") {
    checkmate::assertNumber(alpha_shape, add = arg_check)
    checkmate::assertNumber(alpha_rate, add = arg_check)
  } else {
    checkmate::assertNumber(alpha_sigma, add = arg_check)
    checkmate::assertNumber(alpha_mu, add = arg_check)
  }
  checkmate::assertDouble(
    n_clust,
    upper = num_observations,
    any.missing = FALSE,
    add = arg_check)
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
    output_data, names = "named", min.len = 11, add=arg_check)
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
      "lambda",
      "nu1",
      "nu2"),
    names(output_data),
    .var.name ="output_data required list item names")
  checkmate::assertChoice(output_data$update_type, c("neal", "walker"))
  if (output_data$update_type == "walker") {
    checkmate::assertChoice(
      "weight",
      names(output_data),
      .var.name ="output_data required list item names")
  }
}
