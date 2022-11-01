.check_input_data <- function(
    argument_checker, c14_determinations, c14_uncertainties, calibration_curve){

  checkmate::assertDouble(
    c14_determinations,
    any.missing = FALSE,
    lower = 0,
    null.ok = FALSE,
    typed.missing = FALSE,
    add = argument_checker)
  checkmate::assertDouble(
    c14_uncertainties,
    any.missing = FALSE,
    lower = 0,
    len = length(c14_determinations),
    null.ok = FALSE,
    typed.missing = FALSE,
    add = argument_checker)
  checkmate::assertDataFrame(
    calibration_curve,
    types = "numeric",
    any.missing = FALSE,
    min.cols = 3,
    add = argument_checker)
  checkmate::assertSubset(
    c("calendar_age", "c14_age", "c14_sig"),
    names(calibration_curve),
    .var.name ="calibration_curve required column names",
    add = argument_checker)

}


.check_dpmm_parameters <- function(
    argument_checker,
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
    checkmate::assertDouble(
      A, len = 1, any.missing = FALSE, add = argument_checker)
    checkmate::assertDouble(
      B, len = 1, any.missing = FALSE, add = argument_checker)
    checkmate::assertDouble(
      calendar_ages,
      len = num_observations,
      any.missing = FALSE,
      add = argument_checker)
  }
  checkmate::assertDouble(
    lambda, len = 1, any.missing = FALSE, add = argument_checker)
  checkmate::assertDouble(
    nu1, len = 1, any.missing = FALSE, add = argument_checker)
  checkmate::assertDouble(
    nu2, len = 1, any.missing = FALSE, add = argument_checker)
  assertChoice(alpha_type, c("gamma", "lognorm"))
  if (alpha_type == "gamma") {
    checkmate::assertDouble(
      alpha_shape, len = 1, any.missing = FALSE, add = argument_checker)
    checkmate::assertDouble(
      alpha_rate, len = 1, any.missing = FALSE, add = argument_checker)
  } else {
    checkmate::assertDouble(
      alpha_mu, len = 1, any.missing = FALSE, add = argument_checker)
    checkmate::assertDouble(
      alpha_sigma, len = 1, any.missing = FALSE, add = argument_checker)
  }
  checkmate::assertDouble(
    n_clust,
    upper = num_observations,
    any.missing = FALSE,
    add = argument_checker)
}


.check_slice_parameters <- function(
    argument_checker, slice_width, slice_multiplier){
  checkmate::assertDouble(
    slice_width,
    len = 1,
    lower = 1,
    any.missing = FALSE,
    add = argument_checker)
  checkmate::assertDouble(
    slice_multiplier,
    len = 1,
    lower = 1,
    any.missing = FALSE,
    add = argument_checker)
}


.check_iteration_parameters <- function(argument_checker, n_iter, n_thin){
  checkmate::assertDouble(
    n_iter, len = 1, lower = 1, any.missing = FALSE, add = argument_checker)
  checkmate::assertDouble(
    n_thin, len = 1, lower = 1, any.missing = FALSE, add = argument_checker)
}
