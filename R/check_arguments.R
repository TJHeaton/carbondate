.makeAssertCollection <- function() {
  msgs <- character(0L)
  x <- list(
    push = function(msg) {msgs <<- c(msgs, msg)},
    getMessages = function() msgs,
    isEmpty = function() length(msgs) == 0L)
  return(x)
}

.reportAssertions <- function(collection) {
  if (!collection$isEmpty()) {
    msgs <- collection$getMessages()
    context <- "%i assertions failed:"
    err <- c(sprintf(context, length(msgs)), strwrap(msgs, prefix = " * "))
    stop(simpleError(paste0(err, collapse = "\n"), call = sys.call(1L)))
  }
}

.CheckInteger <- function(arg_check, x, lower = NA, upper = NA) {
  varname <- substitute(x)
  if (!is.numeric(x) || length(x) > 1 || (as.integer(x) - x != 0)) {
    arg_check$push(paste(varname, "must be an integer"))
    return()
  }
  if (!is.na(lower) && x < lower) {
    arg_check$push(paste(varname, "must be more than or equal to", lower))
  }
  if (!is.na(upper) && x > upper) {
    arg_check$push(paste(varname, "must be less than or equal to", upper))
  }
}

.CheckNumber <- function(arg_check, x) {
  if (!is.numeric(x) || length(x) > 1) {
    arg_check$push(paste(substitute(x), "must be a number"))
  }
}

.CheckFlag <- function(arg_check, x) {
  if (!is.logical(x) || length(x) > 1) {
    arg_check$push(paste(substitute(x), "must be a single logical value (TRUE, FALSE OR NA)"))
  }
}

.CheckNumberVector <- function(arg_check, x, min_length = NA, len = NA) {
  varname <- substitute(x)
  if (!is.vector(x)) {
    arg_check$push(paste(varname, "must be a vector"))
  }
  if (!is.numeric(x) || any(is.na(x))) {
    arg_check$push(paste(varname, "must have numeric entries (and not be NA)"))
  }
  if (!is.na(min_length) && length(x) < min_length) {
    arg_check$push(paste(varname, "must have at least", min_length, "elements"))
  }
  if (!is.na(len) && length(x) != len) {
    arg_check$push(paste(varname, "must have exactly", len, "elements"))
  }
}

.CheckCalibrationCurve <- function(arg_check, calibration_curve, F14C_inputs){
  if (!is.data.frame(calibration_curve)) {
    arg_check$push("The calibration curve must be a data frame")
  }
  header_names <- names(calibration_curve)
  if (!("calendar_age_BP" %in% header_names)) {
    arg_check$push("The calibration curve must have the column 'calendar_age_BP'")
  }

  has_f14c <- ("f14c" %in% header_names)
  has_f14c_sig <- ("f14c_sig" %in% header_names)
  has_c14_age <- ("c14_age" %in% header_names)
  has_c14_sig <- ("c14_sig" %in% header_names)

  if (is.na(F14C_inputs)) {
    if (!((has_f14c && has_f14c_sig) || (has_c14_age && has_c14_sig))) {
      arg_check$push(
        "The calibration curve must have the columns ('f14c', 'f14c_sig') and/or the columns ('c14_age', 'c14_sig')")
    }
  } else if (F14C_inputs == TRUE) {
    if (!(has_f14c && has_f14c_sig)) {
      arg_check$push("The calibration curve must have the columns ('f14c', 'f14c_sig')")
    }
  } else if (F14C_inputs == FALSE) {
    if (!(has_c14_age && has_c14_sig)) {
      arg_check$push("The calibration curve must have the columns ('c14_age', 'c14_sig')")
    }
  }

  for (i in ncol(calibration_curve)) {
    if (!is.numeric(calibration_curve[[i]])) {
      arg_check$push("The calibration curve entries must all be numbers")
      break
    }
  }
}


.CheckInputData <- function(arg_check, rc_determinations, rc_sigmas, F14C_inputs){

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
  if (F14C_inputs == TRUE) {
    if (any(rc_determinations > 2) || any(rc_determinations < 0)) {
      warning(
        "You have specified F14C_inputs = TRUE but it looks the rc_determinations may be 14C ages.",
        immediate. = TRUE, call. = FALSE)
    }
  } else {
    if (all(rc_determinations < 2) && all(rc_determinations > 0)) {
      warning(
        "You have specified F14C_inputs = FALSE but it looks the rc_determinations may be F14C concentrations",
        immediate. = TRUE, call. = FALSE)
    }
  }
}


.CheckNBurnAndNEnd <- function(arg_check, n_burn, n_end, n_iter, n_thin) {
  if (!is.na(n_burn)) {
    .CheckInteger(arg_check, n_burn, lower = 0, upper = n_iter - 100 * n_thin)
  } else {
    n_burn <- .SetNBurn(n_burn, n_iter, n_thin) * n_thin
  }
  if (!is.na(n_end)) {
    .CheckInteger(arg_check, n_end, lower = n_burn + 10 * n_thin, upper = n_iter)
  }
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
    reason <- "sensible_initialisation is TRUE"
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
  varname <- deparse(substitute(var))
  if (!is.na(var)) {
    warning_msg <- paste("Provided value of", varname, "was overwritten")
    if (!is.null(reason)) {
      warning_msg <- paste(warning_msg, "since", reason)
    }
    warning(warning_msg)
  }
}


.CheckSliceParameters <- function(
    arg_check, slice_width, slice_multiplier, sensible_initialisation) {
  checkmate::assertNumber(slice_multiplier, lower = 1, add = arg_check)
  if (sensible_initialisation) {
    checkmate::assertNumber(slice_width, lower = 1, na.ok = TRUE, add = arg_check)
  } else {
    checkmate::assertNumber(slice_width, lower = 1, na.ok = FALSE, add = arg_check)
  }
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
  calibration_curve_name <- output_data$input_data$calibration_curve_name
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
    .CheckCalibrationCurve(arg_check, calibration_curve, NA)
  }

}


.CheckMultipleOutputDataConsistent <- function(output_data_list) {
  if (length(output_data_list) == 1) {
    return()
  }
  for (i in 2:length(output_data_list)) {
    first <- output_data_list[[1]]$input_data
    other <- output_data_list[[i]]$input_data
    if (
      !identical(first$rc_determinations, other$rc_determinations)
      || !identical(first$rc_sigmas, other$rc_sigmas)
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


.CheckIntervalWidth <- function(arg_check, interval_width, bespoke_probability) {
  checkmate::assertChoice(
    interval_width, c("1sigma", "2sigma", "bespoke"), add = arg_check)
  if (interval_width == "bespoke") {
    checkmate::assertNumber(
      bespoke_probability, lower = 0, upper = 1, add = arg_check)
  } else {
    if (!is.na(bespoke_probability)) {
      cli::cli_warn(
        paste("You have chosed an interval width of", interval_width),
        paste("The value you have chosen for `bespoke_probability` will
              therefore be ignored"))
    }
  }
}


.CheckCalendarAgeSequence <- function(arg_check, calendar_age_sequence) {
  checkmate::assertNumeric(
    calendar_age_sequence,
    unique = TRUE,
    sorted = TRUE,
    any.missing = FALSE,
    add = arg_check)
}
