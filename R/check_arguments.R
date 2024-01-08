.InitializeErrorList <- function() {
  msgs <- character(0L)
  x <- list(
    push = function(msg) {msgs <<- c(msgs, msg)},
    getMessages = function() msgs,
    isEmpty = function() length(msgs) == 0L)
  return(x)
}

.ReportErrors <- function(error_list) {
  if (!error_list$isEmpty()) {
    msgs <- error_list$getMessages()
    context <- "%i checks failed:"
    err <- c(sprintf(context, length(msgs)), strwrap(msgs, prefix = " * "))
    stop(simpleError(paste0(err, collapse = "\n"), call = sys.call(1L)))
  }
}

.CheckInteger <- function(arg_check, x, lower = NA, upper = NA) {
  varname <- deparse(substitute(x))
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

.CheckNumber <- function(arg_check, x, lower = NA, upper = NA) {
  varname <- deparse(substitute(x))
  if (!is.numeric(x) || length(x) > 1) {
    arg_check$push(paste(varname, "must be a number"))
    return()
  }
  if (!is.na(lower) && x < lower) {
    arg_check$push(paste(varname, "must be more than or equal to", lower))
  }
  if (!is.na(upper) && x > upper) {
    arg_check$push(paste(varname, "must be less than or equal to", upper))
  }
}

.CheckFlag <- function(arg_check, x) {
  varname <- deparse(substitute(x))
  if (!is.logical(x) || length(x) > 1) {
    arg_check$push(paste(varname, "must be a single logical value (TRUE, FALSE OR NA)"))
  }
}


.CheckVector <- function(arg_check, x, min_length = NA, len = NA) {
  varname <- deparse(substitute(x))
  if (!is.vector(x)) {
    arg_check$push(paste(varname, "must be a vector"))
  }
  if (!is.na(min_length) && length(x) < min_length) {
    arg_check$push(paste(varname, "must have at least", min_length, "elements"))
  }
  if (!is.na(len) && length(x) != len) {
    arg_check$push(paste(varname, "must have exactly", len, "elements"))
  }
}


.CheckNumberVector <- function(arg_check, x, min_length = NA, len = NA, lower = NA) {
  varname <- deparse(substitute(x))
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
  if (!is.na(lower) && any(x < lower)) {
    arg_check$push(paste("all entries of", varname, "must be more than or equal to", lower))
  }
}


.CheckIntegerVector <- function(arg_check, x, max_length = NA, len = NA, lower = NA, upper = NA) {
  varname <- deparse(substitute(x))
  if (!is.vector(x)) {
    arg_check$push(paste(varname, "must be a vector"))
  }
  if (!is.numeric(x) || any(is.na(x)) || any(as.integer(x) - x != 0)) {
    arg_check$push(paste(varname, "must contain only integer entries"))
    return()
  }
  if (!is.na(max_length) && length(x) > max_length) {
    arg_check$push(paste(varname, "must have at most", max_length, "elements"))
  }
  if (!is.na(len) && length(x) != len) {
    arg_check$push(paste(varname, "must have exactly", len, "elements"))
  }
  if (!is.na(lower) && any(x < lower)) {
    arg_check$push(paste("all entries of", varname, "must be more than or equal to", lower))
  }
  if (!is.na(upper) && any(x > upper)) {
    arg_check$push(paste("all entries of", varname, "must be less than or equal to", upper))
  }
}

.CheckChoice <- function(arg_check, x, allowed_choices) {
  varname <- deparse(substitute(x))
  if (!(x %in% allowed_choices)) {
    arg_check$push(paste(varname, "must be one of:", paste(allowed_choices, collapse=", ")))
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
  .CheckNumberVector(arg_check, rc_determinations)
  .CheckNumberVector(arg_check, rc_sigmas, lower = 0, len = length(rc_determinations))
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


.CheckDPMMParameters <- function(
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
    .CheckNumber(arg_check, lambda)
    .CheckNumber(arg_check, nu1)
    .CheckNumber(arg_check, nu2)
    .CheckNumber(arg_check, A)
    .CheckNumber(arg_check, B)
    .CheckNumber(arg_check, alpha_shape)
    .CheckNumber(arg_check, alpha_rate)
    .CheckNumber(arg_check, mu_phi)
    .CheckNumberVector(arg_check, calendar_ages, len = num_observations)
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
  .CheckInteger(arg_check, n_clust, upper = num_observations)
}


.WarnIfValueOverwritten <- function(var, reason = NULL) {
  varname <- deparse(substitute(var))
  if (!is.na(var)) {
    warning_msg <- paste("Provided value of", varname, "was overwritten")
    if (!is.null(reason)) {
      warning_msg <- paste(warning_msg, "since", reason)
    }
    warning(warning_msg, immediate. = TRUE)
  }
}


.CheckSliceParameters <- function(arg_check, slice_width, slice_multiplier, sensible_initialisation) {
  .CheckNumber(arg_check, slice_multiplier, lower = 1)
  if (!sensible_initialisation || !is.na(slice_width)) {
    .CheckNumber(arg_check, slice_width, lower = 1)
  }
}


.CheckIterationParameters <- function(arg_check, n_iter, n_thin){
  .CheckInteger(arg_check, n_iter, lower = 10)
  .CheckInteger(arg_check, n_thin, lower = 1, upper = n_iter/10)
}


.CheckOutputData <- function(arg_check, output_data, allowed_update_types) {
  if (!is.list(output_data) || is.null(names(output_data))) {
    arg_check$push("output data must be a named list")
    return()
  }

  list_names <- names(output_data)
  if (!("update_type" %in% list_names)) {
    arg_check$push("output data must contain the item update_type")
  } else {
    .CheckChoice(arg_check, output_data$update_type, allowed_update_types)
  }

  required_list_items <- c("input_data", "input_parameters", "calendar_ages")
  if (output_data$update_type == "Polya Urn") {
    required_list_items <- c(
      required_list_items, "cluster_identifiers", "alpha", "n_clust", "phi", "tau", "mu_phi")
  } else if (output_data$update_type == "Walker") {
    required_list_items <- c(
      required_list_items, "cluster_identifiers", "alpha", "n_clust", "phi", "tau", "mu_phi", "weight")
  } else if (output_data$update_type == "RJPP") {
    required_list_items <- c(required_list_items, "rate_s", "rate_h", "n_internal_changes")
  } else {
    stop(paste("Internal error: unknown update type:", output_data$update_type))
  }

  for (list_item in required_list_items) {
    if (!(list_item %in% list_names)) {
      arg_check$push(paste("output data must contain the item", list_item))
    }
  }
}


.CheckCalibrationCurveFromOutput <- function(arg_check, output_data, calibration_curve) {
  calibration_curve_name <- output_data$input_data$calibration_curve_name
  if (!exists(calibration_curve_name) && is.null(calibration_curve)){
    arg_check$push(
      paste(
        "Calibration curve", calibration_curve_name, "does not exist.",
        "Either ensure a variable with this name exists, or pass the variable in the arguments"
      )
    )
  }
  if (!is.null(calibration_curve)) {
    .CheckCalibrationCurve(arg_check, calibration_curve, NA)
  }
}


.CheckMultipleOutputDataConsistent <- function(arg_check, output_data_list) {
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
      arg_check$push(
        paste0(
          "Output data item [[1]] and item [[", i, "]] are not consistent. ",
          "Ensure all output data given in the list comes from the same ",
          "calibration curve and radiocarbon values, in the same scale."
        )
      )
    }
  }
}


.CheckIntervalWidth <- function(arg_check, interval_width, bespoke_probability) {
  .CheckChoice(arg_check, interval_width, c("1sigma", "2sigma", "bespoke"))
  if (interval_width == "bespoke") {
    .CheckNumber(arg_check, bespoke_probability, lower = 0, upper = 1)
  } else {
    if (!is.na(bespoke_probability)) {
      warning(
        paste0(
          "You have chosed an interval width of ", interval_width, ". The value you have chosen ",
          "for `bespoke_probability` will therefore be ignored"),
        immediate. = TRUE
      )
    }
  }
}


.CheckCalendarAgeSequence <- function(arg_check, calendar_age_sequence) {
  .CheckNumberVector(arg_check, calendar_age_sequence)
  if (length(unique(calendar_age_sequence)) != length(calendar_age_sequence)) {
    arg_check$push("All values in the calendar age sequence must be unique")
  }
  if (is.unsorted(calendar_age_sequence)) {
    arg_check$push("The calendar age sequence must be sorted")
  }
}


.CheckPriorHRateAndPriorHShape <- function(arg_check, prior_h_shape, prior_h_rate) {
  if (is.na(prior_h_shape) && is.na(prior_h_rate)) {
    return()
  }
  if (!is.na(prior_h_shape) && !is.na(prior_h_rate)) {
    .CheckNumber(arg_check, prior_h_rate, lower = 0)
    .CheckNumber(arg_check, prior_h_shape, lower = 0)
  } else {
    arg_check$push(
      "prior_h_shape and prior_h_rate must either both be positive numbers or must both be NA")
  }
}
