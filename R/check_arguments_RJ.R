.CheckRJPPOutputData <- function(arg_check, output_data) {
  if (!is.list(output_data) || is.null(names(output_data))) {
    arg_check$push("output data must be a named list")
    return()
  }
  list_names <- names(output_data)
  if (!("update_type" %in% list_names) || output_data$update_type != "RJPP") {
    arg_check$push("output data must contain the item update_type with value RJPP")
  }
  required_list_items <- c(
    "rate_s",
    "rate_h",
    "calendar_ages",
    "n_internal_changes",
    "input_data",
    "input_parameters"
  )
  for (list_item in required_list_items) {
    if (!(list_item %in% list_names)) {
      arg_check$push(paste("output data must contain the item", list_item))
    }
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