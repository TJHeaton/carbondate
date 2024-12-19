# Code to create output for Pir
library(carbondate)





##### Extra functions to find HPD interval for a set of determinations

# Create HPD interval for a single determination/sample (determined by ident)
FindIndSampleHPDIntervals <- function(ident,
                                      output_data,
                                      interval_width = "2sigma",
                                      bespoke_probability = NA) {

  if(interval_width == "bespoke" && is.na(bespoke_probability)) {
    stop("Chosen bespoke probability but not specified that probability")
  }

  n_iter <- output_data$input_parameters$n_iter
  n_thin <- output_data$input_parameters$n_thin
  n_burn <- floor(n_iter / (2 * n_thin))
  n_end <- floor(n_iter / n_thin) + 1

  rc_determinations <- output_data$input_data$rc_determinations
  rc_sigmas <- output_data$input_data$rc_sigmas
  F14C_inputs <-output_data$input_data$F14C_inputs

  calendar_age_BP <- output_data$calendar_ages[(n_burn+1):n_end, ident]
  rc_age <- rc_determinations[ident]
  rc_sig <- rc_sigmas[ident]

  if(output_data$update_type == "RJPP") {
    bandwidth_selector <- "nrd0" # As all calendar ages are integers
  } else {
    bandwidth_selector <- "SJ" # As continuous
  }

  smoothed_density <- stats::density(calendar_age_BP, bw = bandwidth_selector)
  xrange_BP <- range(calendar_age_BP)

  hpd_probability <- switch(
    interval_width,
    "1sigma" = 0.682,
    "2sigma" = 0.954,
    "bespoke" = bespoke_probability)
  hpd <- carbondate:::FindHPD(smoothed_density$x, smoothed_density$y, hpd_probability)

  merged_hpd_intervals <- c(rbind(hpd$start_ages, hpd$end_ages))
  return(round(merged_hpd_intervals))

}



# Find matrix of intervals for a set of 14C samples/determinations (can change to dataframe is desired)
FindSetSamplesHPDIntervals <- function(ident_set,
                                       output_data,
                                       interval_width = "2sigma",
                                       bespoke_probability = NA) {

  if(interval_width == "bespoke" && is.na(bespoke_probability)) {
    stop("Chosen bespoke probability but not specified that probability")
  }

  n_iter <- output_data$input_parameters$n_iter
  n_thin <- output_data$input_parameters$n_thin
  n_burn <- floor(n_iter / (2 * n_thin))
  n_end <- floor(n_iter / n_thin) + 1

  list_intervals <- lapply(ident_set,
                           FindIndSampleHPDIntervals,
                           output_data = output_data,
                           interval_width = interval_width,
                           bespoke_probability = bespoke_probability)

  # Convert ragged list to a dataframe/matrix with relevant name headings
  max_n_intervals <- max(lengths(list_intervals))
  interval_matrix <- t(sapply(list_intervals, "[", seq(max_n_intervals)))
  colnames(interval_matrix) <- paste(interval_width,
                                     rep(1:(max_n_intervals/2), each = 2),
                                     c("start", "end"),
                                     sep = "_")

  # Create matrix/dataframe with 14C determination information
  sample_info_matrix <- data.frame(sample_ident = ident_set,
                                   radiocarbon_age = output_data$input_data$rc_determinations[ident_set],
                                   radiocarbon_sigma = output_data$input_data$rc_sigmas[ident_set])

  posterior_calendar_ages_ident_set <- output_data$calendar_ages[(n_burn+1):n_end, ident_set]

  # Add columns with posterior mean and median calendar ages
  if(length(ident_set) == 1) { # Deal with case that posterior_calendar_ages_ident_set is a vector
    sample_info_matrix$posterior_calendar_age_mean <- round(mean(posterior_calendar_ages_ident_set),1)
    sample_info_matrix$posterior_calendar_age_median <- round(median(posterior_calendar_ages_ident_set))
  } else {
    sample_info_matrix$posterior_calendar_age_mean <- round(apply(posterior_calendar_ages_ident_set, 2, mean),1)
    sample_info_matrix$posterior_calendar_age_median <- round(apply(posterior_calendar_ages_ident_set, 2, median))
  }

  # Add column specifying determination information
  interval_matrix <- cbind(sample_info_matrix, interval_matrix)

  return(interval_matrix)
}


###########################################################################
####### Run it on the worked example
###########################################################################

# Create example
set.seed(15)

# Set initial values
n_observed <- 40
rc_sigmas <- rep(15, n_observed)

# Create artificial rc_determinations
calendar_age_range <- c(300, 700)
observed_age_range <- c(500, 550)
true_theta <- seq(
  from = observed_age_range[1],
  to = observed_age_range[2],
  length = n_observed)
intcal_mean <- approx(
  x = intcal20$calendar_age_BP,
  y = intcal20$c14_age,
  xout = true_theta)$y
intcal_sd <- approx(
  x = intcal20$calendar_age_BP,
  y = intcal20$c14_sig,
  xout = true_theta)$y
rc_determinations <- rnorm(
  n = n_observed,
  mean = intcal_mean,
  sd = sqrt(rc_sigmas^2 + intcal_sd^2))

# Fit the model
PP_fit_output <- PPcalibrate(
  rc_determinations = rc_determinations,
  rc_sigmas = rc_sigmas,
  calibration_curve = intcal20,
  calendar_age_range = calendar_age_range,
  calendar_grid_resolution = 1,
  n_iter = 1e5,
  n_thin = 10,
  show_progress = TRUE)


# For the 9th determination we get a HPD of 521-501 which we can check later
PlotCalendarAgeDensityIndividualSample(
  9, PP_fit_output, show_hpd_ranges = TRUE, show_unmodelled_density = TRUE)
PlotCalendarAgeDensityIndividualSample(
  38, PP_fit_output, show_hpd_ranges = TRUE, show_unmodelled_density = TRUE)


#########################################################################
####### Testing the HPD functions
#########################################################################

# Now to actually use the function you can create a call such as
FindSetSamplesHPDIntervals(c(9, 2, 19, 38), PP_fit_output)
FindSetSamplesHPDIntervals(9, PP_fit_output)
FindSetSamplesHPDIntervals(38, PP_fit_output)
FindSetSamplesHPDIntervals(1:40, PP_fit_output) # All of first 40 determinations

# Note that there is a bit of an ugly aspect with (e.g. 38) where
# the intervals end and then start the very next year
# this is because we smooth the density first to find the HPD which works
# on a continuous scale in the code

# This also only provides HPD intervals in cal yr BP currently (not cal BC/AD)
# Finally - the intervals as provided currently go from youngest to oldest. Is this ok?

PP_fit_post_mean <- PlotPosteriorMeanRate(PP_fit_output)
PP_fit_post_mean


PlotNumberOfInternalChanges()

