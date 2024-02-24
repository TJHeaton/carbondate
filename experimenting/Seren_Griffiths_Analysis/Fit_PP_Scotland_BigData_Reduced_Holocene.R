set.seed(2)

radiocarbon_age_cutoff <- 10000

Scotland_Data <- read.csv("../SerenGriffithsData/Scotland_14C_Big_Data_191223.csv",
                          header = TRUE)

if(!is.numeric(Scotland_Data$lab_age)) {
  Scotland_Data$lab_age <- as.numeric(Scotland_Data$lab_age)
}
if(!is.numeric(Scotland_Data$lab_error)) {
  Scotland_Data$lab_error <- as.numeric(Scotland_Data$lab_error)
}

# Remove those rows lacking c14 information
remove_obs_id <- which(
  is.na(Scotland_Data$lab_age) |
    is.na(Scotland_Data$lab_error))
if(length(remove_obs_id) != 0) {
  Scotland_Data <- Scotland_Data[-remove_obs_id,]
}

# Remove row with spurious/false precision on c14sig (values less than 2)
remove_spurious_id <- which(Scotland_Data$lab_error < 2)
if(length(remove_spurious_id) != 0) {
  Scotland_Data <- Scotland_Data[-remove_spurious_id,]
}

# Make c14ages positive
Scotland_c14_ages <- as.numeric(Scotland_Data$lab_age)
Scotland_c14_sigs <- as.numeric(Scotland_Data$lab_error)

# Restrict to under radiocarbon age cutoff
recent_Scotland <- which(Scotland_c14_ages < radiocarbon_age_cutoff)
Holocene_Scotland_c14_ages <- Scotland_c14_ages[recent_Scotland]
Holocene_Scotland_c14_sigs <- Scotland_c14_sigs[recent_Scotland]


Holocene_Scotland_Output <- PPcalibrateLargeSets(
  rc_determinations = Holocene_Scotland_c14_ages,
  rc_sigmas = Holocene_Scotland_c14_sigs,
  calibration_curve = intcal20,
  n_iter = 100000,
  prior_n_internal_changepoints_lambda = 10,
  calendar_age_range = c(0, 11650),
  show_progress = TRUE)

Holocene_Scotland_PostMeanRate <- PlotPosteriorMeanRate(
  Holocene_Scotland_Output,
  show_individual_means = FALSE,
  plot_cal_age_scale = "AD",
  denscale = 2)


save.image("../SerenGriffithsData/RWorkspaces/PriorMean10Changes/Holocene_Scotland_Analysis_seed2.RData")


