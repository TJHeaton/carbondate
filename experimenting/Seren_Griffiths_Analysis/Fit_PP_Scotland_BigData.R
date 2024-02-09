set.seed(31)

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


Temp_Scotland_Output <- PPcalibrateLargeSets(
  rc_determinations = Scotland_c14_ages,
  rc_sigmas = Scotland_c14_sigs,
  calibration_curve = intcal20,
  n_iter = 100000,
  show_progress = TRUE)

ScotlandPostMeanRate <- PlotPosteriorMeanRate(Temp_Scotland_Output,
                                              show_individual_means = FALSE,
                                              denscale = 2)

save.image("../SerenGriffithsData/RWorkspaces/ScotlandAnalysis.RData")


