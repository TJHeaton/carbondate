set.seed(16)

England_Data <- read.csv("../SerenGriffithsData/England_14C_Big_Data_191223.csv",
                         header = TRUE)

if(!is.numeric(England_Data$DateResult)) {
  England_Data$DateResult <- as.numeric(England_Data$DateResult)
}
if(!is.numeric(England_Data$DateError)) {
  England_Data$DateError <- as.numeric(England_Data$DateError)
}

# Remove those rows lacking c14 information
remove_obs_id <- which(
  is.na(England_Data$DateResult) |
    is.na(England_Data$DateError))
if(length(remove_obs_id) != 0) {
  England_Data <- England_Data[-remove_obs_id,]
}

# Remove row with spurious/false precision on c14sig (values less than 2)
remove_spurious_id <- which(England_Data$DateError < 2)
if(length(remove_spurious_id) != 0) {
  England_Data <- England_Data[-remove_spurious_id,]
}

# Make c14ages positive
England_c14_ages <- as.numeric(England_Data$DateResult)
England_c14_sigs <- as.numeric(England_Data$DateError)


Temp_England_Output <- PPcalibrateLargeSets(
  rc_determinations = England_c14_ages,
  rc_sigmas = England_c14_sigs,
  calibration_curve = intcal20,
  n_iter = 100000,
  show_progress = TRUE)

EnglandPostMeanRate <- PlotPosteriorMeanRate(Temp_England_Output,
                                             show_individual_means = FALSE,
                                             denscale = 2)

save.image("../SerenGriffithsData/RWorkspaces/EnglandAnalysis.RData")


