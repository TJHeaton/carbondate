set.seed(57)

Wales_Data <- read.csv("../SerenGriffithsData/Wales_14C_Big_Data_191223.csv",
                                  header = TRUE)

if(!is.numeric(Wales_Data$DateResult)) {
Wales_Data$DateResult <- as.numeric(Wales_Data$DateResult)
}
if(!is.numeric(Wales_Data$DateError)) {
  Wales_Data$DateError <- as.numeric(Wales_Data$DateError)
}

# Remove those rows lacking c14 information
remove_obs_id <- which(
  is.na(Wales_Data$DateResult) |
    is.na(Wales_Data$DateError))
if(length(remove_obs_id) != 0) {
  Wales_Data <- Wales_Data[-remove_obs_id,]
}

# Remove row with spurious/false precision on c14sig (values less than 2)
remove_spurious_id <- which(Wales_Data$DateError < 2)
if(length(remove_spurious_id) != 0) {
  Wales_Data <- Wales_Data[-remove_spurious_id,]
}

# Make c14ages positive
Wales_c14_ages <- as.numeric(Wales_Data$DateResult)
Wales_c14_sigs <- as.numeric(Wales_Data$DateError)


Temp_Wales_Output <- PPcalibrateLargeSets(
  rc_determinations = Wales_c14_ages,
  rc_sigmas = Wales_c14_sigs,
  calibration_curve = intcal20,
  n_iter = 100000,
  show_progress = TRUE)

WalesPostMeanRate <- PlotPosteriorMeanRate(Temp_Wales_Output,
                                           show_individual_means = FALSE,
                                           denscale = 2)

save.image("../SerenGriffithsData/RWorkspaces/WalesAnalysis.RData")
