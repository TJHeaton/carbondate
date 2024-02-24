# set.seed(57)
set.seed(12)

radiocarbon_age_cutoff <- 10000

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

recent_Wales <- which(Wales_c14_ages < radiocarbon_age_cutoff)
Holocene_Wales_c14_ages <- Wales_c14_ages[recent_Wales]
Holocene_Wales_c14_sigs <- Wales_c14_sigs[recent_Wales]


Holocene_Wales_Output <- PPcalibrateLargeSets(
  rc_determinations = Holocene_Wales_c14_ages,
  rc_sigmas = Holocene_Wales_c14_sigs,
  calibration_curve = intcal20,
  n_iter = 100000,
  prior_n_internal_changepoints_lambda = 10,
  calendar_age_range = c(0, 11650),
  show_progress = TRUE)

Holocene_Wales_PostMeanRate <- PlotPosteriorMeanRate(
  Holocene_Wales_Output,
  show_individual_means = FALSE,
  plot_cal_age_scale = "AD",
  denscale = 2)

save.image("../SerenGriffithsData/RWorkspaces/PriorMean10Changes/Holocene_Wales_Analysis_ChangeMean10_seed12.RData")
