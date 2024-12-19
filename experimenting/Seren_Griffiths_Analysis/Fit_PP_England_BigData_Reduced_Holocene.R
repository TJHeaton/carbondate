set.seed(27)

radiocarbon_age_cutoff <- 10000

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

# Restrict to under radiocarbon age cutoff
recent_England <- which(England_c14_ages < radiocarbon_age_cutoff)
Holocene_England_c14_ages <- England_c14_ages[recent_England]
Holocene_England_c14_sigs <- England_c14_sigs[recent_England]


Holocene_England_Output <- PPcalibrateLargeSets(
  rc_determinations = Holocene_England_c14_ages,
  rc_sigmas = Holocene_England_c14_sigs,
  calibration_curve = intcal20,
  n_iter = 100000,
  prior_n_internal_changepoints_lambda = 10,
  k_max_internal_changepoints = 40,
  calendar_age_range = c(0, 11650),
  show_progress = TRUE)

Holocene_England_PostMeanRate <- PlotPosteriorMeanRate(
  Holocene_England_Output,
  show_individual_means = FALSE,
  plot_cal_age_scale = "AD",
  denscale = 2)

save.image("../SerenGriffithsData/RWorkspaces/PriorMean10Changes/Holocene_England_Analysis_seed27_kmax_40.RData")


# Test against SPD
TempSPDEngland <- FindSummedProbabilityDistribution(
  c(0, 11650),
  Holocene_England_c14_ages,
  Holocene_England_c14_sigs,
  intcal20, plot_output = TRUE)
