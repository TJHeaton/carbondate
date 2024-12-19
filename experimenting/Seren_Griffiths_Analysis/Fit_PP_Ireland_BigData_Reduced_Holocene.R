#set.seed(9)
set.seed(85)

radiocarbon_age_cutoff <- 10000

Ireland_Data <- read.csv("../SerenGriffithsData/Ireland_14C_Big_Data_191223.csv",
                                  header = TRUE)

if(!is.numeric(Ireland_Data$Conventional.Radiocarbon.age..year.BP...radiocarbon_age.)) {
  Ireland_Data$Conventional.Radiocarbon.age..year.BP...radiocarbon_age. <- as.numeric(Ireland_Data$Conventional.Radiocarbon.age..year.BP...radiocarbon_age.)
}
if(!is.numeric(Ireland_Data$Radiocarbon.age.error..radiocarbon_age_error.)) {
  Ireland_Data$Radiocarbon.age.error..radiocarbon_age_error. <- as.numeric(Ireland_Data$Radiocarbon.age.error..radiocarbon_age_error.)
}

# Remove those rows lacking c14 information
remove_obs_id <- which(
  is.na(Ireland_Data$Conventional.Radiocarbon.age..year.BP...radiocarbon_age.) |
    is.na(Ireland_Data$Radiocarbon.age.error..radiocarbon_age_error.))
if(length(remove_obs_id) != 0) {
  Ireland_Data <- Ireland_Data[-remove_obs_id,]
}

# Remove row with spurious/false precision on c14sig
remove_spurious_id <- which(Ireland_Data$Radiocarbon.age.error..radiocarbon_age_error. < 2)
Ireland_Data <- Ireland_Data[-remove_spurious_id,]


# Make c14ages positive
Ireland_c14_ages <- as.numeric(Ireland_Data$Conventional.Radiocarbon.age..year.BP...radiocarbon_age.)
Ireland_c14_sigs <- as.numeric(Ireland_Data$Radiocarbon.age.error..radiocarbon_age_error.)

# Restrict to under radiocarbon age cutoff
recent_Ireland <- which(Ireland_c14_ages < radiocarbon_age_cutoff)
Holocene_Ireland_c14_ages <- Ireland_c14_ages[recent_Ireland]
Holocene_Ireland_c14_sigs <- Ireland_c14_sigs[recent_Ireland]


Holocene_Ireland_Output <- PPcalibrateLargeSets(
  rc_determinations = Holocene_Ireland_c14_ages,
  rc_sigmas = Holocene_Ireland_c14_sigs,
  calibration_curve = intcal20,
  n_iter = 100000,
  prior_n_internal_changepoints_lambda = 10,
  calendar_age_range = c(0, 11650),
  show_progress = TRUE)

Holocene_Ireland_PostMeanRate <- PlotPosteriorMeanRate(
  Holocene_Ireland_Output,
  show_individual_means = FALSE,
  plot_cal_age_scale = "AD",
  denscale = 2)

save.image("../SerenGriffithsData/RWorkspaces/PriorMean10Changes/Holocene_Ireland_Analysis_seed85.RData")


