set.seed(9)

Republic_Ireland_Data <- read.csv("../SerenGriffithsData/Ireland_14C_Big_Data_191223.csv",
                                  header = TRUE)


if(!is.numeric(Republic_Ireland_Data$Conventional.Radiocarbon.age..year.BP...radiocarbon_age.)) {
  Republic_Ireland_Data$Conventional.Radiocarbon.age..year.BP...radiocarbon_age. <- as.numeric(Republic_Ireland_Data$Conventional.Radiocarbon.age..year.BP...radiocarbon_age.)
}
if(!is.numeric(Republic_Ireland_Data$Radiocarbon.age.error..radiocarbon_age_error.)) {
  Republic_Ireland_Data$Radiocarbon.age.error..radiocarbon_age_error. <- as.numeric(Republic_Ireland_Data$Radiocarbon.age.error..radiocarbon_age_error.)
}


# Remove those rows lacking c14 information
remove_obs_id <- which(
  is.na(Republic_Ireland_Data$Conventional.Radiocarbon.age..year.BP...radiocarbon_age.) |
    is.na(Republic_Ireland_Data$Radiocarbon.age.error..radiocarbon_age_error.))
if(length(remove_obs_id) != 0) {
  Republic_Ireland_Data <- Republic_Ireland_Data[-remove_obs_id,]
}

# Remove row with spurious/false precision on c14sig
remove_spurious_id <- which(Republic_Ireland_Data$Radiocarbon.age.error..radiocarbon_age_error. == 0.486
)
Republic_Ireland_Data <- Republic_Ireland_Data[-remove_spurious_id,]


# Make c14ages positive
Republic_Ireland_c14_ages <- as.numeric(Republic_Ireland_Data$Conventional.Radiocarbon.age..year.BP...radiocarbon_age.)
Republic_Ireland_c14_sigs <- as.numeric(Republic_Ireland_Data$Radiocarbon.age.error..radiocarbon_age_error.)


Temp_RoI_Output <- PPcalibrateLargeSets(
  rc_determinations = Republic_Ireland_c14_ages,
  rc_sigmas = Republic_Ireland_c14_sigs,
  calibration_curve = intcal20,
  n_iter = 100000,
  show_progress = TRUE)

PlotPosteriorMeanRate(Temp_RoI_Output,
                      show_individual_means = FALSE,
                      denscale = 2)

save.image("../SerenGriffithsData/RWorkspaces/IrelandAnalysis.RData")

