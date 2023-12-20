library(devtools)
library(proftools)
load_all()

set.seed(14)

file_pref <- "humans"

Rprof(paste0("experimenting/", file_pref, ".out"))

Test_Output <- PPcalibrate(
  human$c14_age, human$c14_sig, intcal20, n_iter = 5e4, n_thin = 10, calendar_grid_resolution = 10, use_fast = TRUE)

Rprof(NULL)

pd <- readProfileData(paste0("experimenting/", file_pref, ".out"))
flat_profile <- flatProfile(pd)
print(flat_profile)
save(flat_profile, file= paste0("experimenting/", file_pref, ".rda"))

profileCallGraph2Dot(pd, score = "total", filename = paste0("experimenting/", file_pref, ".dot"))


