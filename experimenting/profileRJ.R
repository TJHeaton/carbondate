library(devtools)
library(proftools)
load_all()

set.seed(14)

file_pref <- "humans"
cutoffages <- c(6000, 25000)
data <- read.csv("experimenting/AnimalExtinctions/Humanc14Dates.csv", header = TRUE, na.strings = c("", "greater than"))
keep <- which(data$X14C < cutoffages[2] & data$X14C > cutoffages[1])
data <- data[keep,]
removena <- which(is.na(data$X14C) | is.na(data$X1..Sigma))
if(length(removena) != 0) { # Remove any values where 14C or sigma is NA
  data <- data[-removena,]
}

rc_determinations <- data$X14C
rc_sigmas <- rep(500, length(rc_determinations))# data$X1..Sigma

rm(keep, data, removena)

Rprof(paste0("experimenting/", file_pref, ".out"))

Test_Output <- PPcalibrate(
  rc_determinations, rc_sigmas, intcal20, n_iter = 1e4, n_thin = 10, calendar_grid_resolution = 10, use_fast = TRUE)

Rprof(NULL)

pd <- readProfileData(paste0("experimenting/", file_pref, ".out"))
flat_profile <- flatProfile(pd)
print(flat_profile)
save(flat_profile, file= paste0("experimenting/", file_pref, ".rda"))

profileCallGraph2Dot(pd, score = "total", filename = paste0("experimenting/", file_pref, ".dot"))


