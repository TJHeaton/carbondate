# An analysis of rates for various Pleistocene Animals from Guthrie
# https://www.nature.com/articles/nature04604#Sec3

# We assess RJ MCMC on rates for 14C dates for a range of animals
# We select 14C dates that lie between [6, 25] 14Cyrs BP

##############################################
Species <- "Mammoth" # "Equus" or "Human" or "Mammoth" or "Bison" or "Alces" or "Cervus"

# Main function - you just enter the species and the calibration curve you want (interpolated onto a 5 yearly grid)
cutoffages <- c(6000, 25000)

# Read in the relevant datasets dependent upon species
if(Species == "Human") {
  data <- read.csv("experimenting/AnimalExtinctions/Humanc14Dates.csv", header = TRUE, na.strings = c("", "greater than"))
} else if(Species == "Equus") {
  # Need to force col classes as the spreadsheet has NA in it
  data <- read.csv("experimenting/AnimalExtinctions/Equusc14Dates.csv", header = TRUE, na.strings = c("", "greater than", "NA"),
                   colClasses = c("factor", "factor", "factor", "numeric", "numeric"))
} else if(Species == "Mammoth") {
  data <- read.csv("experimenting/AnimalExtinctions/Mammuthusc14Dates.csv", header = TRUE, na.strings = c("", "greater than"))
} else if(Species == "Alces") { # Moose/Elk
  data <- read.csv("experimenting/AnimalExtinctions/Alcesc14Dates.csv", header = TRUE, na.strings = c("", "greater than"))
} else if(Species == "Cervus") { # Dog
  data <- read.csv("experimenting/AnimalExtinctions/Cervusc14Dates.csv", header = TRUE, na.strings = c("", "greater than"))
} else if(Species == "Bison") { # Bison
  data <- read.csv("experimenting/AnimalExtinctions/Bisonc14Dates.csv", header = TRUE, na.strings = c("", "greater than"))
} else {stop("Unknown Species")}

keep <- which(data$X14C < cutoffages[2] & data$X14C > cutoffages[1])
data <- data[keep,]
removena <- which(is.na(data$X14C) | is.na(data$X1..Sigma))
if(length(removena) != 0) { # Remove any values where 14C or sigma is NA
  data <- data[-removena,]
}

rc_determinations <- data$X14C
rc_sigmas <- data$X1..Sigma

# rm(keep, data, removena)

###### Used if pass to PPcalibrate
calendar_age_range <- c(0,50000)
rate_s <- NA
rate_h <- NA

######
prior_n_internal_changepoints_lambda <- 5
k_max_internal_changepoints <- 30
rescale_factor_rev_jump <- 0.9
default_prior_h_rate <- 0.1
initial_n_internal_changepoints <- 10

n_iter <- 10000
n_thin <- 10
F14C_inputs <- FALSE
use_F14C_space <- TRUE

prior_h_shape <- NA
prior_h_rate <- NA

calendar_age_range <- NA # c(6680, 30000)
calendar_grid_resolution <- 10
grid_extension_factor <- 0.1
show_progress <- TRUE

#set.seed(14)
set.seed(19)

Test_Output <- PPcalibrate(
  rc_determinations = rc_determinations,
  rc_sigmas = rc_sigmas,
  calibration_curve = intcal20,
  F14C_inputs = F14C_inputs,
  n_iter = n_iter,
  n_thin = n_thin,
  use_F14C_space = use_F14C_space,
  prior_h_shape = prior_h_shape,
  prior_h_rate = prior_h_rate,
  show_progress = show_progress,
  calendar_age_range = calendar_age_range,
  calendar_grid_resolution = calendar_grid_resolution,
  prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
  k_max_internal_changepoints = k_max_internal_changepoints,
  rescale_factor_rev_jump = rescale_factor_rev_jump,
  default_prior_h_rate = default_prior_h_rate,
  initial_n_internal_changepoints = initial_n_internal_changepoints,
  grid_extension_factor = grid_extension_factor
  )

PlotPosteriorMeanRate(Test_Output)
PlotPosteriorChangePoints(Test_Output)
PlotPosteriorHeights(Test_Output)
