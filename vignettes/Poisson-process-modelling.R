## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(carbondate)
set.seed(15)

## ----example, out.width="100%", fig.width=10, fig.height=8--------------------
# Set initial values 
n_observed <- 40
n_iter <- 100000
n_thin <- 10
rc_sigmas <- rep(15, n_observed)
calendar_grid_resolution <- 1

# Create artificial rc_determinations
calendar_age_range <- c(400, 700) 
observed_age_range <- c(500, 550) 
true_theta <- seq(from = observed_age_range[1],
                  to = observed_age_range[2],
                  length = n_observed)
intcal_mean <- approx(x = intcal20$calendar_age_BP,
                      y = intcal20$c14_age,
                      xout = true_theta)$y
intcal_sd <- approx(x = intcal20$calendar_age_BP,
                    y = intcal20$c14_sig,
                    xout = true_theta)$y
rc_determinations <- rnorm(n = n_observed,
                                  mean = intcal_mean,
                                  sd = sqrt(rc_sigmas^2 + intcal_sd^2))

# Fit the model
PP_fit_output <- PPcalibrate(
    rc_determinations = rc_determinations,
    rc_sigmas = rc_sigmas,
    calibration_curve = intcal20,
    calendar_age_range = calendar_age_range,
    calendar_grid_resolution = calendar_grid_resolution,
    n_iter = n_iter,
    n_thin = n_thin,
    show_progress = FALSE)

## ----meanrate, out.width="100%", fig.width=10, fig.height=8-------------------
PlotPosteriorMeanRate(PP_fit_output)

## ----changepoint number, out.width="100%", fig.width=10, fig.height=8---------
PlotNumberOfInternalChanges(PP_fit_output)

## ----changepoint locations, out.width="100%", fig.width=10, fig.height=8------
PlotPosteriorChangePoints(PP_fit_output)

## ----changepoint rates, out.width="100%", fig.width=10, fig.height=8----------
PlotPosteriorHeights(PP_fit_output)

