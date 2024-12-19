# Choose fixed value for maximum ylim for plotted posterior rate
#  by setting max_rate_scale, e.g., max_rate_scale <- 15

# Combine Ireland Output
load("../SerenGriffithsData/RWorkspaces/PriorMean10Changes/Holocene_Ireland_Analysis_seed51.RData")
Ireland_run_seed_51 <- Holocene_Republic_Ireland_Output

load("../SerenGriffithsData/RWorkspaces/PriorMean10Changes/Holocene_Ireland_Analysis_seed58.RData")
Ireland_run_seed_58 <- Holocene_Republic_Ireland_Output

load("../SerenGriffithsData/RWorkspaces/PriorMean10Changes/Holocene_Ireland_Analysis_seed85.RData")
Ireland_run_seed_85 <- Holocene_Republic_Ireland_Output

output_data_list <- list(Ireland_run_seed_51,
                         Ireland_run_seed_58,
                         Ireland_run_seed_85)

Holocene_Ireland_Output <- CombineRuns(output_data_list)



pdf("../SerenGriffithsData/OutputPlots/PriorMean_10_Internal_Changes/Combined/Ireland_PP_Rate_Output_Holocene_Combined.pdf",
    width = 12,
    height = 8)
Holocene_Ireland_PostMeanRate <- PlotPosteriorMeanRate(Holocene_Ireland_Output,
                                                                show_individual_means = FALSE,
                                                                denscale = 2,
                                                                n_posterior_samples = 15000,
                                                                plot_cal_age_scale = "AD",
                                                                max_rate_scale = max_rate_scale)



par(new = TRUE,
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 4, 2) + 0.1,
    las = 1)
xlim <- 1950 - rev(range(Holocene_Ireland_PostMeanRate$calendar_age_BP))
ylim <- c(0, 1)
plot(
  NULL,
  NULL,
  type = "n",
  ylim = ylim,
  xlim = xlim,
  axes = FALSE,
  xlab = NA,
  ylab = NA,
  xaxs = "i",
  yaxs = "i")

xticks_minor <- seq(-11000, 2000, by = 100)
axis(1, xticks_minor, labels = FALSE, tcl = -0.2)
xticks_major <- seq(-11000, 2000, by = 1000)
axis(1, xticks_major, labels = FALSE, tcl = -0.5)

Time_Periods <- read.csv("../SerenGriffithsData/TimePeriods/Ireland_Time_Periods.csv",
                         header = TRUE)
Period_Start <- Time_Periods$Island.of.Ireland.start.dates..BC.AD
Period_End <- Time_Periods$Island.of.Ireland.end.dates..BC.AD

n_period <- length(Period_End)

library(RColorBrewer)
cols <- brewer.pal(n_period, "Set3")

add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)
cols <- add.alpha(cols, alpha = 0.4)

for(i in 1:n_period) {
  polygon(x = c(rep(Period_Start[i], 2), rep(Period_End[i], 2)),
          y = c(0, 0.82, 0.82, 0), border = NA,
          col = cols[i])
}

mtext("Island of Ireland", side = 3, cex = 1.2)

# Replot the posterior mean rate over the top
par(new = TRUE,
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 4, 2) + 0.1,
    las = 1)
xlim <- 1950 - rev(range(Holocene_Ireland_PostMeanRate$calendar_age_BP))
ylim <- c(0, max_rate_scale)
plot(
  NULL,
  NULL,
  type = "n",
  ylim = ylim,
  xlim = xlim,
  axes = FALSE,
  xlab = NA,
  ylab = NA,
  xaxs = "i",
  yaxs = "i")

lines(x = 1950 - Holocene_Ireland_PostMeanRate$calendar_age_BP,
      y = Holocene_Ireland_PostMeanRate$rate_mean,
      col = "purple", type = "l")
lines(x = 1950 - Holocene_Ireland_PostMeanRate$calendar_age_BP,
      y = Holocene_Ireland_PostMeanRate$rate_ci_lower,
      col = "purple", type = "l", lty = 2)
lines(x = 1950 - Holocene_Ireland_PostMeanRate$calendar_age_BP,
      y = Holocene_Ireland_PostMeanRate$rate_ci_upper,
      col = "purple", type = "l", lty = 2)



dev.off()

# Also write the output of the posterior rate to a file
output_filename <- paste("../SerenGriffithsData/OutputPlots/PriorMean_10_Internal_Changes/Combined/",
                         "Ireland",
                         "EstimatedPosteriorRateCombined.csv", sep = "")
write.csv(Holocene_Ireland_PostMeanRate, output_filename, row.names = FALSE)
