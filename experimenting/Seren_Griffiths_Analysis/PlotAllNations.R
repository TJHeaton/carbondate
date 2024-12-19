# Plot all nations on a comparable plotting scale

nation_data <- data.frame(
  nation = c("England", "Scotland", "Wales", "Ireland"),
  n_samples = c(10435, 6471 , 2571, 4851)
)

# Choose a maximum rate scale which is proportional to number of samples
max_rate_scale_England <- 8
n_samples_England <- nation_data[nation_data$nation == "England", "n_samples"]

nation_data$max_rate_scale <- max_rate_scale_England * nation_data$n_samples/n_samples_England

# Run each of combined plots with appropriate max rate plotting scale
source("experimenting/Seren_Griffiths_Analysis/CombineIndependentRuns.R")

# England run
max_rate_scale <- 1.5 * nation_data[1, "max_rate_scale"]
source("experimenting/Seren_Griffiths_Analysis/CombineEnglandRuns.R")

# Scotland run
max_rate_scale <- 1.5 * nation_data[2, "max_rate_scale"]
source("experimenting/Seren_Griffiths_Analysis/CombineScotlandRuns.R")

# Wales run
max_rate_scale <- 1.5 * nation_data[3, "max_rate_scale"]
source("experimenting/Seren_Griffiths_Analysis/CombineWalesRuns.R")

# Ireland run
max_rate_scale <- 1.5 * nation_data[4, "max_rate_scale"]
source("experimenting/Seren_Griffiths_Analysis/CombineIrelandRuns.R")



##### Now create a combined plot with everyting on it


all_nations_postrate <- list()
all_nations_postrate[[1]] <- Holocene_England_PostMeanRate
all_nations_postrate[[2]] <- Holocene_Scotland_PostMeanRate
all_nations_postrate[[3]] <- Holocene_Wales_PostMeanRate
all_nations_postrate[[4]] <- Holocene_Republic_Ireland_PostMeanRate

library(RColorBrewer)
nation_data$plot_cols <- brewer.pal(nrow(nation_data), "Set1")

pdf("../SerenGriffithsData/OutputPlots/PriorMean_10_Internal_Changes/Combined/AllNationsPosteriorMeanRates.pdf",
    width = 12,
    height = 8)
# Now plot all posterior rate curves on the same plot
par(mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 2, 4, 2) + 0.1,
    las = 0)
xlim <- 1950 - rev(range(Holocene_Wales_PostMeanRate$calendar_age_BP))
ylim <- c(0, max_rate_scale)
plot(
  NULL,
  NULL,
  type = "n",
  ylim = ylim,
  xlim = xlim,
  axes = FALSE,
  xlab = "Calendar Age (cal AD)",
  ylab = NA,
  xaxs = "i",
  yaxs = "i")
mtext("Posterior Sample Occurrence Rate (Relative)", side = 2, line = 1)
# Add ticks
box()
axis(1)
xticks_minor <- seq(-11000, 2000, by = 100)
axis(1, xticks_minor, labels = FALSE, tcl = -0.2)
xticks_major <- seq(-11000, 2000, by = 1000)
axis(1, xticks_major, labels = FALSE, tcl = -0.5)




# Overlay various nations
for(nation_id in 1:nrow(nation_data)) {

  par(new = TRUE,
      mgp = c(3, 0.7, 0),
      xaxs = "i",
      yaxs = "i",
      mar = c(5, 2, 4, 2) + 0.1)
  xlim <- 1950 - rev(range(all_nations_postrate[[nation_id]]$calendar_age_BP))
  ylim <- c(0, nation_data[nation_id, "max_rate_scale"])
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



  lines(x = 1950 - all_nations_postrate[[nation_id]]$calendar_age_BP,
        y = all_nations_postrate[[nation_id]]$rate_mean,
        col = nation_data$plot_cols[nation_id],
        lwd = 2,
        type = "l")
  # lines(x = 1950 - all_nations_postrate[[nation_id]]$calendar_age_BP,
  #       y = all_nations_postrate[[nation_id]]$rate_ci_lower,
  #       col = nation_data$plot_cols[nation_id],
  #       lty = 2,
  #       type = "l")
  # lines(x = 1950 - all_nations_postrate[[nation_id]]$calendar_age_BP,
  #       y = all_nations_postrate[[nation_id]]$rate_ci_upper,
  #       col = nation_data$plot_cols[nation_id],
  #       lty = 2,
  #       type = "l")
}

leg.txt <- paste(nation_data$nation, " (n = ", nation_data$n_samples, ")", sep = "")

legend("topleft",
       legend = leg.txt,
       lty = 1,
       col = nation_data$plot_cols)

dev.off()

