load("../SerenGriffithsData/RWorkspaces/PriorMean10Changes/Holocene_England_Analysis_seed95_kmax_40.RData")


pdf("../SerenGriffithsData/OutputPlots/PriorMean_10_Internal_Changes/England_PP_Rate_Output_Holocene_seed55.pdf",
    width = 12,
    height = 8)


Holocene_England_PostMeanRate <- PlotPosteriorMeanRate(Holocene_England_Output,
                                                       show_individual_means = FALSE,
                                                       denscale = 2,
                                                       plot_cal_age_scale = "AD")



par(new = TRUE,
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 4, 2) + 0.1,
    las = 1)

xlim <- 1950 - rev(range(Holocene_England_PostMeanRate$calendar_age_BP))
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

Time_Periods <- read.csv("../SerenGriffithsData/TimePeriods/England_Time_Periods.csv",
                         header = TRUE)
Period_Start <- Time_Periods$England.start.dates..BC.AD
Period_End <- Time_Periods$England.end.dates..BC.AD

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

mtext("England", side = 3, cex = 1.2)


dev.off()



