.KullbackLeiblerDivergenceFinal <- function(output_data_1, output_data_2) {
  min_calendar_age = min(output_data_1$calendar_ages, output_data_2$calendar_ages)
  max_calendar_age = max(output_data_1$calendar_ages, output_data_2$calendar_ages)
  calendar_age_sequence = seq(min_calendar_age, max_calendar_age, length.out=1001)

  pred_dens = PlotPredictiveCalendarAgeDensity(list(output_data_1, output_data_2), 5000)

  P = pred_dens[[1]]$density_mean
  Q = pred_dens[[2]]$density_mean

  divergence = .KLD(P, Q)
  title(main = "", sub = paste("KLD = ", signif(divergence, digits=4)))
  return (divergence)
}


.KLD <- function(P, Q) {
  P = P/sum(P)
  Q = Q/sum(Q)
  return (sum(P * log(P / Q)))
}
