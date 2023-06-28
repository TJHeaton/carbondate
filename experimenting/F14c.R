set.seed(7)
walker_output_f14c = WalkerBivarDirichlet(
  kerr$c14_age, kerr$c14_sig, F14C_inputs = FALSE, intcal20, n_iter=1e5, n_thin=10, use_F14C_space = TRUE)
PlotPredictiveCalendarAgeDensity(walker_output_f14c, 5000, show_SPD = TRUE)
PlotPredictiveCalendarAgeDensity(walker_output_f14c, 5000, show_SPD = TRUE, plot_14C_age = FALSE)
PlotCalendarAgeDensityIndividualSample(10, walker_output_f14c)
PlotCalendarAgeDensityIndividualSample(10, walker_output_f14c, plot_14C_age = FALSE)

set.seed(7)
walker_output_c14BP = WalkerBivarDirichlet(
  kerr$c14_age, kerr$c14_sig, FALSE, intcal20, n_iter=1e5, n_thin=10, use_F14C_space = FALSE)
PlotPredictiveCalendarAgeDensity(walker_output_c14BP, 5000, show_SPD = TRUE)
PlotPredictiveCalendarAgeDensity(walker_output_c14BP, 5000, show_SPD = TRUE, plot_14C_age = FALSE)
PlotCalendarAgeDensityIndividualSample(10, walker_output_c14BP)
PlotCalendarAgeDensityIndividualSample(10, walker_output_c14BP, plot_14C_age = FALSE)

set.seed(7)
walker_output_f14c_inputs = WalkerBivarDirichlet(
  kerr$f14c, kerr$f14c_sig, TRUE, intcal20, n_iter=1e5, n_thin=10, use_F14C_space = TRUE)
PlotPredictiveCalendarAgeDensity(walker_output_f14c_inputs, 5000, show_SPD = TRUE)
PlotPredictiveCalendarAgeDensity(walker_output_f14c_inputs, 5000, show_SPD = TRUE, plot_14C_age = FALSE)
PlotCalendarAgeDensityIndividualSample(10, walker_output_f14c_inputs)
PlotCalendarAgeDensityIndividualSample(10, walker_output_f14c_inputs, plot_14C_age = FALSE)

