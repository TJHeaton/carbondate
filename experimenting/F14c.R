set.seed(7)
walker_output_f14c = PolyaUrnBivarDirichlet(
  kerr$c14_ages, kerr$c14_sig, F14C_inputs = FALSE, intcal20, n_iter=1e5, n_thin=10, use_F14C_space = TRUE)
PlotPredictiveCalendarAgeDensity(walker_output_f14c, 5000)

set.seed(7)
walker_output_c14BP = PolyaUrnBivarDirichlet(
  kerr$c14_ages, kerr$c14_sig, FALSE, intcal20, n_iter=1e5, n_thin=10, use_F14C_space = FALSE)
PlotPredictiveCalendarAgeDensity(walker_output_c14BP, 5000)

set.seed(7)
walker_output_f14c_inputs = PolyaUrnBivarDirichlet(
  kerr$f14c, kerr$f14c_sig, TRUE, intcal20, n_iter=1e5, n_thin=10, use_F14C_space = TRUE)
PlotPredictiveCalendarAgeDensity(walker_output_f14c_inputs, 5000)

PlotCalendarAgeDensityIndividualSample(10, walker_output_c14BP)
PlotCalendarAgeDensityIndividualSample(10, walker_output_f14c_inputs)


KullbackLeiblerDivergence(
  WalkerBivarDirichlet(kerr$c14_ages, kerr$c14_sig, FALSE, intcal20, n_iter=1e4, n_thin=10),
  WalkerBivarDirichlet(kerr$c14_ages, kerr$c14_sig, FALSE, intcal20, n_iter=1e4, n_thin=10))
