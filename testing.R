library(devtools)
load_all()
measurements <- read.csv('../custom calcurve/8836_modern.csv', header = TRUE, sep=",")
HOBS2022 <- read.csv('../custom calcurve/HOBS2022.14c', header = FALSE, skip=2, sep="")[, c(1, 4, 5)]
names(HOBS2022) <- c("calendar_age", "c14_age", "c14_sig")

curve_1 = HOBS2022
curve_2 = HOBS2022
curve_3 = HOBS2022

curve_1$calendar_age = 2050 - curve_1$calendar_age
curve_2$calendar_age = 2020 - curve_2$calendar_age

interp_1 = InterpolateCalibrationCurve(NA, curve_1)
interp_2 = InterpolateCalibrationCurve(NA, curve_2)
interp_3 = InterpolateCalibrationCurve(NA, curve_3)

# set.seed(7)
#
# walker_output = WalkerBivarDirichlet(
#   measurements$F14C, measurements$F14C_sd, curve_1, 1e5, 10, slice_width = 70, slice_multiplier = 10, n_clust = 8)

set.seed(7)

walker_output_correct = WalkerBivarDirichlet(
  measurements$F14C, measurements$F14C_sd, curve_2, 1e5, 10, slice_width = 70, slice_multiplier = 10, n_clust = 8)
