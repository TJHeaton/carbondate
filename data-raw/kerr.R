library(usethis)

raw_data <- read.csv("data-raw/kerr2014sss_sup.csv", header = FALSE, sep = ",")
kerr = data.frame(c14_ages = raw_data[, 3], c14_sig = raw_data[, 4])

kerr$f14c = exp(-kerr$c14_ages / 8033)
kerr$f14c_sig = kerr$f14c * kerr$c14_sig / 8033

usethis::use_data(kerr, overwrite = TRUE)
