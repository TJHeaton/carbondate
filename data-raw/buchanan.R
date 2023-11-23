library(usethis)

raw_data <- read.csv("data-raw/buchanan2008pde.csv", header = FALSE, sep = ",")
buchanan <- data.frame(c14_age = raw_data[, 3], c14_sig = raw_data[, 4])

buchanan$f14c <- exp(-buchanan$c14_age / 8033)
buchanan$f14c_sig <- buchanan$f14c * buchanan$c14_sig / 8033

usethis::use_data(buchanan, overwrite = TRUE)
