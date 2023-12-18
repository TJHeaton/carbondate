library(usethis)

intcal13 <- read.csv("data-raw/intcal13.14c", header = FALSE, skip=11)[,1:3]
names(intcal13) <- c(
  "calendar_age_BP",
  "c14_age",
  "c14_sig")

# add columns for F14C
intcal13$f14c <- exp(-intcal13$c14_age / 8033)
intcal13$f14c_sig <- intcal13$f14c * intcal13$c14_sig / 8033

usethis::use_data(intcal13, overwrite = TRUE)
