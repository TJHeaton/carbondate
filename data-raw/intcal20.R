library(usethis)

intcal20 <- read.csv("data-raw/intcal20.14c", header = FALSE, skip=11)[,1:3]
names(intcal20) <- c(
  "calendar_age_BP",
  "c14_age",
  "c14_sig")

# add columns for F14C
intcal20$f14c <- exp(-intcal20$c14_age / 8033)
intcal20$f14c_sig <- intcal20$f14c * intcal20$c14_sig / 8033

usethis::use_data(intcal20, overwrite = TRUE)
