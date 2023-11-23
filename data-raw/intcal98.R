library(usethis)

intcal98 <- read.csv("data-raw/intcal98.14c", header = FALSE, skip=16, sep="")[, c(1, 4, 5)]
names(intcal98) <- c(
  "calendar_age_BP",
  "c14_age",
  "c14_sig")

# convert to calBP
intcal98$calendar_age_BP <- 1950 - intcal98$calendar_age_BP

# add columns for F14C
intcal98$f14c <- exp(-intcal98$c14_age / 8033)
intcal98$f14c_sig <- intcal98$f14c * intcal98$c14_sig / 8033

usethis::use_data(intcal98, overwrite = TRUE)
