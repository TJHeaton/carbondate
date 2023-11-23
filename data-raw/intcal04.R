library(usethis)

intcal04 <- read.csv("data-raw/intcal04.14c", header = FALSE, skip=11)[,1:3]
names(intcal04) <- c(
  "calendar_age_BP",
  "c14_age",
  "c14_sig")

# Remove the last row with a fictitious point
intcal04 <- head(intcal04, -1)

# add columns for F14C
intcal04$f14c <- exp(-intcal04$c14_age / 8033)
intcal04$f14c_sig <- intcal04$f14c * intcal04$c14_sig / 8033

usethis::use_data(intcal04, overwrite = TRUE)

