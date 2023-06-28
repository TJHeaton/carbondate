library(usethis)

intcal09 <- read.csv("data-raw/intcal09.14c", header = FALSE, skip=11)[,1:3]
names(intcal09) <- c(
  "calendar_age_BP",
  "c14_age",
  "c14_sig")

# Remove the last row with a fictitious point
intcal09 = head(intcal09, -1)

# add columns for F14C
intcal09$f14c = exp(-intcal09$c14_age / 8033)
intcal09$f14c_sig = intcal09$f14c * intcal09$c14_sig / 8033

usethis::use_data(intcal09, overwrite = TRUE)
