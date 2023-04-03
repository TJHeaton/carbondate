library(usethis)

intcal98 <- read.csv("data-raw/intcal98.14c", header = FALSE, skip=16, sep="")[, c(1, 4, 5)]
names(intcal98) <- c(
  "calendar_age",
  "c14_age",
  "c14_sig")

# convert to calBP
intcal98$calendar_age = 1950 - intcal98$calendar_age

usethis::use_data(intcal98, overwrite = TRUE)
