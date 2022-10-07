
library(usethis)

intcal20 <- read.table("data-raw/intcal20.14c", sep = ",", header = FALSE, skip = 11)
names(intcal20) <- c(
  "calendar_age",
  "c14_age",
  "c14_sig",
  "Delta14C",
  "DeltaSigma")

usethis::use_data(intcal20, overwrite = TRUE)
