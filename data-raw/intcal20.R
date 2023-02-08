
library(usethis)

intcal20 <- read.csv("data-raw/intcal20.14c", header = FALSE, skip=11)[,1:3]
names(intcal20) <- c(
  "calendar_age",
  "c14_age",
  "c14_sig")

usethis::use_data(intcal20, overwrite = TRUE)
