intcal13 <- read.csv("data-raw/intcal13.14c", header = FALSE, skip=11)[,1:3]
names(intcal13) <- c(
  "calendar_age",
  "c14_age",
  "c14_sig")

usethis::use_data(intcal13, overwrite = TRUE)
