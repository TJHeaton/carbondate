intcal09 <- read.csv("data-raw/intcal09.14c", header = FALSE, skip=11)[,1:3]
names(intcal09) <- c(
  "calendar_age",
  "c14_age",
  "c14_sig")

# Remove the last row with a fictitious point
intcal09 = head(intcal09, -1)

usethis::use_data(intcal09, overwrite = TRUE)
