intcal04 <- read.csv("data-raw/intcal04.14c", header = FALSE, skip=11)[,1:3]
names(intcal04) <- c(
  "calendar_age",
  "c14_age",
  "c14_sig")

# Remove the last row with a fictitious point
intcal04 = head(intcal04, -1)
usethis::use_data(intcal04, overwrite = TRUE)

