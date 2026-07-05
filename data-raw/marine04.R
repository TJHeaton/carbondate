library(usethis)

marine04 <- read.csv("data-raw/marine04.14c", header = FALSE, skip=11)[,1:3]
names(marine04) <- c(
  "calendar_age_BP",
  "c14_age",
  "c14_sig")

# add columns for F14C
marine04$f14c <- exp(-marine04$c14_age / 8033)
marine04$f14c_sig <- marine04$f14c * marine04$c14_sig / 8033

usethis::use_data(marine04, overwrite = TRUE)
