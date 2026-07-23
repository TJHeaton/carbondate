library(usethis)

marine09 <- read.csv("data-raw/marine09.14c", header = FALSE, skip=11)[,1:3]
names(marine09) <- c(
  "calendar_age_BP",
  "c14_age",
  "c14_sig")

# add columns for F14C
marine09$f14c <- exp(-marine09$c14_age / 8033)
marine09$f14c_sig <- marine09$f14c * marine09$c14_sig / 8033

usethis::use_data(marine09, overwrite = TRUE)
