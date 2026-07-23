library(usethis)

marine20 <- read.csv("data-raw/marine20.14c", header = FALSE, skip=11)[,1:3]
names(marine20) <- c(
  "calendar_age_BP",
  "c14_age",
  "c14_sig")

# add columns for F14C
marine20$f14c <- exp(-marine20$c14_age / 8033)
marine20$f14c_sig <- marine20$f14c * marine20$c14_sig / 8033

usethis::use_data(marine20, overwrite = TRUE)
