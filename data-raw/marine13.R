library(usethis)

marine13 <- read.csv("data-raw/marine13.14c", header = FALSE, skip=11)[,1:3]
names(marine13) <- c(
  "calendar_age_BP",
  "c14_age",
  "c14_sig")

# add columns for F14C
marine13$f14c <- exp(-marine13$c14_age / 8033)
marine13$f14c_sig <- marine13$f14c * marine13$c14_sig / 8033

usethis::use_data(marine13, overwrite = TRUE)

