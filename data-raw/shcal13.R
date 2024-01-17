## code to prepare `shcal13` dataset goes here

shcal13 <- read.csv("data-raw/shcal13.14c", header = FALSE, skip=11)[,1:3]
names(shcal13) <- c(
  "calendar_age_BP",
  "c14_age",
  "c14_sig")

# add columns for F14C
shcal13$f14c <- exp(-shcal13$c14_age / 8033)
shcal13$f14c_sig <- shcal13$f14c * shcal13$c14_sig / 8033

usethis::use_data(shcal13, overwrite = TRUE)

