## code to prepare `shcal04` dataset goes here

shcal04 <- read.csv("data-raw/shcal04.14c", header = FALSE, skip=11)[,1:3]
names(shcal04) <- c(
  "calendar_age_BP",
  "c14_age",
  "c14_sig")

# add columns for F14C
shcal04$f14c <- exp(-shcal04$c14_age / 8033)
shcal04$f14c_sig <- shcal04$f14c * shcal04$c14_sig / 8033

usethis::use_data(shcal04, overwrite = TRUE)


