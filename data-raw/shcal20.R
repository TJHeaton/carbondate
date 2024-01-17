## code to prepare `shcal20` dataset goes here

shcal20 <- read.csv("data-raw/shcal20.14c", header = FALSE, skip=11)[,1:3]
names(shcal20) <- c(
  "calendar_age_BP",
  "c14_age",
  "c14_sig")

# add columns for F14C
shcal20$f14c <- exp(-shcal20$c14_age / 8033)
shcal20$f14c_sig <- shcal20$f14c * shcal20$c14_sig / 8033

usethis::use_data(shcal20, overwrite = TRUE)

