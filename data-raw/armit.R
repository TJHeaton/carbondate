library(usethis)

raw_data <- read.csv("data-raw/armit2014rcc_sd01.csv", header = TRUE, sep = ",")

# Remove the two observations with missing xsig values
remove <- which(is.na(raw_data$error))
raw_data <- raw_data[-remove, ]

armit = data.frame(c14_ages = raw_data$X14C.age, c14_sig = raw_data$error)

armit$f14c = exp(-armit$c14_ages / 8033)
armit$f14c_sig = armit$f14c * armit$c14_sig / 8033

usethis::use_data(armit, overwrite = TRUE)
