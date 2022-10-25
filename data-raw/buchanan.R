library(usethis)

raw_data <- read.csv("data-raw/buchanan2008pde.csv", header = FALSE, sep = ",")
buchanan = data.frame(c14_ages = raw_data[, 3], c14_sig = raw_data[, 4])

usethis::use_data(buchanan, overwrite = TRUE)
