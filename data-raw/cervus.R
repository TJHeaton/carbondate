# 14C dates for various Pleistocene Animals from Guthrie
# https://www.nature.com/articles/nature04604#Sec3

raw_data <- read.csv("data-raw/Cervusc14Dates.csv",
                     header = TRUE,
                     na.strings = c("", "greater than"))

source("data-raw/cutoff_ages.R")

keep <- which(raw_data$X14C < cutoff_ages[2] & raw_data$X14C > cutoff_ages[1])
raw_data <- raw_data[keep,]
removena <- which(is.na(raw_data$X14C) | is.na(raw_data$X1..Sigma))
if(length(removena) != 0) { # Remove any values where 14C or sigma is NA
  raw_data <- raw_data[-removena,]
}

cervus <- data.frame(lab_code = raw_data$X14C.Lab.No.,
                    site_code = raw_data$Site.Code,
                    location = raw_data$Locality,
                    c14_age = raw_data$X14C,
                    c14_sig = raw_data$X1..Sigma)

cervus$f14c <- exp(-cervus$c14_age / 8033)
cervus$f14c_sig <- cervus$f14c * cervus$c14_sig / 8033

usethis::use_data(cervus, overwrite = TRUE)
