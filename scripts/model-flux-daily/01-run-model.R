# Univariate first, then multivariate model of fluxes
# on PD, soil moisture, light, temp, and VPD
# or some combo thereof

library(rjags)
load.module('dic')
library(dplyr)

# Read in daily fluxes
flux <- readr::read_csv("data_raw/US-CdM daily.csv") %>%
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j"))) %>%
  relocate(date)

