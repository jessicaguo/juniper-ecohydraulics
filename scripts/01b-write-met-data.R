#### Write out flux and soil data for EWP ####
library(tidyverse)
library(janitor)

raw <- read_csv("data_raw/Flux-tower-data.csv")
names(raw) <- raw[1,]
flux <- raw[-(1:3),] |> 
  clean_names()


raw2 <- read_csv("data_raw/Other-tower-data.csv")
names(raw2) <- raw2[1,]
flux_other <- raw2[-(1:3),] |> 
  clean_names()
