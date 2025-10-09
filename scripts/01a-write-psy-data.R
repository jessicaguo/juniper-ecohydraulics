# This script uses 'suncalc' to calculate the dailiy water potential stats
# Predawn = mean(two hours prior to sunset)
# Midday = mean(two hours after solar noon)
# Daily summarized env variables are also calculated


source("juniper-ecohydraulics-main/source/round.POSIXct.R")
library(suncalc)
library(readr)
library(dplyr)
library(ggplot2)
library(plantecophys)
library(RColorBrewer)
library(tidyverse)

# Obtain sunrise and solar noon
# List 3 half-hours prior to sunrise_round and solarnoon_round
sun_dt <- getSunlightTimes(date = seq(as.Date("2021-01-01"), 
                                      as.Date("2021-12-31"),
                                      by = "day"),
                           lat = 	37.5241,
                           lon = -109.7471,
                           keep = c("sunrise", "solarNoon"),
                           tz = "America/Denver") %>%
  mutate(sunrise_0 = round.POSIXct(sunrise, "30 mins"),
         solarnoon_0 = round.POSIXct(solarNoon, "30 mins")) %>%
  select(-lat, -lon, -sunrise, -solarNoon) %>%
  mutate(sunrise_1 = sunrise_0 - 1*30*60,
         sunrise_2 = sunrise_0 - 2*30*60,
         sunrise_3 = sunrise_0 - 3*30*60,
         solarnoon_1 = solarnoon_0 - 1*30*60,
         solarnoon_2 = solarnoon_0 - 2*30*60,
         solarnoon_3 = solarnoon_0 - 3*30*60) %>%
  tidyr::pivot_longer(-date, 
                      names_to = c("type", "timing"),
                      names_pattern = "(.*)_(.)",
                      values_to = "time") %>%
  arrange(date, type, timing)
rm(round.POSIXct)

# Psychrometer (cleaned)
psy <- read_csv(file = "juniper-ecohydraulics-main/data_cleaned/psy_hourly.csv",
                locale = locale(tz = "America/Denver")) %>%
  mutate(Tree = as.factor(Tree),
         Logger = as.factor(Logger))

# Join together, summarize by tree, logger, date, type
psy_daily <- inner_join(psy, sun_dt, by = c("dt" = "time")) %>%
  group_by(Tree, Logger, date, type) %>%
  summarize(WP = mean(psy),
            WP.SD = sd(psy),
            n = n()) %>%
  mutate(type = case_when(type == "solarnoon" ~ "MD", 
                          type == "sunrise" ~ "PD"),
         doy = lubridate::yday(date))

# Making more psy dfs
psy_30_branch <- psy |> 
  mutate(Logger = case_when(Logger == 1 ~ "a",
                            Logger == 2 ~ "b"),
         branch = paste0(Tree, Logger),
         Logger = case_when(Logger == "a" ~ 1,
                            Logger == "b" ~ 2))

psy_30_shrub <- psy |> 
  group_by(Tree, dt) |> 
  summarize(WP = mean(psy),
            WP_sd = sd(psy),
            WP_n = n()) |> 
  ungroup()

psy_30_site <- psy_30_branch |> 
  group_by(dt) |> 
  summarize(WP = mean(psy),
            WP_sd = sd(psy),
            WP_n = n()) |> 
  ungroup()

psy_daily_branch <- psy_daily |> 
  mutate(Logger = case_when(Logger == 1 ~ "a",
                            Logger == 2 ~ "b"),
         branch = paste0(Tree, Logger),
         Logger = case_when(Logger == "a" ~ 1,
                            Logger == "b" ~ 2)) |> 
  pivot_wider(names_from = type,
              values_from = c(WP, WP.SD, n)) |> 
  rename(MD = WP_MD, PD = WP_PD, MD_sd = WP.SD_MD, PD_sd = WP.SD_PD, MD_n = n_MD, PD_n = n_PD)

psy_daily_shrub <- psy_daily_branch |> 
  group_by(Tree, date) |> 
  summarize(PD_m = mean(PD, na.rm = TRUE),
            PD_sd = sd(PD, na.rm = TRUE),
            PD_n = n(),
            MD_m = mean(MD, na.rm = TRUE),
            MD_sd = sd(MD, na.rm = TRUE),
            MD_n = n()) |> 
  ungroup()

psy_daily_site <- psy_daily_shrub |> 
  group_by(date) |> 
  summarize(PD = mean(PD_m, na.rm = TRUE),
            PD_sd = sd(PD_m, na.rm = TRUE),
            PD_n = n(),
            MD = mean(MD_m, na.rm = TRUE),
            MD_sd = sd(MD_m, na.rm = TRUE),
            MD_n = n()) |> 
  ungroup()

# Write out datasets
write_csv(psy_30_branch, "juniper-ecohydraulics-main/data_cleaned/psy_30_branch.csv")
write_csv(psy_30_shrub, "juniper-ecohydraulics-main/data_cleaned/psy_30_shrub.csv")
write_csv(psy_30_site, "juniper-ecohydraulics-main/data_cleaned/psy_30_site.csv")
write_csv(psy_daily_branch, "juniper-ecohydraulics-main/data_cleaned/psy_daily_branch.csv")
write_csv(psy_daily_shrub, "juniper-ecohydraulics-main/data_cleaned/psy_daily_shrub.csv")
write_csv(psy_daily_site, "juniper-ecohydraulics-main/data_cleaned/psy_daily_site.csv")

# Met data
col_names <- names(read_csv("juniper-ecohydraulics-main/data_raw/Other-tower-data.csv",
                            skip = 1, n_max = 0))

met <- read_csv("juniper-ecohydraulics-main/data_raw/Other-tower-data.csv",
                skip = 4, 
                col_names = col_names, 
                na = c("-9999", ""),
                col_types = cols(TIMESTAMP = "c",
                                 .default = "d")) %>%
  mutate(dt = as.POSIXct(TIMESTAMP, format = "%m/%d/%Y %H:%M", 
                         tz = "America/Denver"),
         date = as.Date(dt),
         doy = lubridate::yday(dt),
         VPD_Avg = RHtoVPD(RH_Avg, AirTemp_Avg)) %>%
  select(date, doy, dt, AirTemp_Avg, RH_Avg, VPD_Avg, 
         Precip_Tot, PAR_In_Avg, contains("VWC")) %>%
  distinct()

# Summarize to daily
met_daily <- met %>%
  group_by(date) %>%
  summarize(n = n(),
            VPD_mean = mean(VPD_Avg, na.rm = TRUE),
            VPD_max = max(VPD_Avg, na.rm = TRUE),
            T_min = min(AirTemp_Avg, na.rm = TRUE),
            T_mean = mean(AirTemp_Avg, na.rm = TRUE),
            T_max = max(AirTemp_Avg, na.rm = TRUE),
            Precip = sum(Precip_Tot, na.rm = TRUE),
            PAR_In = sum(PAR_In_Avg, na.rm = TRUE),
            VWC_5cm = mean(VWC_5cm_Avg, na.rm = TRUE),
            VWC_10cm = mean(VWC_10cm_Avg, na.rm = TRUE),
            VWC_20cm = mean(VWC_20cm_Avg, na.rm = TRUE),
            VWC_50cm = mean(VWC_50cm_Avg, na.rm = TRUE),
            VWC_100cm = mean(VWC_100cm_Avg, na.rm = TRUE)) %>%
  arrange(date)


write_csv(met_daily, "juniper-ecohydraulics-main/data_cleaned/atmdaily.csv")
write_csv(met, "juniper-ecohydraulics-main/data_cleaned/vpd_ppt_30.csv")
