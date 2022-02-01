# This script uses 'suncalc' to calculate the dailiy water potential stats
# Predawn = mean(two hours prior to sunset)
# Midday = mean(two hours after solar noon)
# Daily summarized env variables are also calculated


source("source/round.POSIXct.R")
library(suncalc)
library(readr)
library(dplyr)
library(ggplot2)
library(plantecophys)
library(RColorBrewer)

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


# Psychrometer (cleaned)
psy <- read_csv(file = "data_cleaned/psy_hourly.csv",
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
                          type == "sunrise" ~ "PD"))

# Quick plotting
psy_daily %>%
  select(Tree, Logger, date, type, WP) %>%
  tidyr::pivot_wider(names_from = type, 
                     values_from = WP) %>%
  mutate(month = lubridate::month(date)) %>%
  ggplot(aes(x = PD, y = MD, col = month)) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_point() +
  # facet_wrap(~Tree) +
  theme_bw()

# Met data
col_names <- names(read_csv("data_raw/Other-tower-data.csv",
                            skip = 1, n_max = 0))

met <- read_csv("data_raw/Other-tower-data.csv",
                skip = 4, 
                col_names = col_names, 
                na = c("-9999")) %>%
  mutate(dt = as.POSIXct(TIMESTAMP, format = "%m/%d/%Y %H:%M", 
                         tz = "America/Denver"),
         date = lubridate::date(dt),
         VPD_Avg = RHtoVPD(RH_Avg, AirTemp_Avg)) %>%
  select(date, dt, AirTemp_Avg, RH_Avg, VPD_Avg, 
         Precip_Tot, contains("VWC")) %>%
  filter(date >= as.Date("2021-01-01"))

# Summarize to daily
met_daily <- met %>%
  group_by(`date`) %>%
  summarize(VPD_mean = mean(VPD_Avg, na.rm = TRUE),
            VPD_max = max(VPD_Avg, na.rm = TRUE),
            T_min = min(AirTemp_Avg, na.rm = TRUE),
            T_mean = mean(AirTemp_Avg, na.rm = TRUE),
            T_max = max(AirTemp_Avg, na.rm = TRUE),
            Precip = sum(Precip_Tot, na.rm = TRUE),
            VWC_5cm = mean(VWC_5cm_Avg, na.rm = TRUE),
            VWC_10cm = mean(VWC_10cm_Avg, na.rm = TRUE),
            VWC_20cm = mean(VWC_20cm_Avg, na.rm = TRUE),
            VWC_50cm = mean(VWC_50cm_Avg, na.rm = TRUE),
            VWC_100cm = mean(VWC_100cm_Avg, na.rm = TRUE),
            n = n())

# Quick plots
met_daily %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = VPD_mean, color = "mean")) +
  geom_point(aes(y = VPD_max, color = "max")) +
  geom_bar(aes(y = Precip/6, color = "precip"),
           stat = "identity") +
  scale_y_continuous("VPD (kPa)") +
  theme_bw()

met_daily %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = T_min, color = "min")) +
  geom_point(aes(y = T_mean, color = "mean")) +
  geom_point(aes(y = T_max, color = "max")) +
  geom_bar(aes(y = Precip, color = "precip"),
           stat = "identity") +
  scale_y_continuous("Air Temp (C)") +
  theme_bw()

met_daily %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = VWC_5cm, color = "5")) +
  geom_point(aes(y = VWC_10cm, color = "10")) +
  geom_point(aes(y = VWC_20cm, color = "20")) +
  geom_point(aes(y = VWC_50cm, color = "50")) +
  geom_point(aes(y = VWC_100cm, color = "100")) +
  geom_bar(aes(y = Precip/2, color = "precip"),
           stat = "identity") +
  scale_y_continuous("VWC (%)") +
  scale_color_brewer(palette = "Spectral",
                     breaks = c("5", "10", "20", "50", "100")) +
  theme_bw()

# Write out datasets
save(met_daily, file = "data_cleaned/met_daily.Rdata")
save(psy_daily, file = "data_cleaned/psy_daily.Rdata")
