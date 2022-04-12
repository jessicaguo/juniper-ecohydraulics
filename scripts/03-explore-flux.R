# Explore flux tower lags with VPD and soil moisture

library(readr)
library(dplyr)
library(ggplot2)

# Load daily met data, 
load(file = "data_cleaned/met_daily.Rdata")

# Import flux data
cnames <- read_csv("data_raw/Flux-tower-data.csv",
                     skip = 1,
                     n_max = 0)
flux <- read_csv("data_raw/Flux-tower-data.csv",
                 skip = 4,
                 col_names = colnames(cnames),
                 na = c("", "NA", "INF", "NAN")) %>%
  mutate(TIMESTAMP = as.POSIXct(TIMESTAMP, 
                                format = "%m/%d/%Y %H:%M",
                                tz = "America/Denver"),
         date = as.Date(TIMESTAMP,
                        tz = "America/Denver")) %>%
  relocate(date) %>%
  distinct() # Remove accidentally duplicated rows

# Summarize to daily - sum for ET
daily_flux <- flux %>%
  filter(!is.na(ET)) %>%
  group_by(date) %>%
  summarize(n = n(),
            ET = sum(ET))

ggplot(daily_flux, aes(x = date, y = ET)) +
  geom_point(aes(color = n))

# Which flux values to keep? Depends on number of daily samples
hist(daily_flux$n,
     breaks = 48)

daily_flux %>%
  filter(n > 38) %>%
  nrow()

# Set threshold: at least 80% of 48 values present
flux_daily <- daily_flux %>%
  filter(n > 38)

# Which met values to keep? Remove lack of date
met_daily <- met_daily %>%
  filter(!is.na(date))

# Match and plot data
range(flux_daily$date)
range(met_daily$date)

flux_daily <- flux_daily %>%
  left_join(met_daily, by = c("date")) %>%
  rename(n.flux = n.x,
         n.met = n.y)

ggplot(flux_daily) +
  geom_point(aes(x = date, y = ET, col = "ET")) +
  geom_point(aes(x = date, y = VPD_max, col = "VPD_max")) +
  geom_point(aes(x = date, y = VWC_10cm, col = "VWC 10 cm")) +
  geom_point(aes(x = date, y = VWC_50cm, col = "VWC 50 cm")) +
  scale_y_continuous(expression(paste("ET (mm ", day^-1, ")"))) +
  scale_color_manual(values = c("cornflowerblue",
                                "coral",
                                "sienna",
                                "brown"))

ggplot(flux_daily, aes(x = VPD_max, y = ET)) +
  geom_point(aes(col = VWC_10cm)) 
# VPD drives ET much more strongly when shallow soil is wet


ggplot(flux_daily, aes(x = VPD_max, y = ET)) +
  geom_point(aes(col = VWC_50cm)) 
# VPD drives ET much more variables with deep soil is not wet
