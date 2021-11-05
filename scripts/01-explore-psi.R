# Explore and plot water potentials
library(tidyverse)
source("source/round.POSIXct.R")

# read in data
# psychrometer dataset
psy <- read.csv("data_raw/Merged-psychrometer-data.csv",  encoding="UTF-8") %>%
  rename(chamber_temp = 5,
         dT = 6,
         wet_bulb = 7,
         psy = 8,
         correction_dT = 12,
         correction = 13,
         batt_volt = 14,
         int_batt_temp = 15,
         ext_power = 16,
         ext_power_volt = 17,
         ext_power_current = 18,
         comment = 19) %>%
  mutate(Date = as.POSIXct(Date, format = "%d/%m/%Y", tz = "America/Denver"),
         dt = as.POSIXct(paste0(Date, Time), tz = "America/Denver"))

# read in maintenance days
maint <- read_csv("data_raw/Maintenance-times.csv")

# manual dataset
man <- read_csv("data_raw/Pressure-chamber-data.csv") %>%
  rename(dt = 1,
         WP_bar = 4,
         WP_mpa = 5,
         Psy = 6) %>%
  mutate(dt = as.POSIXct(dt, format = "%m/%d/%Y %H:%M", tz = "America/Denver"),
         round_dt = round.POSIXct(dt, units = "30 mins")) %>%
  left_join(psy, by = c("round_dt" = "dt", 
                        "Tree" = "Tree", 
                        "Logger" = "Logger")) %>%
  mutate(hour = as.numeric(difftime(dt, Date, units = "hours")),
         Date = as.Date(dt))

m1 <- lm(psy ~ WP_mpa, data = man)
summary(m1)

ggplot(man, aes(x = WP_mpa, y = psy)) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_abline(slope = summary(m1)$coef[2, 1],
              intercept = summary(m1)$coef[1,1],
              lty = 2) +
  geom_point(aes(col = hour)) +
  scale_color_gradient(low = "purple", high = "limegreen") +
  # facet_wrap(~Date) +
  scale_x_continuous(expression(paste("Pressure chamber ",  Psi[leaf], " (MPa)"))) +
  scale_y_continuous(expression(paste("Psychrometer ", Psi[stem], " (MPa)"))) +
  coord_equal() + 
  theme_bw()
