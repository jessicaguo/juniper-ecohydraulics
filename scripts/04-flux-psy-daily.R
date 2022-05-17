# Calculate averaged daily pd and md
# Plot with fluxes

library(readr)
library(dplyr)
library(ggplot2)

# Load previuosly cleaned PD and MD
load("data_cleaned/psy_daily.Rdata")

# Summarize to single site level Psy
psy_daily_site <- psy_daily %>%
  group_by(date, type) %>%
  summarize(WP_mean = mean(WP),
            WP_sd = sd(WP),
            n = n())

# Read in daily fluxes
flux <- read_csv("data_raw/US-CdM daily.csv") %>%
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j"))) %>%
  relocate(date)

# Plot whole record of C fluxes, highlight psychrometer period
ggplot(flux, aes(x = date)) +
  # geom_point(aes(y = ET, color = "ET")) +
  geom_rect(aes(ymin = -Inf, ymax = Inf,
                xmin = min(psy_daily$date), 
                xmax = max(psy_daily$date)),
            fill = "gray") +
  geom_point(aes(y = NEE_F, color = "NEE")) +
  geom_point(aes(y = GPP_F, color = "GPP")) +
  geom_point(aes(y = RE_F, color = "RE")) +
  geom_bar(aes(y = P/50, color = "ppt"), stat = "identity") +
  theme_bw()

# Plot whole record of ET flux, highlight psychrometer period
ggplot(flux, aes(x = date)) +
  geom_rect(aes(ymin = -Inf, ymax = Inf,
                xmin = min(psy_daily$date), 
                xmax = max(psy_daily$date)),
            fill = "gray") +
  geom_point(aes(y = ET, color = "ET")) +
  geom_bar(aes(y = P/5, color = "ppt"), stat = "identity") +
  theme_bw()

# Subset to the 2021 period with psychrometer instrumentation
flux_2021 <- flux %>%
  filter(date >= min(psy_daily$date),
         date <= max(psy_daily$date)) 

# Plot C fluxes with ppt and mean PD
ggplot() +
  geom_rect(aes(ymin = -Inf, ymax = Inf,
                xmin = min(psy_daily$date), 
                xmax = max(psy_daily$date)),
            fill = "gray") +
  geom_point(data = flux_2021, aes(x = date,
                                   y = NEE_F, 
                                   color = "NEE")) +
  geom_point(data = flux_2021, aes(x = date, 
                                   y = GPP_F, 
                                   color = "GPP")) +
  geom_point(data = flux_2021, aes(x = date, 
                                   y = RE_F, 
                                   color = "RE")) +
  geom_bar(data = flux_2021, aes(x = date, 
                                 y = P/50, 
                                 color = "ppt"), stat = "identity") +
  geom_point(data = filter(psy_daily_site, type == "PD"),
             aes(x = date,  y = WP_mean/50, col = "WP_PD")) + 
  geom_errorbar(data = filter(psy_daily_site, type == "PD"),
             aes(x = date,  
                 ymin = (WP_mean - WP_sd)/50, 
                 ymax = (WP_mean + WP_sd)/50, 
                 col = "WP_PD"),
             width = 0,
             alpha = 0.5) +
  theme_bw()

ggplot() +
  geom_rect(aes(ymin = -Inf, ymax = Inf,
                xmin = min(psy_daily$date), 
                xmax = max(psy_daily$date)),
            fill = "gray") +
  geom_point(data = flux_2021, aes(x = date,
                                   y = ET, 
                                   color = "ET")) +
  geom_bar(data = flux_2021, aes(x = date, 
                                 y = P/5, 
                                 color = "ppt"), stat = "identity") +
  geom_point(data = filter(psy_daily_site, type == "PD"),
             aes(x = date,  y = WP_mean, col = "WP_PD")) + 
  geom_errorbar(data = filter(psy_daily_site, type == "PD"),
                aes(x = date,  
                    ymin = (WP_mean - WP_sd), 
                    ymax = (WP_mean + WP_sd), 
                    col = "WP_PD"),
                width = 0,
                alpha = 0.5) +
  theme_bw()

# Join site-averaged Psy data to daily flux record
psy_wide <- psy_daily_site %>%
  tidyr::pivot_wider(1:3, names_from = type, values_from = WP_mean)

flux_2021_psy <- flux_2021 %>%
  left_join(psy_wide, by = "date")

flux_2021_psy %>%
  ggplot(aes(x = PD, y = ET)) +
  geom_point(aes(color = date))

flux_2021_psy %>%
  ggplot(aes(x = PD, y = NEE_F)) +
  geom_point(aes(color = date))

flux_2021_psy %>%
  ggplot(aes(x = PD, y = GPP_F)) +
  geom_point(aes(color = date))

# Add NA's to temporal record of psy_daily_site
dates <- data.frame(date = seq(min(psy_daily$date), max(psy_daily$date),
                                       by = "day"))

psy_daily_site <- dates %>%
  left_join(psy_daily_site) # 8 missing days

# Save out psy_daily_site
save(psy_daily_site, file = "data_cleaned/psy_daily_site.Rdata")
