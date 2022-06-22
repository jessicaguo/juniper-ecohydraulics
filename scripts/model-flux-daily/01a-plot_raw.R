library(dplyr)
library(ggplot2)
library(cowplot)
library(ggh4x)

# Plot GPP time series
load("data_cleaned/psy_daily_site_gapfilled.Rdata")
# Read in daily fluxes & clip time range
flux <- readr::read_csv("data_raw/US-CdM daily.csv") %>%
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j")),
         GPP = case_when(GPP_F > 0 ~ GPP_F)) %>%
  relocate(date) %>%
  filter(date >= min(psy_daily_site_gapfilled$date), 
         date <= max(psy_daily_site_gapfilled$date))

fig2 <- ggplot() +
  geom_point(data = filter(flux, GPP_F > 0),
             aes(x = date, y = GPP_F)) +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 months",
               guide = "axis_minor") +
  scale_y_continuous(expression(paste("GPP (mol ", CO[2], " ", m^-2, d^-1, ")"))) +
  scale_color_hp_d(option = "Sprout") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.07, 0.92),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, "cm"))

ggsave(filename = "scripts/model-flux-daily/figs/fig_2.png",
       plot = fig2,
       width = 8, height = 3,
       units = "in")

