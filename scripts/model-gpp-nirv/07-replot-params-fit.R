# Test a scatterplot of GPP resids vs. PD - is there a step function that's appropriate?
# Take labels of each pulse and code the MD vs. PD figure by pulse?

#### For revision
# Replot main fig to be all 3 params for panel a
# Regression params for panel b
# and fit (combined with original data) for panel 

# Make plots of parameters and observed vs. predicted

library(tidyverse)
library(coda)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggh4x)
library(harrypotter)

# Load coda, coda.rep, input data
load("scripts/model-gpp-nirv/coda/coda-resid-pd.Rdata")
load("scripts/model-gpp-nirv/coda/codarep-resid-pd.Rdata")
load("data_cleaned/psy_daily_site_gapfilled.Rdata")

load("scripts/model-gpp-nirv/model-main-resid.Rdata") # resid_df

# Load met data
load("data_cleaned/met_daily.Rdata")
met_in <- met_daily %>%
  filter(date >= as.Date("2021-01-01")) %>%
  mutate(Dmax = scale(VPD_max),
         VWC5 = scale(VWC_5cm),
         VWC10 = scale(VWC_10cm),
         VWC50 = scale(VWC_50cm),
         PAR = scale(PAR_In)) 

# Add seasons
ppt <- met_in %>%
  filter(date >= as.Date("2021-07-01"),
         date <= as.Date("2021-09-30")) %>%
  select(date, Precip) %>%
  mutate(cPrecip = cumsum(Precip))

tot <- sum(ppt$Precip)

ppt <- ppt %>%
  mutate(season = case_when(cPrecip < 0.1 * tot ~ "premonsoon",
                            cPrecip > 0.1 * tot ~ "monsoon"))
monsoon_st <- ppt$date[min(which(ppt$season == "monsoon"))]
monsoon_en <- ppt$date[max(which(ppt$season == "monsoon"))]

# Read in NIRv & scale
sat_dat <- read_csv("data_cleaned/US-Cdm_NIRv.csv") %>%
  rename(date = Date) %>%
  select(-1) %>%
  mutate(NDVI = scale(NDVI),
         NIRv = scale(NIRv),
         smoothed_NIRv = scale(smoothed_NIRv))


flux <- readr::read_csv("data_raw/US-CdM daily.csv") %>%
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j")),
         GPP = case_when(GPP_F > 0 ~ GPP_F)) %>% # should we only use positive values?
  relocate(date) %>%
  filter(date >= min(psy_daily_site_gapfilled$date), 
         date <= max(psy_daily_site_gapfilled$date)) %>%
  mutate(season = case_when(date < monsoon_st ~ "premonsoon",
                            date >= monsoon_st & date <= monsoon_en ~ "monsoon",
                            date > monsoon_en ~ "fall")) %>%
  left_join(sat_dat) %>%
  left_join(psy_daily_site_gapfilled %>%
              filter(type == "PD") %>%
              select(-n)) |> 
  mutate(WP_kalman_scaled = scale(WP_kalman)) |> 
  left_join(met_in)


# Combine psy and flux to plot parameters
# Add residuals of simple model, convert to long format
derived <- cbind.data.frame(flux, resid = as.vector(scale(resid_df$pred.mean))) |> 
  # mutate(lab = "bar(Psi[PD])") |> 
  select(date, WP_kalman_scaled, VWC10, Dmax) |> 
  pivot_longer(-date, names_to = "variable", 
               values_to = "value") |> 
  mutate(lab = case_when(variable == "WP_kalman_scaled" ~ "bar(Psi[PD])",
                         variable == "Dmax" ~ "D",
                         variable == "VWC10" ~ "W[10]^ant"))

resids <- cbind.data.frame(flux, resid = as.vector(scale(resid_df$pred.mean))) |> 
  select(date, resid) |> 
  pivot_longer(-date, names_to = "variable", 
               values_to = "value") 

# Plot raw
parse.labels <- function(x) parse(text = x)

ggplot() +
  geom_line(data = resids, 
            aes(x = date, y = value)) +
  geom_point(data = derived,
             aes(x = date, y = value, color = lab)) +
  scale_y_continuous(expression(paste("Scaled variables"))) +
  scale_color_hp_d(option = "Always",
                   labels = parse.labels,
                   begin = 0.5) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.1, 0.1),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, "cm")) +
  labs(color = "covariate")



# For supplemental fig
# Replot to be models 1 and 2 parameters for panel a
# Lag effects for panel b
# and fit (combined with original data) for panel c