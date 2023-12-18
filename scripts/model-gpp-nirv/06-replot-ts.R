library(tidyverse)
library(ggplot2)
library(cowplot)

load("data_cleaned/psy_daily_site_gapfilled.Rdata")

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
  left_join(met_in)

main <- flux %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = smoothed_NIRv/10, color = "NIRv")) +
  geom_point(aes(y = PAR/10, color = "PAR")) +
  geom_line(aes(y = GPP, color = "GPP"),
            color = "black",
            size = 0.75) +
  scale_y_continuous(expression(paste("GPP (mol ", CO[2], " ", m^-2, d^-1, ")")), 
                     sec.axis = sec_axis(~.*10,
                                         expression(paste("scaled PAR | NIRv")))) +
  scale_color_manual(values = c("#76653B", "#517C8D"),
                   limits = c("PAR", "NIRv")) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.1, 0.1),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, "cm")) +
  labs(color = "covariate")

ggsave("scripts/model-gpp-nirv/figs/ppt_1.png",
       plot = main,
       width = 8,
       height = 3)


### second fig
load("scripts/model-gpp-nirv/model-main-resid.Rdata") # resid_df

flux2 <- readr::read_csv("data_raw/US-CdM daily.csv") %>%
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
              select(-n)) 

# Combine psy and flux to plot parameters
# Add residuals of simple model
derived <- cbind.data.frame(flux2, resid = scale(resid_df$pred.mean)) %>%
  mutate(lab = "bar(Psi[PD])")


# Plot raw
parse.labels <- function(x) parse(text = x)

resid_pd <- derived %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = WP_mean+2, color = lab)) +
  geom_line(aes(y = resid),
            size = 0.75) +
  scale_y_continuous(expression(paste("Scaled resid")), 
                     sec.axis = sec_axis(~.-2,
                                         expression(paste(bar(Psi[PD]), " (MPa)")))) +
  scale_color_manual(values = "#92D050",
                   labels = parse.labels) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.1, 0.1),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, "cm")) +
  labs(color = "covariate")

ggsave("scripts/model-gpp-nirv/figs/ppt_2.png",
       plot = resid_pd,
       width = 5,
       height = 2)
