# Plot the reponses for each pulse
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggh4x)
library(harrypotter)

# Plot env and Psy time series
load("scripts/model-pd-md/met_in.Rdata")
load("scripts/model-pd-md/psy_in.Rdata")

# Calculate mean and SE of all trees/branches
SE <- function(x) {
  sd(x, na.rm = TRUE) / sum(!is.na(x))
}
mean_psy <- psy_in %>%
  group_by(date) %>%
  summarize(MD_mean = mean(MD, na.rm = TRUE),
            MD_sd = sd(MD, na.rm = TRUE),
            MD_se = SE(MD),
            PD_mean = mean(PD, na.rm = TRUE),
            PD_sd = sd(PD, na.rm = TRUE),
            PD_se = SE(PD),
            n = n())

# Read in daily fluxes & clip time range
flux <- readr::read_csv("data_raw/US-CdM daily.csv") %>%
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j")),
         GPP = case_when(GPP_F > 0 ~ GPP_F)) %>%
  relocate(date) %>%
  filter(date >= min(psy_in$date), 
         date <= max(psy_in$date))

# Load predicted sig timeseries
load(file = "scripts/model-pd-md/products/param_pred.Rdata")


##### Find pulses #####
met_in %>% 
  filter(date >=  min(psy_in$date),
         date <=  max(psy_in$date)) %>%
  relocate(Precip, .after = date) %>%
  arrange(desc(Precip))

# Monsoon pulse dates are 
# 07-27-2021 (38.6)
# 08-19-2021 (28.7 for 2 consecutive days)
# 09-01-2021 (17.3)

dates <- as.Date(c("2021-07-27", "2021-08-19", "2021-09-01"))

pulse_list <- list()

for (i in 1:length(dates)) {
  dur <- 10 # pulse duration
  
  d <- seq(dates[i]-2, dates[i] + dur, by = 1)
  
  pulse_list[[i]] <- data.frame(date = d,
                                pulse = rep(i, length(d)),
                                W10 = met_in |> 
                                  filter(date %in% d) |> 
                                  select(VWC_10cm) |> 
                                  pull(),
                                PD = mean_psy |> 
                                  filter(date %in% d) |> 
                                  select(PD_mean) |> 
                                  pull(),
                                PD_sd = mean_psy |> 
                                  filter(date %in% d) |> 
                                  select(PD_sd) |> 
                                  pull(),
                                MD = mean_psy |> 
                                  filter(date %in% d) |> 
                                  select(MD_mean) |> 
                                  pull(),
                                MD_sd = mean_psy |> 
                                  filter(date %in% d) |> 
                                  select(MD_sd) |> 
                                  pull(),
                                GPP = flux |> 
                                  filter(date %in% d) |> 
                                  select(GPP_F) |> 
                                  pull(),
                                sigma = param_pred |> 
                                  filter(date %in% d) |> 
                                  select(sigma.mean) |> 
                                  pull(),
                                sigma.l = param_pred |> 
                                  filter(date %in% d) |> 
                                  select(sigma.lower) |> 
                                  pull(),
                                sigma.u = param_pred |> 
                                  filter(date %in% d) |> 
                                  select(sigma.upper) |> 
                                  pull()
                                )
}
pulse_df <- do.call(rbind, pulse_list)

pulse_df %>%
  select(-PD_sd, -MD, -MD_sd, -sigma.l, -sigma.u) %>%
  tidyr::pivot_longer(-1:-2, 
                      names_to = "variables",
                      values_to = "value") %>%
  mutate(variables = factor(variables, 
                            levels = c("W10", "PD", "GPP", "sigma"))) %>%
ggplot(aes(x = date)) +
  geom_line(aes(y = value)) +
  facet_grid(rows = vars(variables),
             cols = vars(pulse), 
             scales = "free") +
  theme_bw()

