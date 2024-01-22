# Test a scatterplot of GPP resids vs. PD - is there a step function that's appropriate?
# Take labels of each pulse and code the MD vs. PD figure by pulse?

#### For revision
# Replot main fig to be all 3 params for panel a
# All Regression params for panel b
# and all fit (combined with original data) for panel c

# Make plots of parameters and observed vs. predicted

library(tidyverse)
library(coda)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggh4x)
library(harrypotter)

# Load input data

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
                         variable == "VWC10" ~ "W[10]"))

resids <- cbind.data.frame(flux, resid = as.vector(scale(resid_df$pred.mean))) |> 
  select(date, resid) |> 
  pivot_longer(-date, names_to = "variable", 
               values_to = "value") 

# Plot raw
parse.labels <- function(x) parse(text = x)

fig7a <- ggplot() +
  geom_line(data = resids, 
            aes(x = date, y = value)) +
  geom_point(data = derived,
             aes(x = date, y = value, color = lab)) +
  scale_y_continuous(expression(paste("Scaled residuals")),
                     sec.axis = sec_axis(~.*1, "Scaled covariates")) +
  scale_color_manual(values = c("#7B9948", "#00332A", "#7D5535"),
                   labels = parse.labels) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.05, 0.15),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, "cm")) 

#### Second panel: parameters ####
# Load in code from model data
load("scripts/model-gpp-nirv/coda/coda-resid-env.Rdata")
load("scripts/model-gpp-nirv/coda/coda-resid-vwc.Rdata")
load("scripts/model-gpp-nirv/coda/coda-resid-pd.Rdata")

# Plot parameters
param_sum_1 <- tidyMCMC(coda.out2,
                      conf.int = TRUE,
                      conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)

param_sum_2 <- tidyMCMC(coda.out3,
                          conf.int = TRUE,
                          conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)

param_sum_3 <- tidyMCMC(coda.out4,
                        conf.int = TRUE,
                        conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)

# Extract Dsums

Dsums <- param_sum_1 |> 
  filter(grepl("^Dsum", term)) |> 
  bind_rows(param_sum_2 |>  filter((grepl("^Dsum", term)))) |> 
  bind_rows(param_sum_3 |>  filter((grepl("^Dsum", term)))) |> 
  mutate(model = c("Env", "Soil~only", "Psi~only") |> 
           factor(levels = c("Env", "Soil~only", "Psi~only")),
         post_means = round(pred.mean, 0),
         pD = c(6.15, 3.56, 2.98),
         lab2 = paste("italic(D*infinity)%~~%", post_means),
         lab1 = paste("italic(pD)==", pD))

# Combine sets of parameters
B <- param_sum_1 |> 
  filter(grepl("^B", term)) |> 
  bind_rows(param_sum_2 |>  filter((grepl("^B", term)))) |> 
  bind_rows(param_sum_3 |>  filter((grepl("^B", term)))) |> 
  filter(term != "B[1]") |> 
  mutate(model = c(rep("Env", 3), "Soil~only", "Psi~only") |> 
           factor(levels = c("Env", "Soil~only", "Psi~only")),
         Covariate = c("D", "W[10]^ant", "D%*%W[10]^ant", "W[10]^ant", "Psi[PD]") |> 
           factor(levels = c("D", "W[10]^ant", "D%*%W[10]^ant", "Psi[PD]")),
         Sig = case_when(pred.lower*pred.upper > 0 ~ 1))

pos_neg <- B %>%
  filter(Sig == 1)

fig7b <-
  ggplot() +
  geom_hline(yintercept = 0,
             linewidth = 0.8,
             color = "gray70") +
  geom_errorbar(data = B, 
                aes(x = Covariate,
                    ymin = pred.lower,
                    ymax = pred.upper),
                width = 0) +
  geom_point(data = B, 
             aes(x = Covariate, 
                 y = pred.mean)) +
  geom_point(data = pos_neg,
             aes(x = Covariate, y = max(pred.upper) + 0.1),
             pch = 8,
             col = "gray70",
             stroke = 1) +
  geom_text(data = Dsums,
            aes(x = 0.5, y = -0.17,
                label = lab1),
            parse = TRUE,
            hjust = 0,
            vjust = 0,
            size = 4) +
    geom_text(data = Dsums,
              aes(x = 0.5, y = -0.25,
                  label = lab2),
              parse = TRUE,
              hjust = 0,
              vjust = 0,
              size = 4) +
  scale_y_continuous("Covariate effects") +
  scale_x_discrete(labels = scales::parse_format()) +
  facet_grid(~model, scales = "free_x", space = "free_x",
             labeller = label_parsed) +
  theme_bw(base_size = 10) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        panel.grid = element_blank(),
        axis.title.x = element_blank())


# Load in coda rep from model
load("scripts/model-gpp-nirv/coda/codarep-main.Rdata")
load("scripts/model-gpp-nirv/coda/codarep-resid-env.Rdata")
load("scripts/model-gpp-nirv/coda/codarep-resid-vwc.Rdata")
load("scripts/model-gpp-nirv/coda/codarep-resid-pd.Rdata")


# Summarize replicated output
# Residuals from main model
load("scripts/model-gpp-nirv/model-main-resid.Rdata") 
res_mean <- mean(resid_df$pred.mean)
res_sd <- sd(resid_df$pred.mean)

# coda of replicated gpp
main_mat <- rbind(as.matrix(coda.rep1[[1]]),
                  as.matrix(coda.rep1[[2]]),
                  as.matrix(coda.rep1[[3]]))

# coda of replicated resid, then unscale
resid_mat_1 <- (rbind(as.matrix(coda.rep2[[1]]),
                      as.matrix(coda.rep2[[2]]),
                      as.matrix(coda.rep2[[3]])) +
                  res_mean)*res_sd
resid_mat_2 <- (rbind(as.matrix(coda.rep3[[1]]),
                      as.matrix(coda.rep3[[2]]),
                      as.matrix(coda.rep3[[3]])) +
                  res_mean)*res_sd
resid_mat_3 <- (rbind(as.matrix(coda.rep4[[1]]),
                       as.matrix(coda.rep4[[2]]),
                       as.matrix(coda.rep4[[3]])) +
                   res_mean)*res_sd


#  Total predictions from both versions of model
tot_mat_1 <- main_mat + resid_mat_1
tot_mat_2 <- main_mat + resid_mat_2
tot_mat_3 <- main_mat + resid_mat_3

tot_df_1 <- data.frame(pred.mean = apply(tot_mat_1, 2, FUN = mean),
                         pred.lower = apply(tot_mat_1, 2, FUN = quantile, probs = 0.025),
                         pred.upper = apply(tot_mat_1, 2, FUN = quantile, probs = 0.975))
tot_df_2 <- data.frame(pred.mean = apply(tot_mat_2, 2, FUN = mean),
                        pred.lower = apply(tot_mat_2, 2, FUN = quantile, probs = 0.025),
                        pred.upper = apply(tot_mat_2, 2, FUN = quantile, probs = 0.975))
tot_df_3 <- data.frame(pred.mean = apply(tot_mat_3, 2, FUN = mean),
                       pred.lower = apply(tot_mat_3, 2, FUN = quantile, probs = 0.025),
                       pred.upper = apply(tot_mat_3, 2, FUN = quantile, probs = 0.975))

# Empirical data to compare
load("data_cleaned/psy_daily_site_gapfilled.Rdata")
flux <- readr::read_csv("data_raw/US-CdM daily.csv") %>%
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j")),
         GPP = case_when(GPP_F > 0 ~ GPP_F)) %>% # should we only use positive values?
  relocate(date) %>%
  filter(date >= min(psy_daily_site_gapfilled$date), 
         date <= max(psy_daily_site_gapfilled$date))


# Combine to plot df

mod_comp_1 <- cbind.data.frame(flux, tot_df_1) |> 
  select(date, GPP, pred.mean, pred.lower, pred.upper) |> 
  mutate(model = "Env")

mod_comp_2 <- cbind.data.frame(flux, tot_df_2) |> 
  select(date, GPP, pred.mean, pred.lower, pred.upper) |> 
  mutate(model = "Soil~only")

mod_comp_3 <- cbind.data.frame(flux, tot_df_3) |> 
  select(date, GPP, pred.mean, pred.lower, pred.upper) |> 
  mutate(model = "Psi~only")

mod_comp <- bind_rows(mod_comp_1, mod_comp_2, mod_comp_3) |> 
  mutate(model = factor(model, levels = c("Env", "Soil~only", "Psi~only")))
  

# Create annotation/abline layer
sm1 <- summary(lm(pred.mean~GPP, data = mod_comp_1))
sm2 <- summary(lm(pred.mean~GPP, data = mod_comp_2))
sm3 <- summary(lm(pred.mean~GPP, data = mod_comp_3))
sum_df <- data.frame(model = c("Env", "Soil~only", "Psi~only"),
                     slope = c(coef(sm1)[2,1], coef(sm2)[2,1], coef(sm3)[2,1]),
                     int = c(coef(sm1)[1,1], coef(sm2)[1,1], coef(sm3)[1,1]),
                     R2 = c(round(sm1$adj.r.squared, 2),
                            round(sm2$adj.r.squared, 2),
                            sprintf('%.2f', sm3$adj.r.squared)) |> 
                       as.character()) |> 
  mutate(lab = paste("italic(R^2)==", R2),
         lab2 = paste("slope==", round(slope, 2)), 
         model = factor(model, levels = c("Env", "Soil~only", "Psi~only")))
                     
fig7c <-
  mod_comp %>%
  ggplot(aes(x = GPP, y = pred.mean)) +
  geom_abline(intercept = 0, slope = 1, col = "black",
              size = 1) +
  geom_abline(data = sum_df,
              aes(intercept = int,
                  slope = slope),
              col = "black",
              lty = 2) +
  geom_errorbar(aes(ymin = pred.lower, ymax = pred.upper),
                alpha = 0.25) +
  geom_point() +
  geom_text(data = sum_df,
            aes(x = 0.3, y = 0,
            label = lab),
            parse = TRUE,
            hjust = 1,
            vjust = 0,
            size = 4) +
  geom_text(data = sum_df,
            aes(x = 0.3, y = -0.04,
                label = lab2),
            parse = TRUE,
            hjust = 1,
            vjust = 0,
            size = 4) +
  facet_grid(~model,
             labeller = label_parsed) +
  scale_x_continuous(expression(paste("Observed GPP")),
                     limits = c(min(mod_comp$GPP,mod_comp$pred.lower,  na.rm = TRUE), 
                                max(mod_comp$GPP,mod_comp$pred.upper, na.rm = TRUE))) +
  scale_y_continuous(expression(paste("Predicted GPP")),
                     limits = c(min(mod_comp$GPP,mod_comp$pred.lower,  na.rm = TRUE), 
                                max(mod_comp$GPP,mod_comp$pred.upper, na.rm = TRUE))) +
  theme_bw(base_size = 10) +
  # coord_fixed() +
  # coord_fixed(xlim=c(min(pred$GPP,pred$pred.lower,  na.rm = TRUE), 
  #                    max(pred$GPP,pred$pred.upper, na.rm = TRUE)),
  #             ylim=c(min(pred$GPP,pred$pred.lower,  na.rm = TRUE), 
  #                    max(pred$GPP,pred$pred.upper,  na.rm = TRUE))) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))


plot_grid(fig7a, fig7b, fig7c, 
          align = "v",axis = "l",
          labels = "auto",
          nrow = 3)

ggsave("scripts/model-gpp-nirv/figs/fig_8_new.png",
       width = 7, height = 7,
       units = "in")



# For supplemental fig
# For completeness sake, create antecedent plots for models 1 and 2

wB_1 <- param_sum_1 %>%
  filter(grepl("wB", term)) %>%
  tidyr::separate(term, 
                  into = c("Parameter", "ts")) %>%
  mutate(Timestep = case_when(Parameter == "wB"  & ts == 1 ~ "0-2",
                              Parameter == "wB"  & ts == 2 ~ "3-5",
                              Parameter == "wB"  & ts == 3 ~ "6-8",
                              Parameter == "wB"  & ts == 4 ~ "9-11",
                              Parameter == "wB"  & ts == 5 ~ "12-14",
                              Parameter == "wB"  & ts == 6 ~ "15-17",
                              Parameter == "wB"  & ts == 7 ~ "18-20",),
         Timestep = factor(Timestep, levels = c("0-2", "3-5", "6-8", "9-11",
                                                "12-14", "15-17", "18-20")),
         Covariate = case_when(Parameter == "wB" ~ "W[10]^ant"),
         model = "Env")

wB_2 <- param_sum_2 %>%
  filter(grepl("wB", term)) %>%
  tidyr::separate(term, 
                  into = c("Parameter", "ts")) %>%
  mutate(Timestep = case_when(Parameter == "wB"  & ts == 1 ~ "0-2",
                              Parameter == "wB"  & ts == 2 ~ "3-5",
                              Parameter == "wB"  & ts == 3 ~ "6-8",
                              Parameter == "wB"  & ts == 4 ~ "9-11",
                              Parameter == "wB"  & ts == 5 ~ "12-14",
                              Parameter == "wB"  & ts == 6 ~ "15-17",
                              Parameter == "wB"  & ts == 7 ~ "18-20",),
         Timestep = factor(Timestep, levels = c("0-2", "3-5", "6-8", "9-11",
                                                "12-14", "15-17", "18-20")),
         Covariate = case_when(Parameter == "wB" ~ "W[10]^ant"),
         model = "Soil~only")

wB <- wB_1 |> 
  bind_rows(wB_2) |> 
  mutate(model = factor(model, levels = c("Env", "Soil~only")))

ggplot() +
  geom_hline(yintercept = 1/7,
             size = 0.8, 
             color = "gray70") +
  geom_errorbar(data = wB,
                aes(x = Timestep,
                    ymin = pred.lower,
                    ymax = pred.upper),
                width = 0) +
  geom_point(data = wB,
             aes(x = Timestep, y = pred.mean)) +
  scale_y_continuous("Antecedent weights") +
  scale_x_discrete("Timesteps (days)") +
  facet_grid(~model, labeller = label_parsed) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 12, colour = "black"),
        panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))


ggsave("scripts/model-gpp-nirv/figs/fig_S2_new.png",
       width = 8, height = 3,
       units = "in")
