# compare R2 plot for main and main + resid-pd

library(coda)
library(tidyverse)
library(ggplot2)
library(cowplot)

load("scripts/model-gpp-nirv/coda/codarep-main.Rdata")
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


# Summarize replicated output
coda_sum <- tidyMCMC(coda.rep1,
                     conf.int = TRUE,
                     conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)


# Check model fit
pred <- cbind.data.frame(flux, coda_sum)

m1 <- lm(pred.mean ~ GPP, data = pred)
sm <- summary(m1) # R2 = 0.7334, slope = 0.7316

fig_main <- pred %>%
  ggplot(aes(x = GPP, y = pred.mean)) +
  geom_abline(intercept = 0, slope = 1, col = "black",
              size = 1) +
  geom_abline(intercept = coef(sm)[1,1], 
              slope = coef(sm)[2,1], 
              col = "black",
              lty = 2) +
  geom_errorbar(aes(ymin = pred.lower, ymax = pred.upper),
                alpha = 0.25) +
  geom_point() +
  geom_text(x = 0.28, y = 0, 
            label = "italic(R^2)==0.733",
            parse = TRUE,
            hjust = 1,
            vjust = 0,
            size = 4) +
  scale_x_continuous(expression(paste("Observed GPP")),
                     limits = c(min(mod_comp$GPP,mod_comp$pred.lower,  na.rm = TRUE), 
                                max(mod_comp$GPP,mod_comp$pred.upper, na.rm = TRUE))) +
  scale_y_continuous(expression(paste("Predicted GPP")),
                     limits = c(min(mod_comp$GPP,mod_comp$pred.lower,  na.rm = TRUE), 
                                max(mod_comp$GPP,mod_comp$pred.upper, na.rm = TRUE))) +
  theme_bw(base_size = 14) +
  # coord_fixed(xlim=c(min(pred$GPP,pred$pred.lower,  na.rm = TRUE), 
  #                    max(pred$GPP,pred$pred.upper, na.rm = TRUE)),
  #             ylim=c(min(pred$GPP,pred$pred.lower,  na.rm = TRUE), 
  #                    max(pred$GPP,pred$pred.upper,  na.rm = TRUE))) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

fig_main

#### Add in resid-pd

# Residuals from main model
load("scripts/model-gpp-nirv/model-main-resid.Rdata") 
res_mean <- mean(resid_df$pred.mean)
res_sd <- sd(resid_df$pred.mean)

# coda of replicated gpp
load("scripts/model-gpp-nirv/coda/codarep-resid-pd.Rdata")
main_mat <- rbind(as.matrix(coda.rep1[[1]]),
                  as.matrix(coda.rep1[[2]]),
                  as.matrix(coda.rep1[[3]]))

# coda of replicated resid, then unscaled
load("scripts/model-gpp-nirv/coda/codarep-main.Rdata")
resid_pd_mat <- (rbind(as.matrix(coda.rep4[[1]]),
                       as.matrix(coda.rep4[[2]]),
                       as.matrix(coda.rep4[[3]])) +
                   res_mean)*res_sd


tot_mat <- main_mat + resid_pd_mat

tot_sum_df <- data.frame(pred.mean = apply(tot_mat, 2, FUN = mean),
                         pred.lower = apply(tot_mat, 2, FUN = quantile, probs = 0.025),
                         pred.upper = apply(tot_mat, 2, FUN = quantile, probs = 0.975))


mod_comp <- cbind(flux, tot_sum_df) %>%
  select(date, GPP, season, pred.mean, pred.lower, pred.upper) 

m2 <- lm(pred.mean ~ GPP, data = mod_comp)
sm2 <- summary(m2) # R2 = 0.7997, slope = 0.8886
round(sm2$adj.r.squared, 3)
# alternate version of GPP instead of resid 
fig_main_pd <- mod_comp %>%
  ggplot(aes(x = GPP, y = pred.mean)) +
  geom_abline(intercept = 0, slope = 1, col = "black",
              size = 1) +
  geom_abline(intercept = coef(sm2)[1,1], 
              slope = coef(sm2)[2,1], 
              col = "black",
              lty = 2) +
  geom_errorbar(aes(ymin = pred.lower, ymax = pred.upper),
                alpha = 0.25) +
  geom_point() +
  geom_text(x = 0.28, y = 0, 
            label = "italic(R^2)==0.800",
            parse = TRUE,
            hjust = 1,
            vjust = 0,
            size = 4) +
  scale_x_continuous(expression(paste("Observed GPP")),
                     limits = c(min(mod_comp$GPP,mod_comp$pred.lower,  na.rm = TRUE), 
                                max(mod_comp$GPP,mod_comp$pred.upper, na.rm = TRUE))) +
  scale_y_continuous(expression(paste("Predicted GPP")),
                     limits = c(min(mod_comp$GPP,mod_comp$pred.lower,  na.rm = TRUE), 
                                max(mod_comp$GPP,mod_comp$pred.upper, na.rm = TRUE))) +
  theme_bw(base_size = 14) +
  # coord_fixed(xlim=c(min(pred$GPP,pred$pred.lower,  na.rm = TRUE), 
  #                    max(pred$GPP,pred$pred.upper, na.rm = TRUE)),
  #             ylim=c(min(pred$GPP,pred$pred.lower,  na.rm = TRUE), 
  #                    max(pred$GPP,pred$pred.upper,  na.rm = TRUE))) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

fig_main_pd

fig_r2 <- plot_grid(fig_main, fig_main_pd,
                      ncol = 2,
                      labels = c("a", "b"),
                      rel_widths = c(1, 1),
                      label_size = 12)

ggsave(filename = "scripts/model-gpp-nirv/figs/fig_r2.png",
       plot = fig_r2,
       width = 6, height = 3,
       units = "in")

