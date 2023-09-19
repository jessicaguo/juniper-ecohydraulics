# Make plots of parameters and observed vs. predicted

library(coda)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggh4x)
library(harrypotter)

# Load coda, coda.rep, input data
load("scripts/model-gpp-nirv/coda/coda-resid-env.Rdata")
load("scripts/model-gpp-nirv/coda/codarep-resid-env.Rdata")
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
  left_join(met_in)

# Combine psy and flux to plot parameters
# Add residuals of simple model
derived <- cbind.data.frame(flux, resid = scale(resid_df$pred.mean))


# Plot parameters
param_sum <- tidyMCMC(coda.out2,
                      conf.int = TRUE,
                      conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)

# Plot linear model params
B <- param_sum %>%
  filter(grepl("^B", term)) %>%
  tidyr::separate(term, 
                  into = c("Parameter", "Covariate")) %>%
  mutate(Covariate = case_when(Covariate == 1 ~ "Intercept",
                               Covariate == 2 ~ "W[10]^ant",
                               Covariate == 3 ~ "D",
                               Covariate == 4 ~ "D%.%W[10]^ant"),
         Covariate = factor(Covariate, levels = c("Intercept",
                                                  "D",
                                                  "W[10]^ant",
                                                  "D%.%W[10]^ant")),
         Panel = case_when(Covariate == "Intercept" ~ "Intercept",
                           Covariate != "Intercept" ~ "Effects"),
         Panel = factor(Panel, levels = c("Intercept", "Effects")),
         Sig = case_when(Panel == "Effects" & pred.lower*pred.upper > 0 ~ 1)) 

pos_neg <- B %>%
  filter(Sig == 1)

dummy <- data.frame(Panel = c("Effects", "Intercept"),
                    intercept = c(0, NA)) %>%
  mutate(Panel = factor(Panel, levels = c("Intercept", "Effects")))
# labs <- c(
#           expression(D^ant),
#           expression(W[10]^ant),
#           expression(D^ant %*% W[10]^ant))

figS2b <- ggplot() +
  geom_hline(data = dummy, 
             aes(yintercept = intercept),
             linewidth = 0.8, 
             color = "gray70") +
  geom_errorbar(data = filter(B, Covariate != "Intercept"), 
                aes(x = Covariate,
                    ymin = pred.lower,
                    ymax = pred.upper),
                width = 0) +
  geom_point(data = filter(B, Covariate != "Intercept"), 
             aes(x = Covariate, 
                 y = pred.mean)) +
  geom_point(data = pos_neg,
             aes(x = Covariate, y = pred.upper + 0.025),
             pch = 8,
             col = "gray70",
             stroke = 1) +
  scale_y_continuous("Covariate effects") +
  scale_x_discrete(labels = scales::parse_format()) +
  # facet_grid2(cols = vars(Panel),
  #             scales = "free",
  #             independent = "y",
  #             space = "free_x",
  #             labeller = label_parsed) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 16, face = "bold"),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

figS2b

# Plot weights
wB <- param_sum %>%
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
         Covariate = case_when(Parameter == "wB" ~ "W[10]^ant"))


figS2c <- ggplot() +
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
  theme_bw(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text.y = element_text(size = 16, colour = "black"),
        panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
figS2c

# Summarize replicated output
coda_sum <- tidyMCMC(coda.rep2,
                     conf.int = TRUE,
                     conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)


# Check model fit
pred <- cbind.data.frame(derived, coda_sum)

m1 <- lm(pred.mean ~ resid, data = pred)
sm <- summary(m1) # R2 = 0.1163, slope = 0.1165

figS2d <- pred %>%
  ggplot(aes(x = resid, y = pred.mean)) +
  geom_abline(intercept = 0, slope = 1, col = "black",
              size = 1) +
  geom_abline(intercept = coef(sm)[1,1], 
              slope = coef(sm)[2,1], 
              col = "black",
              lty = 2) +
  geom_errorbar(aes(ymin = pred.lower, ymax = pred.upper),
                alpha = 0.25) +
  geom_point() +
  geom_text(x = 2.5, y = -2.75, 
            label = "italic(R^2)==0.116",
            parse = TRUE,
            hjust = 1,
            vjust = 0,
            size = 4) +
  scale_x_continuous(expression(paste("Observed resid (mol ", CO[2], " ", m^-2, d^-1, ")")),
                     limits = c(min(pred$GPP,pred$pred.lower,  na.rm = TRUE), 
                                max(pred$GPP,pred$pred.upper, na.rm = TRUE))) +
  scale_y_continuous(expression(paste("Predicted resid (mol ", CO[2], " ", m^-2, d^-1, ")")),
                     limits = c(min(pred$GPP,pred$pred.lower,  na.rm = TRUE), 
                                max(pred$GPP,pred$pred.upper, na.rm = TRUE))) +
  theme_bw(base_size = 10) +
  # coord_fixed(xlim=c(min(pred$GPP,pred$pred.lower,  na.rm = TRUE), 
  #                    max(pred$GPP,pred$pred.upper, na.rm = TRUE)),
  #             ylim=c(min(pred$GPP,pred$pred.lower,  na.rm = TRUE), 
  #                    max(pred$GPP,pred$pred.upper,  na.rm = TRUE))) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

figS2d

figS2_bot <- plot_grid(figS2c, figS2d,
                      ncol = 2,
                      labels = c("b", "c"),
                      rel_widths = c(1.75, 1),
                      label_size = 12)

figS2 <- plot_grid(figS2b, figS2_bot,
                  ncol = 1,
                  labels = c("a", ""),
                  # rel_widths = c(1.5, 1),
                  label_size = 12)

figS2

ggsave(filename = "scripts/model-gpp-nirv/figs/fig_S2.png",
       plot = figS2,
       width = 8, height = 6,
       units = "in")

