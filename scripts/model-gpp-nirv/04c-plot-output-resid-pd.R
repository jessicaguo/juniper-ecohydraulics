# Make plots of parameters and observed vs. predicted

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
              select(-n)) 

# Combine psy and flux to plot parameters
# Add residuals of simple model
derived <- cbind.data.frame(flux, resid = scale(resid_df$pred.mean)) %>%
  mutate(lab = "Psi[PD]")


# Plot raw
parse.labels <- function(x) parse(text = x)

fig8a <- derived %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = WP_mean+2, color = lab)) +
  geom_line(aes(y = resid),
            size = 0.75) +
  scale_y_continuous(expression(paste("Scaled resid")), 
                     sec.axis = sec_axis(~.-2,
                                         expression(paste(Psi[PD], " (MPa)")))) +
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

fig8a
# Plot parameters
param_sum <- tidyMCMC(coda.out4,
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
                               Covariate == 2 ~ "Psi[PD]"),
         Covariate = factor(Covariate, levels = c("Intercept",
                                                  "Psi[PD]")),
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

fig8b <- ggplot() +
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
  theme_bw(base_size = 10) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 16, face = "bold"),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

fig8b

# Summarize replicated output
coda_sum <- tidyMCMC(coda.rep4,
                     conf.int = TRUE,
                     conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)


# Check model fit
pred <- cbind.data.frame(derived, coda_sum)

m1 <- lm(pred.mean ~ resid, data = pred)
sm <- summary(m1) # R2 = 0.1999, slope = 0.2049

fig8c <- pred %>%
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
  geom_text(x = 2.5, y = -2.35, 
            label = "italic(R^2)==0.199",
            parse = TRUE,
            hjust = 1,
            vjust = 0,
            size = 4) +
  scale_x_continuous(expression(paste("Observed resid")),
                     limits = c(min(pred$GPP,pred$pred.lower,  na.rm = TRUE), 
                                max(pred$GPP,pred$pred.upper, na.rm = TRUE))) +
  scale_y_continuous(expression(paste("Predicted resid")),
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

fig8c


fig8_bot <- plot_grid(fig8b, fig8c,
                     ncol = 2,
                     labels = c("b", "c"),
                     rel_widths = c(1, 1.25),
                     label_size = 12)

fig8 <- plot_grid(fig8a, fig8_bot,
                  ncol = 1,
                  labels = c("a", ""),
                  # rel_widths = c(1.5, 1),
                  label_size = 12)

fig8

ggsave(filename = "scripts/model-gpp-nirv/figs/fig_8.png",
       plot = fig8,
       width = 8, height = 6,
       units = "in")

