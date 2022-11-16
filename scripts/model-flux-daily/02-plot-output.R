# Make plots of parameters and observed vs. predicted

library(coda)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggh4x)
library(harrypotter)

# Load coda, coda.rep, input data
load("scripts/model-flux-daily/coda/coda.Rdata")
load("scripts/model-flux-daily/coda/codarep.Rdata")
load("data_cleaned/psy_daily_site_gapfilled.Rdata")
flux <- readr::read_csv("data_raw/US-CdM daily.csv") %>%
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j")),
         GPP = case_when(GPP_F > 0 ~ GPP_F)) %>%
  relocate(date) %>%
  filter(date >= min(psy_daily_site_gapfilled$date), 
         date <= max(psy_daily_site_gapfilled$date))

# Plot parameters
param_sum <- tidyMCMC(coda.out,
                      conf.int = TRUE,
                      conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)

# Plot B
B <- param_sum %>%
  filter(grepl("^B", term)) %>%
  tidyr::separate(term, 
                  into = c("Parameter", "Covariate")) %>%
  mutate(Covariate = case_when(Covariate == 1 ~ "Intercept",
                               Covariate == 2 ~ "D^ant",
                               Covariate == 3 ~ "W[10]^ant",
                               Covariate == 4 ~ "D^ant%.%W[10]^ant"),
         Covariate = factor(Covariate, levels = c("Intercept",
                                                  "D^ant",
                                                  "W[10]^ant",
                                                  "D^ant%.%W[10]^ant")),
         Panel = case_when(Covariate == "Intercept" ~ "Intercept",
                            Covariate != "Intercept" ~ "Effects"),
         Panel = factor(Panel, levels = c("Intercept", "Effects")),
         Sig = case_when(Panel == "Effects" & pred.lower*pred.upper > 0 ~ 1)) 

pos <- B %>%
  filter(Sig == 1 & pred.mean > 0)

dummy <- data.frame(Panel = c("Effects", "Intercept"),
                    intercept = c(0, NA)) %>%
  mutate(Panel = factor(Panel, levels = c("Intercept", "Effects")))
# labs <- c(
#           expression(D^ant),
#           expression(W[10]^ant),
#           expression(D^ant %*% W[10]^ant))

fig5a <- ggplot() +
  geom_hline(data = dummy, 
             aes(yintercept = intercept),
             size = 0.8, 
             color = "gray70") +
  geom_errorbar(data = filter(B, Covariate != "Intercept"), 
                aes(x = Covariate,
                    ymin = pred.lower,
                    ymax = pred.upper),
                width = 0) +
  geom_point(data = filter(B, Covariate != "Intercept"), 
             aes(x = Covariate, 
                 y = pred.mean)) +
  geom_point(data = pos,
             aes(x = Covariate, y = max(pred.upper + 0.05)),
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
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 16, face = "bold"),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

fig5a

# Plot weights
wAB <- param_sum %>%
  filter(grepl("wA", term) | grepl("wB", term)) %>%
  tidyr::separate(term, 
                  into = c("Parameter", "ts")) %>%
  mutate(Timestep = case_when(Parameter == "wA" ~ as.character(as.numeric(ts) - 1),
                              Parameter == "wB"  & ts == 1 ~ "0-2",
                              Parameter == "wB"  & ts == 2 ~ "3-5",
                              Parameter == "wB"  & ts == 3 ~ "6-8",
                              Parameter == "wB"  & ts == 4 ~ "9-11",
                              Parameter == "wB"  & ts == 5 ~ "12-14",
                              Parameter == "wB"  & ts == 6 ~ "15-17",
                              Parameter == "wB"  & ts == 7 ~ "18-20",),
         Timestep = factor(Timestep, levels = c("0", "1", "2", "3", "4", 
                                                "0-2", "3-5", "6-8", "9-11",
                                                "12-14", "15-17", "18-20")),
         Covariate = case_when(Parameter == "wA" ~ "D^ant",
                               Parameter == "wB" ~ "W[10]^ant"))

dummy2 <- data.frame(Covariate = c("D^ant", "W[10]^ant"),
                    intercept = c(1/5, 1/7)) 

fig5b <- ggplot() +
  geom_hline(data = dummy2, 
             aes(yintercept = intercept),
             size = 0.8, 
             color = "gray70") +
  geom_errorbar(data = wAB,
                aes(x = Timestep,
                    ymin = pred.lower,
                    ymax = pred.upper),
                width = 0) +
  geom_point(data = wAB,
             aes(x = Timestep, y = pred.mean)) +
  scale_y_continuous("Antecedent weights") +
  scale_x_discrete("Timesteps (days)") +
  facet_wrap(~Covariate, ncol = 2,
             scales = "free_x",
             labeller = label_parsed) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text.y = element_text(size = 16, colour = "black"),
        panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

fig5 <- plot_grid(fig5a, fig5b,
                  ncol = 1,
                  labels = "auto",
                  label_size = 12)

fig5
ggsave(filename = "scripts/model-flux-daily/figs/fig_5.png",
       plot = fig5,
       width = 8, height = 6,
       units = "in")

# Summarize replicated output
coda_sum <- tidyMCMC(coda.rep,
                     conf.int = TRUE,
                     conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)


# Check model fit
pred <- cbind.data.frame(flux, coda_sum)

m1 <- lm(pred.mean ~ GPP, data = pred)
sm <- summary(m1) # R2 = 0.741, slope = 0.762

figS2 <- pred %>%
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
  geom_text(x = 0.4, y = 0, 
            label = "italic(R^2)==0.741",
            parse = TRUE,
            hjust = 1,
            vjust = 0,
            size = 5) +
  scale_x_continuous(expression(paste("Observed GPP (mol ", CO[2], " ", m^-2, d^-1, ")"))) +
  scale_y_continuous(expression(paste("Predicted GPP (mol ", CO[2], " ", m^-2, d^-1, ")"))) +
  theme_bw(base_size = 14) +
  coord_fixed(xlim=c(min(pred$GPP,pred$pred.lower,  na.rm = TRUE), 
                     max(pred$GPP,pred$pred.upper, na.rm = TRUE)),
              ylim=c(min(pred$GPP,pred$pred.lower,  na.rm = TRUE), 
                     max(pred$GPP,pred$pred.upper,  na.rm = TRUE))) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))


ggsave(filename = "scripts/model-flux-daily/figs/fig_S2.png",
       plot = figS2,
       width = 5, height = 5,
       units = "in")
