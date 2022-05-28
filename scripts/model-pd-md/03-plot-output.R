# Make plots of parameters and observed vs. predicted

library(coda)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(harrypotter)

# Load coda, coda.rep, input data
load("scripts/model-pd-md/coda/coda.Rdata")
load("scripts/model-pd-md/coda/codarep.Rdata")
load("scripts/model-pd-md/psy_in.Rdata")

# Limit psy_in to both observations
psy_in <- psy_in %>%
  filter(!is.na(PD), !is.na(MD))


# Plot parameters
param_sum <- tidyMCMC(coda.out,
                      conf.int = TRUE,
                      conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)

# Plot Astar and Bstar
AB <- param_sum %>%
  filter(grepl("Astar", term) | grepl("Bstar", term)) %>%
  tidyr::separate(term, 
                  into = c("Parameter", "Covariate")) %>%
  mutate(Parameter = case_when(Parameter == "Astar" ~ "lambda",
                               Parameter == "Bstar" ~ "sigma"),
         Covariate = case_when(Covariate == 1 ~ "Intercept",
                               Covariate == 2 ~ "Dant",
                               Covariate == 3 ~ "W10ant",
                               Covariate == 4 ~ "W10ant*Dant"),
         Panel = case_when(Covariate == "Intercept" ~ "Intercept",
                            Covariate != "Intercept" ~ "Effects"),
         Panel = factor(Panel, levels = c("Intercept", "Effects")),
         Parameter = factor(Parameter, levels = c("sigma", "lambda")),
         Sig = case_when(Panel == "Effects" & pred.lower*pred.upper > 0 ~ 1)) 

pos <- AB %>%
  filter(Sig == 1 & pred.mean > 0)
neg <- AB %>%
  filter(Sig == 1 & pred.mean < 0)

dummy <- data.frame(Panel = c(rep("Effects", 2), "Intercept"),
                    Parameter = c("sigma", "lambda", "sigma"),
                    intercept = c(0, 0, 1)) %>%
  mutate(Panel = factor(Panel, levels = c("Intercept", "Effects")),
         Parameter = factor(Parameter, levels = c("sigma", "lambda")))
labs <- c(
          expression(D^ant),
          expression(W[10]^ant),
          expression(D^ant * W[10]^ant))

ggplot() +
  geom_hline(data = dummy, 
             aes(yintercept = intercept),
             size = 1) +
  geom_errorbar(data = AB, 
                aes(x = Covariate,
                    ymin = pred.lower,
                    ymax = pred.upper),
                width = 0) +
  geom_point(data = AB, 
             aes(x = Covariate, 
                 y = pred.mean)) +
  geom_point(data = pos,
             aes(x = Covariate, y = max(pred.upper + 0.02)),
             pch = 8,
             col = "coral") +
  geom_point(data = neg,
             aes(x = Covariate, y = min(pred.lower - 0.02)),
             pch = 8,
             col = "medium purple") +
  scale_y_continuous("Posterior mean") +
  scale_x_discrete(labels = labs) +
  facet_grid2(rows = vars(Parameter),
              cols = vars(Panel),
              scales = "free",
              independent = "y",
              space = "free_x",
              labeller = label_parsed) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 16),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"))

# Summarize replicated output
coda_sum <- tidyMCMC(coda.rep,
                     conf.int = TRUE,
                     conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)


# Check model fit
pred <- cbind.data.frame(psy_in, coda_sum)

m1 <- lm(pred.mean ~ MD, data = pred)
sm <- summary(m1) # R2 = 0.8985; w/ RE R2 = 0.9198

pred %>%
  ggplot(aes(x = MD, y =pred.mean)) +
  geom_abline(intercept = 0, slope = 1, col = "black",
              size = 1) +
  geom_abline(intercept = coef(sm)[1,1], 
              slope = coef(sm)[2,1], 
              col = "black",
              lty = 2) +
  geom_errorbar(aes(ymin = pred.lower, ymax = pred.upper,
                    color = Tree),
                alpha = 0.25) +
  geom_point(aes(color = Tree)) +
  scale_x_continuous("Observed", breaks = seq(-9, 0, 3)) +
  scale_y_continuous("Predicted", breaks = seq(-9, 0, 3)) +
  scale_colour_hp_d(option = "LunaLovegood", name = "Tree") +
  theme_bw(base_size = 14) +
  coord_equal() +
  theme(panel.grid = element_blank())
