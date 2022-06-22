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
# labs <- c(
#           expression(D^ant),
#           expression(W[10]^ant),
#           expression(D^ant %*% W[10]^ant))

fig3 <- ggplot() +
  geom_hline(data = dummy, 
             aes(yintercept = intercept),
             size = 1, 
             color = "gray") +
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
             col = "coral",
             stroke = 1.25) +
  geom_point(data = neg,
             aes(x = Covariate, y = min(pred.lower - 0.05)),
             pch = 8,
             col = "medium purple", 
             stroke = 1.25) +
  scale_y_continuous("Posterior mean") +
  scale_x_discrete(labels = scales::parse_format()) +
  facet_grid2(rows = vars(Parameter),
              cols = vars(Panel),
              scales = "free",
              independent = "y",
              space = "free_x",
              labeller = label_parsed) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 16, face = "bold"),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

ggsave(filename = "scripts/model-pd-md/figs/fig_3.png",
       plot = fig3,
       width = 8, height = 3,
       units = "in")

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

fig4 <- ggplot() +
  geom_hline(data = dummy2, 
             aes(yintercept = intercept),
             size = 1, 
             color = "gray") +
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

ggsave(filename = "scripts/model-pd-md/figs/fig_4.png",
       plot = fig4,
       width = 8, height = 3,
       units = "in")

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

fig5 <- pred %>%
  ggplot(aes(x = MD, y = pred.mean)) +
  geom_abline(intercept = 0, slope = 1, col = "black",
              size = 1) +
  geom_abline(intercept = coef(sm)[1,1], 
              slope = coef(sm)[2,1], 
              col = "black",
              lty = 2) +
  geom_errorbar(aes(ymin = pred.lower, ymax = pred.upper,
                    col = Tree),
                alpha = 0.15) +
  geom_point(aes(col = Tree)) +
  scale_x_continuous(expression(paste("Observed ", Psi[MD])), breaks = seq(-9, 0, 3)) +
  scale_y_continuous(expression(paste("Predicted ", Psi[MD])), breaks = seq(-9, 0, 3)) +
  scale_colour_hp_d(option = "LunaLovegood", name = "Tree") +
  theme_bw(base_size = 14) +
  coord_fixed(xlim=c(min(pred$MD,pred$pred.lower,  na.rm = TRUE), 
                     max(pred$MD,pred$pred.upper, na.rm = TRUE)),
              ylim=c(min(pred$MD,pred$pred.lower,  na.rm = TRUE), 
                     max(pred$MD,pred$pred.upper,  na.rm = TRUE))) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))


ggsave(filename = "scripts/model-pd-md/figs/fig_5.png",
       plot = fig5,
       width = 4, height = 4,
       units = "in")
