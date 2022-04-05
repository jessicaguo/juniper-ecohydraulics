### Create plots of parameters and model fit

library(coda)
library(broom.mixed)
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(ggh4x)
library(data.table)

# Load input data
load("scripts/daily-model/psy_in.Rdata")

# Load posterior chains
load("scripts/daily-model/coda/coda.Rdata")
load("scripts/daily-model/coda/codarep.Rdata")

# Summarize parameters
coda.param <- tidyMCMC(coda.out,
                       conf.int = TRUE,
                       conf.method = "HPDinterval",
                       drop.pars = NA)

# Plotting tree and population level means
tree.param <- coda.param %>%
  filter(grepl("mu\\.a\\[", term)) %>%
  mutate(Param = as.factor(str_match(term, "mu\\.a\\[(.)")[,2]),
         Tree = as.factor(str_match(term, "mu\\.a\\[[0-9]{1}\\,(.)")[,2]),
         labs = case_when(Param == "1" ~ "Intercept",
                          Param == "2" ~ "D^ant",
                          Param == "3" ~ "W[5]^ant",
                          Param == "4" ~ "W[50]^ant",
                          Param == "5" ~ "D^ant %.% W[5]^ant",
                          Param == "6" ~ "D^ant %.% W[50]^ant",
                          Param == "7" ~ "W[5]^ant %.% W[50]^ant"),
         labs = factor(labs, levels = c("Intercept",
                                        "D^ant", "W[5]^ant", "W[50]^ant",
                                        "D^ant %.% W[5]^ant",
                                        "D^ant %.% W[50]^ant",
                                        "W[5]^ant %.% W[50]^ant")))
pop.param <- coda.param %>%
  filter(grepl("mu\\.alpha", term)) %>%
  mutate(Param = as.factor(str_match(term, "mu\\.alpha\\[(.)")[,2]),
         labs = case_when(Param == "1" ~ "Intercept",
                          Param == "2" ~ "D^ant",
                          Param == "3" ~ "W[5]^ant",
                          Param == "4" ~ "W[50]^ant",
                          Param == "5" ~ "D^ant %.% W[5]^ant",
                          Param == "6" ~ "D^ant %.% W[50]^ant",
                          Param == "7" ~ "W[5]^ant %.% W[50]^ant"),
         labs = factor(labs, levels = c("Intercept",
                                        "D^ant", "W[5]^ant", "W[50]^ant",
                                        "D^ant %.% W[5]^ant",
                                        "D^ant %.% W[50]^ant",
                                        "W[5]^ant %.% W[50]^ant")))

ggplot() +
  geom_hline(yintercept = 0,
             lty = 2) +
  geom_errorbar(data = tree.param,
                aes(x = fct_rev(labs),
                    ymin = conf.low,
                    ymax = conf.high,
                    col = Tree),
                width = 0,
                alpha = 0.5,
                position = position_dodge(width = .5)) +
  geom_point(data = tree.param,
             aes(x = fct_rev(labs), 
                 y = estimate,
                 col = Tree),
             position = position_dodge(width = .5)) +
  geom_errorbar(data = pop.param,
                aes(x = fct_rev(labs),
                    ymin = conf.low,
                    ymax = conf.high),
                width = 0) +
  geom_point(data = pop.param,
             aes(x = fct_rev(labs),
                 y = estimate),
             size = 3) +
  scale_x_discrete("Mean", 
                   labels = parse(text = rev(levels(pop.param$labs)))) +
  scale_y_continuous("Posterior mean") +
  coord_flip() +
  theme_bw(base_size = 14) 

# Plotting tree and population level sigs
tree.sigs <- coda.param %>%
  filter(grepl("sig\\.a\\[", term)) %>%
  mutate(Param = as.factor(str_match(term, "sig\\.a\\[(.)")[,2]),
         Tree = as.factor(str_match(term, "sig\\.a\\[[0-9]{1}\\,(.)")[,2]),
         labs = case_when(Param == "1" ~ "Intercept",
                          Param == "2" ~ "D^ant",
                          Param == "3" ~ "W[5]^ant",
                          Param == "4" ~ "W[50]^ant",
                          Param == "5" ~ "D^ant %.% W[5]^ant",
                          Param == "6" ~ "D^ant %.% W[50]^ant",
                          Param == "7" ~ "W[5]^ant %.% W[50]^ant"),
         labs = factor(labs, levels = c("Intercept",
                                        "D^ant", "W[5]^ant", "W[50]^ant",
                                        "D^ant %.% W[5]^ant",
                                        "D^ant %.% W[50]^ant",
                                        "W[5]^ant %.% W[50]^ant")))
pop.sigs <- coda.param %>%
  filter(grepl("sig\\.alpha", term)) %>%
  mutate(Param = as.factor(str_match(term, "sig\\.alpha\\[(.)")[,2]),
         labs = case_when(Param == "1" ~ "Intercept",
                          Param == "2" ~ "D^ant",
                          Param == "3" ~ "W[5]^ant",
                          Param == "4" ~ "W[50]^ant",
                          Param == "5" ~ "D^ant %.% W[5]^ant",
                          Param == "6" ~ "D^ant %.% W[50]^ant",
                          Param == "7" ~ "W[5]^ant %.% W[50]^ant"),
         labs = factor(labs, levels = c("Intercept",
                                        "D^ant", "W[5]^ant", "W[50]^ant",
                                        "D^ant %.% W[5]^ant",
                                        "D^ant %.% W[50]^ant",
                                        "W[5]^ant %.% W[50]^ant")))

ggplot() +
  geom_hline(yintercept = 0,
             lty = 2) +
  geom_errorbar(data = tree.sigs,
                aes(x = fct_rev(labs),
                    ymin = conf.low,
                    ymax = conf.high,
                    col = Tree),
                width = 0,
                alpha = 0.5,
                position = position_dodge(width = .5)) +
  geom_point(data = tree.sigs,
             aes(x = fct_rev(labs), 
                 y = estimate,
                 col = Tree),
             position = position_dodge(width = .5)) +
  geom_errorbar(data = pop.sigs,
                aes(x = fct_rev(labs),
                    ymin = conf.low,
                    ymax = conf.high),
                width = 0) +
  geom_point(data = pop.sigs,
             aes(x = fct_rev(labs),
                 y = estimate),
             size = 3) +
  scale_x_discrete("Standard deviation",
                   labels = parse(text = rev(levels(pop.param$labs)))) +
  scale_y_continuous("Posterior mean") +
  coord_flip() +
  theme_bw(base_size = 14)

# Weights (population level)

weight.param <- coda.param %>%
  filter(grepl("^w[ABC]{1}\\[", term)) %>%
  mutate(Param = str_match(term, "^w(.)\\[")[,2],
         step.num = str_match(term, "^w[ABC]\\[(.)")[,2],
         labs = case_when(Param == "A" ~ "D^ant",
                          Param == "B" ~ "W[10]^ant",
                          Param == "C" ~ "W[50]^ant"),
         step.size = case_when(Param %in% c("A", "B") ~ step.num,
                               Param == "C" & step.num == "1" ~ "8-12",
                               Param == "C" & step.num == "2" ~ "13-17",
                               Param == "C" & step.num == "3" ~ "18-22",
                               Param == "C" & step.num == "4" ~ "23-27",
                               Param == "C" & step.num == "5" ~ "28-32",
                               Param == "C" & step.num == "6" ~ "33-37",
                               Param == "C" & step.num == "7" ~ "38-42",
                               Param == "C" & step.num == "8" ~ "43-47",
                               Param == "C" & step.num == "9" ~ "48-52",
                               Param == "C" & step.num == "10" ~ "53-57"))

# scales that vary by facet
scales <- list(scale_x_discrete(),
               scale_x_discrete(),
               scale_x_discrete(labels = c("8-12", "13-17", "18-22",
                                           "23-27", "28-32", "33-37",
                                           "38-42", "43-47", "48-52", 
                                           "53-57")))
ggplot(weight.param, 
       aes(x = step.num,
           y = estimate)) +
  geom_hline(yintercept = 1/7,
             lty = 2) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0) +
  geom_point() +
  labs(x = "Days past", y = "Posterior mean") +
  facet_wrap(~labs, 
             scales = "free_x",
             labeller = label_parsed,
             ncol = 1) +
  facetted_pos_scales(x = scales) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank())


# Summarize replicated output
coda.sum <- tidyMCMC(coda.rep,
                     conf.int = TRUE,
                     conf.method = "HPDinterval") %>%
  rename(pd.mean = estimate,
         pd.lower = conf.low,
         pd.upper = conf.high)

# Check model fit
pred <- cbind.data.frame(psy_in, coda.sum) %>%
  filter(!is.na(PD))
m1 <- lm(pd.mean ~ PD, data = pred)
summary(m1) # R2 = 0.9599

# Fit by tree
pred2 <- data.table(pred)
fit_bytree <- pred2[, .(m = list(summary(lm(pd.mean ~ I(PD))))), 
                      by = Tree]
getr <- function(x) {x$adj.r.squared}
getslope <- function(x) {x$coef[2,1]}
getint <- function(x) {x$coef[1,1]}
fit_tree <- data.frame(Tree = unique(pred$Tree),
                       r2 = round(unlist(lapply(fit_bytree$m, FUN = getr)), 3),
                       slope = round(unlist(lapply(fit_bytree$m, FUN = getslope)), 3),
                       int = round(unlist(lapply(fit_bytree$m, FUN = getint)), 3))


# Plot model fit
ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(data = fit_tree,
              aes(slope = slope,
                  intercept = int,
                  color = Tree),
              lty = 2) +
  geom_errorbar(data = pred,
                aes(x = PD, 
                    ymin = pd.lower, 
                    ymax = pd.upper,
                    color = Tree),
                alpha = 0.25) +
  geom_point(data = pred,
             aes(x = PD, 
                 y = pd.mean,
                 color = Tree)) +
  geom_text(data = fit_tree, 
            aes(label = r2), 
            x = Inf, y = -Inf, hjust = 1.1, vjust = -0.2,
            size = 4) +
  scale_x_continuous(expression(paste("Observed ", Psi[PD]))) +
  scale_y_continuous(expression(paste("Predicted ", Psi[PD]))) +
  theme_bw(base_size = 12) +
  facet_wrap(~Tree, ncol = 4) +
  theme(strip.background = element_blank()) +
  guides(color = "none")
