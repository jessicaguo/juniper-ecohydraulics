# Calculate
library(coda)
library(dplyr)
library(ggplot2)
library(harrypotter)
library(ggh4x)
library(cowplot)

# Load env and Psy time series
# load("data_cleaned/met_daily.Rdata")
load("scripts/model-pd-md/met_in.Rdata")
load("scripts/model-pd-md/psy_in.Rdata")

# Load coda
load("scripts/model-pd-md/coda/coda.Rdata")
Bstars <- rbind(coda.out[[1]], coda.out[[2]], coda.out[[3]])[, 13:16]
str(Bstars)
colnames(Bstars)

Astars <- rbind(coda.out[[1]], coda.out[[2]], coda.out[[3]])[, 5:8]
str(Astars)
colnames(Astars)

wA <- rbind(coda.out[[1]], coda.out[[2]], coda.out[[3]])[, 51:55]
colnames(wA)
wB <- rbind(coda.out[[1]], coda.out[[2]], coda.out[[3]])[, 56:62]
colnames(wB)

# Create antecedent variables
doy <- unique(psy_in$doy)
nlagA = 5
nlagB = 7
pA = 1
pB = 3

Dant <- matrix(NA, nrow = length(doy), ncol = 1)
W10ant <- matrix(NA, nrow = length(doy), ncol = 1)
DTemp <- matrix(NA, nrow = length(doy), ncol = nlagA)
w10Temp <- matrix(NA, nrow = length(doy), ncol = nlagB)

# Functions for an apply
pred_Dant <- function(weights, met) {
  Dant <- numeric(length = length(doy))
  DTemp <- matrix(NA, nrow = length(doy), ncol = nlagA)

  for(i in 1:length(doy)) {
    for(k in 1:nlagA){
      DTemp[i,k] <- mean(met$Dmax[(doy[i]-k*pA+1):(doy[i]-k*pA+pA)])*weights[k]
    }
    Dant[i] <- sum(DTemp[i,])
  }
  return(Dant)
}

# pred_Dant(weights = wA[3,], met = met_in)

pred_W10ant <- function(weights, met) {
  W10ant <- numeric(length = length(doy))
  w10Temp <- matrix(NA, nrow = length(doy), ncol = nlagB)
  for(i in 1:length(doy)) {
    for(k in 1:nlagB){
    w10Temp[i,k] <- mean(met$VWC10[(doy[i]-k*pB+1):(doy[i]-k*pB+pB)])*weights[k]
    }
    W10ant[i] <- sum(w10Temp[i,])
  }
  return(W10ant)
}

# pred_W10ant(weights = wB[1,], met = met_in)

# Trim input
met_input <- met_in%>%
  select(date, Dmax, VWC10)

# Apply to all 9000 iterations pf wA and wB, matrix output
Dant_mat <- apply(wA, MARGIN = 1, FUN = pred_Dant, met = met_input)
W10ant_mat <- apply(wB, MARGIN = 1, FUN = pred_W10ant, met = met_input)

# Create posterior summary
met_pred <- data.frame(doy = doy,
                       date = unique(psy_in$date))
met_pred$Dant.mean <- apply(Dant_mat, 1, FUN = mean)
met_pred$Dant.lower <- apply(Dant_mat, 1, FUN = quantile, probs = 0.025)
met_pred$Dant.upper <- apply(Dant_mat, 1, FUN = quantile, probs = 0.975)

met_pred$W10ant.mean <- apply(W10ant_mat, 1, FUN = mean)
met_pred$W10ant.lower <- apply(W10ant_mat, 1, FUN = quantile, probs = 0.025)
met_pred$W10ant.upper <- apply(W10ant_mat, 1, FUN = quantile, probs = 0.975)

# Save/load as needed
# save(met_pred, file = "scripts/model-pd-md/products/met_pred.Rdata")
load(file = "scripts/model-pd-md/products/met_pred.Rdata")

# Use Dant.mean and W10ant.mean to predict timeseries of sigma or lambda
pred_param <- function(betas, covariates) {
  param <- numeric(length = nrow(covariates))
  for(i in 1:nrow(covariates)) {
    param[i] <- betas[1] + betas[2] * covariates$Dant.mean[i] + 
             betas[3] * covariates$W10ant.mean[i] + 
             betas[4] * covariates$Dant.mean[i] * covariates$W10ant.mean[i]
  }
  return(param)
}


# pred_param(Bstars[1, ], met_pred)

# Apply to all 9000 iterations of Bstars, matrix output
sigma_mat <- apply(Bstars, MARGIN = 1, FUN = pred_param, covariates = met_pred)
dim(sigma_mat)

lambda_mat <- apply(Astars, MARGIN = 1, FUN = pred_param, covariates = met_pred)
dim(lambda_mat)

param_pred <- met_pred 

param_pred$sigma.mean <- apply(sigma_mat, 1, FUN = mean)
param_pred$sigma.lower <- apply(sigma_mat, 1, FUN = quantile, probs = 0.025)
param_pred$sigma.upper <- apply(sigma_mat, 1, FUN = quantile, probs = 0.975)

param_pred$lambda.mean <- apply(lambda_mat, 1, FUN = mean)
param_pred$lambda.lower <- apply(lambda_mat, 1, FUN = quantile, probs = 0.025)
param_pred$lambda.upper <- apply(lambda_mat, 1, FUN = quantile, probs = 0.975)


param_pred <- param_pred %>%
  mutate(strategy = case_when(sigma.lower > 1 ~ "extreme~anisohydry~(sigma > 1)",
                              sigma.upper < 1 ~ "isohydry~(sigma < 1)",
                              sigma.lower < 1 & sigma.upper > 1 ~ "anisohydry~(sigma %~~% 1)"),
         strategy = factor(strategy, levels = c("isohydry~(sigma < 1)", 
                                                "anisohydry~(sigma %~~% 1)", 
                                                "extreme~anisohydry~(sigma > 1)")))

# save(param_pred, file = "scripts/model-pd-md/products/param_pred.Rdata")
load(file = "scripts/model-pd-md/products/param_pred.Rdata")

fig7a <- ggplot(param_pred, aes(x = date)) +
  geom_hline(yintercept = 1, 
             color = "gray",
             size = 1) +
  geom_errorbar(aes(ymin = sigma.lower, ymax = sigma.upper,
                    color = strategy),
                width = 0,
                alpha = 0.5) +
  geom_point(aes(y = sigma.mean, color = strategy)) +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 months",
               guide = "axis_minor") +
  scale_y_continuous(expression(sigma)) +
  scale_color_hp_d(option = "HermioneGranger", direction = -1,
                   labels = scales::label_parse()) +
                   # labels = scales::label_parse("isohydry~(sigma < 1)",
                   #            "anisohydry~(sigma %~~% 1)",
                   #            "extreme~anisohydry~(sigma > 1)")) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(size = 16, face = 'bold'),
        axis.title.x = element_blank(),
        ggh4x.axis.ticks.length.minor = rel(1),
        legend.title = element_blank(),
        legend.position = c(0.15, 0.88),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, "cm")) +
  guides(color = guide_legend(override.aes = list(linetype = c(0, 0, 0))))

fig7b <- ggplot(param_pred, aes(x = date)) +
  # geom_hline(yintercept = 1, 
  #            color = "gray",
  #            size = 1) +
  geom_errorbar(aes(ymin = lambda.lower, ymax = lambda.upper),
                width = 0,
                alpha = 0.25) +
  geom_point(aes(y = lambda.mean)) +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 months",
               guide = "axis_minor") +
  scale_y_continuous(expression(lambda)) +
  # scale_color_hp_d(option = "HermioneGranger", direction = -1) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'),
        axis.title.x = element_blank(),
        ggh4x.axis.ticks.length.minor = rel(1)) +
  guides(color = guide_legend(override.aes = list(linetype = c(0, 0, 0))))

fig7 <- plot_grid(fig7a, fig7b, ncol = 1,
                  align = "v")

# ggsave(filename = "scripts/model-pd-md/figs/fig_7.png",
#        plot = fig7,
#        width = 8, height = 6,
#        units = "in")

# ggsave(filename = "scripts/model-pd-md/figs/fig_7a.png",
#        plot = fig7a,
#        width = 8, height = 3,
#        units = "in")

