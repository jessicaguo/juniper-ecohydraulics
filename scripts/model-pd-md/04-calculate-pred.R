# Calculate
library(coda)
library(dplyr)
library(ggplot2)
library(harrypotter)
library(ggh4x)

# Load env and Psy time series
load("data_cleaned/met_daily.Rdata")
load("scripts/model-pd-md/met_in.Rdata")
load("scripts/model-pd-md/psy_in.Rdata")

# Load coda
load("scripts/model-pd-md/coda/coda.Rdata")
Bstars <- rbind(coda.out[[1]], coda.out[[2]], coda.out[[3]])[, 13:16]
str(Bstars)
colnames(Bstars)

# Mean and sd from whole dataset, used to reverse the center/scaling
Dmax_m <- mean(met_daily$VPD_max)
Dmax_sd <- sd(met_daily$VPD_max)
W10_m <- mean(met_daily$VWC_10cm)
W10_sd <- sd(met_daily$VWC_10cm)

# Functions to center/uncenter
convert_dmax <- function(x, center) {
  if(center == FALSE) {
    return(x*Dmax_sd + Dmax_m)
  } else {
    return((x - Dmax_m)/Dmax_sd)
  }
}

convert_w10 <- function(x, center) {
  if(center == FALSE) {
    return(x*W10_sd + W10_m)
  } else {
    return((x - W10_m)/W10_sd)
  }
}

met_lim <- met_in %>%
  filter(date >= min(psy_in$date), date <= max(psy_in$date))

range(met_lim$VPD_max)
range(met_lim$Dmax)
range(met_lim$VWC_10cm)
range(met_lim$VWC10)

in_dmax <- convert_dmax(seq(0, 7, by = 0.05), center = TRUE)
in_w10 <- convert_w10(seq(4, 17, by = 0.1), center = TRUE)

df <- expand.grid(in_dmax, in_w10)
str(df)

# Function for an apply
pred_sigma <- function(input, vec) {
  return(vec[1] + vec[2] * input[1] + vec[3] * input[2] + 
           vec[4] * input[1] * input[2])
}

pred_sigma(df[1,], Bstars[1,])


out <- matrix(NA, nrow = nrow(df),  ncol = nrow(Bstars))
for(b in 1:nrow(Bstars)){
  out[,b] <- apply(df, MARGIN = 1, FUN = pred_sigma, vec = Bstars[b,])
}
# save(out, file = "scripts/model-pd-md/products/sigma_mat.Rdata")
load(file = "scripts/model-pd-md/products/sigma_mat.Rdata")

df$pred.mean <- apply(out, 1, FUN = mean)
df$pred.lower <- apply(out, 1, FUN = quantile, probs = 0.025)
df$pred.upper <- apply(out, 1, FUN = quantile, probs = 0.975)
colnames(df)[1:2] <- c("Dmax", "W10")
head(df)

df$Dmax_orig <- convert_dmax(df$Dmax, center = FALSE)
df$W10_orig <- convert_w10(df$W10, center = FALSE)
head(df)

hydry_pred <- df %>%
  mutate(strategy = case_when(pred.lower > 1 ~ "extreme anisohydry",
                              pred.upper < 1 ~ "isohydry",
                              pred.lower < 1 & pred.upper > 1 ~ "anisohydry"),
         strategy = factor(strategy, levels = c("isohydry", 
                                                "anisohydry", 
                                                "extreme anisohydry")))

# save(hydry_pred, file = "scripts/model-pd-md/products/hydry_pred.Rdata")
load(file = "scripts/model-pd-md/products/hydry_pred.Rdata")

fig6 <- ggplot(hydry_pred, aes(x = Dmax_orig, y = W10_orig)) +
  geom_tile(aes(fill = strategy),
            alpha = 0.8) +
  scale_x_continuous(expression(paste(D[max]^ant, " (kPa)")), 
                     guide = "axis_minor") +
  scale_y_continuous(expression(paste(W[10]^ant, " (%)")), 
                     guide = "axis_minor") +
  scale_fill_hp_d(option = "HermioneGranger", direction = -1) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'),
        ggh4x.axis.ticks.length.minor = rel(1),
        legend.title = element_blank()) 

# ggsave(filename = "scripts/model-pd-md/figs/fig_6.jpg",
#        plot = fig6, 
#        width = 6, height = 3,
#        units = "in")

min(df$W10)
min(df$W10_orig)
