# compare sigma to flux
library(dplyr)
library(ggplot2)
library(harrypotter)

# Load predicted sig timeseries
load(file = "scripts/model-pd-md/products/param_pred.Rdata")

# Read in daily fluxes & clip time range
flux <- readr::read_csv("data_raw/US-CdM daily.csv") %>%
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j"))) %>%
  relocate(date) 

pred <- param_pred %>%
  left_join(flux[, 1:7], by = "date")

fig8a <- pred %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = sigma.mean, color = "sigma")) +
  geom_point(aes(y = GPP_F*5, color = "GPP")) +
  scale_y_continuous(expression(sigma),
                     sec.axis = sec_axis(~./5, name = "GPP")) +
  scale_color_manual(values = c("forestgreen", "sienna3")) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y.left = element_text(color = "sienna3"),
        axis.title.y.right = element_text(color = "forestgreen"),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'),
        ggh4x.axis.ticks.length.minor = rel(1),
        legend.title = element_blank()) +
  guides(color = "none")

m1 <- lm(GPP_F ~ VPD_F_1_1_1 * SWC_1_2_1, data = pred)
summary(m1)

m2 <- lm(GPP_F ~ sigma.mean, data = pred)
summary(m2)

m3 <- lm(GPP_F ~ VPD_F_1_1_1, data = pred)
summary(m3)


cor(pred$VPD_F_1_1_1, pred$SWC_1_2_1, use = "complete.obs")
cor(pred$VPD_F_1_1_1, pred$sigma.mean, use = "complete.obs")
cor(pred$GPP_F, pred$sigma.mean, use = "complete.obs")
cor(pred$GPP_F, pred$VPD_F_1_1_1, use = "complete.obs")
# Correlation test, due to variation in x
ct <- cor.test(pred$sigma.mean, pred$GPP_F,
         method = "pearson")
ct$estimate

# Plot p-value and R2
fit <- data.frame(value = c(round(ct$estimate, 3), 0.001),
                  type = c("R", "P"),
                  symbol = c("==", "<"),
                  lat = c(1.35, 1.35),
                  lon = c(0.015, -0.015)) %>%
  mutate(lab = paste0("italic(", type, ") ", symbol, " ", value))

fig8b <- ggplot() +
  geom_errorbarh(data = pred, 
                 aes(xmin = sigma.lower,
                     xmax = sigma.upper,
                     y = GPP_F,
                     color = date),
                 alpha = 0.5) +
  geom_point(data = pred,
             aes(x= sigma.mean, 
                 y = GPP_F,
                 color = date)) +
  geom_abline(intercept = coef(sm)[1,1], 
              slope = coef((sm))[2,1], 
              col = "black",
              lty = 2) +
  geom_text(data = fit,  
            aes(x = lat, y = lon,
                label = lab), 
            parse = TRUE,
            hjust = 1) +
  # scale_color_gradient(low = "coral", high = "mediumpurple",
  #                      breaks = seq(as.Date("2021-06-01"), 
  #                                   as.Date("2021-11-01"),
  #                                   by = "month")) +
  scale_x_continuous(expression(sigma)) +
  scale_y_continuous("GPP") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'),
        ggh4x.axis.ticks.length.minor = rel(1),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, "cm"))

fig8 <- plot_grid(fig8a, fig8b, rel_widths = c(2, 1.5))

ggsave(filename = "scripts/model-pd-md/figs/fig_8.png",
       plot = fig8,
       width = 8, height = 2,
       units = "in")
