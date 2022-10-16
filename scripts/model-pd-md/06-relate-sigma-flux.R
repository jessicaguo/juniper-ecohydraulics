# compare sigma to flux
library(dplyr)
library(ggplot2)
library(harrypotter)
library(cowplot)
library(ggh4x)

# Load predicted sig timeseries
load(file = "scripts/model-pd-md/products/param_pred.Rdata")

# Read in daily fluxes & clip time range
flux <- readr::read_csv("data_raw/US-CdM daily.csv") %>%
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j"))) %>%
  relocate(date) 

pred <- param_pred %>%
  left_join(flux[, 1:7], by = "date")

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
        legend.position = c(0.15, 0.86),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, "cm")) +
  guides(color = guide_legend(override.aes = list(linetype = c(0, 0, 0))))

fig7b <- pred %>%
  filter(GPP_F > 0) %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = GPP_F, color = "GPP")) +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 months",
               guide = "axis_minor") +
  scale_y_continuous("GPP") +
  scale_color_manual(values = c("forestgreen")) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 16, face = 'bold'), 
        ggh4x.axis.ticks.length.minor = rel(1),
        legend.title = element_blank()) +
  guides(color = "none")

fig7 <- plot_grid(fig7a, fig7b, ncol = 1,
                  align = "v")
# ggsave(filename = "scripts/model-pd-md/figs/fig_7_gpp.png",
#        plot = fig7,
#        width = 8, height = 5,
#        units = "in")

# Remove GPP < 0
pred2 <- filter(pred, GPP_F > 0)

fig8a <- pred2 %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = sigma.mean, color = "sigma")) +
  geom_point(aes(y = GPP_F*5, color = "GPP")) +
  scale_y_continuous(expression(sigma),
                     sec.axis = sec_axis(~./5, name = "GPP")) +
  scale_color_manual(values = c("forestgreen", "black")) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y.left = element_text(color = "black", face = "bold",
                                         size = 16),
        axis.title.y.right = element_text(color = "forestgreen"),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        ggh4x.axis.ticks.length.minor = rel(1),
        legend.title = element_blank()) +
  guides(color = "none")

# m1 <- lm(GPP_F ~ VPD_F_1_1_1 * SWC_1_2_1, data = pred)
# summary(m1)
# 
m2 <- lm(GPP_F ~ sigma.mean, data = pred2)
sm <- summary(m2)
# 
# m3 <- lm(GPP_F ~ VPD_F_1_1_1, data = pred)
# summary(m3)


# cor(pred$VPD_F_1_1_1, pred$SWC_1_2_1, use = "complete.obs")
# cor(pred$VPD_F_1_1_1, pred$sigma.mean, use = "complete.obs")
# cor(pred$GPP_F, pred$sigma.mean, use = "complete.obs")
# cor(pred$GPP_F, pred$VPD_F_1_1_1, use = "complete.obs")
# Correlation test, due to variation in x

ct <- cor.test(pred2$sigma.mean, pred2$GPP_F,
         method = "pearson")
ct$estimate

# Plot p-value and R2
fit <- data.frame(value = c(round(ct$estimate, 3), 0.001),
                  type = c("r", "p"),
                  symbol = c("==", "<"),
                  lat = c(1.35, 1.35),
                  lon = c(0.015, -0.015)) %>%
  mutate(lab = paste0("italic(", type, ") ", symbol, " ", value))

fig8b <- ggplot() +
  geom_errorbarh(data = pred2, 
                 aes(xmin = sigma.lower,
                     xmax = sigma.upper,
                     y = GPP_F,
                     col = strategy),
                 alpha = 0.25) +
  geom_point(data = pred2,
             aes(x= sigma.mean, 
                 y = GPP_F,
                 col = strategy)) +
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
  scale_color_hp_d(option = "HermioneGranger", direction = -1) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title.x = element_text(face = 'bold', size = 16),
        ggh4x.axis.ticks.length.minor = rel(1),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, "cm")) +
  guides(color = "none")

fig8 <- plot_grid(fig8a, fig8b, rel_widths = c(2, 1.5),
                  align = "h")

# ggsave(filename = "scripts/model-pd-md/figs/fig_8.png",
#        plot = fig8,
#        width = 8, height = 3,
#        units = "in")

# ggsave(filename = "scripts/model-pd-md/figs/fig_8b.png",
#        plot = fig8b,
#        width = 4, height = 4,
#        units = "in")
