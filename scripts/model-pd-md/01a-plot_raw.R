library(dplyr)
library(ggplot2)
library(cowplot)
library(ggh4x)
library(harrypotter)

# Plot env and Psy time series
load("scripts/model-pd-md/met_in.Rdata")
load("scripts/model-pd-md/psy_in.Rdata")
man <- readr::read_csv("data_raw/Pressure-chamber-data.csv") %>%
  mutate(Timestamp = as.POSIXct(Timestamp, format = "%m/%d/%Y %H:%M", 
                                tz = "America/Denver"),
         date = as.Date(Timestamp))

# m1 <- lm(Psych ~ Psi, data = man)
# summary(m1)

# Calculate rectangles for 3 sub-seasons of growing season
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
rect <- data.frame(season = c("premonsoon", "monsoon", "fall"),
                   xmin = c(min(psy_in$date), monsoon_st, as.Date("2021-10-01")),
                   xmax = c(monsoon_st, monsoon_en + 1, max(psy_in$date))) %>%
  mutate(mid = as.Date(rowMeans(cbind(xmin, xmax)), origin = "1970-01-01"))



#  All together in single plot
rect <- data.frame(xmin = min(man$date, na.rm = TRUE),
                   xmax = max(psy_in$date))

met_ppt <- met_in %>%
  mutate(ppt = case_when(Precip > 0 ~ Precip)) 

fig1 <- ggplot() +
  geom_rect(data = rect, aes(ymin = -Inf, ymax = Inf,
                             xmin = xmin, xmax = xmax),
            alpha = 0.1) +
  geom_point(data = met_ppt,
             aes(x = date, y = VPD_max, color = "VPD")) +
  geom_point(data = met_ppt,
             aes(x = date, y = VWC_5cm, color = "5 cm")) +
  geom_point(data = met_ppt,
             aes(x = date, y = VWC_10cm, color = "10 cm")) +
  geom_bar(data = met_ppt,
           aes(x = date, y = ppt/3), stat = "identity") +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 months",
               guide = "axis_minor") +
  scale_y_continuous("VWC (%) | VPD (kPa)", 
                     sec.axis = sec_axis(~.*3, name = "Precipitation (mm)")) +
  scale_color_hp_d(option = "Mischief") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.07, 0.92),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, "cm"),
        ggh4x.axis.ticks.length.minor = rel(1))

# ggsave(filename = "scripts/model-pd-md/figs/fig_1.png",
#        plot = fig1,
#        width = 8, height = 3,
#        units = "in")


# Alternatively, separate VPD + precip from SWC
# Limit to period of psy instrumentation
met_ppt <- met_in %>%
  mutate(ppt = case_when(Precip > 0 ~ Precip)) %>%
  filter(date >= min(psy_in$date) &
           date <= max(psy_in$date))

fig1a <- ggplot() +
  geom_rect(data = rect,
          aes(xmin = xmin, xmax = xmax,
              ymin = -Inf, ymax = Inf,
              fill = season),
          alpha = 0.25) + 
  geom_text(data = rect,
          aes(x = mid, y = 7,
              label = season)) +
  geom_line(data = met_ppt,
             aes(x = date, y = VPD_max, color = "D"),
            size = 1.25) +
  geom_line(data = met_ppt,
             aes(x = date, y = PAR_In/100 -5, color = "PAR"),
            size = 1.25) +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 months",
               guide = "axis_minor") +
  scale_y_continuous("D (kPa)", 
                     sec.axis = sec_axis(~(.+5)*100,
                                         expression(paste("PAR (mmol ",
                                                          m^-2, s^-1, ")")))) +
  scale_color_hp_d(option = "Hufflepuff") +
  scale_fill_manual(values = c("gray90", "gray70", "gray90")) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.066, 0.12),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, "cm"),
        ggh4x.axis.ticks.length.minor = rel(1)) +
  guides(fill = "none")
fig1a
  
fig1b <- ggplot() +
  geom_rect(data = rect,
            aes(xmin = xmin, xmax = xmax,
                ymin = -Inf, ymax = Inf,
                fill = season),
            alpha = 0.25) + 
  geom_text(data = rect,
            aes(x = mid, y = 18,
                label = season)) +
  geom_line(data = met_ppt,
             aes(x = date, y = VWC_5cm, color = "5 cm"),
            size = 1.25) +
  geom_line(data = met_ppt,
             aes(x = date, y = VWC_10cm, color = "10 cm"),
            size = 1.25) +
  geom_bar(data = met_ppt,
           aes(x = date, y = ppt*0.4), stat = "identity",
           fill = "black") +
  geom_point(data = met_ppt,
             aes(x = date, y = T_max*0.4, color = "Temp"),
             size = 1.25) +
  geom_text(data = data.frame(x = min(rect$xmin),
                              y = c(met_ppt$VWC_5cm[1], met_ppt$VWC_10cm[1], met_ppt$T_max[4]*0.4+1),
                              lab = c("5~cm", "10~cm", "T[air]")),
            aes(x = x, y = y, label = lab),
            parse = TRUE,
            vjust = -0.5,
            hjust = 0) +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 months",
               guide = "axis_minor") +
  scale_y_continuous("VWC (%)",
                     sec.axis = sec_axis(~./0.4,
                                         expression(paste("Precip (mm) | ", T[air], " (°C)")))) +
  scale_color_hp_d(option = "Mischief",
                   limits = c("5 cm", "10 cm", "Temp")) +
  scale_fill_manual(values = c("gray90", "gray70", "gray90")) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.066, 0.85),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, "cm"),
        ggh4x.axis.ticks.length.minor = rel(1)) +
  guides(fill = "none",
         color = "none")
fig1b 

fig1_all <- plot_grid(fig1a, fig1b, 
                      ncol = 1,
                      align = "v",
                      labels = "auto",
                      label_size = 12)
ggsave(filename = "scripts/model-pd-md/figs/fig_1_sep.png",
       plot = fig1_all,
       width = 8, height = 5,
       units = "in")

#### Plot water potentials
# Summarize PD and MD across trees, including SE
SE <- function(x) {
  sd(x, na.rm = TRUE) / sum(!is.na(x))
}
mean_psy <- psy_in %>%
  group_by(date) %>%
  summarize(MD_mean = mean(MD, na.rm = TRUE),
            MD_sd = sd(MD, na.rm = TRUE),
            MD_se = SE(MD),
            PD_mean = mean(PD, na.rm = TRUE),
            PD_sd = sd(PD, na.rm = TRUE),
            PD_se = SE(PD),
            n = n())

# Manual measurements
mean_man <- man %>%
  filter(!is.na(Psi)) %>%
  group_by(date) %>%
  summarize(Psi_mean = mean(Psi, na.rm = TRUE),
            Psi_sd = sd(Psi, na.rm = TRUE),
            Psi_se = SE(Psi),
            n = n())

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
rect <- data.frame(season = c("premonsoon", "monsoon", "fall"),
                   xmin = c(min(mean_man$date), monsoon_st, as.Date("2021-10-01")),
                   xmax = c(monsoon_st, monsoon_en + 1, max(psy_in$date))) %>%
  mutate(mid = as.Date(rowMeans(cbind(xmin, xmax)), origin = "1970-01-01"))


fig2a <- ggplot() +
  geom_rect(data = rect,
            aes(xmin = xmin, xmax = xmax,
                ymin = -Inf, ymax = Inf,
                fill = season),
            alpha = 0.25) +
  geom_text(data = rect,
            aes(x = mid, y = -7,
                label = season)) +
  geom_errorbar(data = mean_psy,
                aes(x = date, 
                    ymin = PD_mean - PD_sd,
                    ymax = PD_mean + PD_sd, color = "PD"),
                alpha = 0.5,
                width = 0) +
  geom_errorbar(data = mean_psy,
                aes(x = date, 
                    ymin = MD_mean - MD_sd,
                    ymax = MD_mean + MD_sd, color = "MD"),
                alpha = 0.5,
                width = 0) +
  geom_point(data = mean_psy,
             aes(x = date, 
                 y = PD_mean, 
                 color = "PD")) +
  geom_point(data = mean_psy,
             aes(x = date, 
                 y = MD_mean,
                 color = "MD")) +
  geom_errorbar(data = mean_man,
                aes(x = date, 
                    ymin = Psi_mean - Psi_sd,
                    ymax = Psi_mean + Psi_sd, color = "chamber"),
                alpha = 0.5,
                width = 0) +
  geom_point(data = mean_man,
             aes(x = date, 
                 y = Psi_mean,
                 color = "chamber"),
             size = 2, shape = 15) +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 months",
               guide = "axis_minor") +
  scale_y_continuous(expression(paste(Psi, " (MPa)"))) +
  scale_color_hp_d(option = "Sprout") +
  scale_fill_manual(values = c("gray90", "gray70", "gray90")) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.07, 0.9),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, "cm"),
        ggh4x.axis.ticks.length.minor = rel(1)) +
  guides(color = guide_legend(override.aes = list(shape = c(15, 16, 16),
                                                  linetype = c(0, 0, 0))),
         fill = "none")

# Read in daily fluxes & clip time range
load("data_cleaned/psy_daily_site_gapfilled.Rdata")
flux <- readr::read_csv("data_raw/US-CdM daily.csv") %>%
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j")),
         GPP = case_when(GPP_F > 0 ~ GPP_F)) %>%
  relocate(date) %>%
  filter(date >= min(psy_daily_site_gapfilled$date), 
         date <= max(psy_daily_site_gapfilled$date))

fig2b <- ggplot() +
  geom_rect(data = rect,
            aes(xmin = xmin, xmax = xmax,
                ymin = -Inf, ymax = Inf,
                fill = season),
            alpha = 0.25) +
  geom_text(data = rect,
            aes(x = mid, y = 0,
                label = season)) +
  geom_point(data = filter(flux, GPP_F > 0),
             aes(x = date, y = GPP_F)) +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 months",
               guide = "axis_minor") +
  scale_y_continuous(expression(paste("GPP (mol ", CO[2], " ", m^-2, d^-1, ")"))) +
  scale_fill_manual(values = c("gray90", "gray70", "gray90")) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank()) +
  guides(fill = "none")

fig2 <- plot_grid(fig2a, fig2b,
          ncol = 1,
          align = "v",
          labels = "auto",
          label_size = 12)

ggsave(filename = "scripts/model-pd-md/figs/fig_2.png",
       plot = fig2,
       width = 8, height = 5,
       units = "in")


# Join to label each psy with season and plot in MD vs. PD space
mean_psy2 <- mean_psy %>%
  left_join(ppt, by = "date") %>%
  select(-Precip) %>%
  mutate(season = case_when(is.na(cPrecip) & date <= as.Date("2021-06-30") ~ "premonsoon",
                            season == "premonsoon" ~ "premonsoon",
                            season == "monsoon" ~ "monsoon",
                            is.na(cPrecip) & date >= as.Date("2021-10-01") ~ "fall"),
         season = factor(season, levels = c("premonsoon", "monsoon", "fall"))) %>%
  left_join(met_in, by = "date")

# VPD as color axis
ggplot(mean_psy2, aes(x = PD_mean, y = MD_mean)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_errorbar(aes(ymin = MD_mean - MD_se, ymax = MD_mean + MD_se,
                    col = VPD_max),
                width = 0,
                alpha = 0.25) +
  geom_errorbarh(aes(xmin = PD_mean - PD_se, xmax = PD_mean + PD_se,
                     col = VPD_max),
                 width = 0,
                 alpha = 0.25) +
  geom_point(aes(col = VPD_max)) +
  facet_wrap(~season) +
  coord_equal() +
  scale_color_hp(option = "Gryffindor", direction = -1,
                 name = bquote(D[max])) +
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank())

fig3 <- ggplot() +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(data = mean_psy2,
             aes(x = PD_mean, y = MD_mean,
                 col = VWC_10cm)) +
  facet_wrap(~season) +
  scale_x_continuous(expression(paste(Psi[PD], " (MPa)")), 
                     limits = c(-5, 0),
                     breaks = c(-4, -2, 0)) +
  scale_y_continuous(expression(paste(Psi[MD], " (MPa)")),
                     limits = c(-5, 0),
                     breaks = c(-4, -2, 0)) +
  coord_equal() +
  scale_color_hp(option = "Sprout", direction = 1,
                 name = bquote(W["10"])) +
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
fig3

ggsave(filename = "scripts/model-pd-md/figs/fig_3.png",
       plot = fig3,
       width = 8, height = 3,
       units = "in")
