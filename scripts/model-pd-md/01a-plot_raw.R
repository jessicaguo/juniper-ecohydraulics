library(dplyr)
library(ggplot2)
library(cowplot)
library(ggh4x)

# Plot env and Psy time series
load("scripts/model-pd-md/met_in.Rdata")
load("scripts/model-pd-md/psy_in.Rdata")
man <- readr::read_csv("data_raw/Pressure-chamber-data.csv") %>%
  mutate(Timestamp = as.POSIXct(Timestamp, format = "%m/%d/%Y %H:%M", 
                                tz = "America/Denver"),
         date = as.Date(Timestamp))

# m1 <- lm(Psych ~ Psi, data = man)
# summary(m1)


rect <- data.frame(xmin = min(man$date, na.rm = TRUE),
                   xmax = max(psy_in$date))

met_ppt <- met_in %>%
  mutate(ppt = case_when(Precip > 0 ~ Precip)) 

ggplot() +
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
        legend.position = c(0.1, 0.92),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, "cm"),
        ggh4x.axis.ticks.length.minor = rel(1))


# Summarize PD and MD across trees
mean_psy <- psy_in %>%
  group_by(date) %>%
  summarize(MD_mean = mean(MD, na.rm = TRUE),
            MD_sd = sd(MD, na.rm = TRUE),
            PD_mean = mean(PD, na.rm = TRUE),
            PD_sd = sd(PD, na.rm = TRUE),
            n = n())

mean_man <- man %>%
  filter(!is.na(Psi)) %>%
  group_by(date) %>%
  summarize(Psi_mean = mean(Psi, na.rm = TRUE),
            Psi_sd = sd(Psi, na.rm = TRUE),
            n = n())
ggplot() +
  geom_errorbar(data = mean_psy,
                aes(x = date, 
                    ymin = PD_mean - PD_sd,
                    ymax = PD_mean + PD_sd, color = "PD"),
                alpha = 0.5,
                width = 0) +
  geom_point(data = mean_psy,
             aes(x = date, 
                 y = PD_mean, 
                 color = "PD")) +
  geom_errorbar(data = mean_psy,
                aes(x = date, 
                    ymin = MD_mean - MD_sd,
                    ymax = MD_mean + MD_sd, color = "MD"),
                alpha = 0.5,
                width = 0) +
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
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.1, 0.92),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, "cm"),
        ggh4x.axis.ticks.length.minor = rel(1))



