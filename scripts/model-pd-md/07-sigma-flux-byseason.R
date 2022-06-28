#  Use the 10% of JAS precip to determine monsoon start day
library(dplyr)
library(ggplot2)
library(ggh4x)
library(harrypotter)

# Load met data
load("data_cleaned/met_daily.Rdata") 

ppt <- met_daily %>%
  filter(date >= as.Date("2021-07-01"),
         date <= as.Date("2021-09-30")) %>%
  select(date, Precip) %>%
  mutate(cPrecip = cumsum(Precip))

tot <- sum(ppt$Precip)

ppt <- ppt %>%
  mutate(season = case_when(cPrecip < 0.1 * tot ~ "premonsoon",
                            cPrecip > 0.1 * tot ~ "monsoon"))

ggplot(ppt, aes(x = date, y = Precip)) +
  geom_bar(stat = "identity",
           aes(color = season))

### Plot sigma and flux by season
# Load predicted sig timeseries
load(file = "scripts/model-pd-md/products/param_pred.Rdata")

# Read in daily fluxes & clip time range
flux <- readr::read_csv("data_raw/US-CdM daily.csv") %>%
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j"))) %>%
  relocate(date) 

pred <- param_pred %>%
  left_join(flux[, 1:7], by = "date") %>%
  left_join(met_daily, by = "date") %>%
  left_join(ppt, by = c("date", "Precip")) %>%
  mutate(season = case_when(is.na(cPrecip) & date <= as.Date("2021-06-30") ~ "premonsoon",
                            season == "premonsoon" ~ "premonsoon",
                            season == "monsoon" ~ "monsoon",
                            is.na(cPrecip) & date >= as.Date("2021-10-01") ~ "fall"),
         season = factor(season, levels = c("premonsoon", "monsoon", "fall")),
         GPP_F = case_when(GPP_F < 0 ~ NA_real_,
                         GPP_F >= 0 ~ GPP_F))


ggplot(pred, aes(x = date)) +
  geom_hline(yintercept = 1, 
             color = "gray",
             size = 1) +
  geom_errorbar(aes(ymin = sigma.lower, ymax = sigma.upper,
                    color = strategy),
                width = 0,
                alpha = 0.5) +
  geom_point(aes(y = sigma.mean, color = strategy)) +
  geom_line(aes(y = GPP_F * 4 + 0.25)) +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 months",
               guide = "axis_minor") +
  scale_y_continuous(expression(sigma),
                     sec.axis = sec_axis(~./4 - 0.25/4 , name = "GPP")) +
  scale_color_hp_d(option = "HermioneGranger", direction = -1,
                   labels = scales::label_parse()) +
  facet_wrap(~season, scale = "free_x", ncol = 3) +
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


ggplot() +
  geom_errorbarh(data = pred, 
                 aes(xmin = sigma.lower,
                     xmax = sigma.upper,
                     y = GPP_F,
                     col = strategy),
                 alpha = 0.25) +
  geom_point(data = pred,
             aes(x= sigma.mean, 
                 y = GPP_F,
                 col = strategy)) +
  facet_wrap(~season, scale = "free_x") +
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

# Extract sigma and GPP from monsoon, run correlation tests

# Shift GPP back by 1 day at a time, see what the correlation becomes
mon <- pred %>%
  filter(season == "monsoon") %>%
  mutate(sigma_0 = sigma.mean,
         sigma_1 = lead(sigma.mean, n = 1),
         sigma_2 = lead(sigma.mean, n = 2),
         sigma_3 = lead(sigma.mean, n = 3),
         sigma_4 = lead(sigma.mean, n = 4),
         sigma_5 = lead(sigma.mean, n = 5))

out <- data.frame(lead = 0:5,
                  cor = NA,
                  cor.lower = NA,
                  cor.upper = NA,
                  p = NA)
for(i in 1:6) {
  sigma_ind <- which(colnames(mon) == paste0("sigma_", i-1))
  m <- cor.test(mon[, sigma_ind], mon$GPP_F)
  out$cor[i] <- m$estimate
  out$cor.lower[i] <- m$conf.int[1]
  out$cor.upper[i] <- m$conf.int[2]
  out$p[i] <- m$p.value
}

ggplot(out, aes(x = lead, y = cor)) +
  geom_errorbar(aes(ymin = cor.lower,
                    ymax = cor.upper),
                width = 0,
                alpha = 0.35) +
  geom_point() +
  scale_x_continuous("Offset (days)") +
  scale_y_continuous(bquote(italic(R))) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12))
