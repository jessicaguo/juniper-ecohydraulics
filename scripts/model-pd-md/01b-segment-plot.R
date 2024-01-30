# Conduct breakpoint analysis of Martinez-Vilalta figure during monsoon
# And update to new fig 3 with trend lines
library(tidyverse)
library(segmented)
library(harrypotter)


# Plot env and Psy time series
load("scripts/model-pd-md/met_in.Rdata")
load("scripts/model-pd-md/psy_in.Rdata")

# Summarize psy to mean by shrub
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

# Add seasons
ppt <- met_in %>%
  filter(date >= as.Date("2021-07-01"),
         date <= as.Date("2021-09-30")) %>%
  dplyr::select(date, Precip) %>%
  mutate(cPrecip = cumsum(Precip))

tot <- sum(ppt$Precip)

ppt <- ppt %>%
  mutate(season = case_when(cPrecip < 0.1 * tot ~ "premonsoon",
                            cPrecip > 0.1 * tot ~ "monsoon"))
monsoon_st <- ppt$date[min(which(ppt$season == "monsoon"))]
monsoon_en <- ppt$date[max(which(ppt$season == "monsoon"))]

# Join mean psychrometer and met/precip data
mean_psy2 <- mean_psy %>%
  left_join(ppt, by = "date") %>%
  dplyr::select(-Precip) %>%
  mutate(season = case_when(is.na(cPrecip) & date <= as.Date("2021-06-30") ~ "premonsoon",
                            season == "premonsoon" ~ "premonsoon",
                            season == "monsoon" ~ "monsoon",
                            is.na(cPrecip) & date >= as.Date("2021-10-01") ~ "fall"),
         season = factor(season, levels = c("premonsoon", "monsoon", "fall"))) %>%
  left_join(met_in, by = "date") |> 
  drop_na(PD_mean) # remove since won't be plotted/predicted anyway

#### Test lm for each season ####
ggplot(mean_psy2) +
  geom_point(aes(x = PD_mean, y = MD_mean, color = PD_mean)) +
  facet_wrap(~season)

m1 <- lm(MD_mean ~ PD_mean, data = filter(mean_psy2, season == "premonsoon"))
m2 <- lm(MD_mean ~ PD_mean, data = filter(mean_psy2, season == "monsoon"))
m3 <- lm(MD_mean ~ PD_mean, data = filter(mean_psy2, season == "fall"))
summary(m1)
summary(m2)
summary(m3)


# testing whether data have breakpoints
m.seg2 <- selgmented(m2,
                    seg.Z = ~PD_mean,
                    Kmax = 2,
                    type = "score",
                    alpha = 0.11)

summary(m.seg2)
m.seg2$psi
slope(m.seg2)

my.fitted <- fitted(m.seg2)
my.model <- data.frame(PD_mean = filter(mean_psy2, season == "monsoon") |> pull(PD_mean), MD_mean = my.fitted)

ggplot(my.model, aes(x = PD_mean, y = MD_mean)) + 
  geom_line() 

# Add fitted points to dataset
mean_psy3 <- mean_psy2 |> 
  mutate(fitted = c(fitted(m1), fitted(m.seg2), fitted(m3)))
predict(m.seg2, newdata = data.frame(PD_mean = seq(-5, 0, by = 0.1)))

# Create df for geom_abline
fit_df <- data.frame(season = c("premonsoon", "monsoon", "monsoon", "fall"),
                     slope = c(coef(m1)[[2]], slope(m.seg2)$PD_mean[1,1], slope(m.seg2)$PD_mean[2,1], coef(m3)[[2]]),
                     int = c(coef(m1)[[1]], intercept(m.seg2)$PD_mean[1], intercept(m.seg2)$PD_mean[2], coef(m3)[[1]]))
# Create df for geom_line
fitted_df <- data.frame(PD_mean = rep(seq(-5, 0, by = 0.1), 3)) |> 
  mutate(season = rep(c("premonsoon", "monsoon", "fall"), each = 51),
         fitted = c(predict(m1, newdata = data.frame(PD_mean = seq(-5, 0, by = 0.1))),
                    predict(m.seg2, newdata = data.frame(PD_mean = seq(-5, 0, by = 0.1))),
                    predict(m3, newdata = data.frame(PD_mean = seq(-5, 0, by = 0.1)))))
                    
# Plot
fig3 <- ggplot() +
  geom_abline(slope = 1, intercept = 0, color = "gray40") +
  geom_point(data = mean_psy3,
             aes(x = PD_mean, y = MD_mean,
                 col = VWC_10cm)) +
  geom_line(data = fitted_df,
            aes(x = PD_mean, y = fitted),
            lty = 2) +
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


ggsave(filename = "scripts/model-pd-md/figs/fig_3_new.png",
       plot = fig3,
       width = 8, height = 3,
       units = "in")

