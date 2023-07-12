load("scripts/model-pd-md/met_in.Rdata")
load("scripts/model-pd-md/psy_in.Rdata")


### Test WP againts VWC
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
            n = n()) %>%
  left_join(select(met_in, date, starts_with("VWC")))

cor(mean_psy$MD_mean, mean_psy$VWC_5cm)
cor(mean_psy$MD_mean, mean_psy$VWC_10cm)
cor(mean_psy$MD_mean, mean_psy$VWC_20cm)
cor(mean_psy$MD_mean, mean_psy$VWC_50cm)
ggplot(mean_psy, aes(y = MD_mean)) +
  geom_point(aes(x = VWC_5cm, color = "5 cm")) +
  geom_point(aes(x = VWC_10cm, color = "10 cm")) +
  geom_point(aes(x = VWC_20cm, color = "20 cm"))

cor(mean_psy$PD_mean, mean_psy$VWC_5cm, use = "complete.obs")
cor(mean_psy$PD_mean, mean_psy$VWC_10cm, use = "complete.obs")
cor(mean_psy$PD_mean, mean_psy$VWC_20cm, use = "complete.obs")
cor(mean_psy$PD_mean, mean_psy$VWC_50cm, use = "complete.obs")
ggplot(mean_psy, aes(y = PD_mean)) +
  geom_point(aes(x = VWC_5cm, color = "5 cm")) +
  geom_smooth(aes(x = VWC_5cm, color = "5 cm"),
              method = "lm", se = FALSE) +
  geom_point(aes(x = VWC_10cm, color = "10 cm")) +
  geom_smooth(aes(x = VWC_10cm, color = "10 cm"),
              method = "lm", se = FALSE) +
  geom_point(aes(x = VWC_20cm, color = "20 cm")) +
  geom_smooth(aes(x = VWC_20cm, color = "20 cm"),
              method = "lm", se = FALSE)

met_in %>%
  ggplot() + 
  geom_point(aes(x = date, y = PAR_In))
