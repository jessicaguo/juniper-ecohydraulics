# Use imputation package to gapfill daily PD and MD for flux analysis

library(imputeTS)

# Load site-summarized daily psy data
load("data_cleaned/psy_daily_site.Rdata")
psy_wide <- psy_daily_site %>%
  select(-WP_sd, -n) %>%
  tidyr::pivot_wider(names_from = "type", 
                     values_from = "WP_mean") %>%
  select(-"NA")

# Check on NA's
ggplot_na_distribution(psy_wide$PD)
statsNA(psy_wide$MD)

# na_kalman and na_seadec recommended for complext time series with trend and seasonality
imp_kalman <- na_kalman(psy_wide,
                        model = "auto.arima")
ggplot_na_imputations(x_with_na = psy_wide$PD,
                      x_with_imputations = imp_kalman$PD)
ggplot_na_imputations(x_with_na = psy_wide$MD,
                      x_with_imputations = imp_kalman$MD)

imp_seadec <- na_seadec(psy_wide,
                        algorithm = "interpolation",
                        find_frequency = TRUE)
ggplot_na_imputations(x_with_na = psy_wide$PD,
                      x_with_imputations = imp_seadec$PD)

# Use na_kalman with auto.arima
psy_kalman <- imp_kalman %>%
  tidyr::pivot_longer(-date, 
                      names_to = "type", 
                      values_to = "WP_kalman")
psy_daily_site_gapfilled <- psy_kalman %>%
  left_join(psy_daily_site)

save(psy_daily_site_gapfilled, file = "data_cleaned/psy_daily_site_gapfilled.Rdata")
