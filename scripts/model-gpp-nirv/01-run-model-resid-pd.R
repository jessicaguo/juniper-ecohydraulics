# Univariate model of  with GPP ~ APAR * LUE
# where LUE ~ SWC, VPD or some combo/antecedent thereof
# GPP daily in mol CO2 m^-2 d^-1
# Can remove all points <= 0 (Should we? Not according to Dario Papale)

library(rjags)
load.module('dic')
library(dplyr)
library(readr)
library(ggplot2)
library(mcmcplots)
library(postjags)
library(broom.mixed)

# Load gapfilled water potentials
load("data_cleaned/psy_daily_site_gapfilled.Rdata")

# Load met data
load("data_cleaned/met_daily.Rdata")
met_in <- met_daily %>%
  filter(date >= as.Date("2021-01-01")) %>%
  mutate(Dmax = scale(VPD_max),
         VWC5 = scale(VWC_5cm),
         VWC10 = scale(VWC_10cm),
         VWC50 = scale(VWC_50cm),
         PAR = scale(PAR_In)) 

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

# Read in NIRv
sat_dat <- read_csv("data_cleaned/US-Cdm_NIRv.csv") %>%
  rename(date = Date) %>%
  select(-1) %>%
  mutate(NDVI = scale(NDVI),
         NIRv = scale(NIRv),
         smoothed_NIRv = scale(smoothed_NIRv))

sat_dat %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = NDVI/5, color = "NDVI")) +
  geom_point(aes(y = NIRv, color = "NIRv")) +
  geom_point(aes(y = smoothed_NIRv, color = "smoothed NIRv"))



# Read in daily fluxes, clip time range, assign seasons
flux <- readr::read_csv("data_raw/US-CdM daily.csv") %>%
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j")))%>%
         # GPP = case_when(GPP_F > 0 ~ GPP_F)) %>% # should we only use positive values?
  relocate(date) %>%
  filter(date >= min(psy_daily_site_gapfilled$date), 
         date <= max(psy_daily_site_gapfilled$date)) %>%
  mutate(season = case_when(date < monsoon_st ~ "premonsoon",
                            date >= monsoon_st & date <= monsoon_en ~ "monsoon",
                            date > monsoon_en ~ "fall")) %>%
  left_join(sat_dat) %>%
  left_join(psy_daily_site_gapfilled %>%
                                filter(type == "PD") %>%
                                select(-n)) 

# Simple linear model of GPP_F with smoothed_NIRv
m <- lm(GPP_F ~ smoothed_NIRv, data = flux)
summary(m)

# Combine psy and flux to plot parameters
# Add residuals of simple model
derived <- flux %>%
  left_join(met_in) %>%
  left_join(psy_daily_site_gapfilled %>%
              filter(type == "PD") %>%
              select(-n)) %>%
  mutate(resid_gpp_nirv = scale(m$residuals))

derived %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = GPP_F*10, color = "GPP"),
             color = "black") +
  # geom_point(aes(y = smoothed_NIRv, color = "smooth NIRv")) +
  geom_point(aes(y = NDVI, color = "NDVI")) +
  geom_point(aes(y = VWC10, color = "VWC10")) +
  geom_point(aes(y = Dmax, color = "Dmax")) 

derived %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = resid_gpp_nirv, color = "resid_gpp_nirv")) +
  # geom_point(aes(y = VWC10, color = "VWC10")) +
  # geom_point(aes(y = Dmax, color = "Dmax")) +
  geom_point(aes(y = WP_mean, color = "PD")) 

# Create list of data inputs for model
dat_list1 <- list(GPP = flux$GPP,
                  NIRv = as.vector(flux$smoothed_NIRv),
                  N = nrow(flux))


#### PART 2: resid ~ PD ###
load("scripts/model-gpp-nirv/model-main-resid.Rdata") # resid_df

dat_list3 <- list(
  resid = as.vector(scale(resid_df$pred.mean)),
  PD = as.vector(scale(flux$WP_kalman)),
  N = nrow(flux),
  NParam = 2)

# Function to generate initials
init3 <- function() {
  list(B = rnorm(dat_list3$NParam, 0, 10),
       tau = runif(1, 0, 1)
  )
}

inits_list3 <- list(init3(), init3(), init3())

# Initialize model
jm3 <- jags.model("scripts/model-gpp-nirv/model-resid-pd.jags",
                 data = dat_list3,
                 inits = saved_state3[[2]],
                 n.chains = 3)

update(jm3, n.iter = 10000)

# Monitor
params3 <- c("deviance", "Dsum", "R2_resid",
            "B",
            "tau", "sig")

coda.out3 <- coda.samples(jm3,
                         variable.names = params3,
                         n.iter = 3000,
                         n.thin = 50)

# save(coda.out1, file = "scripts/model-gpp-nirv/coda/coda-env.Rdata")

# Inspect chains visually
mcmcplot(coda.out3, parms = c("deviance", "Dsum", "R2_resid",
                              "B",
                              "tau", "sig"))
caterplot(coda.out3, regex = "^B\\[", reorder = FALSE)

# Check convergence diagnostic
gel3 <- gelman.diag(coda.out3, multivariate = FALSE)
gel3$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("Dsum", rowname) |grepl("deviance", rowname) | grepl("R2", rowname))

gel3$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("sig", rowname))

gel3$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^B", rowname))


# Save state
final <- initfind(coda.out3, OpenBUGS = FALSE)
final[[1]]
saved_state3 <- removevars(final, variables = c(2:3, 6, 8:9))
saved_state3[[1]]

save(saved_state3, file = "scripts/model-gpp-nirv/inits/saved_state-resid-pd.Rdata")


# Run replicated data
coda.rep3 <- coda.samples(jm3,
                         variable.names = "resid.rep",
                         n.iter = 3000,
                         n.thin = 50)

# Save out
# save(coda.rep, file = "scripts/model-gpp-nirv/coda/codarep-env.Rdata")

coda_sum3 <- tidyMCMC(coda.rep3,
                      conf.int = TRUE,
                      conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)

pred3 <- cbind.data.frame(resid_df %>%
                            mutate(resid_mean = scale(pred.mean)) %>%
                            select(resid_mean), coda_sum3)

m3 <- lm(pred.mean ~ resid_mean, data = pred3)
sm3 <- summary(m3) # R2 = 0.1449

ggplot(pred3) +
  geom_abline(intercept = 0, slope = 1, col = "black",
              size = 1) +
  geom_abline(intercept = coef(sm2)[1,1], 
              slope = coef(sm2)[2,1], 
              col = "black",
              lty = 2) +
  geom_point(aes(x = resid_mean,
                 y = pred.mean)) 


