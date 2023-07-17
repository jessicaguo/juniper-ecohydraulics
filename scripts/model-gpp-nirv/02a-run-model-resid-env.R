# Univariate model of  with GPP ~ fAPAR * PAR * LUE
# where LUE ~ SWC, VPD or some combo/antecedent thereof
# GPP daily in mol CO2 m^-2 d^-1
# Residual model of LUE
# Linear function of Dmax (concurrent) and VWC10 (antecedent)

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

# Read in NIRv & scale
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
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j")),
         GPP = case_when(GPP_F > 0 ~ GPP_F)) %>% # should we only use positive values?
  relocate(date) %>%
  filter(date >= min(psy_daily_site_gapfilled$date), 
         date <= max(psy_daily_site_gapfilled$date)) %>%
  mutate(season = case_when(date < monsoon_st ~ "premonsoon",
                            date >= monsoon_st & date <= monsoon_en ~ "monsoon",
                            date > monsoon_en ~ "fall")) %>%
  left_join(sat_dat) %>%
  left_join(met_in)

# Load residuals from model-main.jags
load("scripts/model-gpp-nirv/model-main-resid.Rdata") # resid_df

# Combine psy and flux to plot parameters
# Add residuals of simple model
derived <- cbind.data.frame(flux, resid = scale(resid_df$pred.mean))


derived %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = resid)) +
  geom_point(aes(y = Dmax, color = "Dmax")) +
  geom_point(aes(y = VWC10, color = "VWC10")) +
  scale_color_manual(values = c("VWC10" = "forestgreen",
                                "Dmax" = "royalblue")) +
  theme_bw() +
  labs(color = "covariate")


#### PART 2, verson with concurrent Dmax and antecedent VWC10 ####

dat_list2 <- list(
  resid = as.vector(scale(resid_df$pred.mean)),
  Dmax = as.vector(met_in$Dmax), # Already scaled
  W10 = as.vector(met_in$VWC10), # Already scaled
  N = nrow(flux),
  doy = flux$DOY,
  nlagA = 1,
  nlagB = 7,
  pA = 1,
  pB = 3,
  # weights, of length nlag
  alphaA = rep(1, 1),
  alphaB = rep(1, 7),
  NParam = 4)

# Function to generate initials
init2 <- function() {
  list(B = rnorm(dat_list2$NParam, 0, 10),
       tau = runif(1, 0, 1),
       deltaA = runif(dat_list2$nlagA, 0, 1),
       deltaB = runif(dat_list2$nlagB, 0, 1)
  )
}

inits_list2 <- list(init2(), init2(), init2())

load(file = "scripts/model-gpp-nirv/inits/saved_state-resid-env.Rdata")

# Initialize model
jm2 <- jags.model("scripts/model-gpp-nirv/model-resid-env.jags",
                 data = dat_list2,
                 inits = saved_state2[[2]],
                 n.chains = 3)

# update(jm2, n.iter = 10000)
dic.samples(jm2, n.iter = 3000, thin = 50)

# Monitor
params2 <- c("deviance", "Dsum", "R2_resid",
            "B",
            "tau", "sig",
            "wA", "wB",
            "deltaA", "deltaB")

coda.out2 <- coda.samples(jm2,
                         variable.names = params2,
                         n.iter = 3000,
                         n.thin = 50)

save(coda.out2, file = "scripts/model-gpp-nirv/coda/coda-resid-env.Rdata")
load(file = "scripts/model-gpp-nirv/coda/coda-resid-env.Rdata")

# Inspect chains visually
mcmcplot(coda.out2, parms = c("deviance", "Dsum", "R2_resid",
                              "B",
                              "tau", "sig",
                              "wA", "wB"))
caterplot(coda.out2, regex = "^B\\[", reorder = FALSE)
caterplot(coda.out2, regex = "^wA", reorder = FALSE)
caterplot(coda.out2, regex = "^wB\\[", reorder = FALSE)

# Check convergence diagnostic
gel2 <- gelman.diag(coda.out2, multivariate = FALSE)
gel2$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("Dsum", rowname) |grepl("deviance", rowname) | grepl("R2", rowname))

gel2$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("sig", rowname))

gel2$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^B", rowname))

gel2$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^wB", rowname))

# Posterior mean of Dsum and R2
coda2 <- tidyMCMC(coda.out2,
                  conf.int = TRUE,
                  conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)
coda2 %>%
  filter(grepl("Dsum", term) |
           grepl("R2_resid", term))

# Save state
# final <- initfind(coda.out2, OpenBUGS = FALSE)
# final[[1]]
# saved_state2 <- removevars(final, variables = c(2:3, 6, 8:9))
# saved_state2[[1]]
# 
# save(saved_state2, file = "scripts/model-gpp-nirv/inits/saved_state-resid-env.Rdata")


# Run replicated data
coda.rep2 <- coda.samples(jm2,
                         variable.names = "resid.rep",
                         n.iter = 3000,
                         n.thin = 50)

# Save out
save(coda.rep2, file = "scripts/model-gpp-nirv/coda/codarep-resid-env.Rdata")

coda_sum2 <- tidyMCMC(coda.rep2,
                      conf.int = TRUE,
                      conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)

pred2 <- cbind.data.frame(resid_df %>%
                            mutate(resid_mean = scale(pred.mean)) %>%
                            select(resid_mean), coda_sum2)

m2 <- lm(pred.mean ~ resid_mean, data = pred2)
sm2 <- summary(m2) # R2 = 0.1163

ggplot(pred2) +
  geom_abline(intercept = 0, slope = 1, col = "black",
              linewidth = 1) +
  geom_abline(intercept = coef(sm2)[1,1], 
              slope = coef(sm2)[2,1], 
              col = "black",
              lty = 2) +
  geom_point(aes(x = resid_mean,
                 y = pred.mean)) 


