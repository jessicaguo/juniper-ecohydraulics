# Univariate model of  with GPP ~ fAPAR * PAR * LUE
# where LUE ~ SWC, VPD or some combo/antecedent thereof
# GPP daily in mol CO2 m^-2 d^-1
# Main model of GPP ~ fAPAR * PAR, save out residuals for the second model

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

# Plot GPP timeseries with smoothed_NIRv and PAR
flux %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = GPP*10, color = "GPP"),
             color = "black") +
  geom_point(aes(y = smoothed_NIRv, color = "smooth NIRv")) +
  geom_point(aes(y = PAR, color = "PAR")) +
  theme_bw() +
  labs(color = "covariate")


# Create list of data inputs for model
dat_list1 <- list(GPP = flux$GPP,
                  NIRv = as.vector(flux$smoothed_NIRv), # Already scaled
                  PAR = as.vector(flux$PAR), # Already scaled
                  N = nrow(flux),
                  NParam = 4)


# Function to generate initials
init1 <- function() {
  list(A = rnorm(dat_list1$NParam, 0, 10),
       tauGPP = runif(1, 0, 1)
  )
}

inits_list1 <- list(init1(), init1(), init1())

# Alternatively, load saved state
load("scripts/model-gpp-nirv/inits/saved_state-main.Rdata")

# Initialize model
jm1 <- jags.model("scripts/model-gpp-nirv/model-main.jags",
                 data = dat_list1,
                 inits = saved_state1[[2]],
                 n.chains = 3)

# update(jm1, n.iter = 10000)
dic.samples(jm1, n.iter = 3000, thin = 50)

# Monitor
params1 <- c("deviance", "DsumGPP", "R2_GPP",
            "A", "resid",
            "tauGPP", "sigGPP")

coda.out1 <- coda.samples(jm1,
                         variable.names = params1,
                         n.iter = 3000,
                         n.thin = 50)

save(coda.out1, file = "scripts/model-gpp-nirv/coda/coda-main.Rdata")
# Inspect chains visually
mcmcplot(coda.out1, parms = c("deviance", "DsumGPP", "R2_GPP",
                             "A",
                             "tauGPP", "sigGPP"))
caterplot(coda.out1, regex = "^A\\[", reorder = FALSE)

# Check convergence diagnostic
gel1 <- gelman.diag(coda.out1, multivariate = FALSE)
gel1$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("Dsum", rowname) |grepl("deviance", rowname) | grepl("R2", rowname))

gel1$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("sig", rowname))

gel1$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^A", rowname))

# Posterior mean of Dsum and R2
coda1 <- tidyMCMC(coda.out1,
                  conf.int = TRUE,
                  conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)
coda1 %>%
  filter(grepl("Dsum", term) |
           grepl("R2_GPP", term))

# Save state
# final <- initfind(coda.out1, OpenBUGS = FALSE)
# final[[1]]
# saved_state1 <- removevars(final, variables = c(2:5))
# saved_state1[[1]]
# 
# save(saved_state1, file = "scripts/model-gpp-nirv/inits/saved_state-main.Rdata")


# Run replicated data
coda.rep1 <- coda.samples(jm1,
                         variable.names = "GPP.rep",
                         n.iter = 3000,
                         n.thin = 50)
# Save out
save(coda.rep1, file = "scripts/model-gpp-nirv/coda/codarep-main.Rdata")


coda_sum1 <- tidyMCMC(coda.rep1,
                     conf.int = TRUE,
                     conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)

pred1 <- cbind.data.frame(flux, coda_sum1)

m1 <- lm(pred.mean ~ GPP, data = pred1)
sm1 <- summary(m1) # R2 = 0.7308

ggplot(pred1) +
  geom_abline(intercept = 0, slope = 1, col = "black",
              size = 1) +
  geom_abline(intercept = coef(sm1)[1,1], 
              slope = coef(sm1)[2,1], 
              col = "black",
              lty = 2) +
  geom_point(aes(x = GPP,
                 y = pred.mean)) +
  coord_fixed()

# Process posterior means of residuals
resid_df <- tidyMCMC(coda.out1,
                      conf.int = TRUE,
                      conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high) %>%
  filter(grepl("resid", term))

save(resid_df, file = "scripts/model-gpp-nirv/model-main-resid.Rdata")


