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
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j"))) %>%
         # GPP = case_when(GPP_F > 0 ~ GPP_F)) %>% # should we only use positive values?
  relocate(date) %>%
  filter(date >= min(psy_daily_site_gapfilled$date), 
         date <= max(psy_daily_site_gapfilled$date)) %>%
  mutate(season = case_when(date < monsoon_st ~ "premonsoon",
                            date >= monsoon_st & date <= monsoon_en ~ "monsoon",
                            date > monsoon_en ~ "fall")) %>%
  left_join(sat_dat)

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
  geom_point(aes(y = Dmax, color = "Dmax")) +
  geom_point(aes(y = WP_mean, color = "PD")) 

derived %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = resid_gpp_nirv, color = "resid_gpp_nirv")) +
  geom_point(aes(y = VWC10, color = "VWC10")) +
  geom_point(aes(y = Dmax, color = "Dmax")) +
  geom_point(aes(y = WP_mean, color = "PD")) 

# Create list of data inputs for model
dat_list <- list(GPP = flux$GPP_F,
                 NIRv = as.vector(flux$smoothed_NIRv),
                 W10 = as.vector(met_in$VWC10),
                 Dmax = as.vector(met_in$Dmax),
                 N = nrow(flux),
                 doy = flux$DOY,
                 nlagA = 1,
                 nlagB = 7,
                 pA = 1,
                 pB = 3,
                 # weights, of length nlag
                 alphaA = rep(1, 1),
                 alphaB = rep(1, 7),
                 NParam = 7)

# Function to generate initials
init <- function() {
  list(B = rnorm(dat_list$NParam, 0, 10),
       tau = runif(1, 0, 1),
       deltaA = runif(dat_list$nlagA, 0, 1),
       deltaB = runif(dat_list$nlagB, 0, 1)
  )
}

inits_list <- list(init(), init(), init())

which(colnames(coda.out[[1]]) == "deviance")
mean(coda.out[[1]][,35])
mean(coda.out[[2]][,35])
mean(coda.out[[3]][,35])

reinits <- list(saved_state[[2]][[2]],
                saved_state[[2]][[1]],
                saved_state[[2]][[2]])

# Alternatively, load saved state
# load("scripts/model-gpp-nirv/inits/saved_state-env.Rdata")

# Initialize model
jm <- jags.model("scripts/model-gpp-nirv/model-env.jags",
                 data = dat_list,
                 inits = inits_list,
                 n.chains = 3)

update(jm, n.iter = 10000)

# Monitor
params <- c("deviance", "Dsum", "R2",
            "B",
            "tau", "sig",
            "wA", "wB", 
            "deltaA", "deltaB")

coda.out <- coda.samples(jm,
                         variable.names = params,
                         n.iter = 3000,
                         n.thin = 50)

# save(coda.out, file = "scripts/model-gpp-nirv/coda/coda-env.Rdata")

# Check deviance and pD
# dic.samples(jm, n.iter = 3000)


# Inspect chains visually
mcmcplot(coda.out, parms = c("deviance", "Dsum", "R2",
                             "B", "sig",
                             "wA", "wB"))
caterplot(coda.out, regex = "^B\\[", reorder = FALSE)
caterplot(coda.out, parms = "wA", reorder = FALSE)
caterplot(coda.out, parms = "wB", reorder = FALSE)

# Check convergence diagnostic
gel <- gelman.diag(coda.out, multivariate = FALSE)
gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("Dsum", rowname) |grepl("deviance", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("sig", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^B", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^wA", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^wB", rowname))


# Save state
final <- initfind(coda.out, OpenBUGS = FALSE)
final[[1]]
saved_state <- removevars(final, variables = c(2:3, 6, 8:9))
saved_state[[1]]

save(saved_state, file = "scripts/model-gpp-nirv/inits/saved_state-env.Rdata")

# Run replicated data
coda.rep <- coda.samples(jm,
                         variable.names = "GPP_resid.rep",
                         n.iter = 3000,
                         n.thin = 50)

# Save out
# save(coda.rep, file = "scripts/model-gpp-nirv/coda/codarep-env.Rdata")


