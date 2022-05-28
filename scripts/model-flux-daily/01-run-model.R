# Univariate first with ET, then multivariate model of fluxes
# on PD, soil moisture, light, temp, and VPD
# or some combo thereof

library(rjags)
load.module('dic')
library(dplyr)
library(ggplot2)
library(mcmcplots)
library(postjags)

# Load gapfilled water potentials
doy22 <- data.frame(date = seq(as.Date("2021-01-01"), 
                               as.Date("2021-12-31"),
                               by = "day"))
load("data_cleaned/psy_daily_site_gapfilled.Rdata")
psy <- doy22 %>%
  left_join(psy_daily_site_gapfilled) %>%
  tidyr::pivot_wider(1:3, 
                     names_from = type, 
                     values_from = WP_kalman) %>%
  select(-`NA`) %>%
  mutate(PD_scaled = scale(PD),
         MD_scaled = scale(MD))
str(psy)
range(psy$date)

# Load met data
load("data_cleaned/met_daily.Rdata")
met_in <- met_daily %>%
  mutate(Dmax = scale(VPD_max),
         VWC5 = scale(VWC_5cm),
         VWC10 = scale(VWC_10cm),
         VWC50 = scale(VWC_50cm)) %>%
  filter(date >= as.Date("2021-01-01"))

# Read in daily fluxes & clip time range
flux <- readr::read_csv("data_raw/US-CdM daily.csv") %>%
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j"))) %>%
  relocate(date) %>%
  filter(date >= min(psy_daily_site_gapfilled$date) + 14, # to allow for lags of up to 1 week
         date <= max(psy_daily_site_gapfilled$date))
range(flux$date)

# Combine psy and flux to see parameters
foo <- flux %>%
  left_join(psy) %>%
  left_join(met_in)

foo %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = ET, color = "ET")) +
  geom_point(aes(y = Dmax, color = "Dmax")) +
  geom_point(aes(y = PD, color = "PD")) +
  geom_point(aes(y = VWC10, color = "vWC10"))

foo %>%
  ggplot(aes(x = PD, y = ET)) +
  geom_point(aes(color = VWC5))

foo %>%
  ggplot(aes(x = VWC5, y = ET)) +
  geom_point(aes(color = Dmax))

cor(foo$PD, foo$Dmax)
cor(foo$Dmax, foo$VWC5)
cor(foo$PD_scaled, foo$VWC5)
cor(foo$ET, foo$VWC5, use = "complete.obs")
cor(foo$ET, foo$MD, use = "complete.obs")

# Create list of data inputs for model
dat_list <- list(ET = flux$ET,
                PD = as.vector(psy$PD),
                Dmax = as.vector(met_in$Dmax),
                VWC5 = as.vector(met_in$VWC5),
                N = nrow(flux),
                doy = flux$DOY,
                nlagA = 7,
                nlagB = 7,
                nlagC = 7,
                pA = 1,
                pB = 1,
                pC = 3,
                # weights, of length nlag
                alphaA = rep(1, 7), 
                alphaB = rep(1, 7),
                alphaC = rep(1, 7),
                NParam = 7)

# Function to generate initials
init <- function() {
  list(B = rnorm(dat_list$NParam, 0, 10),
       tau = runif(1, 0, 1),
       deltaA = runif(dat_list$nlagA, 0, 1),
       deltaB = runif(dat_list$nlagB, 0, 1),
       deltaC = runif(dat_list$nlagC, 0, 1)
  )
}

inits_list <- list(init(), init(), init())

# Initialize model
jm <- jags.model("scripts/model-flux-daily/model.jags",
                 data = dat_list,
                 inits = inits_list,
                 n.chains = 3)

update(jm, n.iter = 10000)

# Monitor
params <- c("deviance", "Dsum",
            "B",
            "tau", "sig",
            "wA", "wB", "wC", 
            "deltaA", "deltaB", "deltaC")
coda.out <- coda.samples(jm,
                         variable.names = params,
                         n.iter = 3000,
                         n.thin = 15)

# save(coda.out, file = "scripts/model-flux-daily/coda/coda.Rdata")

# Check deviance and pD
# dic.samples(jm, n.iter = 3000)


# Inspect chains visually
mcmcplot(coda.out, parms = c("deviance", "Dsum", 
                             "B", "sig",
                             "wA", "wB", "wC"))
caterplot(coda.out, parms = "wA", reorder = FALSE)
caterplot(coda.out, parms = "wB", reorder = FALSE)
caterplot(coda.out, parms = "wC", reorder = FALSE)

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

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^wC", rowname))

# Save state
final <- initfind(coda.out, OpenBUGS = FALSE)
final[[1]]
saved_state <- removevars(final, variables = c(2, 6, 8:10))
saved_state[[1]]

save(saved_state, file = "scripts/model-flux-daily/inits/saved_state.Rdata")

