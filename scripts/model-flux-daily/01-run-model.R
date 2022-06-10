# Univariate first with GPP, then multivariate model of fluxes
# on PD, soil moisture, light, temp, and VPD
# or some combo thereof
# GPP daily in mol CO2 m^-2 d^-1
# Can remove <0 or make zero
# Average daily PPFD between 9 am and 5 pm

library(rjags)
load.module('dic')
library(dplyr)
library(ggplot2)
library(mcmcplots)
library(postjags)
library(broom.mixed)

# Load gapfilled water potentials
load("data_cleaned/psy_daily_site_gapfilled.Rdata")

# Load met data
load("data_cleaned/met_daily.Rdata")
met_in <- met_daily %>%
  mutate(Dmax = scale(VPD_max),
         VWC5 = scale(VWC_5cm),
         VWC10 = scale(VWC_10cm),
         VWC50 = scale(VWC_50cm),
         PAR = scale(PAR_In)) %>%
  filter(date >= as.Date("2021-01-01"))

# Read in daily fluxes & clip time range
flux <- readr::read_csv("data_raw/US-CdM daily.csv") %>%
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j"))) %>%
  relocate(date) %>%
  filter(date >= min(psy_daily_site_gapfilled$date), 
         date <= max(psy_daily_site_gapfilled$date))
range(flux$date)

# Combine psy and flux to see parameters
foo <- flux %>%
  left_join(met_in)

foo %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = GPP_F*10, color = "GPP")) +
  geom_point(aes(y = VWC10, color = "VWC10")) +
  geom_point(aes(y = Dmax, color = "Dmax"))

cor(foo$GPP_F, foo$Dmax)
cor(foo$GPP_F, foo$VWC10)
cor(foo$VWC10, foo$Dmax)

# Create list of data inputs for model
dat_list <- list(GPP = flux$GPP_F,
                 W10 = as.vector(met_in$VWC10),
                 Dmax = as.vector(met_in$Dmax),
                 N = nrow(flux),
                 doy = flux$DOY,
                 nlagA = 5,
                 nlagB = 7,
                 pA = 1,
                 pB = 3,
                 # weights, of length nlag
                 alphaA = rep(1, 5),
                 alphaB = rep(1, 7),
                 NParam = 4)

# Function to generate initials
init <- function() {
  list(B = rnorm(dat_list$NParam, 0, 10),
       tau = runif(1, 0, 1),
       deltaA = runif(dat_list$nlagA, 0, 1),
       deltaB = runif(dat_list$nlagB, 0, 1)
       # deltaC = runif(dat_list$nlagC, 0, 1)
  )
}

inits_list <- list(init(), init(), init())

# Alternatively, load saved state
load("scripts/model-flux-daily/inits/saved_state.Rdata")

# Initialize model
jm <- jags.model("scripts/model-flux-daily/model.jags",
                 data = dat_list,
                 inits = saved_state[[2]],
                 n.chains = 3)

# update(jm, n.iter = 10000)

# Monitor
params <- c("deviance", "Dsum",
            "B",
            "tau", "sig",
            "wA", "wB", 
            "deltaA", "deltaB")
            # "wC", "deltaC")
coda.out <- coda.samples(jm,
                         variable.names = params,
                         n.iter = 3000,
                         n.thin = 50)

# save(coda.out, file = "scripts/model-flux-daily/coda/coda.Rdata")

# Check deviance and pD
# dic.samples(jm, n.iter = 3000)


# Inspect chains visually
mcmcplot(coda.out, parms = c("deviance", "Dsum", 
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
saved_state <- removevars(final, variables = c(2, 5, 7:8))
saved_state[[1]]

# save(saved_state, file = "scripts/model-flux-daily/inits/saved_state.Rdata")

# Run replicated data
coda.rep <- coda.samples(jm,
                         variable.names = "GPP.rep",
                         n.iter = 3000,
                         n.thin = 50)

# Save out
# save(coda.rep, file = "scripts/model-flux-daily/coda/codarep.Rdata")


