# Univariate model of  with GPP ~ fAPAR * PAR * LUE
# where LUE ~ SWC, VPD or some combo/antecedent thereof
# GPP daily in mol CO2 m^-2 d^-1
# Residual model of LUE
# Linear function of predawn WP (concurrent) 

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
  left_join(psy_daily_site_gapfilled %>%
                                filter(type == "PD") %>%
                                select(-n)) %>%
  left_join(met_in)


# Load residuals from model-main.jags
load("scripts/model-gpp-nirv/model-main-resid.Rdata") # resid_df

# Combine psy and flux to plot parameters
# Add residuals of simple model
derived <- cbind.data.frame(flux, resid = scale(resid_df$pred.mean))


derived %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = resid)) +
  geom_point(aes(y = WP_mean +2 , color = "predawn WP")) +
  geom_point(aes(y = Dmax, color = "Dmax")) +
  scale_y_continuous(expression(paste("Scaled resid")), 
                     sec.axis = sec_axis(~.-2,
                                         expression(paste(bar(Psi[PD]), " (MPa)")))) +
  theme_bw() +
  labs(color = "covariate")

#### PART 2, verson with concurrent predawn WP ####

dat_list5 <- list(
  resid = as.vector(scale(resid_df$pred.mean)),
  PD = as.vector(scale(flux$WP_kalman)), # scaled, gap-filled predawns
  Dmax = as.vector(flux$Dmax), # Already scaled
  N = nrow(flux),
  NParam = 4)

# Function to generate initials
init5 <- function() {
  list(B = rnorm(dat_list5$NParam, 0, 10),
       tau = runif(1, 0, 1)
  )
}

inits_list5 <- list(init5(), init5(), init5())

load(file = "scripts/model-gpp-nirv/inits/saved_state-resid-pd-D.Rdata")

# Initialize model
jm5 <- jags.model("scripts/model-gpp-nirv/model-resid-pd-D.jags",
                 data = dat_list5,
                 inits = saved_state5[[2]],
                 n.chains = 3)

# update(jm5, n.iter = 10000)
dic.samples(jm5, n.iter = 3000, thin = 50)

# Monitor
params5 <- c("deviance", "Dsum", "R2_resid",
            "B",
            "tau", "sig")

coda.out5 <- coda.samples(jm5,
                         variable.names = params5,
                         n.iter = 3000,
                         n.thin = 50)

save(coda.out5, file = "scripts/model-gpp-nirv/coda/coda-resid-pd-D.Rdata")
load(file = "scripts/model-gpp-nirv/coda/coda-resid-pd-D.Rdata")

# Inspect chains visually
mcmcplot(coda.out5, parms = c("deviance", "Dsum", "R2_resid",
                              "B",
                              "tau", "sig"))
caterplot(coda.out5, regex = "^B\\[", reorder = FALSE)

# Check convergence diagnostic
gel5 <- gelman.diag(coda.out5, multivariate = FALSE)
gel5$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("Dsum", rowname) |grepl("deviance", rowname) | grepl("R2", rowname))

gel5$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("sig", rowname))

gel5$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^B", rowname))

# Posterior mean of Dsum and R2
coda5 <- tidyMCMC(coda.out5,
                      conf.int = TRUE,
                      conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)
coda5 %>%
  filter(grepl("Dsum", term) |
           grepl("R2_resid", term))


# # Save state
# final <- initfind(coda.out5, OpenBUGS = FALSE)
# final[[1]]
# saved_state5 <- removevars(final, variables = c(2:4))
# saved_state5[[1]]
# 
# save(saved_state5, file = "scripts/model-gpp-nirv/inits/saved_state-resid-pd-D.Rdata")


# Run replicated data
coda.rep5 <- coda.samples(jm5,
                         variable.names = "resid.rep",
                         n.iter = 3000,
                         n.thin = 50)

# Save out
save(coda.rep5, file = "scripts/model-gpp-nirv/coda/codarep-resid-pd-D.Rdata")
load(file = "scripts/model-gpp-nirv/coda/codarep-resid-pd-D.Rdata")

coda_sum5 <- tidyMCMC(coda.rep5,
                      conf.int = TRUE,
                      conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)

pred5 <- cbind.data.frame(resid_df %>%
                            mutate(resid_mean = scale(pred.mean)) %>%
                            select(resid_mean), coda_sum5)

m5 <- lm(pred.mean ~ resid_mean, data = pred5)
sm5 <- summary(m5) # R2 = 0.2027

ggplot(pred5) +
  geom_abline(intercept = 0, slope = 1, col = "black",
              linewidth = 1) +
  geom_abline(intercept = coef(sm5)[1,1], 
              slope = coef(sm5)[2,1], 
              col = "black",
              lty = 2) +
  geom_point(aes(x = resid_mean,
                 y = pred.mean)) +
  geom_errorbar(aes(x = resid_mean,
                    ymin = pred.lower,
                    ymax = pred.upper),
                width = 0,
                alpha = 0.2) +
  scale_y_continuous("Predicted residuals") +
  scale_x_continuous("Observed residuals") +
  theme_bw(base_size = 14)


