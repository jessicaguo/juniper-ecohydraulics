# Now multivariate model of fluxes NEE, GPP, and ET
# on PD, soil moisture, light, temp, and VPD
# or some combo thereof
# Because the fluxes are partially missing, use OpenBUGS
# Must also setwd() to model-flux-daily
setwd("scripts/model-flux-daily/")

library(R2OpenBUGS)
library(dplyr)
library(ggplot2)
library(mcmcplots)
library(postjags)

# Load gapfilled water potentials
doy22 <- data.frame(date = seq(as.Date("2021-01-01"), 
                               as.Date("2021-12-31"),
                               by = "day"))
load("../../data_cleaned/psy_daily_site_gapfilled.Rdata")
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
load("../../data_cleaned/met_daily.Rdata")
met_in <- met_daily %>%
  mutate(Dmax = scale(VPD_max),
         VWC5 = scale(VWC_5cm),
         VWC10 = scale(VWC_10cm),
         VWC50 = scale(VWC_50cm)) %>%
  filter(date >= as.Date("2021-01-01"))

# Read in met_flux & clip to 2021
met_flux <- readr::read_csv("../../data_raw/US-CdM daily.csv") %>%
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j"))) %>%
  relocate(date) %>%
  filter(date >= as.Date("2021-01-01")) 
  
# Read in daily fluxes & clip time range
flux <- readr::read_csv("../../data_raw/US-CdM daily.csv") %>%
  mutate(date = as.Date(as.POSIXct(paste0(Year, DOY), format = "%Y%j"))) %>%
  relocate(date) %>%
  filter(date >= min(psy_daily_site_gapfilled$date) + 5, # to allow for lags of up to 1 week
         date <= max(psy_daily_site_gapfilled$date))
range(flux$date)

flux_mat <- cbind(flux$NEE_F,
                  flux$GPP_F,
                  flux$ET)
# Combine psy and flux to see parameters
foo <- flux %>%
  left_join(psy) %>%
  left_join(met_in)

summary(lm(ET ~ (PPFD_IN+PD+VWC5+Dmax)^2, data = foo))

foo %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = GPP_F, color = "GPP")) +
  geom_point(aes(y = Dmax, color = "Dmax")) +
  geom_point(aes(y = PD, color = "PD")) +
  geom_point(aes(y = VWC10, color = "vWC10")) +
  geom_point(aes(y = PPFD_IN/1500, color = "PPFD"))

foo %>%
  ggplot(aes(x = PD, y = ET)) +
  geom_point(aes(color = VWC5))

foo %>%
  ggplot(aes(x = VWC5, y = ET)) +
  geom_point(aes(color = Dmax))

foo %>%
  ggplot(aes(x = PPFD_IN, y = NEE_F)) +
  geom_point(aes(color = VWC5))

cor(foo$PPFD_IN, foo$VWC10)
cor(foo$PD, foo$Dmax)
cor(foo$Dmax, foo$VWC5)
cor(foo$PD_scaled, foo$VWC5)
cor(foo$ET, foo$VWC5, use = "complete.obs")
cor(foo$ET, foo$MD, use = "complete.obs")

# Create list of data inputs for model
dat_list <- list(flux = as.matrix(flux_mat),
                PD = as.vector(psy$PD),
                Dmax = as.vector(met_in$Dmax),
                VWC5 = as.vector(met_in$VWC5),
                PPFD = as.vector(met_flux$PPFD_IN),
                N = nrow(flux_mat),
                doy = flux$DOY,
                nlagA = 10,
                nlagB = 5,
                nlagC = 10,
                nlagD = 5,
                pA = 1,
                pB = 1,
                pC = 3,
                pD = 1,
                # weights, of length nlag
                alphaA = matrix(rep(1, 10*3), ncol = 3), 
                alphaB = matrix(rep(1, 5*3), ncol = 3),
                alphaC = matrix(rep(1, 10*3), ncol = 3),
                alphaD = matrix(rep(1, 5*3), ncol = 3),
                NParam = 11,
                R = diag(x = 1, nrow = 3, ncol = 3))

# Function to initialize precision matrix, used below in inits() function
f_mat <- as.matrix(flux_mat)
omega.gen <- function(x){
  noise <- rnorm(n = nrow(f_mat)*ncol(f_mat), mean = 0, sd = 10)
  nois.mat <- matrix(noise, ncol = ncol(f_mat))
  return(solve(cov(f_mat + nois.mat, use = "complete.obs")))
}


# Function to generate initials
init <- function() {
  list(B = matrix(rnorm(dat_list$NParam * 3, 0, 10), ncol = 3),
       omega = round(omega.gen(), 4),
       deltaA = matrix(runif(dat_list$nlagA*3, 0, 1), ncol = 3),
       deltaB = matrix(runif(dat_list$nlagB*3, 0, 1), ncol = 3),
       deltaC = matrix(runif(dat_list$nlagC*3, 0, 1), ncol = 3),
       deltaD = matrix(runif(dat_list$nlagD*3, 0, 1), ncol = 3)
  )
}

inits_list <- list(init(), init(), init(), init())

# Or create initials from lowest deviance chains
load("inits/saved_state_mv.Rdata")
new_inits <- list(saved_state[[2]][[1]],
                  saved_state[[2]][[4]],
                  saved_state[[2]][[4]],
                  saved_state[[2]][[1]])

# Compile and adapt BUGS model
params <- c("deviance", "Dsum",
            "B",
            "Sig", "Rho", "omega", 
            "wA", "wB", "wC", "wD", 
            "deltaA", "deltaB", "deltaC", "deltaD")

start<-proc.time()
bm <- bugs(model.file = "model_mv.txt",
           data = dat_list,
           inits = inits_list,
           n.iter = 4000,
           n.burnin = 1000,
           n.chains = 4,
           parameters.to.save = params,
           codaPkg = TRUE,
           debug = T)
end<-proc.time()
(end-start)/60

# Change to coda object
coda_out <- read.bugs(bm)

if(!dir.exists("coda")){
  dir.create('coda', showWarnings = FALSE)
}
save(coda_out, file = "coda/coda_mv.Rdata")

# Inspect chains visually
mcmcplot(coda_out, parms = c("deviance", "Dsum", 
                             "B", "Sig", "Rho",
                             "wA", "wB", "wC", "wD"))
which(colnames(coda_out[[1]]) == "deviance")
mean(coda_out[[1]][,118])
mean(coda_out[[2]][,118])
mean(coda_out[[3]][,118])
mean(coda_out[[4]][,118])

caterplot(coda_out, regex = "^B\\[[0-9]\\,1\\]", reorder = FALSE)
caterplot(coda_out, regex = "^B\\[[0-9]\\,2\\]", reorder = FALSE)
caterplot(coda_out, regex = "^B\\[[0-9]\\,3\\]", reorder = FALSE)

caterplot(coda_out, regex = "wA\\[[0-9]\\,1\\]", reorder = FALSE)
caterplot(coda_out, regex = "wA\\[[0-9]\\,2\\]", reorder = FALSE)
caterplot(coda_out, regex = "wA\\[[0-9]\\,3\\]", reorder = FALSE)
caterplot(coda_out, regex = "wB\\[[0-9]\\,1\\]", reorder = FALSE)
caterplot(coda_out, regex = "wB\\[[0-9]\\,2\\]", reorder = FALSE)
caterplot(coda_out, regex = "wB\\[[0-9]\\,3\\]", reorder = FALSE)
caterplot(coda_out, regex = "wC\\[[0-9]\\,1\\]", reorder = FALSE)
caterplot(coda_out, regex = "wC\\[[0-9]\\,2\\]", reorder = FALSE)
caterplot(coda_out, regex = "wC\\[[0-9]\\,3\\]", reorder = FALSE)

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
final <- initfind(coda_out, OpenBUGS = TRUE)
final[[1]]
saved_state <- removevars(final, variables = c(2:4, 9:11))
saved_state[[1]]

save(saved_state, file = "inits/saved_state_mv.Rdata")

