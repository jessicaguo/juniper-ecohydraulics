# Run model script

library(dplyr)
library(rjags)
load.module("dic")
library(mcmcplots)
library(postjags)
library(broom.mixed)
library(ggplot2)

load("scripts/daily-model/met_in.Rdata")
load("scripts/daily-model/psy_in.Rdata")
load("scripts/daily-model/branch_in.Rdata")

# Center and scale env vars
Dmax <- scale(met_in$VPD_max)
VWC5 <- scale(met_in$VWC_5cm)
VWC50 <- scale(met_in$VWC_50cm)

# Estimates of sd
psy_in %>%
  group_by(Tree) %>%
  summarize(meanPD = mean(PD, na.rm = T)) %>%
  summarize(sd(meanPD))

psy_in %>%
  group_by(branchID) %>%
  summarize(Tree = unique(Tree),
            meanPD = mean(PD, na.rm = T)) %>%
  group_by(Tree) %>%
  summarize(sd.tree = sd(meanPD, na.rm = T))
# Prepare data list
dat_list <- list(N = nrow(psy_in),
                 pd = psy_in$PD,
                 branch = psy_in$branchID,
                 # time step number
                 nlagA = 7,
                 nlagB = 7,
                 nlagC = 7,
                 # time step size,
                 pA = 1,
                 pB = 1,
                 pC = 3,
                 # weights, of length nlag
                 alphaA = rep(1, 7), 
                 alphaB = rep(1, 7),
                 alphaC = rep(1, 7),
                 doy = psy_in$doy,
                 Dmax = as.vector(Dmax),
                 VWC5 = as.vector(VWC5),
                 VWC50 = as.vector(VWC50),
                 NBranch = nrow(branch_in),
                 NParam = 7,
                 tree = branch_in$tree,
                 NTree = length(unique(branch_in$tree)),
                 Sa = matrix(rep(4, 7*7), nrow = 7), # NParam x NTree matrix
                 Salpha = rep(2, 7) # vector of NParam
                 )

# Function to generate initials
init <- function() {
  list(mu.alpha = c(runif(1, -12, 0),
                    rnorm(dat_list$NParam-1, 0, 10)),
       tau = runif(1, 0, 1),
       tau.eps.alpha = runif(dat_list$NParam, 0, 1),
       tau.eps.a = matrix(runif(dat_list$NTree * dat_list$NParam, 0, 1),
                          nrow = dat_list$NParam),
       deltaA = runif(dat_list$nlagA, 0, 1),
       deltaB = runif(dat_list$nlagB, 0, 1),
       deltaC = runif(dat_list$nlagC, 0, 1)
  )
}

inits_list <- list(init(), init(), init())

# Initialize model
jm <- jags.model("scripts/daily-model/modelc.jags",
                 data = dat_list,
                 inits = saved_state[[2]],
                 n.chains = 3)

update(jm, n.iter = 10000)

# Monitor
params <- c("deviance", "Dsum",
            "mu.alpha", "sig.alpha",
            "mu.a", "sig.a",
            "a", "tau", "sig",
            "tau.eps.alpha", "tau.eps.a",
            "wA", "wB", "wC",
            "deltaA", "deltaB", "deltaC"
            )
coda.out <- coda.samples(jm,
                         variable.names = params,
                         n.iter = 15000,
                         n.thin = 5)

# Inspect chains visually
mcmcplot(coda.out, parms = c("deviance", "Dsum", "mu.alpha", 
                             "sig.alpha", "sig",
                             "wA", "wB", "wC"))
caterplot(coda.out, parms = "sig.alpha", reorder = FALSE)
caterplot(coda.out, parms = "sig.a", reorder = FALSE)
caterplot(coda.out, parms = "mu.alpha", reorder = FALSE)
caterplot(coda.out, regex = "mu.a\\[1", reorder = FALSE)
caterplot(coda.out, regex = "mu.a\\[2", reorder = FALSE)
caterplot(coda.out, regex = "mu.a\\[3", reorder = FALSE)
caterplot(coda.out, regex = "mu.a\\[4", reorder = FALSE)
caterplot(coda.out, regex = "mu.a\\[5", reorder = FALSE)
caterplot(coda.out, regex = "mu.a\\[6", reorder = FALSE)
caterplot(coda.out, regex = "mu.a\\[7", reorder = FALSE)
caterplot(coda.out, parms = "wA", reorder = FALSE)
caterplot(coda.out, parms = "wB", reorder = FALSE)
caterplot(coda.out, parms = "wC", reorder = FALSE)

# Check convergence diagnostic
gel <- gelman.diag(coda.out, multivariate = FALSE)
gel$psrf %>%
  data.frame() %>%
  head(1)

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("sig", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("mu", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("w", rowname))

# Save state
final <- initfind(coda.out, OpenBUGS = FALSE)
final[[1]]
saved_state <- removevars(final, variables = c(1:2, 6, 8:10, 14:16))
saved_state[[1]]

# Run model for replicated data
coda.rep <- coda.samples(jm, 
                         variable.names = "pd.rep",
                         n.iter = 15000,
                         n.thin = 5)

# Summarize replicated output
coda.sum <- tidyMCMC(coda.rep,
                     conf.int = TRUE,
                     conf.method = "HPDinterval") %>%
  rename(pd.mean = estimate,
         pd.lower = conf.low,
         pd.upper = conf.high)

# Check model fit
pred <- cbind.data.frame(psy_in, coda.sum)

m1 <- lm(pd.mean ~ PD, data = pred)
summary(m1) # R2 = 0.8291; R2 = 0.959

pred %>%
  filter(!is.na(PD)) %>%
  ggplot(aes(x = PD, y = pd.mean)) +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  geom_errorbar(aes(ymin = pd.lower, ymax = pd.upper,
                    color = as.factor(branch_fixed)),
                alpha = 0.25) +
  geom_point(aes(color = as.factor(branch_fixed))) +
  scale_x_continuous("Observed") +
  scale_y_continuous("Predicted") +
  theme_bw(base_size = 12) +
  facet_wrap(~Tree) +
  coord_equal()


pred %>%
  filter(Tree == 1) %>%
  ggplot(aes(x = PD, y = pd.mean)) +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  geom_errorbar(aes(ymin = pd.lower, ymax = pd.upper),
                alpha = 0.5) +
  geom_point() +
  geom_smooth(method = "lm", 
              se = FALSE) +
  scale_x_continuous("Observed") +
  scale_y_continuous("Predicted") +
  theme_bw(base_size = 12) +
  facet_wrap(~branch_fixed,
             nrow = 2)

pred %>%
  filter(Tree == 1) %>%
  ggplot() +
  geom_pointrange(aes(x = VPD_max, y = pd.mean,
                      ymin = pd.lower, ymax = pd.upper,
                      color = "Predicted")) +
  geom_point(aes(x = VPD_max, y = PD,
                 color = "Observed")) +
  scale_x_continuous("VPD max") +
  scale_y_continuous("PD") +
  theme_bw(base_size = 12) +
  facet_wrap(~branch_fixed, 
             nrow = 2) +
  guides(color = "none")

pred %>%
  filter(Tree == 1) %>%
  ggplot() +
  geom_pointrange(aes(x = VPD_max, y = pd.mean,
                      ymin = pd.lower, ymax = pd.upper,
                      color = "Predicted")) +
  geom_point(aes(x = VPD_max, y = PD,
                 color = "Observed")) +
  scale_x_continuous("VPD max") +
  scale_y_continuous("PD") +
  theme_bw(base_size = 12) +
  facet_wrap(~branch_fixed, 
             nrow = 2) +
  guides(color = "none")

pred %>%
  filter(Tree == 1) %>%
  ggplot() +
  geom_pointrange(aes(x = date, y = pd.mean,
                      ymin = pd.lower, ymax = pd.upper,
                      color = "Predicted")) +
  geom_point(aes(x = date, y = PD,
                 color = "Observed")) +
  scale_x_date("Date") +
  scale_y_continuous("PD") +
  theme_bw(base_size = 12) +
  facet_wrap(~branch_fixed, 
             nrow = 2,
             scales = "free") +
  guides(color = "none")
