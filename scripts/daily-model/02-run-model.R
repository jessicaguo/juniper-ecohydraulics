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
  list(mu.alpha = rnorm(dat_list$NParam, 0, 10),
       tau = runif(1, 0, 1),
       tau.eps.alpha = runif(dat_list$NParam, 0, 1),
       tau.eps.a = matrix(runif(dat_list$NTree * dat_list$NParam, 0, 1),
                          nrow = dat_list$NParam)
  )
}

inits_list <- list(init(), init(), init())

# Initialize model
jm <- jags.model("scripts/daily-model/modelb.jags",
                 data = dat_list,
                 inits = inits_list,
                 n.chains = 3)

update(jm, n.iter = 10000)

# Monitor
params <- c("deviance", "Dsum",
            "mu.alpha", "sig.alpha",
            "mu.a", "sig.a",
            "a", "tau", "sig",
            "tau.eps.alpha", "tau.eps.a"
            )
coda.out <- coda.samples(jm,
                         variable.names = params,
                         n.iter = 15000,
                         n.thin = 5)

# Inspect chains visually
mcmcplot(coda.out, parms = c("deviance", "Dsum", "mu.alpha", 
                             "sig.alpha", "sig"))
caterplot(coda.out, parms = "sig.alpha", reorder = FALSE)
caterplot(coda.out, parms = "sig.a", reorder = FALSE)
caterplot(coda.out, parms = "mu.alpha", reorder = FALSE)
caterplot(coda.out, regex = "mu.a\\[2", reorder = FALSE)
caterplot(coda.out, regex = "mu.a\\[3", reorder = FALSE)
caterplot(coda.out, regex = "mu.a\\[4", reorder = FALSE)
caterplot(coda.out, regex = "mu.a\\[5", reorder = FALSE)
caterplot(coda.out, regex = "mu.a\\[6", reorder = FALSE)
caterplot(coda.out, regex = "mu.a\\[7", reorder = FALSE)

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

# Save state
final <- initfind(coda.out, OpenBUGS = FALSE)
final[[1]]
saved_state <- removevars(final, variables = c(1:3, 5:7))
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
summary(m1) # R2 = 0.8291

pred %>%
  filter(!is.na(PD)) %>%
  ggplot(aes(x = PD, y = pd.mean)) +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  geom_errorbar(aes(ymin = pd.lower, ymax = pd.upper),
                alpha = 0.5) +
  geom_point() +
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
