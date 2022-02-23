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

m <- lm(PD ~ (scale(VPD_max) + scale(VWC_5cm) + scale(VWC_50cm))^2, 
        data = psy_in)
m <- lm(PD ~ (VPD_max + VWC_5cm + VWC_50cm)^2, 
        data = psy_in)
m <- lm(PD ~ (scale(VPD_max, center = F) + 
                scale(VWC_5cm, center = F) + 
                scale(VWC_50cm, center = F))^2, 
        data = psy_in)
summary(m)

# Plot with scale vars
ggplot(psy_in, aes(x = scale(VPD_max), y = PD)) +
  geom_vline(xintercept = 0) +
  geom_point(aes(color = as.factor(branch_fixed))) +
  geom_smooth(aes(color = as.factor(branch_fixed)),
             method = lm,
             se = FALSE) +
  scale_y_continuous(expression(paste(Psi[PD], " (MPa)")),
                     limits = c(-10, 0),
                     breaks = seq(-9, 0, 3)) +
  scale_x_continuous(expression(paste(D[max], " (kPa)"))) +
  facet_wrap(~Tree, nrow = 1) +
  guides(color = "none") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank())

ggplot(psy_in, aes(x = scale(VWC_5cm), y = PD)) +
  geom_vline(xintercept = 0) +
  geom_point(aes(color = as.factor(branch_fixed))) +
  geom_smooth(aes(color = as.factor(branch_fixed)),
              method = lm,
              se = FALSE) +
  scale_y_continuous(expression(paste(Psi[PD], " (MPa)")),
                     limits = c(-10, 0),
                     breaks = seq(-9, 0, 3)) +
  scale_x_continuous(expression(paste(VWC[5]))) +
  facet_wrap(~Tree, nrow = 1) +
  guides(color = "none") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank())

ggplot(psy_in, aes(x = scale(VWC_50cm), y = PD)) +
  geom_vline(xintercept = 0) +
  geom_point(aes(color = as.factor(branch_fixed))) +
  geom_smooth(aes(color = as.factor(branch_fixed)),
              method = lm,
              se = FALSE) +
  scale_y_continuous(expression(paste(Psi[PD], " (MPa)")),
                     limits = c(-10, 0),
                     breaks = seq(-9, 0, 3)) +
  scale_x_continuous(expression(paste(VWC[50]))) +
  facet_wrap(~Tree, nrow = 1) +
  guides(color = "none") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank())

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
                 # Dmax = met_in$VPD_max,
                 # VWC5 = met_in$VWC_5cm,
                 # VWC50 = met_in$VWC_50cm,
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

# Alternative, start from saved state
load("scripts/daily-model/inits/saved_state.Rdata")

# Initialize model
jm <- jags.model("scripts/daily-model/modelb.jags",
                 data = dat_list,
                 inits = list(saved_state[[2]][[2]],
                              saved_state[[2]][[1]],
                              saved_state[[2]][[2]]),
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
                         n.iter = 1000,
                         n.thin = 20)

# save(coda.out, file = "scripts/daily-model/coda/coda.Rdata")

# Inspect chains visually
mcmcplot(coda.out, parms = c("deviance", "Dsum", "mu.alpha", 
                             "sig.alpha", "sig",
                             "wA", "wB", "wC"))
caterplot(coda.out, parms = "sig.alpha", reorder = FALSE)
caterplot(coda.out, regex = "sig.a\\[1", reorder = FALSE)
caterplot(coda.out, regex = "sig.a\\[2", reorder = FALSE)
caterplot(coda.out, regex = "sig.a\\[3", reorder = FALSE)
caterplot(coda.out, regex = "sig.a\\[4", reorder = FALSE)
caterplot(coda.out, regex = "sig.a\\[5", reorder = FALSE)
caterplot(coda.out, regex = "sig.a\\[6", reorder = FALSE)
caterplot(coda.out, regex = "sig.a\\[7", reorder = FALSE)
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

mean(coda.out[[1]][,1])
mean(coda.out[[2]][,1])
mean(coda.out[[3]][,1])

save(saved_state, file = "scripts/daily-model/inits/saved_state.Rdata")


# Run model for replicated data
coda.rep <- coda.samples(jm, 
                         variable.names = "pd.rep",
                         n.iter = 1000,
                         n.thin = 20)

save(coda.rep, file = "scripts/daily-model/coda/codarep.Rdata")


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
  filter(Tree == 3) %>%
  ggplot(aes(x = PD, y = pd.mean)) +
  geom_vline(xintercept = 0, col = "blue", lwd = 2) +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  geom_errorbar(aes(ymin = pd.lower, ymax = pd.upper),
                alpha = 0.5) +
  geom_point() +
  geom_smooth(method = "lm", 
              se = FALSE,
              fullrange = TRUE) +
  scale_x_continuous("Observed") +
  scale_y_continuous("Predicted") +
  theme_bw(base_size = 12) +
  facet_wrap(~branch_fixed,
             nrow = 2,
             scales = "free")

# pred %>%
#   filter(Tree == 1) %>%
#   ggplot() +
#   geom_pointrange(aes(x = VPD_max, y = pd.mean,
#                       ymin = pd.lower, ymax = pd.upper,
#                       color = "Predicted")) +
#   geom_point(aes(x = VPD_max, y = PD,
#                  color = "Observed")) +
#   scale_x_continuous("VPD max") +
#   scale_y_continuous("PD") +
#   theme_bw(base_size = 12) +
#   facet_wrap(~branch_fixed, 
#              nrow = 2) +
#   guides(color = "none")
# 
# pred %>%
#   filter(Tree == 1) %>%
#   ggplot() +
#   geom_pointrange(aes(x = VPD_max, y = pd.mean,
#                       ymin = pd.lower, ymax = pd.upper,
#                       color = "Predicted")) +
#   geom_point(aes(x = VPD_max, y = PD,
#                  color = "Observed")) +
#   scale_x_continuous("VPD max") +
#   scale_y_continuous("PD") +
#   theme_bw(base_size = 12) +
#   facet_wrap(~branch_fixed, 
#              nrow = 2) +
#   guides(color = "none")
# 
# pred %>%
#   filter(Tree == 1) %>%
#   ggplot() +
#   geom_pointrange(aes(x = date, y = pd.mean,
#                       ymin = pd.lower, ymax = pd.upper,
#                       color = "Predicted")) +
#   geom_point(aes(x = date, y = PD,
#                  color = "Observed")) +
#   scale_x_date("Date") +
#   scale_y_continuous("PD") +
#   theme_bw(base_size = 12) +
#   facet_wrap(~branch_fixed, 
#              nrow = 2,
#              scales = "free") +
#   guides(color = "none")
