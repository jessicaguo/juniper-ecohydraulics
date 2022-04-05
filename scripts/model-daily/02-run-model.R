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

# Van Genuchten function for these soils
# vg <- function(x) {
#   theta.s <- -12.211
#   theta.r <- 91.395
#   alpha <- 33.772
#   n <- 1.0205
#   sat <- (1/(1 + ((alpha*x)^n)))^(1 - (1/n))
#   wp <- theta.r + (theta.s - theta.r)*sat
#   return(wp)
# }
# 
# # add SWP to met data
# met_in <- met_in %>%
#   mutate(SWP_10cm = vg(VWC_10cm),
#          SWP_50cm = vg(VWC_50cm))

# Limit psy_in to PD only
psy_in <- psy_in %>%
  filter(!is.na(PD))

# Dmax, VWC5, VWC10, VWC50 are already scaled to entire record

# Plotting for each tree
for(i in 1:7) {
  ggplot() +
    geom_point(data = met_in,
               aes(x = date, y = Dmax, 
                   col = "Dmax")) +
    geom_point(data = met_in,
               aes(x = date, y = VWC10, 
                   col = "10 cm")) +
    geom_point(data = met_in,
               aes(x = date, y = VWC50, 
                   col = "50 cm")) +
    geom_bar(data = met_in,
             aes(x = date, y = Precip/10,
                 col = "ppt"),
             stat = "identity") +
    geom_point(data = filter(psy_in, Tree == i),
               aes(x = date, y = PD, shape = Logger,
                   col = "PD")) +
    scale_y_continuous("Variables", limits = c(-9, 4)) +
    theme_bw(base_size = 12)
  
  ggsave(filename = paste0("scripts/daily-model/plots/ts_met_pd_", i, ".png"),
         height = 4,
         width = 8, 
         units = "in")
}

# Plot with scaled vars
psy2 <- psy_in %>%
  left_join(met_in)
ggplot(psy2, aes(x = Dmax, y = PD)) +
  geom_vline(xintercept = 0) +
  geom_point(aes(color = as.factor(branch_fixed),
                 shape = Logger)) +
  geom_smooth(aes(color = as.factor(branch_fixed)),
              method = lm,
              se = FALSE) +
  scale_y_continuous(expression(paste(Psi[PD], " (MPa)"))) +
  scale_x_continuous(expression(paste("Scaled ", D[max], " (kPa)"))) +
  facet_grid(rows = vars(Tree), 
             cols = vars(Logger),
             scales = "free") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_blank())
ggsave(filename = "scripts/daily-model/plots/PD_Dmax.png",
       height = 8,
       width = 6, 
       units = "in")

ggplot(psy2, aes(x = VWC10, y = PD)) +
  geom_vline(xintercept = 0) +
  geom_point(aes(color = as.factor(branch_fixed),
                 shape = Logger)) +
  geom_smooth(aes(color = as.factor(branch_fixed)),
              method = lm,
              se = FALSE) +
  scale_y_continuous(expression(paste(Psi[PD], " (MPa)"))) +
  scale_x_continuous(expression(paste("Scaled ", VWC[10]))) +
  facet_grid(rows = vars(Tree), 
             cols = vars(Logger),
             scales = "free") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_blank())
ggsave(filename = "scripts/daily-model/plots/PD_VWC10.png",
       height = 8,
       width = 6, 
       units = "in")

ggplot(psy2, aes(x = VWC50, y = PD)) +
  geom_vline(xintercept = 0) +
  geom_point(aes(color = as.factor(branch_fixed),
                 shape = Logger)) +
  geom_smooth(aes(color = as.factor(branch_fixed)),
              method = lm,
              se = FALSE) +
  scale_y_continuous(expression(paste(Psi[PD], " (MPa)"))) +
  scale_x_continuous(expression(paste("Scaled ", VWC[50]))) +
  facet_grid(rows = vars(Tree), 
             cols = vars(Logger),
             scales = "free") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_blank())
ggsave(filename = "scripts/daily-model/plots/PD_VWC50.png",
       height = 8,
       width = 6, 
       units = "in")

ggplot(psy2, aes(x = VWC10, y = PD)) +
  geom_vline(xintercept = 0) +
  geom_point(aes(color = as.factor(Logger))) +
  geom_smooth(aes(color = as.factor(Logger)),
              method = lm,
              se = FALSE) +
  scale_y_continuous(expression(paste(Psi[PD], " (MPa)"))) +
  scale_x_continuous(expression(paste("Scaled ", D[max], " (kPa)"))) +
  facet_grid(rows = vars(Tree), 
             cols = vars(Logger),
             scales = "free") +
  theme_bw(base_size = 12) +
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
                 logger = psy_in$Logger,
                 # time step number
                 nlagA = 7,
                 nlagB = 7,
                 nlagC = 7,
                 # time step size (single value if modelb, vector of values if modeld),
                 pA = 1,
                 pB = 1,
                 pC = 5,
                 # weights, of length nlag
                 alphaA = rep(1, 7), 
                 alphaB = rep(1, 7),
                 alphaC = rep(1, 7),
                 doy = psy_in$doy,
                 Dmax = as.vector(met_in$Dmax),
                 VWC10 = as.vector(met_in$VWC10),
                 VWC50 = as.vector(met_in$VWC50),
                 NBranch = nrow(branch_in),
                 NLogger = length(unique(psy_in$Logger)),
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
                 inits = inits_list,
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
                         n.iter = 3000,
                         n.thin = 5)

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
caterplot(coda.out, regex = "^a\\[2", reorder = FALSE)
caterplot(coda.out, parms = "wA", reorder = FALSE)
caterplot(coda.out, parms = "wB", reorder = FALSE)
caterplot(coda.out, parms = "wC", reorder = FALSE)

foo <- psy_in %>%
  ungroup() %>%
  count(branchID)

psy_in %>%
  ungroup() %>%
  group_by(Tree, Logger) %>%
  summarize(min_PD = min(PD, na.rm = TRUE))

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
which(colnames(coda.out[[1]]) == "deviance")

mean(coda.out[[1]][,1])
mean(coda.out[[2]][,1])
mean(coda.out[[3]][,1])

mean(coda.out[[1]][,453])
mean(coda.out[[2]][,453])
mean(coda.out[[3]][,453])

save(saved_state, file = "scripts/daily-model/inits/saved_state.Rdata")


# Run model for replicated data
coda.rep <- coda.samples(jm, 
                         variable.names = "pd.rep",
                         n.iter = 3000,
                         n.thin = 15)

# save(coda.rep, file = "scripts/daily-model/coda/codarep.Rdata")


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
