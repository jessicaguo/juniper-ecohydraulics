# Run model script

library(dplyr)
library(rjags)
load.module("dic")
library(mcmcplots)
library(postjags)
library(broom.mixed)
library(ggplot2)

load("scripts/model-pd-md/met_in.Rdata")
load("scripts/model-pd-md/psy_in.Rdata")
load("scripts/model-pd-md/branch_in.Rdata")

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

# Limit psy_in to both observations
psy_in <- psy_in %>%
  filter(!is.na(PD), !is.na(MD))

# Plotting sigma for each tree
psy_in %>%
  ggplot(aes(x = PD, y = MD)) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(color = VPD_max)) +
  facet_wrap(~Tree) +
  scale_color_gradient(low = "cornflowerblue",
                       high = "coral") +
  theme_bw()


# Prepare data list
dat_list <- list(N = nrow(psy_in),
                 md = psy_in$MD,
                 pd = psy_in$PD,
                 # branch = psy_in$branchID,
                 # logger = psy_in$Logger,
                 # time step number
                 tree = as.numeric(psy_in$Tree),
                 NTree = length(unique(psy_in$Tree)),
                 nlagA = 5,
                 nlagB = 7,
                 # nlagC = 7,
                 # time step size (single value if modelb, vector of values if modeld),
                 pA = 1,
                 pB = 3,
                 # pC = 3,
                 # weights, of length nlag
                 alphaA = rep(1, 10), 
                 alphaB = rep(1, 10),
                 alphaC = rep(1, 10),
                 doy = psy_in$doy,
                 Dmax = as.vector(met_in$Dmax),
                 VWC10 = as.vector(met_in$VWC10),
                 VWC50 = as.vector(met_in$VWC50),
                 # NBranch = nrow(branch_in),
                 # NLogger = length(unique(psy_in$Logger)),
                 NParam = 4
                 # tree = branch_in$tree,
                 # NTree = length(unique(branch_in$tree)),
                 # Sa = matrix(rep(4, 7*7), nrow = 7), # NParam x NTree matrix
                 # Salpha = rep(2, 7) # vector of NParam
)

# Function to generate initials
init <- function() {
  list(A = rnorm(dat_list$NParam, 0, 10),
       B = rnorm(dat_list$NParam, 0, 10),
       tau = runif(1, 0, 1),
       sig.eps.sig = runif(1, 0, 1),
       sig.eps.lam = runif(1, 0, 1),
       deltaA = runif(dat_list$nlagA, 0, 1),
       deltaB = runif(dat_list$nlagB, 0, 1)
       # deltaC = runif(dat_list$nlagC, 0, 1)
  )
}

inits_list <- list(init(), init(), init())

# Alternative, start from saved state
load("scripts/model-pd-md/inits/saved_state.Rdata")

# updated_inits <- saved_state[[2]]
# for(i in 1:3){
#   updated_inits[[i]]$deltaA <- updated_inits[[i]]$deltaA[1:5]
#   updated_inits[[i]]$deltaB <- updated_inits[[i]]$deltaB[1:7]
#   # updated_inits[[i]]$deltaB <- c(updated_inits[[i]]$deltaB, runif(3, 0, 1))
# }

# Initialize model
jm <- jags.model("scripts/model-pd-md/model.jags",
                 data = dat_list,
                 inits = saved_state[[2]],
                 n.chains = 3)

update(jm, n.iter = 5000)

# Monitor
params <- c("deviance", "Dsum",
            "A", "B",
            "Astar", "Bstar",
            "tau", "sig.eps.sig", "sig.eps.lam",
            "wA", "wB", 
            "Sigs",
            "deltaA", "deltaB", 
            "Estar.sig", "Estar.lam"
)
coda.out <- coda.samples(jm,
                         variable.names = params,
                         n.iter = 3000,
                         n.thin = 50)

# save(coda.out, file = "scripts/model-pd-md/coda/coda.Rdata")
load(file = "scripts/model-pd-md/coda/coda.Rdata")

# Inspect chains visually
mcmcplot(coda.out, parms = c("deviance", "Dsum", 
                             "Astar", "Bstar", "Sigs",
                             "wA", "wB",
                             "Estar.sig", "Estar.lam"))
caterplot(coda.out, regex = "^Astar\\[", reorder = FALSE)
caterplot(coda.out, regex = "^Bstar\\[", reorder = FALSE)
caterplot(coda.out, parms = "wA", reorder = FALSE)
caterplot(coda.out, parms = "wB", reorder = FALSE)
caterplot(coda.out, parms = "Estar.sig", reorder = FALSE)
caterplot(coda.out, parms = "Estar.lam", reorder = FALSE)

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
  filter(grepl("^Astar", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^Bstar", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^wA", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^wB", rowname))

# Save state
# final <- initfind(coda.out, OpenBUGS = FALSE)
# final[[1]]
# saved_state <- removevars(final, variables = c(2, 4:8, 14:15))
# saved_state[[1]]
# 
# ind <- which(colnames(coda.out[[1]]) == "deviance")
# mean(coda.out[[1]][,ind])
# mean(coda.out[[2]][,ind])
# mean(coda.out[[3]][,ind])
# 
# if(!dir.exists("scripts/model-pd-md/inits")) {
#   dir.create("scripts/model-pd-md/inits")
# }
# save(saved_state, file = "scripts/model-pd-md/inits/saved_state.Rdata")


# Run model for replicated data, time series of sigma and lambda
coda.rep <- coda.samples(jm, 
                         variable.names = c("md.rep"),
                         n.iter = 3000,
                         n.thin = 50)

save(coda.rep, file = "scripts/model-pd-md/coda/codarep.Rdata")
# load(file = "scripts/model-pd-md/coda/codarep.Rdata")

# Summarize replicated output
coda_sum <- tidyMCMC(coda.rep,
                     conf.int = TRUE,
                     conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)


# Check model fit
pred <- cbind.data.frame(psy_in, coda_sum)

m1 <- lm(pred.mean ~ MD, data = pred)
summary(m1) # R2 = 0.8985; w/ RE R2 = 0.9198

pred %>%
  ggplot(aes(x = MD, y =pred.mean)) +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  geom_errorbar(aes(ymin = pred.lower, ymax = pred.upper),
                alpha = 0.25) +
  geom_point() +
  scale_x_continuous("Observed") +
  scale_y_continuous("Predicted") +
  theme_bw(base_size = 12) +
  facet_wrap(~Tree) +
  coord_equal()

