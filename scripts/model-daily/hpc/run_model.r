#!/usr/bin/env Rscript

library(dplyr)
library(rjags)
load.module("dic")

# Load data
load("met_in.Rdata")
load("psy_in.Rdata")
load("branch_in.Rdata")
# Load saved_state
load("saved_state.Rdata")

# Center and scale env vars
Dmax <- scale(met_in$VPD_max)
VWC5 <- scale(met_in$VWC_5cm)
VWC50 <- scale(met_in$VWC_50cm)

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



# Initialize model
jm <- jags.model("modelb.jags",
                 data = dat_list,
                 inits = list(saved_state[[2]][[1]],
                              saved_state[[2]][[1]],
                              saved_state[[2]][[3]]),
                 n.chains = 3)

update(jm, n.iter = 100000)

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

save(coda.out, file = "coda.Rdata")


# Run model for replicated data
coda.rep <- coda.samples(jm, 
                         variable.names = "pd.rep",
                         n.iter = 1000,
                         n.thin = 20)

save(coda.rep, file = "codarep.Rdata")
