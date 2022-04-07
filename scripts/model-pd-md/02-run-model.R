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
                 nlagA = 7,
                 nlagB = 7,
                 # nlagC = 7,
                 # time step size (single value if modelb, vector of values if modeld),
                 pA = 1,
                 pB = 1,
                 # pC = 3,
                 # weights, of length nlag
                 alphaA = rep(1, 7), 
                 alphaB = rep(1, 7),
                 # alphaC = rep(1, 7),
                 doy = psy_in$doy,
                 Dmax = as.vector(met_in$Dmax),
                 VWC10 = as.vector(met_in$VWC10),
                 # VWC50 = as.vector(met_in$VWC50),
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
       deltaA = runif(dat_list$nlagA, 0, 1),
       deltaB = runif(dat_list$nlagB, 0, 1)
       # deltaC = runif(dat_list$nlagC, 0, 1)
  )
}

inits_list <- list(init(), init(), init())

# Alternative, start from saved state
load("scripts/model-pd-md/inits/saved_state.Rdata")

# Initialize model
jm <- jags.model("scripts/model-pd-md/model.jags",
                 data = dat_list,
                 inits = inits_list,
                 n.chains = 3)

update(jm, n.iter = 10000)

# Monitor
params <- c("deviance", "Dsum",
            "A", "B",
            "tau", "sig",
            "wA", "wB", #"wC",
            "deltaA", "deltaB"#, "deltaC"
)
coda.out <- coda.samples(jm,
                         variable.names = params,
                         n.iter = 3000,
                         n.thin = 5)

# save(coda.out, file = "scripts/model-pd-md/coda/coda.Rdata")

# Inspect chains visually
mcmcplot(coda.out, parms = c("deviance", "Dsum", 
                             "A", "B", "sig",
                             "wA", "wB"))
caterplot(coda.out, regex = "^A\\[", reorder = FALSE)
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
  filter(grepl("^A", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^B", rowname))

gel$psrf %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(grepl("^w", rowname))

# Save state
final <- initfind(coda.out, OpenBUGS = FALSE)
final[[1]]
saved_state <- removevars(final, variables = c(3, 6, 8:9))
saved_state[[1]]
ind <- which(colnames(coda.out[[1]]) == "deviance")

mean(coda.out[[1]][,ind])
mean(coda.out[[2]][,ind])
mean(coda.out[[3]][,ind])

if(!dir.exists("scripts/model-pd-md/inits")) {
  dir.create("scripts/model-pd-md/inits")
}
save(saved_state, file = "scripts/model-pd-md/inits/saved_state.Rdata")


# Run model for replicated data, time series of sigma and lambda
coda.rep <- coda.samples(jm, 
                         variable.names = c("md.rep", "sigma", "lambda"),
                         n.iter = 3000,
                         n.thin = 15)

# save(coda.rep, file = "scripts/model-pd-md/coda/codarep.Rdata")


# Summarize replicated output
coda_sum <- tidyMCMC(coda.rep,
                     conf.int = TRUE,
                     conf.method = "HPDinterval") %>%
  rename(pred.mean = estimate,
         pred.lower = conf.low,
         pred.upper = conf.high)

rep_sum <- coda_sum %>%
  filter(grepl("md.rep", term)) %>%
  rename(md.mean = pred.mean,
         md.lower = pred.lower,
         md.upper = pred.upper)

sig_sum <- coda_sum %>%
  filter(grepl("sigma", term)) %>%
  rename(sig.mean = pred.mean,
         sig.lower = pred.lower,
         sig.upper = pred.upper) %>%
  select(-1, -3)

lam_sum <- coda_sum %>%
  filter(grepl("lambda", term)) %>%
  rename(lam.mean = pred.mean,
         lam.lower = pred.lower,
         lam.upper = pred.upper) %>%
  select(-1, -3)

# Check model fit
pred <- cbind.data.frame(psy_in, rep_sum)

m1 <- lm(md.mean ~ MD, data = pred)
summary(m1) # R2 = 0.8985

pred %>%
  ggplot(aes(x = MD, y = md.mean)) +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  geom_errorbar(aes(ymin = md.lower, ymax = md.upper),
                alpha = 0.25) +
  geom_point() +
  scale_x_continuous("Observed") +
  scale_y_continuous("Predicted") +
  theme_bw(base_size = 12) +
  facet_wrap(~Tree) +
  coord_equal()


# Match and plot sigma
pred_params <- cbind.data.frame(psy_in, sig_sum, lam_sum)

fig_a <- ggplot(pred_params, aes(x = date,
                     y = sig.mean)) +
  geom_errorbar(aes(ymin = sig.lower,
                   ymax = sig.upper)) +
  geom_point(aes(color = Tree)) +
  guides(color = "none")

fig_b <- ggplot(pred_params, aes(x = date,
                        y = lam.mean)) +
  geom_errorbar(aes(ymin = lam.lower,
                    ymax = lam.upper)) +
  geom_point(aes(color = Tree)) +
  guides(color = "none")


fig_c <- filter(met_in, 
       date >= min(psy_in$date),
       date <= max(psy_in$date)) %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = VPD_max, color = "max D")) +
  geom_point(aes(y = VWC_10cm, color = "VWC 10 cm")) +
  guides(color = "none")


cowplot::plot_grid(fig_a, fig_b, fig_c, 
                   ncol = 1, align = "v")             
