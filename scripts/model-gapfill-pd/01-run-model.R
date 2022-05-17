# Predict missing predawn values with gapfilling model
# BUGS model requires setting working directory where model file is located
setwd("~/Projects/hot-spots-psi/scripts/model-gapfill-pd")

library(dplyr)
library(ggplot2)
library(R2OpenBUGS)
library(mcmcplots)
library(postjags)

# Load data
load(file = "../../data_cleaned/psy_daily.Rdata")
min(psy_daily$doy)

# Restrict to PD, summarize by tree, widen
pd <- psy_daily %>%
  filter(type == "PD") %>%
  group_by(Tree, doy) %>%
  summarize(WP_mean = mean(WP, na.rm = TRUE)) %>%
  tidyr::pivot_wider(names_from = Tree,
                     names_prefix = "Tree_",
                     values_from = WP_mean) %>%
  mutate(mpd = rowMeans(.[,2:8], na.rm = TRUE),
         n = rowSums(!is.na(.[,2:8]))) %>%
  arrange(doy)

# Check visually
ggplot(pd, aes(x = doy)) +
  geom_point(aes(y = Tree_1, col = "1")) +
  geom_point(aes(y = Tree_2, col = "2")) +
  geom_point(aes(y = Tree_3, col = "3")) +
  geom_point(aes(y = Tree_4, col = "4")) +
  geom_point(aes(y = Tree_5, col = "5")) +
  geom_point(aes(y = Tree_6, col = "6")) +
  geom_point(aes(y = Tree_7, col = "7")) +
  geom_point(aes(y = mpd),
             col = "black",
             size = 2)

# Calculate sd ranges
apply(pd[,2:8], 2, FUN = sd, na.rm = TRUE)

# Prepare data
pd_mat <- as.matrix(pd[,2:8])
dim(pd_mat) # 157 by 7

dat_list <- list(pd = pd_mat,
                 mpd = array(dim = c(nrow(pd), 1), data = pd$mpd),
                 N = nrow(pd))

# Create initial values
init <- function(){
  list(sd_pd = runif(1, 0.5, 5),
       rho_pd = runif(1, 0, 1),
       tau.eps = runif(1, 0, 0.5))
}
inits_list <- list(init(), init(), init())

# Parameters to monitor
params <- c("Dsum","eps", "tau.eps", "sig.eps",
          "sd_pd", "rho_pd",
          "pd.rep")

# Compile model
bugsmod <- bugs(model.file = "model.txt", #model file
                data = dat_list,
                inits = inits_list,
                parameters.to.save = params, #parameters to monitor
                n.iter = 3000, n.chains = 3, n.burnin = 1000, n.thin = 1,
                codaPkg = TRUE, debug = FALSE,
                DIC = FALSE)		

# Convert bugs object  to coda object (when codaPkg = TRUE)
coda <- read.bugs(bugsmod)

# Visualize chains
mcmcplot(coda, parms = c("Dsum","eps", "tau.eps", "sig.eps",
                         "sd_pd", "rho_pd"))

# Save state
newinits <- initfind(coda, OpenBUGS = TRUE)
newinits[[1]]
saved_state <- removevars(initsin = newinits, variables=c(1,2,3,6))
#check both items in list
saved_state[[1]]
saved_state[[2]]#this goes into bugs() to reinitialize
save(saved_state, file = "inits/inits.R") #for local


#check convergence
gel <- gelman.diag(coda, multivariate = F)
str(gel)
gel$psrf[match("Dsum[1]", row.names(gel$psrf)):match("Dsum[12]", row.names(gel$psrf)),]
gel$psrf[match("eps[1]", row.names(gel$psrf)):match("eps[12]", row.names(gel$psrf)),]
gel$psrf[match("rho_pd", row.names(gel$psrf)),]
gel$psrf[match("sd_pd", row.names(gel$psrf)),]
gel$psrf[match("sig.eps", row.names(gel$psrf)),]

