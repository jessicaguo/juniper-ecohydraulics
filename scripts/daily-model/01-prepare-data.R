# Prepare daily water potential and env data for JAGS model
# Identify separate branches per tree for a random effect

library(readr)
library(dplyr)
library(ggplot2)

# Load daily summarized variables
load("data_cleaned/met_daily.Rdata")
load("data_cleaned/psy_daily.Rdata")

# Read in maintenance times
maint <- read_csv("data_raw/Maintenance-times.csv") %>%
  mutate(st_date = as.Date(Start),
         en_date = as.Date(End))

# Separate PD and MD columns
# Label branches by maintenance dates
psy_all <- psy_daily %>%
  select(Tree, Logger, date, type, WP) %>%
  tidyr::pivot_wider(names_from = type,
              values_from = WP) %>%
  left_join(met_daily) %>%
  mutate(period = case_when(date > maint$en_date[1] &
                              date < maint$st_date[2] ~ 1,
                            date > maint$en_date[2] &
                              date < maint$st_date[3] ~ 2,
                            date > maint$en_date[3] &
                              date < maint$st_date[4] ~ 3,
                            date > maint$en_date[4] &
                              date < maint$st_date[5] ~ 4,
                            date > maint$en_date[5] ~ 5))

# Render NA an outlier in Tree 7 period 1
# Where PD dropped below -3.5 on one day only
psy_all[psy_all$Tree == 7 &
        psy_all$PD < -3.5 &
        psy_all$period == 1, ] <- NA

# Check number of observations per period
psy_all %>%
  group_by(Tree, Logger, period) %>%
  summarize(npd = sum(!is.na(MD))) %>%
  print(n = 65)

# Remove tree/logger/period combo with only one observation
psy_all <- psy_all %>%
  filter(!(Tree == 4 & Logger == 1 & period == 3))

# Create association between period within logger
# and branch within tree
branch_id <- psy_all %>%
  group_by(Tree, Logger, period) %>%
  count() %>%
  group_by(Logger, period) %>%
  mutate(branch_fixed = cur_group_id()) %>%
  select(-n) %>%
  group_by(Tree) %>%
  mutate(branchID = cur_group_rows())

# Join branch_fixed to psy_all
psy <- psy_all %>%
  left_join(branch_id)

###### Prepare data for model #######

# Scale/center across all available data
# Use 2021 data and doy as index
# Will center/scale in model script
met_in <- met_daily %>%
  mutate(Dmax = scale(VPD_max),
         VWC5 = scale(VWC_5cm),
         VWC10 = scale(VWC_10cm),
         VWC50 = scale(VWC_50cm)) %>%
  filter(date >= as.Date("2021-01-01"))

psy_in <- psy %>%
  mutate(doy = lubridate::yday(date))

branch_in <- psy %>%
  group_by(branchID) %>%
  summarize(tree = unique(Tree))

save(met_in, file = "scripts/daily-model/met_in.Rdata")
save(psy_in, file = "scripts/daily-model/psy_in.Rdata")
save(branch_in, file = "scripts/daily-model/branch_in.Rdata")
