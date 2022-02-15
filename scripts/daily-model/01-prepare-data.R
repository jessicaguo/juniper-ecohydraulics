# Prepare daily water potential and env data for JAGS model
# Identify separate branches per tree for a random effect

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

# Use 2021 data and doy as index
# Will center/scale in model script
met_in <- met_daily %>%
  filter(date >= as.Date("2021-01-01"))

psy_in <- psy %>%
  mutate(doy = lubridate::yday(date))

branch_in <- psy %>%
  group_by(branchID) %>%
  summarize(tree = unique(Tree))

save(met_in, file = "scripts/daily-model/met_in.Rdata")
save(psy_in, file = "scripts/daily-model/psy_in.Rdata")
save(branch_in, file = "scripts/daily-model/branch_in.Rdata")
