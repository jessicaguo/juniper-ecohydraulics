# Make figure for IOS-IEP proposal - between branch correlations
library(tidyverse)

load("data_cleaned/psy_daily.Rdata")
man <- read_csv("data_raw/Pressure-chamber-data.csv")

ggplot(man) +
  geom_point(aes(x = Timestamp, y = Psi,
                 color = as.factor(Tree),
                 shape = as.factor(Logger)))

psy_both <- psy_daily |> 
  pivot_wider(names_from = Logger, values_from = c(WP, WP.SD)) |> 
  select(-n, -doy) |> 
  filter(!is.na(sum(WP_1 +  WP_2)))

psy_both |> 
  filter(Tree %in% c(3,7),
         type == "PD",
         date >= as.Date("2021-07-01"),
         date <= as.Date("2021-09-15")) |> 
  ggplot(aes(x = date)) +
geom_point(aes(y = WP_1,
               color = "branch 1")) +
  geom_point(aes(y  = WP_2,
                 color = "branch 2")) +
  facet_wrap(~Tree) +
  scale_y_continuous(expression(paste(Psi[PD], " (MPa)"))) +
  scale_color_brewer(type = "qual", palette = 7) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.85, 0.25),
        legend.title = element_blank(),
        axis.title.x = element_blank())

ggsave("plots/psy_branch_partial.png",
       height = 2.5,
       width  = 5.5)
