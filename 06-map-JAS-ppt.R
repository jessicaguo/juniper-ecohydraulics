# Calculate and make map of 4 corners area with %JAS precip
# Using PRISM 4 km normals for 1991-2020
# Insert lat/lons for DBG, Wasatch Range, and Bears Ears


# library(rgdal) # rgdal will be retired by end of 2023, transition to sf
library(maps)
library(sf)
library(dplyr)
library(ggplot2)
library(prism)
library(raster)
library(viridis)
library(ggthemes)

#### Add locs ####

loc <- data.frame(name = c("Wasatch Front", "Bears Ears", "Desert Botanical Garden", "Austin"),
                  x = c(-111.9175, -109.7471, -111.9411, -97.743056),
                  y = c(40.575, 37.5241, 33.46472, 30.267222))

#### Select states ####
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) %>%
  filter(ID %in% c("arizona", "new mexico", "colorado", "utah", "texas"))
head(states)
class(states)
projection(states)

ggplot(states) +
  geom_sf()

#### Calculate from PRISM ####
prism_set_dl_dir("data_raw/prism_ppt_1991_2020")
fn <- prism_archive_ls()
length(fn)

# Stack data frames and sum
RS_annual <- pd_stack(fn) # RasterStack object
RS_annual_sum <- calc(RS_annual, fun = sum)
dim(RS_annual_sum)

RS_monsoon <- pd_stack(fn[7:9])
RS_monsoon_sum <- calc(RS_monsoon, fun = sum)
class(RS_monsoon_sum)
dim(RS_monsoon_sum)

RL_monsoon_prop <- RS_monsoon_sum / RS_annual_sum * 100
RL_monsoon_prop
class(RL_monsoon_prop)
dim(RL_monsoon_prop)
projection(RL_monsoon_prop)
# clip by 4 corners states
RL_monsoon_prop_crop <- mask(RL_monsoon_prop, states) # vs. crop, which provides the outermost extent
dim(RL_monsoon_prop_crop)
class(RL_monsoon_prop_crop)
crs(RL_monsoon_prop_crop)

spdf_monsoon_prop <- as(RL_monsoon_prop_crop, "SpatialPixelsDataFrame")
class(spdf_monsoon_prop)
df_monsoon_prop <- as.data.frame(spdf_monsoon_prop) %>%
  rename(`% JAS` = layer)
class(df_monsoon_prop)

fig_base <- ggplot() +
  geom_tile(data = df_monsoon_prop,
              aes(x = x, y = y, fill = `% JAS`)) +
  geom_sf(data = states,
          fill = NA) +
  # geom_point(data = loc, aes(x = x, y = y)) +
  # geom_label(data = loc, aes(x = x, y = y, label = name), 
  #            size = 2.5, 
  #            label.size = NA,
  #            hjust = 0, 
  #            vjust = 1) +
  scale_fill_viridis() +
  theme_map() +
  theme(legend.position = "right",
        legend.key.width = unit(1, "cm"),
        legend.background = element_blank())
ggsave("plots/map_JAS.png",
       plot = fig_base,
       width = 6,
       height = 4, 
       units = "in")

fig_all <- ggplot() +
  geom_tile(data = df_monsoon_prop,
            aes(x = x, y = y, fill = `% JAS`)) +
  geom_sf(data = states,
          fill = NA) +
  geom_point(data = loc, aes(x = x, y = y)) +
  geom_label(data = loc, aes(x = x, y = y, label = name),
             size = 2.5,
             label.size = NA,
             alpha = 0.6,
             hjust = 0,
             vjust = 1) +
  scale_fill_viridis() +
  theme_map() +
  theme(legend.position = "right",
        legend.key.width = unit(1, "cm"),
        legend.background = element_blank())
ggsave("plots/map_labeled.png",
       plot = fig_all,
       width = 6,
       height = 4, 
       units = "in")
