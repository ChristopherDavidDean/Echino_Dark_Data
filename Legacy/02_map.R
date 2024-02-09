## *****************************************************************************
##
## Project: Plunder of dinosaur holotypes
##
## Purpose of script: 'Flight path map'
##
## Author: Emma M. Dunne, 2023
## Copyright (c) E. Dunne, 2023
## Email: emma.dunne@fau.de
##
## Date Created: 2023-08-29
## Last Modified:
##
## *****************************************************************************
##
## Notes:
##   
##
## *****************************************************************************


## Load packages:
library(tidyverse)
library(ggmap)


## Load capitals dataset:
holo_centroids <- read_csv("data/holo_centroids.csv")


## Pear down holotype data:
holo_foreign_trunc <- holo_foreign %>% 
  select(type_country, country_reposited) %>% # select origin andrepo countries
  rename(from = type_country, to = country_reposited)  # rename these columns

## Find the number of times countries match
holo_nodes <- holo_foreign_trunc %>% 
  group_by(from, to) %>% # group when from = to country
  tally() %>% # number of matches
  rename(weight = n) %>% # rename n column to be weight
  ungroup


## Redistribute weights for better plotting:
holo_nodes <- holo_nodes %>%
  mutate(weight2 = case_when(weight == 1 ~ 1,
                             weight == 2 | weight == 3 ~ 2,
                             weight >= 4 & weight <= 7 ~ 3,
                             weight >=11 & weight <= 13 ~ 4,
                             weight == 17 ~ 5,)) 


## Add coordinates from the capitals data:
holo_edges <- inner_join(holo_nodes, holo_centroids, by = c("from" = "country")) %>% 
  rename(x = lon, y = lat) # from country
holo_edges <- inner_join(holo_edges, holo_centroids, by = c("to" = "country")) %>% 
  rename(xend = lon, yend = lat) # to country


## Make plot 
worldmap <- borders("world", colour="#D8C9B2", fill="#D8C9B2")

ggplot(holo_centroids) + worldmap +
  geom_curve(data = holo_edges, aes(x = x, y = y, xend = xend, yend = yend, size = weight2),
             curvature = 0.33, alpha = 0.5, color = "#E07701") +
  scale_size_continuous(guide = FALSE, range = c(0.25, 2)) +
  geom_point(aes(x = lon, y = lat),
             shape = 21, size = 3, fill = 'white',
             color = 'black', stroke = 0.5) +
  geom_text(aes(x = lon, y = lat, label = country), 
            hjust = 0, nudge_x = 2, nudge_y = 0,
            size = 3, color = "black", fontface = "bold") +
  theme_minimal()

#ggsave("plots/flight_map.pdf", width = 17, height = 10)





# ## Add the coordinates back into the holo_foreign data:
# origin_countries <- holo_foreign %>% 
#   select(accepted_name, type_country) %>% 
#   rename(country = type_country) %>% 
#   left_join(holo_capitals, by = "country") %>% 
#   rename(origin = country)
# 
# repo_countries <- holo_foreign %>% 
#   select(accepted_name, country_reposited) %>% 
#   rename(country = country_reposited) %>% 
#   left_join(holo_capitals, by = "country") %>% 
#   rename(country_repo = country,
#          lat_r = lat, lon_r = lon) %>% 
#   select(country_repo, lat_r, lon_r)
# 
# ## Join these:
# holo_locs <- cbind(origin_countries, repo_countries)
# holo_locs <- distinct(holo_locs, accepted_name, .keep_all = TRUE)
# View(holo_locs)
# 
# 
# worldmap <- borders("world", colour="#CEB998", fill="#CEB998")
# 
# flight_map <- ggplot() + worldmap + 
#   geom_curve(data = holo_locs, aes(x = lon, y = lat, xend = lon_r, yend = lat_r), col = "#857457", size = .4) + 
#   geom_text_repel(data=holo_capitals, aes(x = lon, y = lat, label = country), col = "black", size = 2, segment.color = NA) +
#   geom_point(data = holo_capitals, aes(x = lon, y = lat), col = "#970027") + 
#   theme_void()
# flight_map
# 
# ggsave("plots/flight_map.pdf", width = 10, height = 6)