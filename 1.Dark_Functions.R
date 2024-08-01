################################################################################
### MUSEUM `DARK DATA` ILLUMINATES THE EVOLUTIONARY HISTORY OF FOSSIL GROUPS ###
################################################################################

# Christopher D. Dean & Jeffrey R. Thompson
# 2024
# Script written by Christopher D. Dean

################################################################################
#                             FILE 1: FUNCTIONS                                #
################################################################################

################################################################################
# 0. INDEX AND BASIC SETUP
################################################################################

################################################################################

# 1. SETUP
# --- simple.grain
# --- simple.rank
#
# 2. BASIC COMPARISONS
# --- setup.table 
# --- table.plots
# --- country.clean
#
# 3. SPATIAL OCCUPANCY THROUGH TIME
# --- get_grid
# --- get_grid_im
#
# 4. SPATIAL RANGES THROUGH TIME
# --- make.range
# --- lat.range.fun
# --- space.results
# --- org.funct
# --- spatial.tests
# --- geo.range.fun
#
# 5. TEMPORAL RANGE COMPARISON
# --- temp_compare
# --- perc.range
#
# 6. DIVERSITY ANALYSIS
# --- run.div
# --- run.col
# --- setup.div
# --- freq.data
# --- run.iNEXT
#
# 7. SEDIMENT AFFINITY
# --- make.fam.lith
# --- setup.SRA
# --- test.func
# --- plot.SRA
# --- shift_legend
#
# 8. SCIENTIFIC COLONIALISM
# --- setup.col.sec
# --- network.plots
# --- reshape.df
# --- research.prep
# --- continent.stats
# --- top.countries
# --- network.stats
#
# 9. MODELLING COMPARISONS
# --- run.models
# --- ca.binning
# --- model.results
#
################################################################################

# Load Packages
library(tidyr)
library(tibble)
library(raster)
library(latticeExtra)
library(rasterVis)
library(sp)
library(stringr)
library(patchwork)
library(sf)
library(maps)
library(mapview)
library(mapdata)
library(ggplot2)
library(palaeoverse)
library(countrycode)
library(geoscale)
library(vegan)
library(RColorBrewer)
library(divDyn)
library(chronosphere)
library(devtools)
library(maditr)
library(MASS)
library(effects)
library(deeptime)
library(wesanderson)
library(viridis)
library(sjPlot)
library(rmacrostrat)

################################################################################
# 1. SETUP
################################################################################
##### simple.grain #####

# Make grainsize into fine and coarse grained
simple.grain <- function(data.set){
  for (l in 1:nrow(data.set)){
    if(is.na(data.set$Finalised_grainsize[l]) == T){
      data.set$Finalised_grainsize[l] <- ""
    }
    if (data.set$Finalised_grainsize[l]=="Fine Grained"  | data.set$Finalised_grainsize[l] == "Slate" |
        data.set$Finalised_grainsize[l]=="Chert" | data.set$Finalised_grainsize[l] =="Mudstone/Grainstone" |
        data.set$Finalised_grainsize[l]=="Shale" | data.set$Finalised_grainsize[l] == "Wackestone" | data.set$Finalised_grainsize[l] == "Mudstone" | 
        data.set$Finalised_grainsize[l]=="Siltstone" | data.set$Finalised_grainsize[l]== "Mudstone/Wackestone" |
        data.set$Finalised_grainsize[l]=="Micrite" | data.set$Finalised_grainsize[l]== "Mudstone/Siltstone" |
        data.set$Finalised_grainsize[l]=="Claystone" | data.set$Finalised_grainsize[l]=="Floatstone"){
      data.set$Finalised_grainsize_simp[l]<-"Fine Grained"
    }
    else if (data.set$Finalised_grainsize[l]=="Grainstone" | data.set$Finalised_grainsize[l] == "Packstone"| 
             data.set$Finalised_grainsize[l]=="Sandstone" | data.set$Finalised_grainsize[l]=="Reefal" |
             data.set$Finalised_grainsize[l] == "Coquina" | data.set$Finalised_grainsize[l] == "Bindstone" |
             data.set$Finalised_grainsize[l] == "Boundstone"){ 
      data.set$Finalised_grainsize_simp[l]<-"Coarse Grained"
    }
    else{
      data.set$Finalised_grainsize_simp[l] <- NA
    }
  }
  data.set <- data.set
}

##### simple.rank #####

# Make anything including or above 'class' into a single category
simple.rank <- function(data.set){
  for (l in 1:nrow(data.set)){
    if(is.na(data.set$Rank[l]) == T){
      data.set$Rank[l] <- ""
    }
    if (data.set$Rank[l]=="Class"  | data.set$Rank[l]=="Order" |
        data.set$Rank[l]=="Infraorder"){
      data.set$Rank_simp[l]<-"Above Family"
    }
    else if (data.set$Rank[l]=="Family"){ 
      data.set$Rank_simp[l]<-"Family"
    }
    else if (data.set$Rank[l]=="Genus"){ 
      data.set$Rank_simp[l]<-"Genus"
    }
    else if (data.set$Rank[l]=="Species"){ 
      data.set$Rank_simp[l]<-"Species"
    }
    else{
      data.set$Rank_simp[l] <- NA
    }
  }
  data.set <- data.set
}

################################################################################
# 2. BASIC COMPARISONS
################################################################################
##### setup.table #####
# Make table of tallied results to plot, taking specific column.
setup.table <- function(mdata, 
                        pbdata, 
                        pubdata, 
                        alldata, 
                        column, 
                        pivot = T, 
                        useNA = F){
  if(useNA == T){
    temp <- list(as.data.frame(table(mdata[[column]], useNA = "ifany")),
                 as.data.frame(table(pbdata[[column]], useNA = "ifany")),
                 as.data.frame(table(pubdata[[column]], useNA = "ifany")),
                 as.data.frame(table(alldata[[column]], useNA = "ifany")))
  }else{
    temp <- list(as.data.frame(table(mdata[[column]])),
                 as.data.frame(table(pbdata[[column]])),
                 as.data.frame(table(pubdata[[column]])),
                 as.data.frame(table(alldata[[column]])))
  }
  temp <- temp %>%
    purrr::reduce(full_join, by = "Var1") %>%
    dplyr::rename(temp.name = Var1,
                  `Dark Data` = Freq.x,
                  PBDB = Freq.y,
                  Published = Freq.x.x,
                  All = Freq.y.y)
  temp[is.na(temp)] <- 0
  names(temp)[names(temp) == "temp.name"] <- column
  if(pivot == T){
    temp <- tidyr::pivot_longer(temp, cols = c("Dark Data", "PBDB", 
                                               "Published", "All"))
    temp$value <- as.numeric(temp$value)
  }
  assign(deparse(substitute(column)), temp, envir = .GlobalEnv)
}

##### table.plots #####
# Plots results from setup.table, as both count and proportion
table.plots <- function(pivot, x, removeNA = T, colour, labs){
  if(removeNA == T){
    pivot <- na.omit(pivot)
  }
  a <- ggplot(data = pivot, aes(x = name, y = value, fill = !!sym(x))) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(values=colour) +
    guides(fill=guide_legend(title=x)) +
    ylab("Total occurrences") +
    xlab("Dataset") +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "white"), 
          panel.grid.minor = element_line(colour = "white")) +
    guides(alpha = "none")
  b <- ggplot(data = pivot, aes(x = name, y = value, fill = !!sym(x))) +
    geom_bar(stat = 'identity', position = 'fill') +
    scale_fill_manual(values=colour) +
    guides(fill=guide_legend(title=x)) +
    ylab("Proportion of dataset") +
    xlab("Dataset") +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "white"), 
          panel.grid.minor = element_line(colour = "white")) +
    guides(alpha = "none")
  
  ggarrange(a, b, 
            align='hv',
            labels = labs,
            nrow = 1, 
            ncol = 2,
            common.legend = T, 
            legend = "bottom")
}

##### country.clean #####

country.clean <- function(df){
  ##### Setup data for analysis #####
  df[df == "England" | df == "Wales" | df == "Scotland" | 
       df == "Northern Ireland"] <- "United Kingdom"
  df[df == "USA"] <- "United States"
  df[df == "Czechia"] <- "Czech Republic"
  df$Country <- gsub("&", "and", df$Country)
  
  df$aff_code <- countrycode(df$Museum_country, "country.name", "iso3c")
  df$samp_code <- countrycode(df$Country, "country.name", "iso3c")
  return(df)
}

################################################################################
# 3. SPATIAL OCCUPANCY THROUGH TIME
################################################################################

##### get_grid #####

# Creates a raster of chosen resolution, and attaches associated grid cell IDs to occurrences/collections
get_grid <- function(data, res, e, r = "N"){ # data is first output from combine_data (fossil.colls). Res is chosen resolution in degrees
  if (class(r) == "character"){
    r <- raster(res = res, val = 1, ext = e) # Value must be added because extract uses values
    r <<- r
  }
  crs <- r@crs
  xy <- SpatialPointsDataFrame(cbind.data.frame(data$lng, data$lat), data, proj4string = crs)
  Final <- raster::extract(r, xy, sp = TRUE, cellnumbers = TRUE)
  Final <- as.data.frame(Final)
  singletons <- Final %>%
    dplyr::select(cells, Locality) %>%
    dplyr::distinct() %>%
    dplyr::group_by(cells) %>%
    dplyr::summarize(Coll_count = n())
  singletons$Coll_count[singletons$Coll_count == 1] <- "Singleton"
  singletons$Coll_count[singletons$Coll_count != "Singleton"] <- "Non-singleton"
  Final <- left_join(Final, singletons, by = "cells")
  stack_name <- paste(deparse(substitute(data)), ".spat.bin", sep = "")
  assign(stack_name, Final, envir = .GlobalEnv)
}

##### get_grid_im #####
# Set up background info
#countries <- maps::map("world", plot=FALSE, fill = TRUE) # find map to use as backdrop
#states <- maps::map("state", plot = FALSE, fill = TRUE)
#countries <<- maptools::map2SpatialPolygons(countries, IDs = countries$names, proj4string = CRS("+proj=longlat")) # Turn map into spatialpolygons
#states <<- maptools::map2SpatialPolygons(states, IDs = states$names, proj4string = CRS("+proj=longlat")) # Turn map into spatialpolygons
#
## Function for making maps
#get_grid_im <- function(data, res, name, ext){ # Data is first output from combine_data (fossil.colls). Res is chosen resolution in degrees. name is user inputted string related to data inputted, for display on graphs. 
#  xy <- cbind(as.double(data$lng), as.double(data$lat))
#  #xy <- unique(xy)
#  r <- raster::raster(ext = ext, res = res)
#  r <- raster::rasterize(xy, r, fun = 'count')
#  #r[r > 0] <- 1 # Remove if you want values instead of pure presence/absence.
#  countries <- maps::map("world", plot=FALSE, fill = TRUE) # find map to use as backdrop
#  countries <<- maptools::map2SpatialPolygons(countries, IDs = countries$names, proj4string = CRS("+proj=longlat")) # Turn map into spatialpolygons
#  mapTheme <- rasterVis::rasterTheme(region=viridis(8))
#  (r <- rasterVis::levelplot(r, margin=F, par.settings=mapTheme,  main = paste("Total ", (substitute(name)), " per Grid Cell", sep = "")) + #create levelplot for raster
#          #   latticeExtra::layer(sp.polygons(states, col = "white", fill = NA), under = T)  + # Plots state lines
#          latticeExtra::layer(sp.polygons(countries, col = 0, fill = "light grey"), under = T)) # Plots background colour
#}

################################################################################
# 4. SPATIAL RANGES THROUGH TIME
################################################################################
##### make.range #####
# Function for finding ranges of organisms
make.range <- function(occdf, name, method, rank, coords = T){
  if(rank == "Genus"){
    temp.data <- occdf %>%
      dplyr::filter(Rank == "Genus" | Rank == "Species") %>%
      dplyr::filter(is.na(p_lng) == F)
    temp.data <- group_apply(temp.data, group = c("bin_assignment"), 
                             fun = tax_range_space, 
                             name = "Genus", 
                             lng = "p_lng", 
                             lat = "p_lat", 
                             method = method, 
                             coords = coords)
  }
  if(rank == "Family"){
    temp.data <- occdf %>%
      dplyr::filter(Rank == "Family" | Rank == "Genus" | Rank == "Species") %>%
      dplyr::filter(is.na(Family) == F) %>%
      dplyr::filter(is.na(p_lng) == F)
    temp.data <- group_apply(temp.data, group = c("bin_assignment"), 
                             fun = tax_range_space, 
                             name = "Family", 
                             lng = "p_lng", 
                             lat = "p_lat", 
                             method = method, 
                             coords = coords)
  }
  if(method == "con"){
    temp.data <- temp.data[temp.data$area != 0, ]
  }
  if(method == "lat"){
    temp.data <- temp.data[temp.data$range_lat != 0, ]
  }
  temp.data$Data <- name
  assign(paste(deparse(substitute(occdf)), ".range.", method, sep = ""), temp.data, 
         envir = .GlobalEnv)
}

##### lat.range.fun #####
# Function for finding latitudinal ranges
lat.range.fun <- function(pbddata, 
                          musdata, 
                          pubdata, 
                          alldata, 
                          name, rank = "Genus"){
  # Combine ranges from all datasets
  lat.ranges <- rbind(make.range(pbddata, "PBDB", "lat", rank), 
                      make.range(musdata, "Museum", "lat", rank), 
                      make.range(pubdata, "Published", "lat", rank), 
                      make.range(alldata, "All", "lat", rank))
  # lat.ranges$bin_assignment <- as.numeric(lat.ranges$bin_assignment)
  
  assign(paste("lat.ranges.", name, ".", rank, sep = ""), lat.ranges, envir = .GlobalEnv)
  
  lat.results <- data.frame()
  for(i in sort(unique(lat.ranges$bin_assignment))){
    # Get data for bin, and remove any taxa that can't be compared
    temp.data <- lat.ranges %>%
      dplyr::filter(bin_assignment == i) %>%
      dplyr::group_by(taxon) %>%
      dplyr::filter(n() >2)
    # If no taxa remain, skip to next bin
    if(nrow(temp.data) < 1){
      next
    }
    # Get difference between Total range and specific range for each taxon
    temp.data <- temp.data %>%
      left_join(dplyr::filter(., Data == "All"), by = "taxon") %>%
      mutate(Difference = range_lat.y- range_lat.x, 
             Perc = (range_lat.x/range_lat.y)*100, 
             Mean_shift = ((max_lat.y + min_lat.y)/2) - ((max_lat.x + min_lat.x)/2)) %>%
      dplyr::select(taxon, Data.x, Difference, Perc, Mean_shift) %>%
      dplyr::filter(Data.x != "All") %>%
      dplyr::rename(Taxon = taxon, Data = Data.x) %>%
      as.data.frame() %>%
      arrange(Taxon) %>%
      mutate(bin = i)
    # Save results
    lat.results <- rbind(lat.results, temp.data)
  }
  lat.results$Data[lat.results$Data == "Museum"] <- "Dark Data"
  return(lat.results)
}

##### space.results #####
# Function for quickly summarising differences between datasets
space.results <- function(results){
  results %>%
    group_by(Data) %>%
    dplyr::summarize(Mean.diff = mean(Difference), 
                     Median.diff = median(Difference), 
                     Mean.perc = mean(Perc), 
                     Median.perc = median(Perc))
}

##### org.funct #####
org.funct <- function(data, var){
  test.2 <- data %>%
    dplyr::group_by(!!sym(var), taxon) %>%
    dplyr::summarize(count= n()) %>%
    st_cast("POLYGON") %>%
    st_convex_hull() %>%
    na.omit(geometry)
  if(any(st_is_valid(test.2) == F)){
    test.2 <- test.2[-which(st_is_valid(test.2) == F),]
  }
  test.2 <- test.2 %>%
    dplyr::mutate(Ex_Area = st_area(geometry)) %>%
    as.data.frame() %>%
    dplyr::select(taxon, Ex_Area)
  return(test.2)
}

##### spatial.tests #####
# Function that generates results of a variety of statistical tests. 
spatial.tests <- function(geo.results, focal){
  kruskal.test <- geo.results %>%
    kruskal_test(as.formula(paste(focal, " ~ Data", sep = "")))
  kruskal.effect <- geo.results %>%
    kruskal_effsize(as.formula(paste(focal, " ~ Data", sep = "")))
  dunn.test <- geo.results %>%
    dunn_test(as.formula(paste(focal, " ~ Data", sep = "")), 
              p.adjust.method = "bonferroni")
  wilcox.test <- geo.results %>%
    wilcox_test(as.formula(paste(focal, " ~ Data", sep = "")),
                p.adjust.method = "bonferroni")
  
  return(list(kruskal.test = kruskal.test, 
              kruskal.effect = kruskal.effect,
              dunn.test = dunn.test, 
              wilcox.test = wilcox.test))
}

##### geo.range.fun #####
# Calculates results for geographic area for taxa through time
geo.range.fun <- function(pbddata, 
                          musdata, 
                          pubdata, 
                          alldata, 
                          name, 
                          rank = "Genus"){
  # Combine ranges from all datasets
  geo.ranges <- rbind(make.range(pbddata, "PBDB", "con", rank, coords =  T), 
                      make.range(musdata, "Museum", "con", rank, coords =  T), 
                      make.range(pubdata, "Published", "con", rank, coords =  T),
                      make.range(alldata, "All", "con", rank, coords =  T))
  
  assign(paste("geo.ranges.", name, sep = ""), geo.ranges, envir = .GlobalEnv)
  
  geo.results <- data.frame()
  for(i in sort(unique(geo.ranges$bin_assignment))){
    # Get data for bin, and remove any taxa that can't be compared
    temp.data <- geo.ranges %>%
      dplyr::filter(bin_assignment == i) %>%
      group_by(taxon) %>%
      dplyr::filter(n_distinct(Data) > 2)
    
    # If no taxa remain, skip to next bin
    if(nrow(temp.data) < 1){
      next
    }
    temp_xy <- st_as_sf(temp.data, coords = c("p_lng", "p_lat"), crs = 4326)
    temp_xy <- st_transform(temp_xy, crs = c("+proj=eqearth"))
    polys <- temp_xy %>%
      dplyr::group_by(Data, taxon) %>%
      dplyr::summarize(count= n()) %>%
      st_cast("POLYGON") %>%
      st_convex_hull()
    
    # Setup biggest contribution
    test <- temp_xy %>%
      mutate(PBDB = case_when(Data == "Museum" ~ "PBDB",
                              Data == "Published" ~ "PBDB"), 
             Museum = case_when(Data == "PBDB" ~ "Museum",
                                Data == "Published" ~ "Museum"),
             Published = case_when(Data == "Museum" ~ "Published",
                                   Data == "PBDB" ~ "Published")
      )
    
    mus.test <- org.funct(test, "Museum")
    names(mus.test)[names(mus.test) == "Ex_Area"] <- "Museum_Ex_Area"
    PBD.test <- org.funct(test, "PBDB")
    names(PBD.test)[names(PBD.test) == "Ex_Area"] <- "PBDB_Ex_Area"
    pub.test <- org.funct(test, "Published")
    names(pub.test)[names(pub.test) == "Ex_Area"] <- "Published_Ex_Area"
    
    comb.test <- full_join(mus.test, pub.test, by = "taxon")
    comb.test <- full_join(comb.test, PBD.test, by = "taxon")
    
    print(comb.test)
    
    if(any(st_is_valid(polys) == F)){
      polys <- polys[-which(st_is_valid(polys) == F),]
      polys <- polys %>%
        group_by(taxon) %>%
        dplyr::filter(n_distinct(Data) > 2)
    }
    
    # Getting full area
    combined_polys <- polys %>%
      dplyr::group_by(taxon) %>%
      dplyr::summarise() %>%
      st_convex_hull() %>%
      dplyr::mutate(area = st_area(geometry))
    
    # Combining
    polys$area <- st_area(polys)
    temp.data <- polys %>%
      as.data.frame() %>%
      dplyr::full_join(combined_polys, by = "taxon") %>%
      mutate(Area = area.x,
             Difference = area.y - area.x, 
             Perc = (area.x/area.y)*100) %>%
      dplyr::select(taxon, Area, Data, Difference, Perc, count) %>%
      dplyr::filter(Data != "All") %>%
      dplyr::rename(Taxon = taxon) %>%
      as.data.frame() %>%
      arrange(Taxon) %>%
      mutate(bin = i)
    
    # Added area
    added.area <- combined_polys %>%
      as.data.frame() %>%
      dplyr::select(-geometry) %>%
      dplyr::full_join(comb.test, by = "taxon") %>%
      dplyr::mutate(Museum = area - Museum_Ex_Area,
                    Published = area - Published_Ex_Area, 
                    PBDB = area - PBDB_Ex_Area) %>%
      dplyr::select(taxon, Museum, PBDB, Published, area) %>%
      tidyr::pivot_longer(cols = c(Museum, PBDB, Published)) %>%
      dplyr::mutate(Percentage.added = (value/area)*100) %>%
      dplyr::rename(Taxon = taxon, 
                    Data = name,
                    Added.area = value)
    print(added.area)
    
    temp.data$Area <- as.numeric(temp.data$Area)
    temp.data$Difference <- as.numeric(temp.data$Difference)
    temp.data$Perc <- as.numeric(temp.data$Perc)
    temp.data <- left_join(temp.data, added.area, by = c("Taxon", "Data"))
    temp.data$Added.area <- as.numeric(temp.data$Added.area)
    temp.data$Added.area[temp.data$Added.area < 0] <- 0
    geo.results <- rbind(geo.results, temp.data)
    print(i)
  }
  geo.results$Data[geo.results$Data == "Museum"] <- "Dark Data"
  return(geo.results)
}

################################################################################
# 5. TEMPORAL RANGE COMPARISON
################################################################################
##### temp_compare #####
# Function
temp_compare <- function(occdf, bin_filter){
  ranks <- c("Species", "Genus", "Family")
  comp <- data.frame()
  for(i in ranks){
    if(i == "Species"){
      test.all.dat <- all.dat %>%
        dplyr::filter(Rank == "Species") %>%
        dplyr::filter(n_bins < bin_filter)
      test.occdf <- occdf %>%
        dplyr::filter(Rank == "Species") %>%
        dplyr::filter(n_bins < bin_filter)
      a <- tax_range_time(test.all.dat, name = "Combined_name", 
                          by = "name", plot = F)
      b <- tax_range_time(test.occdf, name = "Combined_name", 
                          by = "name", plot = F)
    }
    if(i == "Genus"){
      test.all.dat <- all.dat %>%
        dplyr::filter(Rank == "Genus" | Rank == "Species") %>%
        dplyr::filter(is.na(Genus) == F) %>%
        dplyr::filter(n_bins < bin_filter)
      test.occdf <- occdf %>%
        dplyr::filter(Rank == "Genus" | Rank == "Species") %>%
        dplyr::filter(n_bins < bin_filter)
      a <- tax_range_time(test.all.dat, name = "Genus", by = "name", plot = F)
      b <- tax_range_time(test.occdf, name = "Genus", by = "name", plot = F)
    }
    if(i == "Family"){
      test.all.dat <- all.dat %>%
        dplyr::filter(Rank == "Family" | Rank == "Genus" | Rank == "Species") %>%
        dplyr::filter(is.na(Family) == F) %>%
        dplyr::filter(n_bins < bin_filter) %>%
        dplyr::filter(Family != "NO_FAMILY_SPECIFIED")
      test.occdf <- occdf %>%
        dplyr::filter(Rank == "Family" | Rank == "Genus" | Rank == "Species") %>%
        dplyr::filter(is.na(Family) == F) %>%
        dplyr::filter(n_bins < bin_filter) %>%
        dplyr::filter(Family != "NO_FAMILY_SPECIFIED")
      a <- tax_range_time(test.all.dat, name = "Family", by = "name", plot = F)
      b <- tax_range_time(test.occdf, name = "Family", by = "name", plot = F)
    }
    temp_comp <- inner_join(a, b, by = "taxon")
    temp_comp$range_diff <- temp_comp$range_myr.x - temp_comp$range_myr.y
    temp_comp$bottom_diff <- temp_comp$max_ma.x - temp_comp$max_ma.y
    temp_comp$top_diff <- temp_comp$min_ma.x - temp_comp$min_ma.y
    temp_comp$Group <- i
    temp_comp$Range_perc <- (temp_comp$range_myr.y/temp_comp$range_myr.x)*100
    comp <- rbind(comp, temp_comp)
  }
  comp$dataframe <- deparse(substitute(occdf))
  assign(paste("comp.", deparse(substitute(occdf)), sep = ""), comp, 
         envir = .GlobalEnv)
  comp_stats <- comp %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(Mean_range_diff = mean(range_diff), 
                     Mean_top_diff = mean(top_diff),
                     Mean_bot_diff = mean(bottom_diff),
                     Correct_perc = sum(range_diff == 0)/n(), 
                     Range_perc = mean(Range_perc),
                     Count = n()
    )
  comp_stats$dataframe <- deparse(substitute(occdf))
  assign(paste("comp.stats.", deparse(substitute(occdf)), sep = ""), comp_stats, 
         envir = .GlobalEnv)
}

##### perc.range #####
# Function to find the position of occurrence within the temporal range of the genera
perc.range <- function(occdf, label, rank = "Genus"){
  # Setup occurrence dataframe
  if(rank == "Species"){
    occdf <- occdf %>%
      dplyr::filter(Rank == "Species")
    # Setup reference dataframe (all occurrences)
    all.occdf <- all.dat %>%
      dplyr::filter(Rank == "Species")
  }else if(rank == "Genus"){
    occdf <- occdf %>%
      dplyr::filter(Rank == "Genus" | Rank == "Species")
    # Setup reference dataframe (all occurrences)
    all.occdf <- all.dat %>%
      dplyr::filter(Rank == "Genus" | Rank == "Species")
  }else if (rank == "Family"){
    occdf <- occdf %>%
      dplyr::filter(Rank == "Genus" | Rank == "Species" | Rank == "Family") %>%
      dplyr::filter(Family != "NO_FAMILY_SPECIFIED")
    # Setup reference dataframe (all occurrences)
    all.occdf <- all.dat %>%
      dplyr::filter(Rank == "Genus" | Rank == "Species" | Rank == "Family") %>%
      dplyr::filter(is.na(Family) != T) %>%
      dplyr::filter(Family != "NO_FAMILY_SPECIFIED")
  }
  # Get temporal range of genera using all data as a reference
  a <- tax_range_time(all.occdf, name = rank, by = "name", plot = F)
  # For each row of the supplied dataframe
  occdf$Position <- unlist(lapply(seq(1,nrow(occdf),1), function(x){
    # Find the minimum age of the range for the Genera
    min <- a$min_ma[[which(a$taxon == occdf[[rank]][x])]]
    # Find the maximum age of the range for the Genera
    max <- a$max_ma[[which(a$taxon == occdf[[rank]][x])]]
    # Get the midpoint age of the occurrence
    input <- ((occdf$max_ma[x] + occdf$min_ma[x])/2)
    # Calculate the position of the occurrence within the age range
    return(((input - min) * 100)/(max-min))
  }))
  # Find singleton taxa (no range)
  singletons <- a %>%
    filter(n_occ == 1)
  # Remove singleton taxa
  occdf <- occdf[!occdf[[rank]] %in% singletons$taxon,]
  # Label data
  occdf$Data <- label
  occdf$Rank <- rank
  # Remove redundant data
  occdf <- occdf %>%
    dplyr::select(Data, Position, Rank)
  return(occdf)
}

################################################################################
# 6. DIVERSITY ANALYSIS
################################################################################
##### run.div #####
# Run divDyn, either for raw results for SQS, on a chosen dataset at chosen taxonomic rank.
run.div <- function(dataset, Rank = "Genus", type = "raw"){
  if(Rank == "Genus"){
    temp.data <- dataset %>%
      dplyr::filter(Rank == "Genus" | Rank == "Species")
  }
  if(Rank == "Species"){
    temp.data <- dataset %>%
      dplyr::filter(Rank == "Species")
  }
  if(type == "raw"){
    if(Rank == "Genus"){
      data.div <- divDyn::divDyn(x = temp.data, bin = "bin_assignment", 
                                 tax = "Genus", noNAStart = T)
    }else{
      data.div <- divDyn::divDyn(x = temp.data, bin = "bin_assignment", 
                                 tax = "Combined_name", noNAStart = T)
    }
  }
  if(type == "SQS"){
    quorum_levels <- round(seq(from = 0.3, to = 0.6, by = 0.1), 1)
    data.div <- list()
    for(i in 1:length(quorum_levels)){
      if(Rank == "Genus"){
        data.div.temp <- divDyn::subsample(x = temp.data, 
                                           bin = "bin_assignment", 
                                           tax = "Genus", 
                                           type = "sqs", 
                                           q = quorum_levels[i],
                                           iter = 500, 
                                           noNAStart = T)
        data.div.temp$quorum_level <- quorum_levels[i]
      }else{
        data.div.temp <- divDyn::subsample(x = temp.data, 
                                           bin = "bin_assignment", 
                                           tax = "Combined_name", 
                                           type = "sqs", 
                                           q = quorum_levels[i],
                                           iter = 500, 
                                           noNAStart = T)
        data.div.temp$quorum_level <- quorum_levels[i]
      }
      data.div[[i]] <- data.div.temp
    }
  }
  return(data.div)
}

##### run.col #####
run.col <- function(dataset, Rank = "Genus"){
  if(Rank == "Genus"){
    temp.data <- dataset %>%
      dplyr::filter(Rank == "Genus" | Rank == "Species")
    data.col <- divDyn::binstat(x = temp.data, bin = "bin_assignment",
                                coll = "Locality", 
                                tax = "Genus", noNAStart = T)
  }
  if(Rank == "Species"){
    temp.data <- dataset %>%
      dplyr::filter(Rank == "Species")
    data.col <- divDyn::binstat(x = temp.data, bin = "bin_assignment", 
                                coll = "Locality", 
                                tax = "Combined_name", noNAStart = T)
  }
  return(data.col)
}

##### setup.div #####
# Combines outputs from divDyn to provide a dataframe of diversity values. 
setup.div <- function(pbddata,
                      musdata,
                      pubdata,
                      alldata){
  # Only select relevant info (SIB diversity)
  m.only.div <- musdata %>%
    dplyr::select(bin_assignment, divSIB)
  pb.only.div <- pbddata %>%
    dplyr::select(bin_assignment, divSIB)
  pub.all.div <- pubdata %>%
    dplyr::select(bin_assignment, divSIB)
  all.div <- alldata %>%
    dplyr::select(bin_assignment, divSIB)
  
  # Set up dataset for merging, rename bin column
  div.tot <- stages
  names(div.tot)[names(div.tot) == "bin"]  <- "bin_assignment"
  
  # Merge data
  div.tot <- left_join(div.tot, m.only.div, by = "bin_assignment")
  names(div.tot)[names(div.tot) == "divSIB"] <- "Museum_SIB"
  div.tot <- left_join(div.tot, pb.only.div, by = "bin_assignment")
  names(div.tot)[names(div.tot) == "divSIB"] <- "PBDB_SIB"
  div.tot <- left_join(div.tot, pub.all.div, by = "bin_assignment")
  names(div.tot)[names(div.tot) == "divSIB"] <- "Pub_SIB"
  div.tot <- left_join(div.tot, all.div, by = "bin_assignment")
  names(div.tot)[names(div.tot) == "divSIB"] <- "All_SIB"
  
  # Remove rows with NAs, then remove unneccessary columns
  div.tot <- div.tot[rowSums(is.na(div.tot[,c("Museum_SIB", "PBDB_SIB", 
                                              "Pub_SIB",  "All_SIB")]))!=4,]
  div.tot <- div.tot %>%
    dplyr::mutate(bin_midpoint = (max_ma + min_ma)/2) %>%
    dplyr::select(-c(max_ma, min_ma, abbr, color, name))
  return(div.tot)
}

##### freq.data #####
# Function to generate frequency data for use in iNEXT
freq.data <- function(dataset, Rank = "Genus", bin.dataset, filter_data = T){
  if(Rank == "Genus"){
    temp.data <- dataset %>%
      dplyr::filter(Rank == "Genus" | Rank == "Species")
    names(temp.data)[names(temp.data) == "Genus"] <- "accepted_name"
  }
  if(Rank == "Species"){
    temp.data <- dataset %>%
      dplyr::filter(Rank == "Species")
    names(temp.data)[names(temp.data) == "Combined_name"] <- "accepted_name"
  }
  taxa_data <- subset(temp.data, select=c(accepted_name, Locality, bin_assignment))
  freq_data <- lapply(bin.dataset$bin, function(f) {
    tmp <- taxa_data %>% dplyr::filter(bin_assignment == f) %>% 
      dplyr::count(., accepted_name) %>% 
      dplyr::arrange(desc(n)) %>% 
      tibble::add_row(accepted_name = "n", n = sum(.$n), .before = 1) %>%
      dplyr::select(n)
    freq_raw <- as.numeric(tmp$n)
    freq_raw
  })
  if(deparse(substitute(bin.dataset)) == "series"){
    names(freq_data) <- bin.dataset$bin # give each list element its correct interval name
  }else{
    names(freq_data) <- bin.dataset$name # give each list element its correct interval name
  }
  if(filter_data == T){
    # Set to 2 to avoid issue of singletons
    freq_data <- freq_data[lapply(freq_data,length)>2] 
    # Remove any time periods with only 1 occurrence of each taxa
    #freq_data <- freq_data[names(freq_data) %in% names(Filter(function(vector) all(vector[-1] == 1), freq_data)) == F] 
    # Keep freq_data with above 15 total occurrences (otherwise fails)
    freq_data <- freq_data[lapply(freq_data,function(l)l[[1]]>15) == T]
    freq_data <- freq_data[lapply(freq_data,length)>3] 
  }
  return(freq_data)
}

##### run.iNEXT #####
# Function to run iNEXT, esting estimateD. Uses set quora of 0.3-0.6
run.iNEXT <- function(freq_data, bin.dataset){
  ## And create a new list to store output
  estD_output <- list() 
  ## Create a vector of quorum levels that we want to compute
  ## 0.4 is considered the 'standard', but the fashion now is to plot multiple quorum levels
  quorum_levels <- round(seq(from = 0.3, to = 0.6, by = 0.1), 1)
  
  for(t in 1:1000){
    tryCatch({
      for(i in 1:length(quorum_levels)) {
        # estimateD function in iNEXT, running over each quorum level
        estD_tmp <- iNEXT::estimateD(freq_data, q = c(0,1,2), datatype = "incidence_freq", 
                                     base = "coverage", level = quorum_levels[i])
        # filter to the diversity estimates (order = 1):
        estD_tmp <- dplyr::filter(estD_tmp, Order.q == 0)
        # organise the output:
        estD_tmp$quorum_level <- quorum_levels[i]
        names(estD_tmp)[names(estD_tmp) == "Assemblage"] <- "name"
        
        estD_tmp <- left_join(estD_tmp, bin.dataset, by = "name")
        # add the output to the newly created list
        estD_output[[i]] <- estD_tmp
      }
      print(estD_output)
      if(length(estD_output) > 0){
        print("success")
        break
      }
    }, error = function(e){print(paste("Failed run ", t, sep = "")) 
      print(conditionMessage(e))})
  }
  return(estD_output)
}

################################################################################
# 7. SEDIMENT AFFINITY
################################################################################
##### make.fam.lith #####
make.fam.lith <- function(data){
  fam.lith.data <- data %>%
    dplyr::filter(Rank == "Family"| Rank == "Genus" | Rank == "Species") %>%
    dplyr::filter(Family != "NO_FAMILY_SPECIFIED" & 
                    Family != "Archaeocidaridae or miocidaridae") %>%
    dplyr::filter(!is.na(Family)) %>%
    dplyr::filter(!is.na(Finalised_lith)) %>%
    dplyr::filter(Finalised_lith != "Mixed") %>%
    dplyr::filter(Finalised_lith != "")
  return(fam.lith.data)
}

##### setup.SRA #####
setup.SRA <- function(fam.lith.data, bin.dataset, data.source){
  families <- c("Archaeocidaridae", "Lepidesthidae", "Lepidocentridae", 
                "Palaechinidae", "Proterocidaridae")
  fam_list <- list()
  
  for(f in families){
    print(f)
    fam_df <- fam.lith.data %>%
      dplyr::filter(Family == f)
    num.rows <- nrow(fam_df)
    if(num.rows < 2){
      next
    }
    fam_df <- fam_df %>%
      dplyr::group_by(bin_assignment) %>%
      {as.data.frame(table(.$bin_assignment, .$Finalised_lith))} %>%
      filter(Var2 != "") %>%
      tidyr::pivot_wider(names_from = Var2, 
                         values_from = Freq) %>%
      dplyr::select(Var1, Carbonate, Siliciclastic) %>%
      dplyr::rename(bin = Var1) %>%
      dplyr::mutate(Total = Carbonate + Siliciclastic, 
                    Ac = Carbonate/Total, 
                    As = Siliciclastic/Total)
    
    bin.dataset2 <- bin.dataset %>%
      dplyr::filter(bin %in% fam_df$bin)
    
    fam_df_time <- dplyr::full_join(bin.dataset2, fam_df)
    
    bin.list <- list()
    for(b in fam_df_time$bin){
      print(b)
      bin.list[[match(b, bin.dataset2$bin)]] <- test.func(b, 
                                                          fam_df_time, 
                                                          fam.lith.data = fam.lith.data, 
                                                          f)
    }
    
    sim.aff <- bind_rows(bin.list)
    names(sim.aff)[names(sim.aff) == "Carbonate"] <- "sim.Ac"
    names(sim.aff)[names(sim.aff) == "Siliciclastic"] <- "sim.As"
    sim.aff <- sim.aff %>% dplyr::mutate_all(~replace(., is.nan(.), 0))
    sim.aff[is.na(sim.aff)] <- 0
    fam_df_time <- dplyr::full_join(fam_df_time, sim.aff, by = "bin")
    fam_df_time <- fam_df_time %>%
      dplyr::mutate(Carb_SRA = (Ac - sim.Ac)/SDrndC, 
                    Family = f)
    fam_df_time$Data <- data.source
    fam_list[[match(f, families)]] <- fam_df_time
    names(fam_list)[[match(f, families)]] <- f
  }
  return(fam_list)
}

##### test.func #####
test.func <- function(x, fam_df_time, fam.lith.data, f) {
  bin_total <- fam_df_time %>%
    dplyr::filter(bin == x)
  total.occs <- as.numeric(bin_total[names(bin_total) == "Total"])
  combined <- data.frame(Finalised_lith = c("Carbonate", "Siliciclastic"))
  for(i in 1:1000){
    tmp <- fam.lith.data %>% 
      dplyr::filter(bin_assignment == x) %>%
      dplyr::slice_sample(., n = total.occs, replace = TRUE) %>%
      dplyr::count(., Finalised_lith) %>%
      dplyr::mutate(Affinity = n/total.occs) %>%
      dplyr::select(-n)
    colnames(tmp)[2] <- i
    combined <- dplyr::full_join(combined, tmp, by = "Finalised_lith")
  }
  rownames(combined) <- combined[,1]
  combined[,1] <- NULL
  combined[is.na(combined)] <- 0
  comb.sd <- data.frame(bin = x, 
                        SDrndC = sd(combined[1,], na.rm = T), 
                        SDrndS = sd(combined[2,], na.rm = T))
  combined <- t(as.data.frame(rowMeans(combined, na.rm = T)))
  combined <- as.data.frame(combined)
  combined$bin <- x
  combined <- dplyr::full_join(combined, comb.sd, by = "bin")
  
  return(combined)
}

##### plot.SRA #####
plot.SRA <- function(fam_list, bin.choice, fill_var, 
                     group_var, colour){
  test <- fam_list %>%
    dplyr::select(Carb_SRA, Family, bin, Data) %>%
    dplyr::filter(!is.na(Carb_SRA)) %>%
    dplyr::filter(Carb_SRA != -Inf) %>%
    dplyr::filter(Carb_SRA != Inf)
  test <- left_join(test, bin.choice, by = "bin")
  a <- ggplot(test,   aes_string(x = "mid_ma", y= "Carb_SRA", color = fill_var)) + 
    ## Each quorum level is called individually to be plotted:
    ## Set our line and point sizes (and shapes):
    geom_line() +
    geom_point(size = 2) +
    ## Add our colours, theme, and axes labels:
    scale_colour_manual(values = wes_palette(colour)) +
    scale_x_reverse() +
    labs(x = "Time (Ma)", y = "Standardized Relative Affinity") +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "white"), 
          panel.grid.minor = element_line(colour = "white")) +
    geom_hline(yintercept = 0, 
               linetype = "dashed", 
               color = "light grey") +
    theme(plot.margin = margin(0.1,0.2,0, 1, "cm")) +
    coord_geo(dat = series2, 
              xlim = c(485.4, 251.902),
              ylim = c(-10, 6),
              rot = list(0),
              size = 4,
              pos = list("b"),
              abbrv = T) +
    facet_wrap(group_var, nrow =3)

  if((length(unique(combined.SRA[[group_var]])) %% 2) != 0){
    a <- shift_legend(a)
    a <- ggpubr::as_ggplot(a)
  }
  return(a)
}

##### shift_legend #####
shift_legend <- function(p) {
  # ...
  # to grob
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  
  # establish name of empty panels
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  # example of names:
  #[1] "panel-3-2" "panel-3-3"
  
  # now we just need a simple call to reposition the legend
  reposition_legend(p, 'center', panel=names)
}

################################################################################
# 8. COLONIALISM
################################################################################

# For this section, code was reworked from that provided for the following publication:

# Raja, N.B., Dunne, E.M., Matiwane, A., Khan, T.M., Nätscher, P.S., Ghilardi, A.M. 
# and Chattopadhyay, D., 2022. Colonial history and global economics distort our 
# understanding of deep-time biodiversity. Nature ecology & evolution, 6(2), pp.145-154.

# We thank the authors for providing this code and for enabling this research. 

##### setup.col.sec #####
# Function to run and plot network analysis
setup.col.sec <- function(data){
  ##### Setup data for analysis #####
  df <-  data %>%
    dplyr::select(Locality, Museum_country, Country, lat, lng) %>%
    distinct() %>%
    dplyr::filter(!is.na(Museum_country)) %>%
    dplyr::filter(Museum_country != "") %>%
    dplyr::filter(!is.na(Country))  %>%
    dplyr::rename(samp_country = Country, 
                  aff_country = Museum_country) %>%
    dplyr::mutate(collection_no = 1:nrow(.)) %>%
    dplyr::select(-c(Locality, lat, lng))
  
  df[df == "England" | df == "Wales" | df == "Scotland" | df == "Northern Ireland"] <- "United Kingdom"
  df[df == "USA"] <- "United States"
  df[df == "Czechia"] <- "Czech Republic"
  
  df$aff_code <- countrycode(df$aff_country, "country.name", "iso3c")
  df$samp_code <- countrycode(df$samp_country, "country.name", "iso3c")
  return(df)
}

##### network.plots #####
network.plots <- function(df){
  # each country gets an id
  countries <- unique(c(df$aff_country, df$samp_country))
  countries <- gsub("&", "and", countries)
  countries <- data.frame(id=1:length(countries),
                          label=countries)
  countries$region <- countrycode::countryname(countries$label, 
                                               destination = "continent") #may change
  countries <- countries[order(countries$region),]
  
  regs <- unique(countries$region)
  cols <- palette.colors(length(regs))
  countries$cols <- plyr::mapvalues(countries$region, from=regs, to=cols)
  
  # By region
  df$samp_region <- countrycode::countryname(df$samp_country, destination="region23")
  df$aff_region <- countrycode::countryname(df$aff_country, destination="region23")
  
  df$samp_continent <- countrycode::countryname(df$samp_country, destination="continent")
  df$aff_continent <- countrycode::countryname(df$aff_country, destination="continent")
  
  # each region gets an id
  regions <- unique(c(df$samp_region, df$aff_region))
  regions <- data.frame(id=1:length(regions),
                        label=regions)
  
  ##### Create edge list #####
  edges <- df %>%  
    dplyr::group_by(samp_country, aff_country) %>% 
    dplyr::summarise(weight = n()) %>% 
    dplyr::ungroup()
  
  #add node id
  edges <- edges %>% 
    dplyr::left_join(countries, by = c("aff_country" = "label")) %>% 
    dplyr::rename(from = id)
  
  edges <- edges %>% 
    dplyr::left_join(countries %>%  dplyr::select(id, label), 
                     by = c("samp_country" = "label")) %>% 
    dplyr::rename(to = id)
  
  edges <- na.omit(edges)
  
  ###### Plot #####
  edges <- edges[edges$samp_country != "",] 
  edges2 <- edges[edges$samp_country != edges$aff_country,]
  
  df2 <- countries %>%  dplyr::filter(!is.na(region))
  df2 <- df2 %>% dplyr::left_join(edges %>%  dplyr::select(aff_country, weight), 
                                  by=c("label"="aff_country"))
  
  df2$weight[is.na(df2$weight)] <- 0
  #df2 <- df2[df2$label %in% unique(edges2$samp_country, edges2$aff_country),]
  
  df2 <- df2 %>%  dplyr::group_by(label, region) %>% 
    dplyr::summarise(weight=sum(weight))
  
  df2 <- df2[order(df2$region, df2$weight, decreasing = T),]
  regs <- unique(df2$region)
  
  df2$y <- 0
  df2$x <- 0
  
  seq_max <- c(20,22,21, 22, 18)
  x_mid <- c(5,4,2,3,1)
  f <- c(3,5,8,5,10)
  
  for(i in 1:length(regs)){
    temp <- df2[df2$region == regs[i],]
    df2[df2$region == regs[i],]$y <- seq(seq_max[i],1, length.out = nrow(temp))
    set.seed(42)
    df2[df2$region == regs[i],]$x <- jitter(rep(x_mid[i], nrow(temp)), factor=f[i])
  }
  
  df2[df2$label=="United States",]$y <-  24
  df2[df2$label=="China",]$y <- 22
  df2[df2$label=="Germany",]$y <- 22.8
  
  edges2 <- edges2 %>% 
    dplyr::left_join(df2 %>%  dplyr::select(label, y,x), 
                     by=c("samp_country" = "label")) %>% 
    dplyr::rename(x2 = x, y2=y) %>% 
    dplyr::left_join(df2 %>%  dplyr::select(label, y,x), 
                     by=c("aff_country" = "label")) %>% 
    dplyr::rename(x1=x, y1=y)
  
  p.breaks <- df2 %>% dplyr::group_by(region) %>% 
    dplyr::summarise(x=median(x)) %>% 
    arrange(x)
  
  edges2$region[edges2$weight < 26] <- "Other"
  
  edges2$region <- factor(edges2$region,levels=c(regs, "Other"))
  df2$weight_alpha <- df2$weight
  df2$weight_alpha[df2$label %in% df2$label[order(df2$weight, decreasing = T)][1:25]] <- 9999
  
  curvature = 0.3
  alpha = 0.8
  breaks <- c(50,100,500,1000,
              2000, 5000, 10000)
  lab_size <- 3
  
  palv <- viridis::viridis(5)
  names(palv) <- c("Europe", "Americas", "Oceania", "Asia", "Africa")
  palv["Europe"] <- "#75428f"
  palv["Americas"] <- "#365091"
  palv["Africa"] <- "#fdaa25"
  palv["Asia"] <- "#027218"
  
  p1 <- ggplot() + geom_curve(data=na.omit(edges2[edges2$x2 > edges2$x1,]), 
                              aes(x=x1, y=y1, xend=x2, yend=y2, col=region
                              ), 
                              curvature = curvature, 
                              alpha=alpha,
                              size=0.5) +
    geom_curve(data=na.omit(edges2[edges2$x1 > edges2$x2,]), 
               aes(x=x1, y=y1, xend=x2, yend=y2, col=region), 
               curvature = -curvature, alpha=alpha) +
    geom_point(data=df2, 
               aes(x=x, y=y, 
                   size=weight, col=region,
                   alpha=weight_alpha)) +
    geom_text(data=df2[df2$x < 3,], aes(x=x, y=y, 
                                        label=ifelse(label %in% df2$label[order(df2$weight, decreasing = T)][1:25], label, NA)),
              hjust=1, vjust=0, size=lab_size, fontface="bold", nudge_x = -0.1) + 
    geom_text(data=df2[df2$x > 3,], aes(x=x, y=y, 
                                        label=ifelse(label %in% df2$label[order(df2$weight, decreasing = T)][1:25], label, NA)),
              hjust=0, vjust=0, size=lab_size, nudge_y = 0.2, nudge_x = 0.1, fontface="bold") +
    scale_x_continuous(breaks=p.breaks$x, limits=c(0,6)) +
    scale_size_continuous(breaks=breaks, 
                          range=c(1,15)) +
    scale_color_manual(values=c("Other"="#e5e5e560", 
                                palv)) +
    scale_alpha_continuous(breaks=c(0, 50,100,200,500, 1000,2000, 9999), 
                           range=c(0.3,1)) +
    guides(alpha="none") +
    labs(col="Region", size="Number of \noutgoing nodes") +
    theme_void() 
  return(p1)
}

##### reshape.df #####
# Function for setting up network statistics
reshape.df <- function(x){
  x <- sort(x, decreasing = TRUE)[1:15]
  names <- names(x)
  x <- data.frame(value=x)
  if(rownames(x)[1] == "1"){
    x <- na.omit(x)
    rownames(x) <-  names[complete.cases(names)]
  }
  x$code <- rownames(x)
  x$country <- countrycode(x$code, "iso3c", "country.name")
  return(x)
}

##### research.prep #####
# Function to plot foreign vs. domestic contributions to fossil collections
research.prep <- function(df, data){
  colls_n <- data.frame(table(df$aff_code), stringsAsFactors = FALSE)
  colnames(colls_n) <- c("code", "freq")
  total_colls <- length(unique(df$collection_no))
  colls_n$freq <- colls_n$freq/total_colls
  colls_n <- colls_n[order(colls_n$freq, decreasing = TRUE),]
  topcountries <- colls_n[1:15,]
  
  topcountries <- na.omit(topcountries)
  topcountries$country <- countrycode(topcountries$code, 
                                      origin = "iso3c", 
                                      destination = "country.name")
  
  dat2 <- df[df$aff_code != df$samp_code,]
  dat2 <- unique(dat2[,c("collection_no", "aff_code")])
  
  colls_n2 <- data.frame(table(dat2$aff_code), stringsAsFactors = FALSE)
  colnames(colls_n2) <- c("code", "foreign")
  colls_n2$foreign <- colls_n2$foreign/total_colls
  
  colls_foreign <-colls_n2[colls_n2$code %in% topcountries$code,]
  
  topcountries <- merge(topcountries, colls_n2, all = TRUE)
  topcountries[is.na(topcountries)] <- 0
  topcountries$local <- topcountries$freq - topcountries$foreign
  
  if(nrow(topcountries) > 15){
    topcountries <- topcountries[order(topcountries$freq, decreasing = TRUE),]
    topcountries <- topcountries[1:15,]
  }
  
  topcountries <- topcountries %>%  dplyr::select(local, foreign, country) %>% 
    pivot_longer(cols=c("local", "foreign"), names_to = c("type"), values_to="freq")
  
  topcountries$type <- factor(topcountries$type, levels=c("local", "foreign"))
  topcountries$Data <- data
  
  return(topcountries)
}

##### continent.stats #####
continent.stats <- function(df, data){
  codes <- c("samp_code", "aff_code")
  iso2c_to_continents <- c(CAN = "North America" , USA = "North America")
  
  
  a <- bind_rows(lapply(codes, function(x){
    dat2 <- df[c(x, "collection_no")]
    dat2 <- unique(dat2[,c("collection_no", x)])
    dat2$Continent <- countrycode(dat2[,x], 
                                  origin = "iso3c", 
                                  destination = "continent", 
                                  custom_match = iso2c_to_continents)
    
    dat2 <- data.frame(table(dat2$Continent), stringsAsFactors = FALSE)
    
    colnames(dat2) <- c("Continent", "Total")
    
    dat2$Frequency <- dat2$Total/sum(dat2$Total)
    dat2$Origin <- ifelse(x == "samp_code", "Country", "Museum")
    dat2 <- dat2 %>%
      mutate(Continent = str_replace(Continent, "Americas", "South America"))
    return(dat2)
  }))
  a$Data <- data
  return(a)
}

##### top.countries #####
top.countries <- function(df, code, sum = NA){
  colls_n <- data.frame(table(df[code]), stringsAsFactors = FALSE)
  colnames(colls_n) <- c(code, "freq")
  total_colls <- length(unique(df$collection_no))
  colls_n$freq <- colls_n$freq/total_colls
  colls_n <- colls_n[order(colls_n$freq, decreasing = TRUE),]
  topcountries <- colls_n[1:15,]
  
  topcountries <- na.omit(topcountries)
  if(is.na(sum) == T){
    return(topcountries)
  }else{
    perc <- sum(topcountries$freq[1:sum])*100
    return(list(topcountries[1:sum,], perc))
  }
}

##### network.stats #####
network.stats <- function(df){
  
  pal <- c("#f0ffe9", "#ffe599", "#bbe487", "#4e9755", "#173109")
  
  # Setup
  gr <- df[,c("aff_code", "samp_code")]
  gr <- gr[gr$aff_code != gr$samp_code,] # remove self nodes
  gr <- igraph::graph_from_data_frame(gr, directed = T)
  
  # The degree of a node in a network is the number of connections it has to other 
  # nodes and the degree distribution is the probability distribution of these degrees 
  # over the whole network. 
  gr$degree <- igraph::degree(gr) 
  
  # Detecting the amount of influence a node has over the flow of information in a graph. 
  gr$betweenness <- igraph::betweenness(gr) 
  
  # The more central a node is, the closer it is to all other nodes. 
  gr$closeness <- igraph::closeness(gr) 
  
  # Run function
  degree_top15 <- reshape.df(gr$degree)
  betweenness_top15 <- reshape.df(gr$betweenness)
  closeness_top15 <- reshape.df(gr$closeness)
  
  # Plot Degree
  p3 <- ggplot(degree_top15, aes(x=reorder(country, value), y=value)) +
    geom_bar(stat="identity", fill=pal[4]) + coord_flip(expand=FALSE) +
    labs(y="Degree", x="") +
    ggthemes::theme_hc()
  
  # Plot betweenness
  p4 <- ggplot(betweenness_top15, aes(x=reorder(country, value), y=value)) +
    geom_bar(stat="identity", fill=pal[4]) + coord_flip(expand=FALSE) +
    labs(y="Betweenness", x="") +
    ggthemes::theme_hc()
  
  # Plot closeness
  p5 <- ggplot(closeness_top15, aes(x=reorder(country, value), y=value)) +
    geom_bar(stat="identity", fill=pal[4]) + coord_flip(expand=FALSE) +
    labs(y="Closeness", x="") +
    ggthemes::theme_hc()
  
  # Combine plots
  p6 <- ggarrange(p3, p4, p5, 
                  align='hv',
                  labels = c("A", "B", "C"),
                  nrow = 1, 
                  ncol = 3)
  return(p6)
}

################################################################################
# 9. MODELLING COMPARISONS
################################################################################
##### run.models #####
run.models <- function(dataset, bin.dataset, div = "raw", na = T){
  
  if(na == T){
    # Reduce to N. America
    dataset <- dataset %>%
      dplyr::filter(Continent == "North America")
  } 
  
  ##### DIVERSITY #####
  # Get diversity estimates
  if(div == "raw"){
    model.div <- run.div(dataset, "Genus", type = "raw")
    model.div <- model.div %>%
      dplyr::select(bin_assignment, divSIB) %>%
      dplyr::rename(bin = bin_assignment)
  } 
  if(div == "SQS"){
    model.div <- run.div(dataset, "Genus", type = "SQS")
    # Keep relevant data only
    model.div <- model.div[[4]] %>%
      dplyr::select(bin_assignment, divCSIB) %>%
      dplyr::rename(bin = bin_assignment)
  }
  comb.mod <- full_join(comb.mod, model.div, by = "bin")
  if(div == "raw"){
    comb.mod$divSIB <- log10(comb.mod$divSIB+1) # Add 1 to stop infinity appearing
  }else{
    comb.mod$divSIB <- log10(comb.mod$divCSIB)
  }
  
  comb.mod <- comb.mod %>%
    dplyr::filter(bin_midpoint > 253 & bin_midpoint < 481.55000)
  
  div.model <- comb.mod %>%
    dplyr::filter(!is.na(divSIB)) %>%
    dplyr::filter(!is.na(colls)) %>%
    dplyr::filter(!is.na(hardie)) 
  
  if(na == T){
    mod.form <- as.formula("divSIB ~ mean_sl + carb + sili + u + colls + hardie + wilkinson")
  }else{
    mod.form <- as.formula("divSIB ~ mean_sl + colls + u + hardie + wilkinson")
  }
  
  # Full model
  full.model <- gls(mod.form, data = div.model, method = "ML", 
                    correlation = corARMA(p=1))
  
  dredged <- dredge(full.model, rank = "AICc")
  top_model <- get.models(dredged, subset = 1)[[1]]
  
  # model results
  model.results <- list()
  model.results[[1]] <- summary(top_model)
  if(dredged$delta[2] < 2){
    model.results[[2]] <- summary(model.avg(dredged, subset = delta <= 2))
  }else{
    model.results[[2]] <- summary(model.avg(dredged, subset = delta <= dredged$delta[2]))
  }
  model.results[[3]] <- dredged
  model.results[[4]] <- r.squaredLR(top_model)
  names(model.results) <- c("Top.model", "Mod.Av", "All.mods", "R.squared")
  return(model.results)
}

##### ca.binning #####
ca.binning <- function(data){
  test <- data %>%
    dplyr::filter(max_ma < max(stages$max_ma)) %>%
    dplyr::filter(min_ma > min(stages$min_ma))
  test <- bin_time(occdf =  test, bins = stages, method = 'mid')
  test <- test %>%
    dplyr::group_by(bin_midpoint) %>%
    dplyr::summarize(mean = mean(MgCa))
  return(test)
}

##### model.results #####
model.results <- function(model.list){
  # Getting top model data
  a <- bind_rows(lapply(seq_along(model.list), function(x){
    name <- names(model.list[x])
    test <- as.data.frame(model.list[[x]][["All.mods"]][1,])
    test$Data <- name
    return(test)
  }))
  # Getting top model covariates/signif
  b <- bind_rows(lapply(seq_along(model.list), function(x){
    name <- names(model.list[x])
    test <- as.data.frame(model.list[[x]][["Top.model"]]["tTable"])
    colnames(test) <- c("Estimate", "Std. Error", 
                        "T value", "p value")
    test <- cbind(data.frame(Data = name, Covariate = rownames(test)), test)
    rownames(test) <- NULL
    test$Signif <- ifelse(test$`p value`< 0.05, "*", "")
    return(test)
  }))
  # Getting model average data
  c <- bind_rows(lapply(seq_along(model.list), function(x){
    name <- names(model.list[x])
    test <- as.data.frame(model.list[[x]][["Mod.Av"]]["coefmat.full"])
    colnames(test) <- c("Estimate", "Std. Error", 
                        "Adjusted SE", "Z value", "p value")
    test <- cbind(data.frame(Data = name, Covariate = rownames(test)), test)
    rownames(test) <- NULL
    test$Signif <- ifelse(test$`p value`< 0.05, "*", "")
    return(test)
  }))
  d <- bind_rows(lapply(seq_along(model.list), function(x){
    name <- names(model.list[x])
    test <- as.data.frame(model.list[[x]][["R.squared"]][1])
    test <- cbind(data.frame(Data = name), test)
    colnames(test) <- c("Data", "R.squared")
    return(test)
  }))
  return(list(a, b, c, d))
}