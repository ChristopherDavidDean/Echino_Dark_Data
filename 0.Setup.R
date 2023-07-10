################################################################################
#                                                                              #
# ECHINOID PAPERS DATA SETUP AND FUNCTIONS                                     #
# Code by C. D. Dean, started 29.09.23                                         #
#                                                                              #
################################################################################

################################################################################
# 1. SETUP
################################################################################

# Load packages
library(palaeoverse)
library(countrycode)
library(dplyr)
library(ggplot2)
library(geoscale)
library(vegan)
library(RColorBrewer)
library(divDyn)
library(patchwork)
library(chronosphere)
library(devtools)
library(tidyr)
library(vcd)
library(maditr)
library(MASS)
library(effects)
library(deeptime)
library(wesanderson)
library(viridis)
library(raster)
library(MuMIn)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data
m.dat <- read.csv("Final_Database_for_analysis.csv")
m.dat[m.dat == "?"] <- NA
m.dat[m.dat == ""] <- NA

# Remove occurrences without taph. grade and Triassic occurrences
m.dat <- m.dat %>%
  filter(Max_period != "Triassic") %>%
  filter(Min_period != "Triassic") %>%
  filter(Max_period != "Triassic/Jurassic") %>%
  filter(Min_period != "Triassic/Jurassic") %>%
  filter(Max_period != "Jurassic") %>%
  filter(Min_period != "Jurassic") %>%
  filter(is.na(Preservation_score) == F) 

# Sort stages - Deeptime
data(stages)
data(periods)
colnames(stages)[2] <- "max_ma"
colnames(stages)[3] <- "min_ma"
stages$bin <- 1:nrow(stages)

# Resolve promise
periods

################################################################################
# 2. ASSIGNING AGES
################################################################################

# Separate datasets
period <- m.dat %>%
  dplyr::filter(Age_resolution == "Period" | Age_resolution == "Series")
stage <- m.dat %>%
  dplyr::filter(Age_resolution == "Stage") 
  

# Assign numerical ages to periods
colnames(period)[28] <- "max_ma"
colnames(period)[29] <- "min_ma"
period <- look_up(occdf = period, early_interval = "max_ma", late_interval = "min_ma")

# Assign numerical ages to stages
colnames(stage)[30] <- "max_ma"
colnames(stage)[31] <- "min_ma"
stage <- look_up(occdf = stage, early_interval = "max_ma", late_interval = "min_ma")

# Reorganise column headings
stage <- stage[,-c(30:31)]
period <- period[,-c(30:31)]
colnames(period)[28:29] <- c("Max_period", "Min_period")

# Bind datasets
m.dat <- rbind(period, stage)

# Rename columns
names(m.dat)[names(m.dat) == "interval_max_ma"] <- "max_ma"
names(m.dat)[names(m.dat) == "interval_min_ma"] <- "min_ma"

# Remove NAs (NOTE CAN REMOVE THIS LATER!)
m.dat <- m.dat %>%
  filter(is.na(max_ma) == F)

# Assign occurrences to bins
m.dat <- bin_time(m.dat, bins = stages, method = 'majority')

# Create factors for later
order_ind <- c("Permian", "Pennsylvanian", "Mississippian", "Devonian", "Silurian", "Ordovician")
m.dat$Max_period <- as.factor(m.dat$Max_period)
m.dat$Min_period <- as.factor(m.dat$Min_period)
m.dat$Min_period <- factor(m.dat$Max_period, levels = order_ind)
m.dat$Min_period <- factor(m.dat$Min_period, levels = order_ind)

################################################################################
# 3. PALAEOROTATION
################################################################################

# Filter for specimens which have geographic coords and make numeric
m.dat.rotate <- m.dat %>%
  filter(is.na(lat) != T)
m.dat.rotate$lat <- as.numeric(m.dat.rotate$lat)
m.dat.rotate$lng <- as.numeric(m.dat.rotate$lng)

# Palaeorotate
m.dat.rotate <- palaeorotate(m.dat.rotate, 
             lng = 'lng', 
             lat = 'lat', 
             age = "bin_midpoint",
             model = "PALEOMAP",
             method = "point")

################################################################################
# 4. FUNCTIONS
################################################################################

# Bar plot
make.bar.plot <- function(dataset, fill, legend, colour, flip = FALSE){
  a <- ggplot(dataset) +
    aes(x = Preservation_score, fill = fill) +
    geom_bar() +
    scale_fill_manual(values=wes_palette(colour)) +
    guides(fill=guide_legend(title=legend)) +
    ylab("Frequency") +
    xlab("Preservation Score") +
    theme_bw()
  if(flip == TRUE){
    a <- a + coord_flip()
  }
  return(a)
}

# Proportional bar plot
make.prop.bar.plot <- function(dataset, fill, legend, colour, flip = FALSE){
  a <- ggplot(dataset) +
    aes(x = Preservation_score, fill = fill) +
    geom_bar(position = 'fill') +
    scale_fill_manual(values=wes_palette(colour)) +
    guides(fill=guide_legend(title=legend)) +
    ylab("Proportion of dataset") +
    xlab("Preservation Score") +
    theme_bw() +
    scale_y_continuous(expand = c(0,0), limits = c(0,1))
  if(flip == TRUE){
    a <- a + coord_flip()
  }
  return(a)get_grid_im <- function(data, res, name, ext){ # Data is first output from combine_data (fossil.colls). Res is chosen resolution in degrees. name is user inputted string related to data inputted, for display on graphs. 
    xy <- cbind(as.double(data$lng), as.double(data$lat))
    #xy <- unique(xy)
    r <- raster::raster(ext = ext, res = res)
    r <- raster::rasterize(xy, r, fun = 'count')
    #r[r > 0] <- 1 # Remove if you want values instead of pure presence/absence.
    countries <- maps::map("world", plot=FALSE, fill = TRUE) # find map to use as backdrop
    countries <<- maptools::map2SpatialPolygons(countries, IDs = countries$names, proj4string = CRS("+proj=longlat")) # Turn map into spatialpolygons
    mapTheme <- rasterVis::rasterTheme(region=viridis(8))
    print(rasterVis::levelplot(r, margin=F, par.settings=mapTheme,  main = paste("Total ", (substitute(name)), " per Grid Cell", sep = "")) + #create levelplot for raster
            #   latticeExtra::layer(sp.polygons(states, col = "white", fill = NA), under = T)  + # Plots state lines
            latticeExtra::layer(sp.polygons(countries, col = 0, fill = "light grey"), under = T)) # Plots background colour
    hist(r, breaks = 20,
         main = paste((substitute(name)), " per Grid Cell", sep = ""),
         xlab = "Number of Collections", ylab = "Number of Grid Cells",
         col = "springgreen")
    r <<- r
  }
}

# Set preservation score for Logistic Regression models
set_Pres_score <- function(dataset, level){
  # if score = 5, class as 1, otherwise class as 0.
  dataset$LR_Pres_score <- 0
  if(length(level) > 1){
    for(n in 1:length(level)){
      dataset$LR_Pres_score[dataset$Preservation_score == level[n]] <- 1
    }
  }else{
    # Set values in 'new_column' to 1 where the original_column is equal to 5
    dataset$LR_Pres_score[dataset$Preservation_score == level] <- 1
  }
  dataset$LR_Pres_score <- as.factor(dataset$LR_Pres_score)
  dataset <- dataset
}

# Function to split out and reorganise data by preservation score for correlations
split <- function(data, Var2, order = FALSE, stage = FALSE){
  for(t in Var2){
    temp.data <- filter(data, Var2 == t)
    if(order == TRUE){
      temp.data <- temp.data[match(order_ind, temp.data$Var1),]
    }
    if(stage == TRUE){
      temp_stages <- stages[55:92,]
      temp_stages$bin_midpoint <- (temp_stages$max_ma + temp_stages$min_ma)/2 
      temp.data <- merge(temp_stages, temp.data, by = "bin_midpoint", all = TRUE)
      temp.data$Preservation_score <- t
      temp.data$Freq[is.na(temp.data$Freq)] <- 0
    }
    assign(paste(deparse(substitute(data)), t, sep = "."), temp.data, envir = .GlobalEnv)
  }
}

# Function for making maps
get_grid_im <- function(data, res, name, ext){ # Data is first output from combine_data (fossil.colls). Res is chosen resolution in degrees. name is user inputted string related to data inputted, for display on graphs. 
  xy <- cbind(as.double(data$lng), as.double(data$lat))
  #xy <- unique(xy)
  r <- raster::raster(ext = ext, res = res)
  r <- raster::rasterize(xy, r, fun = 'count')
  #r[r > 0] <- 1 # Remove if you want values instead of pure presence/absence.
  countries <- maps::map("world", plot=FALSE, fill = TRUE) # find map to use as backdrop
  countries <<- maptools::map2SpatialPolygons(countries, IDs = countries$names, proj4string = CRS("+proj=longlat")) # Turn map into spatialpolygons
  mapTheme <- rasterVis::rasterTheme(region=viridis(8))
  print(rasterVis::levelplot(r, margin=F, par.settings=mapTheme,  main = paste("Total ", (substitute(name)), " per Grid Cell", sep = "")) + #create levelplot for raster
          #   latticeExtra::layer(sp.polygons(states, col = "white", fill = NA), under = T)  + # Plots state lines
          latticeExtra::layer(sp.polygons(countries, col = 0, fill = "light grey"), under = T)) # Plots background colour
  hist(r, breaks = 20,
       main = paste((substitute(name)), " per Grid Cell", sep = ""),
       xlab = "Number of Collections", ylab = "Number of Grid Cells",
       col = "springgreen")
  r <<- r
}