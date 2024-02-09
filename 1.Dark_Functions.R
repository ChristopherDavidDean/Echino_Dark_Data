################################################################################
### MUSEUM `DARK DATA` ILLUMINATES THE EVOLUTIONARY HISTORY OF FOSSIL GROUPS ###
################################################################################

# Christopher D. Dean & Jeffrey R. Thompson
# 2024
# Script written by Christopher D. Dean

################################################################################
#                             FILE 1: FUNCTIONS                                #
################################################################################

# function that automatically installs necessary packages that the user is lacking.

#===== iPAK =====
ipak <- function(pkg){ # Function to install packages. Read in character vector of any packages required. 
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
} 

#===== REQUIRED PACKAGES =====

library(tidyr)
library(tibble)
library(raster)
library(plyr)
library(latticeExtra)
library(rasterVis)
library(sp)
library(maptools)
library(stringr)
library(patchwork)
library(sf)
library(rgeos)
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

#==================================== GET_GRID ===========================================

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

#=============================================== GET_GRID_IM ===========================================================

# Set up background info
countries <- maps::map("world", plot=FALSE, fill = TRUE) # find map to use as backdrop
states <- maps::map("state", plot = FALSE, fill = TRUE)
countries <<- maptools::map2SpatialPolygons(countries, IDs = countries$names, proj4string = CRS("+proj=longlat")) # Turn map into spatialpolygons
states <<- maptools::map2SpatialPolygons(states, IDs = states$names, proj4string = CRS("+proj=longlat")) # Turn map into spatialpolygons

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

