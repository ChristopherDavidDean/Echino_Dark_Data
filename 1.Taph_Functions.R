################################################################################
######## TAPHONOMIC CONTROLS ON A MULTI-ELEMENT SKELETAL FOSSIL RECORD #########
################################################################################

# Jeffrey R. Thompson, Christopher D. Dean, Madeline Ford, Timothy A. M. Ewin
# 2024
# Script written by Christopher D. Dean

################################################################################
#                            FILE 1: FUNCTIONS                                 #
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
library(ggpubr)
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
library(sjPlot)

################################################################################
# 2. FUNCTIONS
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
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "white"), 
          panel.grid.minor = element_line(colour = "white"), 
          legend.position = "bottom") +
    guides(alpha = "none")
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
    theme(panel.grid.major = element_line(colour = "white"), 
          panel.grid.minor = element_line(colour = "white")) +
    guides(alpha = "none") +
    scale_y_continuous(expand = c(0,0), limits = c(0,1))
  if(flip == TRUE){
    a <- a + coord_flip()
  }
  return(a)
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