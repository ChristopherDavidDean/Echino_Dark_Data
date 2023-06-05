################################################################################
#                                                                              #
# ECHINOID PAPERS DATA SETUP                                                   #
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
