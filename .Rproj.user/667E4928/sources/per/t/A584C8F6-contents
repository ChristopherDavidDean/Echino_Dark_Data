################################################################################
#                                                                              #
# ECHINOID DARK DATA                                                           #
# C.D. Dean, M. Clapham, T.A.M. Ewin, J. Thompson                              #
# Code by C. D. Dean, started 29.09.23                                         #
#                                                                              #
################################################################################

#================================= SETUP =======================================

# Load packages
library(palaeoverse)
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
library(countrycode)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Create factors for later
order_ind <- c("Permian", "Pennsylvanian", "Mississippian", "Devonian", "Silurian", "Ordovician")

#===== Load and prepare museum dataset =====

# Load Museum dataset
m.dat <- read.csv("Palaeozoic_Occurrences_JRT_Finalizing_3.csv")

# Truncate columns and adjust column names
m.dat <- m.dat[,1:37] 
colnames(m.dat)[17] <- "max_ma"
colnames(m.dat)[18] <- "min_ma"
m.dat[m.dat == "?"] <- NA
m.dat[m.dat == ""] <- NA

# Remove Triassic occurrences
m.dat <- m.dat %>%
  filter(Period != "Triassic") %>%
  filter(Period != "Triassic/Jurassic")

# Sort stages - Deeptime
data(stages)
colnames(stages)[2] <- "max_ma"
colnames(stages)[3] <- "min_ma"
stages$bin <- 1:nrow(stages)

# Assign numerical ages to each specimen
m.dat <- look_up(occdf = m.dat, early_interval = "max_ma", late_interval = "min_ma")

# Rename columns
names(m.dat)[names(m.dat) == "max_ma"] <- "old_max_ma"
names(m.dat)[names(m.dat) == "min_ma"] <- "old_min_ma"
names(m.dat)[names(m.dat) == "interval_max_ma"] <- "max_ma"
names(m.dat)[names(m.dat) == "interval_min_ma"] <- "min_ma"

# Remove occurrences without associated temporal data (in second dataset)
m.dat.NA <- m.dat %>%
  filter(is.na(max_ma))

m.dat <- m.dat %>%
  filter(!is.na(max_ma))

# Assign occurrences to bins (in second dataset)
m.dat <- bin_time(m.dat, bins = stages, method = 'all')

#Combine datasets
m.dat <- dplyr::bind_rows(m.dat, m.dat.NA)

# Order Periods by index
m.dat$Period <- as.factor(m.dat$Period)
m.dat$Period <- factor(m.dat$Period, levels = order_ind)

# Made occurrence based dataset instead of specimen based
m.dat.occ <- m.dat %>%
  dplyr::select(-c(Museum.Number, Original.Order, Preservation.Score, Description, id)) %>%
  dplyr::distinct()

#===== Load and prepare PBDB dataset =====
 
# Load dataset
pb.dat <- read.csv("pbdb_echinoids.csv", skip = 18)

# Assign occurrences to bins
pb.dat <- bin_time(pb.dat, bins = stages, method = 'all')

# Covert country code into full words
pb.dat$Country <- countrycode(sourcevar = pb.dat$cc, "eurostat", "country.name")
pb.dat$Country[is.na(pb.dat$Country) == TRUE] <- "Greece"

# Convert lithologies to broad categories



# Load all Palaeozoic occurrences
pb.all <- read.csv("all_palaeozoic.csv", skip = 18)

#=========================== BASIC COMPARISONS =================================

# Maps
e <<- extent(-180, 180, -90, 90)
get_grid_im(m.dat.occ, 1,  "Museum Echinoids", ext = e)
get_grid_im(pb.dat, 1,  "PBDB Echinoids", ext = e)
get_grid_im(pb.all, 1,  "PBDB", ext = e)


# Per lithology
table(pb.dat$lithology1)

# Formations
table(pb.dat$formation)
table(m.dat$Formation)

# Tax range
tax_range_time(pb.dat, name = "genus", by = "name", plot = TRUE)
m.dat.occ <- m.dat.occ %>%
  filter(!is.na(Genus)) %>%
  filter(!is.na(max_ma))
tax_range_time(m.dat.occ, name = "Genus", by = "name", plot = TRUE)

#=========================== OCCUPANCY COMPARISON ==============================

# Read in all palaeozoic marine occurrenceS (full spatial coverage)
read.csv() 


