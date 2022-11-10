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
order_ind <- c("Permian", "Pennsylvanian", "Mississippian", "Devonian", 
               "Silurian", "Ordovician")

# Load required functions
source("0.Functions_Dark_data_test.R")

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
m.dat <- look_up(occdf = m.dat, 
                 early_interval = "max_ma", 
                 late_interval = "min_ma")

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

# Combine datasets
m.dat <- dplyr::bind_rows(m.dat, m.dat.NA)

# Order Periods by index
m.dat$Period <- as.factor(m.dat$Period)
m.dat$Period <- factor(m.dat$Period, levels = order_ind)

# Made occurrence based dataset instead of specimen based
m.dat.occ <- m.dat %>%
  dplyr::select(-c(Museum.Number, Original.Order, 
                   Preservation.Score, Description, id)) %>%
  dplyr::distinct()

#===== Load and prepare PBDB datasets =====
 
# Load Paleozoic echinoid dataset
pb.dat <- read.csv("pbdb_echinoids.csv", skip = 18)

# Assign occurrences to bins
pb.dat <- bin_time(pb.dat, bins = stages, method = 'all')

# Covert country code into full words
pb.dat$Country <- countrycode(sourcevar = pb.dat$cc, "eurostat", "country.name")
pb.dat$Country[is.na(pb.dat$Country) == TRUE] <- "Greece"
pb.dat$Continent <- countrycode(sourcevar = pb.dat$Country, "country.name", "continent")
pb.dat$Continent[which(pb.dat$Country == "United States" | pb.dat$Country== "Canada")] <- "North America"
pb.dat$Continent[which(pb.dat$Continent == "Americas")] <- "South America"

# Convert lithologies to broad categories


# Load all Palaeozoic occurrences
pb.all <- read.csv("all_palaeozoic.csv", skip = 18)

# Assign occurrences to bins (WARNING: Takes a looooong time)
pb.all <- pb.all %>%
  filter(max_ma < stages[nrow(stages), 2])
pb.all <- bin_time(pb.all, bins = stages, method = 'all')

#=========================== BASIC COMPARISONS =================================

#===== Comparative maps =====

# Set extent
e <<- extent(-180, 180, -90, 90)

# Make maps
get_grid_im(m.dat.occ, 1,  "Museum Echinoids", ext = e)
get_grid_im(pb.dat, 1,  "PBDB Echinoids", ext = e)
get_grid_im(pb.all, 1,  "PBDB", ext = e)

#===== Table comparisons =====

# Per lithology - TO DO
table(pb.dat$lithology1)

# Formations
table(pb.dat$formation)
table(m.dat$Formation)
length(unique(m.dat$Formation))
length(unique(pb.dat$formation))

# Continent/Country
table(pb.dat$Country)
table(m.dat$Country)
country.long <- rbind(cbind(as.data.frame(table(m.dat$Country)), data = c(rep("M", nrow(as.data.frame(table(m.dat$Country)))))), 
      cbind(as.data.frame(table(pb.dat$Country)), data = c(rep("P", nrow(as.data.frame(table(pb.dat$Country)))))))

table(m.dat$Continent)
table(pb.dat$Continent)
cont.long <- rbind(cbind(as.data.frame(table(m.dat$Continent)), data = c("M", "M", "M", "M")),
      cbind(as.data.frame(table(pb.dat$Continent)), data = c("P", "P", "P", "P", "P", "P")))

#===== Taxonomic range =====
# filter genera from m.dat so that datasets match - TO DO

# Calculate temporal range for PBDB genera
tax_range_time(pb.dat, name = "genus", by = "name", plot = TRUE)

# Calculate temporal range for Museum genera
m.dat.occ <- m.dat.occ %>%
  filter(!is.na(Genus)) %>%
  filter(!is.na(max_ma))
tax_range_time(m.dat.occ, name = "Genus", by = "name", plot = TRUE)

#========================== OCCUPANCY COMPARISONS ==============================

# Set resolution and extent
e <<- extent(-180, 180, -90, 90)
res <- 1

# Setup data for get_grid()
m.dat.occ$lat <- as.numeric(m.dat.occ$lat)
m.dat.occ$lng <- as.numeric(m.dat.occ$lng)
m.dat.occ <- m.dat.occ %>%
  filter(!is.na(lat) & !is.na(lng)) %>%
  dplyr::rename(collection_no = Locality)

# Apply get_grid()
get_grid(m.dat.occ, res, e)
get_grid(pb.dat, res, e)
get_grid(pb.all, res, e)

# Find how many unique cells each dataset occupies
length(unique(pb.dat.binned$cells))
length(unique(m.dat.occ.binned$cells))
length(unique(pb.all.binned$cells))

# Compare to see whether the PBDB misses any cells
m.dat.occ.cells <- unique(m.dat.occ.binned$cells)
pb.all.cells <- unique(pb.all.binned$cells)
pb.dat.cells <- unique(pb.dat.binned$cells)

length(setdiff(m.dat.occ.cells, pb.all.cells))
setdiff(pb.dat.cells, pb.all.cells)

#===== Occupancy through time =====
# Loop through stages, providing occupancy estimates for each stage
results <- as.data.frame(matrix(NA, nrow = 0, ncol = 6))
stages$bin_midpoint <- (stages$max_ma+stages$min_ma)/2
for(i in unique(pb.all$bin_assignment)){
  m.cells <- unique(m.dat.occ.binned$cells[which(m.dat.occ.binned$bin_assignment == i)])
  pb.cells <- unique(pb.dat.binned$cells[which(pb.dat.binned$bin_assignment == i)])
  all.cells <- unique(pb.all.binned$cells[which(pb.all.binned$bin_assignment == i)])
  if(length(setdiff(m.cells, all.cells)) > 0){
    add.cells <- length(setdiff(m.cells, all.cells))
    all.cells <- c(m.cells, all.cells)
  }
  results <- rbind(results, 
                       c(i, 
                         stages$bin_midpoint[stages$bin == i], 
                         length(all.cells), 
                         "m.dat",
                         length(m.cells), 
                         (length(m.cells)/length(all.cells))*100), 
                       c(i,
                         stages$bin_midpoint[stages$bin == i], 
                         length(all.cells), 
                         "pb.dat",
                         length(pb.cells), 
                         (length(pb.cells)/length(all.cells))*100
                         )
                       )
}
colnames(results) <- c("bin", "bin_midpoint", "all.cells", "data", "n.cells", "occ")
results$data <- as.factor(results$data)
results$bin_midpoint <- as.numeric(results$bin_midpoint)
results$n.cells <- as.numeric(results$n.cells)
results$occ <- as.numeric(results$occ)
results$all.cells <- as.numeric(results$all.cells)

results <- results %>%
  filter(bin > 55 & bin < 93)

a <- ggplot(results, aes(x=bin_midpoint, y=occ)) + 
  geom_line(aes(colour = data)) +
  scale_x_reverse() +
  theme_classic() +
  ylab("Occupancy (%)") +
  xlab("Time (Ma)")

a <- ggplot(results, aes(x=bin_midpoint, y=all.cells)) + 
  geom_line() +
  scale_x_reverse() +
  theme_classic() +
  ylab("No. of occupied cells") +
  xlab("Time (Ma)")

gggeo_scale(a)

#===== Correlations =====
# Reshape data for tests
wide.results<- dcast(results, bin+bin_midpoint+all.cells~data, value.var = c("n.cells"))
colnames(wide.results)[[4]] <- "m.dat.cells"
colnames(wide.results)[[5]] <- "pb.dat.cells"
wide.results <- cbind(wide.results, dcast(results, bin+bin_midpoint+all.cells~
                                            data, value.var = c("occ"))[,4:5])

# Log data
wide.results$all.cells <- log10(wide.results$all.cells)
wide.results$m.dat.cells <- log10(wide.results$m.dat.cells)
wide.results$pb.dat.cells <- log10(wide.results$pb.dat.cells)
wide.results$m.dat <- log10(wide.results$m.dat)
wide.results$pb.dat <- log10(wide.results$pb.dat)

# Occupancy vs. occupancy
cor.test(wide.results$m.dat, wide.results$pb.dat, method = "spearman")
cor.test(wide.results$pb.dat.cells, wide.results$m.dat.cells, method = "spearman")

# Occupancy vs. Total occ.
cor.test(wide.results$m.dat.cells, wide.results$all.cells, method = "spearman")
cor.test(wide.results$pb.dat.cells, wide.results$all.cells, method = "spearman")

# Occupancy chi squared
chisq.test(wide.results$m.dat.cells, wide.results$all.cells)
chisq.test(wide.results$m.dat.cells, wide.results$pb.dat.cells)

#======================== RANGE SIZE THROUGH TIME ==============================

#======================== DIVERSITY COMPARISON =================================

#======================== TEMPORAL RANGE COMPARISON ============================
