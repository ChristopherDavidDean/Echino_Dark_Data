################################################################################
### MUSEUM `DARK DATA` ILLUMINATES THE EVOLUTIONARY HISTORY OF FOSSIL GROUPS ###
################################################################################

# Christopher D. Dean & Jeffrey R. Thompson
# 2024
# Script written by Christopher D. Dean

################################################################################
#                              FILE 2: SETUP                                   #
################################################################################

################################################################################
# 1. SETUP
################################################################################

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data
m.dat <- read.csv("Specimen_data/Final_Database_for_analysis.csv")
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
names(stages)[names(stages) == "max_age"] <- "max_ma"
names(stages)[names(stages) == "min_age"] <- "min_ma"
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
names(period)[names(period) == "Max_period"] <- "max_ma"
names(period)[names(period) == "Min_period"] <- "min_ma"
period <- look_up(occdf = period, early_interval = "max_ma", late_interval = "min_ma")

# Assign numerical ages to stages
names(stage)[names(stage) == "max_stage"] <- "max_ma"
names(stage)[names(stage) == "min_stage"] <- "min_ma"
stage <- look_up(occdf = stage, early_interval = "max_ma", late_interval = "min_ma")

# Reorganise column headings
stage <- stage[,-c(31:32)]
period <- period[,-c(31:32)]
colnames(period)[29:30] <- c("Max_period", "Min_period")

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

# Load periods and rename columns for binning
data(periods)
colnames(periods)[1] <- "bin"
colnames(periods)[2] <- "max_ma"
colnames(periods)[3] <- "min_ma"


##### PERIOD LEVEL TIME #####

# Bin into periods
m.dat.period <- bin_time(m.dat, bins = periods, method = 'majority')

# Create factors
order_ind <- c("Permian", "Carboniferous", "Devonian", "Silurian", "Ordovician")
m.dat.period$bin_assignment <- as.factor(m.dat.period$bin_assignment)
m.dat.period$bin_assignment <- factor(m.dat.period$bin_assignment, levels = order_ind)

# Load period colour information
df <- palaeoverse::GTS2020
df <- df[155:159,]
myColours <- df$colour

# Assign colours
names(myColours) <- levels(m.dat.period$bin_assignment)
custom_colours <- scale_colour_manual(name = "bin_assignment", values = myColours[1:5])

##### SERIES LEVEL TIME #####

# Load series data
series <- read.csv("Additional_data/series.csv")

# Bin into series
order_ind <- c("Permian", "Pennsylvanian", "Mississippian", "Devonian", "Silurian", "Ordovician")
m.dat.series <- palaeoverse::bin_time(m.dat, bins = series, method = 'majority')
m.dat.series$bin_assignment <- as.factor(m.dat.series$bin_assignment)
m.dat.series$bin_assignment <- factor(m.dat.series$bin_assignment, levels = order_ind)

# Assign colours
myColours <- series$color
names(myColours) <- levels(m.dat.series$bin_assignment)
custom_colours <- scale_colour_manual(name = "bin_assignment", values = myColours[1:6])

# Adjust names for gggeo_scale
series2 <- series %>%
  dplyr::rename(name = bin, 
                max_age = max_ma, 
                min_age = min_ma)

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
# 4. MACROSTRAT SETUP
################################################################################

# Read in data
carb.macro <- read.csv('https://macrostrat.org/api/v2/units?lith_type=carbonate&environ_class=marine&response=long&format=csv', stringsAsFactors = FALSE)
sili.macro <- read.csv('https://macrostrat.org/api/v2/units?lith_type=siliciclastic&environ_class=marine&response=long&format=csv', stringsAsFactors = FALSE)

# Change column names
names(carb.macro)[names(carb.macro) == 't_age'] <- "min_ma"
names(carb.macro)[names(carb.macro) == 'b_age'] <- "max_ma"
names(sili.macro)[names(sili.macro) == 't_age'] <- "min_ma"
names(sili.macro)[names(sili.macro) == 'b_age'] <- "max_ma"

# Filter to relevant data
carb.macro <- carb.macro %>%
  filter(max_ma < 485.41) %>%
  filter(min_ma > 251.901)
sili.macro <- sili.macro %>%
  filter(max_ma < 485.40) %>%
  filter(min_ma > 251.901)

##### STAGE #####

# Bin data
carb.macro <- bin_time(carb.macro, stages, method = "all")
sili.macro <- bin_time(sili.macro, stages, method = "all")

# Count
carb.macro.count  <- carb.macro  %>%
  dplyr::group_by(bin_midpoint) %>%
  dplyr::summarise(count = n(), lith = "carb") 
sili.macro.count <- sili.macro  %>%
  dplyr::group_by(bin_midpoint) %>%
  dplyr::summarise(count = n(), lith = "sili") 

# Col_area
carb.macro.area <- carb.macro  %>%
  dplyr::group_by(bin_midpoint) %>%
  dplyr::summarise(count = sum(col_area), lith = "carb") 
sili.macro.area <- sili.macro  %>%
  dplyr::group_by(bin_midpoint) %>%
  dplyr::summarise(count = sum(col_area), lith = "sili") 

macro.area <- rbind(carb.macro.area, sili.macro.area)
macro.count <- rbind(carb.macro.count, sili.macro.count)

##### PERIOD #####

# Bin data
carb.macro.period <- bin_time(carb.macro, series, method = "majority")
sili.macro.period <- bin_time(sili.macro, series, method = "majority")

# Count
carb.macro.count.period  <- carb.macro.period  %>%
  dplyr::group_by(bin_midpoint) %>%
  dplyr::summarise(count = n(), lith = "carb") 
sili.macro.count.period <- sili.macro.period  %>%
  dplyr::group_by(bin_midpoint) %>%
  dplyr::summarise(count = n(), lith = "sili") 

# Col_area
carb.macro.area.period <- carb.macro.period  %>%
  dplyr::group_by(bin_midpoint) %>%
  dplyr::summarise(count = sum(col_area), lith = "carb") 
sili.macro.area.period <- sili.macro.period  %>%
  dplyr::group_by(bin_midpoint) %>%
  dplyr::summarise(count = sum(col_area), lith = "sili") 

macro.area.period <- rbind(carb.macro.area.period, sili.macro.area.period)
macro.count.period <- rbind(carb.macro.count.period, sili.macro.count.period)


################################################################################
# 5. SETUP FOR ANALYSIS
################################################################################

###########################################
###### CREATE BASIC METADATA DATASETS #####
###########################################

##### LITHOLOGY #####

# Remove data without lithological info
l.m.dat <- m.dat %>%
  filter(is.na(Finalised_lith) == F)

##### GRAIN SIZE #####

# Make dataset for grain size
g.m.dat <- m.dat %>%
  filter(is.na(Finalised_grainsize) == F) 

# Assign specific grain size categories into either "Fine Grained" or "Coarse Grained"
g.m.dat <- simple.grain(g.m.dat)

# Remove remaining specimens without grain size
g.m.dat <- g.m.dat %>%
  filter(is.na(Finalised_grainsize) == F) 

##### FAMILY #####

# Make plot of preservation scores against family, NA values removed
f.m.dat <- m.dat %>%
  filter(is.na(Family) == F) %>%
  filter(Family != 'Triadotiaridae') %>%
  filter(Family != 'Cravenechinidae') %>%
  filter(Family != 'Archaeocidaridae or miocidaridae')

################################################
###### BEGIN SETUP OF COMPARATIVE DATASETS #####
################################################

# Set extent and resolution
e <- extent(-180, 180, -90, 90)
res <- 1

# Add midpoint to stages
stages$bin_midpoint <- (stages$max_ma+stages$min_ma)/2

# Made occurrence based dataset instead of specimen based
m.dat.occ <- m.dat.rotate %>%
  dplyr::select(Rank, Family, Genus, Species, Museum_code, Museum_country,  
                Publication, Museum_only, In_PBDB, Locality, Continent, Country, 
                lat, lng, Latlng_Precision, Group, Formation, Member, Age_resolution, 
                Max_period, Min_period, early_stage, late_stage, max_ma, min_ma, 
                interval_mid_ma, n_bins, bin_assignment, bin_midpoint, 
                overlap_percentage, Finalised_lith, Finalised_grainsize, 
                p_lng, p_lat) %>%
  dplyr::distinct()

### Split into museum and publication datasets ###
# Museum dataset
m.only.dat <- m.dat.occ %>%
  filter(Museum_only == "1")
# Publication dataset 
m.p.dat <- m.dat.occ %>%
  filter(Museum_only =="0")

###########################################
###### LOAD AND PREPARE PBDB DATASETS #####
###########################################

##### All Echinoids #####

# Load Paleozoic echinoid dataset
pb.dat <- read.csv("Specimen_data/pbdb_data_21092023.csv", skip = 19)

# Assign occurrences to bins
pb.dat <- bin_time(pb.dat, bins = stages, method = 'majority')

# Palaeorotate
pb.dat <- palaeorotate(pb.dat, 
                       lng = 'lng', 
                       lat = 'lat', 
                       age = "bin_midpoint",
                       model = "PALEOMAP",
                       method = "point")

# Covert country code into full words and sort continents
pb.dat$Country <- countrycode(sourcevar = pb.dat$cc, "eurostat", "country.name")
pb.dat$Country[is.na(pb.dat$Country) == TRUE] <- "Greece"
pb.dat$Continent <- countrycode(sourcevar = pb.dat$Country, 
                                "country.name", "continent")
pb.dat$Continent[which(pb.dat$Country == "United States" | 
                         pb.dat$Country== "Canada")] <- "North America"
pb.dat$Continent[which(pb.dat$Continent == "Americas")] <- "South America"

##### All occurrences #####

# Load all palaeozoic occurrences
pb.all <- read.csv("Specimen_data/all_palaeozoic_binned.csv")
# Load all palaeozoic echinodermata occurrences
pb.echino <- read.csv("Specimen_data/Echinodermata.csv", skip = 19)
pb.echino <- pb.echino %>%
  dplyr::filter(max_ma < 485.4000) %>%
  dplyr::filter(min_ma > 251.2000)
pb.echino <- bin_time(pb.echino, bins = stages, method = 'majority')

# Assign occurrences to bins (WARNING: Takes a looooong time)
#pb.all <- read.csv("Specimen_data/all_palaeozoic.csv", skip = 18)
#pb.all <- pb.all %>%
#  filter(max_ma < stages[nrow(stages), 2])
#pb.all <- bin_time(pb.all, bins = stages, method = 'majority')
#write.csv(pb.all, file = "Specimen_data/all_palaeozoic_binned.csv")

###########################################
###### FINALISE DATASETS FOR ANALYSIS #####
###########################################

# Museum dataset
m.only.dat <- m.only.dat %>%
  mutate(Combined_name = case_when(
    Rank == "Species" ~ paste(Genus, Species, sep = " "), 
    Rank == "Genus" ~ paste(Genus))) %>%
  mutate(In_PBDB = NA)

# PBDB only dataset
pb.only.dat <- pb.dat %>%
  dplyr::filter(In_our_dataset == "No") %>%
  dplyr::rename(Rank = accepted_rank,
                Combined_name = accepted_name,
                Family = family,
                Genus = genus,
                Species = species_name, 
                Museum_code = museum,
                Publication = primary_reference,
                Description = occurrence_comments,
                Locality = collection_no,
                Latlng_Precision = latlng_precision,
                Group = stratgroup,
                Formation = formation,
                Member = member,
                Formation_Notes = stratcomments,
                Finalised_lith = Lithology, 
                Finalised_grainsize = Grain_size) %>%
  dplyr::select(Rank, Combined_name, Family, Genus, Species, Museum_code,
                Museum_country, Publication, Description, Locality, Continent,  
                Country, lat, lng, Latlng_Precision, Group, Formation, 
                Member, Formation_Notes, max_ma, min_ma, p_lng, p_lat,
                Finalised_grainsize, Finalised_lith, In_our_dataset, 
                id, n_bins, bin_assignment, overlap_percentage) 

pb.only.dat$Age_resolution <- NA
pb.only.dat$Max_period <- NA
pb.only.dat$Min_period <- NA
pb.only.dat$early_stage <- NA
pb.only.dat$late_stage <- NA
pb.only.dat$interval_mid_ma <- (pb.only.dat$max_ma + pb.only.dat$min_ma)/2
pb.only.dat$In_PBDB <- "Yes"
pb.only.dat$Museum_only <- 0

# Published record dataset
m.p.dat <- m.p.dat %>%
  dplyr::mutate(Combined_name = case_when(
    Rank == "Species" ~ paste(Genus, Species, sep = " "), 
    Rank == "Genus" ~ paste(Genus))) %>%
  dplyr::select(-c(bin_midpoint))
pb.only.dat$Locality <- as.character(pb.only.dat$Locality)
pb.only.dat <- pb.only.dat %>%
  dplyr::select(-c(Description, Formation_Notes, In_our_dataset, id))
pub.all.dat <- rbind(pb.only.dat, m.p.dat)

# Complete dataset
all.dat <- m.only.dat %>%
  dplyr::select(-c(bin_midpoint))
all.dat <- rbind(all.dat, pub.all.dat)

# Fix grainsize with function
m.only.dat <- simple.grain(m.only.dat)
pb.only.dat <- simple.grain(pb.only.dat)
pub.all.dat <- simple.grain(pub.all.dat)
all.dat <- simple.grain(all.dat)

# Bin according to series level data
pb.only.series <- bin_time(pb.only.dat, bins = series, method = 'majority')
m.only.series <- bin_time(m.only.dat, bins = series, method = 'majority')
pub.all.series <- bin_time(pub.all.dat, bins = series, method = 'majority')
all.series <- bin_time(all.dat, bins = series, method = 'majority')

####################
###### SUMMARY #####
####################

# Museum only
head(m.only.dat)
# PBDB only
head(pb.only.dat)
# All published record
head(pub.all.dat)
# All data
head(all.dat)

combined_data <- list(m.only.dat, pb.only.dat, pub.all.dat, all.dat)
