################################################################################
######## TAPHONOMIC CONTROLS ON A MULTI-ELEMENT SKELETAL FOSSIL RECORD #########
################################################################################

# Jeffrey R. Thompson, Christopher D. Dean, Madeline Ford, Timothy A. M. Ewin
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
m.dat <- bin_time(occdf = m.dat, bins = stages, method = 'majority')

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
m.dat.period <- bin_time(occdf = m.dat, bins = periods, method = 'majority')

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
m.dat.series <- palaeoverse::bin_time(occdf = m.dat, bins = series, method = 'majority')
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
carb.macro <- bin_time(occdf = carb.macro, bins = stages, method = "all")
sili.macro <- bin_time(occdf = sili.macro, bins = stages, method = "all")

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
carb.macro.period <- bin_time(occdf = carb.macro, bins = series, method = "majority")
sili.macro.period <- bin_time(occdf = sili.macro, bins = series, method = "majority")

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