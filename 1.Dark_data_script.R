################################################################################
#                                                                              #
# ECHINOID DARK DATA                                                           #
# C.D. Dean, M. Clapham, T.A.M. Ewin, J. Thompson                              #
# Code by C. D. Dean, started 29.09.23                                         #
#                                                                              #
################################################################################

################################################################################
# 1. SETUP
################################################################################

############################################
###### LOAD AND PREPARE MUSEUM DATASET #####
############################################

# Load setup file
source("0.Setup.R")

# Load required functions
source("0.Functions_Dark_data_test.R")

# Made occurrence based dataset instead of specimen based
m.dat.occ <- m.dat %>%
  dplyr::select(Rank, Family, Genus, Species, Publication, Museum_only, Locality,
                Continent, Country, lat, lng, Latlng_Precision, Group, Formation, 
                Member, Age_resolution, Age, Max_period, Min_period, early_stage,
                late_stage, max_ma, min_ma, interval_mid_ma, n_bins, bin_assignment, 
                bin_midpoint, overlap_percentage, Finalised_lith, Finalised_grainsize) %>%
  dplyr::distinct()

# Split into museum and publication datasets
m.only.dat <- m.dat.occ %>%
  filter(Museum_only == "1")
m.p.dat <- m.dat.occ %>%
  filter(Museum_only =="0")

###########################################
###### LOAD AND PREPARE PBDB DATASETS #####
###########################################
 
##### All Echinoids #####

# Load Paleozoic echinoid dataset
pb.dat <- read.csv("pbdb_data_09052023.csv", skip = 19)

# Assign occurrences to bins
pb.dat <- bin_time(pb.dat, bins = stages, method = 'majority')

# Covert country code into full words and sort continents
pb.dat$Country <- countrycode(sourcevar = pb.dat$cc, "eurostat", "country.name")
pb.dat$Country[is.na(pb.dat$Country) == TRUE] <- "Greece"
pb.dat$Continent <- countrycode(sourcevar = pb.dat$Country, "country.name", "continent")
pb.dat$Continent[which(pb.dat$Country == "United States" | pb.dat$Country== "Canada")] <- "North America"
pb.dat$Continent[which(pb.dat$Continent == "Americas")] <- "South America"

##### All occurrences #####

# Load all palaeozoic occurrences
pb.all <- read.csv("all_palaeozoic_binned.csv")

# Assign occurrences to bins (WARNING: Takes a looooong time)
#pb.all <- read.csv("all_palaeozoic.csv", skip = 18)
#pb.all <- pb.all %>%
#  filter(max_ma < stages[nrow(stages), 2])
#pb.all <- bin_time(pb.all, bins = stages, method = 'majority')
#write.csv(pb.all, file = "all_palaeozoic_binned.csv")

#####################################
###### PREPARE RELATED DATASETS #####
#####################################

##### CARB/SILIC COLLECTIONS: MACROSTRAT #####

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
  filter(max_ma < 470) %>%
  filter(min_ma > 252)
sili.macro <- sili.macro %>%
  filter(max_ma < 470) %>%
  filter(min_ma > 252)

# Bin data
carb.macro <- bin_time(carb.macro, stages, method = "all")
sili.macro <- bin_time(sili.macro, stages, method = "all")

carb.macro  <- carb.macro  %>%
  dplyr::group_by(bin_midpoint) %>%
  dplyr::summarise(count = n(), lith = "carb") 
sili.macro <- sili.macro  %>%
  dplyr::group_by(bin_midpoint) %>%
  dplyr::summarise(count = n(), lith = "sili") 

carb.macro <- carb.macro  %>%
  dplyr::group_by(bin_midpoint) %>%
  dplyr::summarise(count = sum(col_area), lith = "carb") 
sili.macro <- sili.macro  %>%
  dplyr::group_by(bin_midpoint) %>%
  dplyr::summarise(count = sum(col_area), lith = "sili") 

macro <- rbind(carb.macro, sili.macro)

a <- ggplot(macro, aes(x=bin_midpoint, y=count)) + 
  geom_line(aes(colour = lith)) +
  scale_x_reverse() +
  theme_bw() +
  ylab("Jaccard Similarity") +
  xlab("Time (Ma)") +
  geom_point(aes(color = lith))
a

################################################################################
# 2. COMPARISONS BETWEEN OUR RECORDS AND PBDB
################################################################################

# Tests to run -
# How does taxonomy differ? total numbers, percentage 'incorrect'. Percentage at
# higher taxonomic grade than possible (e.g. species rather than genus)
# Any other major differences?

##### TEMPORAL UNCERTAINTY #####

# Make temporal range column
pb.dat$temp_range <- pb.dat$max_ma - pb.dat$min_ma
m.only.dat$temp_range <- m.only.dat$max_ma - m.only.dat$min_ma

temp.compare <- data.frame(group = "Museum", vals = m.only.dat$temp_range)
temp.compare <- rbind(temp.compare, data.frame(group = "PBDB", vals = pb.dat$temp_range))

# Plot
temp.compare %>%
  ggplot( aes(x=group, y=vals, fill=group)) +
  geom_boxplot() 

# Mann Whitney U test
wilcox.test(vals ~ group, data = temp.compare)

##### DATA COMPLETENESS #####
# Lithology
# Formation
# Age
# Location precision

################################################################################
# 3. BASIC COMPARISONS
################################################################################

##### PRESERVATION SCORE #####

# Format into table and run Chi Squared test
test <- chisq.test(table(m.dat$Preservation_score, m.dat$Museum_only))
test
test$expected # compare against what the test would have expected

# Mosaic plot
mosaic(~ Preservation_score + Museum_only,
       direction = c("v", "h"),
       data = m.dat,
       labeling_args = list(just_labels = "right", 
                            set_varnames = c(Preservation_score = "Preservation Score", 
                                             Finalised_grainsize = "Grain size")),
       shade = TRUE
)

##### COMPARATIVE MAPS #####
# Set extent
e <<- extent(-180, 180, -90, 90)

# Make maps
get_grid_im(m.dat.occ, 1,  "Museum Echinoids", ext = e)
get_grid_im(pb.dat, 1,  "PBDB Echinoids", ext = e)
get_grid_im(pb.all, 1,  "PBDB", ext = e)

# Difference map (think about this more)
xy <- cbind(as.double(m.dat.occ$lng), as.double(m.dat.occ$lat))
xy <- unique(xy)
r <- raster::raster(ext = e, res = 1)
r.m <- raster::rasterize(xy, r, fun = 'count')

xy <- cbind(as.double(pb.dat$lng), as.double(pb.dat$lat))
xy <- unique(xy)
r <- raster::raster(ext = e, res = 1)
r.p <- raster::rasterize(xy, r, fun = 'count')

countries <- maps::map("world", plot=FALSE, fill = TRUE) # find map to use as backdrop
countries <<- maptools::map2SpatialPolygons(countries, IDs = countries$names, proj4string = CRS("+proj=longlat")) # Turn map into spatialpolygons
mapTheme <- rasterVis::rasterTheme(region=viridis(8))
print(rasterVis::levelplot((r.m - r.p), margin=F, par.settings=mapTheme,  main = paste("Total ", (substitute(name)), " per Grid Cell", sep = "")) + #create levelplot for raster
        #   latticeExtra::layer(sp.polygons(states, col = "white", fill = NA), under = T)  + # Plots state lines
        latticeExtra::layer(sp.polygons(countries, col = 0, fill = "light grey"), under = T)) # Plots background colour

##### TABLE COMPARISONS AND GRAPHS #####

# Formations - how many are shared/different?
table(pb.dat$formation)
table(m.dat$Formation)
length(unique(m.dat$Formation))
length(unique(m.only.dat$Formation))
length(unique(m.p.dat$Formation))
length(unique(pb.dat$formation))

# Country - How many shared/different?
table(pb.dat$Country)
table(m.dat$Country)
country.long <- rbind(cbind(as.data.frame(table(m.dat$Country)), data = c(rep("M", nrow(as.data.frame(table(m.dat$Country)))))), 
      cbind(as.data.frame(table(pb.dat$Country)), data = c(rep("P", nrow(as.data.frame(table(pb.dat$Country)))))))

# Continent - % difference through time, % of that dataset through time.
table(m.dat$Continent)
table(pb.dat$Continent)
cont.long <- rbind(cbind(as.data.frame(table(m.dat$Continent)), data = c("M", "M", "M", "M")),
      cbind(as.data.frame(table(pb.dat$Continent)), data = c("P", "P", "P", "P", "P", "P")))

# Per lithology - % difference through time, % of dataset through time.
table(pb.dat$Lithology)
table(m.dat.occ$Finalised_lith)

################################################################################
# 4. SPATIAL COMPARISONS
################################################################################

##### OCCUPANCY COMPARISON #####

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
gggeo_scale(a)

a <- ggplot(results, aes(x=bin_midpoint, y=all.cells)) + 
  geom_line() +
  geom_line(aes(x = bin_midpoint, y = n.cells, color = data)) +
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
chisq.test(wide.results$all.cells, wide.results$pb.dat.cells)
chisq.test(wide.results$m.dat.cells, wide.results$pb.dat.cells)

##### SPATIAL RANGE COMPARISON #####

# Set up museum data
m.dat.occ.gen <- m.dat.occ %>%
  filter(is.na(Genus) == F)

# Calculate spatial range for PBDB genera
p.range <- tax_range_space(pb.dat, name = "genus", lng = "lng", lat = "lat", method = "con")

# Calculate spatial range for Museum genera
m.range <- tax_range_space(m.dat.occ.gen, name = "Genus", lng = "lng", lat = "lat", method = "con")

# Comparison tests/difference between genera
p.range <- p.range %>%
  filter(taxon %in% intersect(unique(m.range$taxon), unique(p.range$taxon))) %>%
  dplyr::select(taxon, area) %>%
  distinct()
m.range <- m.range %>%
  filter(taxon %in% intersect(unique(m.range$taxon), unique(p.range$taxon))) %>%
  dplyr::select(taxon, area) %>%
  distinct()

chisq.test(p.range$area, m.range$area)

#==== Comparison between families
# Set up museum data
m.dat.occ.fam <- m.dat.occ %>%
  filter(is.na(Family) == F)

# Calculate spatial range for PBDB genera
p.range.fam <- tax_range_space(pb.dat, name = "family", lng = "lng", lat = "lat", method = "con")

# Calculate spatial range for Museum genera
m.range.fam <- tax_range_space(m.dat.occ.fam, name = "Family", lng = "lng", lat = "lat", method = "con")

# Comparison tests/difference between genera
p.range.fam <- p.range.fam %>%
  filter(taxon %in% intersect(unique(m.range.fam$taxon), unique(p.range.fam$taxon))) %>%
  dplyr::select(taxon, area) %>%
  distinct() %>%
  mutate(data = "PBDB")
m.range.fam <- m.range.fam %>%
  filter(taxon %in% intersect(unique(m.range.fam$taxon), unique(p.range.fam$taxon))) %>%
  dplyr::select(taxon, area) %>%
  distinct() %>%
  mutate(data = "Museum")

chisq.test(m.range.fam$area, p.range.fam$area)
m.range.fam$area - p.range.fam$area

#===== Comparison through time =====
# Apply group apply with tax_range_space
pb.range <- palaeoverse::group_apply(pb.dat,
                                 group = c("bin_assignment"),
                                 fun = tax_range_space,
                                 name = "genus",
                                 lng = "lng",
                                 lat = "lat",
                                 method = "con")

m.range <- palaeoverse::group_apply(m.dat.occ.gen,
                                 group = c("bin_assignment"),
                                 fun = tax_range_space,
                                 name = "Genus",
                                 lng = "lng",
                                 lat = "lat",
                                 method = "con")
test <- pb.range %>%
  dplyr::select(taxon, area, bin_assignment) %>%
  distinct() %>%
  filter(area != 0, taxon != "") %>%
  mutate(data = "PBDB")

test2 <- m.range %>%
  dplyr::select(taxon, area, bin_assignment) %>%
  distinct() %>%
  filter(area != 0) %>%
  mutate(data = "Museum")

range.res <- bind_rows(filter(test2, taxon %in% intersect(unique(test$taxon), unique(test2$taxon))),
          filter(test, taxon %in% intersect(unique(test$taxon), unique(test2$taxon))))

range <- (dcast(range.res, taxon + bin_assignment ~ data, value.var = "area"))
range <- range[order(bin_assignment),]
range$abs.difference <- log10(range$Museum - range$PBDB)
range$perc.difference <- (range$PBDB/range$Museum)

testy <- range %>%
  filter(taxon == "Archaeocidaris")

plot(testy$bin_assignment, testy$abs.difference)

################################################################################
# 5. TEMPORAL COMPARISONS
################################################################################

##### TEMPORAL RANGE COMPARISON #####

# Calculate temporal range for PBDB genera
tax_range_time(pb.dat, name = "genus", by = "name", plot = TRUE)

# Calculate temporal range for Museum genera
m.dat.occ <- m.dat.occ %>%
  filter(!is.na(Genus)) %>%
  filter(!is.na(max_ma))
tax_range_time(m.dat.occ, name = "Genus", by = "name", plot = TRUE)

# Comparison tests

##### DIVERSITY COMPARISON #####
# Species
test <- divDyn(m.dat.occ, tax = "Species", bin = "bin_assignment")
test <- subsample(m.dat.occ,iter=100, q=0.6, tax="Species", 
                  bin="bin_assignment", coll = 'collection_no', output="dist", 
                  type="sqs", duplicates = TRUE, useFailed = TRUE)
# Genera

# Families

# Spatial subsampling?
# Difference in amount of data available... might impact feasibility of the analysis itself.

################################################################################
# 6. MODELLING COMPARISONS
################################################################################

# Get climate/sampling variables
table(pb.all.binned$lithdescript)
