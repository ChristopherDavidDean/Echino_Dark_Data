################################################################################
######## TAPHONOMIC CONTROLS ON A MULTI-ELEMENT SKELETAL FOSSIL RECORD #########
################################################################################

# Jeffrey R. Thompson, Christopher D. Dean, Madeline Ford, Timothy A. M. Ewin
# 2024
# Script written by Christopher D. Dean

################################################################################
#                            FILE 3: MAIN SCRIPT                               #
################################################################################

################################################################################
# 1. SETUP
################################################################################

# Load functions file
source("1.Taph_Functions.R")

# Load setup file
source("2.Taph_Setup.R")

################################################################################
# 2. BAR PLOTS AND CHI-SQUARED TESTS
################################################################################

##########################
##### TAXONOMIC RANK #####
##########################

# Make bar plot
a <- make.bar.plot(m.dat, m.dat$Rank, "Taxonomic Rank", colour = "Zissou1", FALSE)

# Make proportional bar plot
b <- make.prop.bar.plot(m.dat, m.dat$Rank, "Taxonomic Rank", colour = "Zissou1", FALSE)

(p1 <- ggarrange(a, b, 
                align='hv',
                labels = c("A", "B"),
                nrow = 1, 
                ncol = 2, 
                common.legend = T, 
                legend = "bottom"))

# Format into table and run Chi Squared test
tax.res <- chisq.test(table(m.dat$Preservation_score, m.dat$Rank))
tax.res
tax.res$expected # compare against what the test would have expected

# Mosaic plot
vcd::mosaic(~ Preservation_score + Rank,
       direction = c("v", "h"),
       data = m.dat,
       labeling_args = list(
         set_varnames = c(Preservation_score = "Preservation Score", 
                          Rank = "Taxonomic Rank")),
       shade = TRUE
)

#####################
##### LITHOLOGY #####
#####################

# Make bar plot
make.bar.plot(l.m.dat, l.m.dat$Finalised_lith, "Lithology", colour = "Zissou1", FALSE)

# Make proportional bar plot
make.prop.bar.plot(l.m.dat, l.m.dat$Finalised_lith, "Lithology", colour = "Zissou1", FALSE)

# Format into table and run Chi Squared test
res.lith <- chisq.test(table(l.m.dat$Preservation_score, l.m.dat$Finalised_lith))
res.lith
res.lith$expected # compare against what the test would have expected

# Mosaic plot
vcd::mosaic(~ Preservation_score + Finalised_lith,
       direction = c("v", "h"),
       data = l.m.dat,
       labeling_args = list(
         set_varnames = c(Preservation_score = "Preservation Score", 
                          Finalised_lith = "Lithology")),
       shade = TRUE
)

######################
##### GRAIN SIZE #####
######################

# Make bar plot
make.bar.plot(g.m.dat, g.m.dat$Finalised_grainsize_simp, "Grain Size", 
              colour = "Darjeeling2", flip = FALSE)

# Make proportional bar plot
make.prop.bar.plot(g.m.dat, g.m.dat$Finalised_grainsize_simp, "Grain Size", 
                   colour = "Darjeeling2", flip = FALSE)

# Format into table and run Chi Squared test
res.grain <- chisq.test(table(g.m.dat$Preservation_score, g.m.dat$Finalised_grainsize))
res.grain
res.grain$expected # compare against what the test would have expected

# Mosaic plot
vcd::mosaic(~ Preservation_score + Finalised_grainsize_simp,
       direction = c("v", "h"),
       data = g.m.dat,
       labeling_args = list(just_labels = "right", 
                            set_varnames = c(Preservation_score = "Preservation Score", 
                                             Finalised_grainsize_simp = "Grain size")),
       shade = TRUE
)

##################
##### FAMILY #####
##################

ggplot(f.m.dat) +
  aes(x = Preservation_score, fill = Family) +
  geom_bar() +
  ylab("Frequency") +
  xlab("Preservation Score") +
  scale_fill_manual(values=c(wes_palette("Royal1"), 
                             wes_palette("Zissou1"), 
                             wes_palette("Royal2"))) +
  theme_bw() 

ggplot(f.m.dat) +
  aes(x = Preservation_score, fill = Family) +
  geom_bar(position = 'fill') +
  ylab("Frequency") +
  ylab("Proportion of total") +
  xlab("Preservation Score") +
  scale_fill_manual(values=c(wes_palette("Royal1"), 
                             wes_palette("Zissou1"), 
                             wes_palette("Royal2"))) +
  theme_bw() +
  scale_y_continuous(expand = c(0,0), limits = c(0,1))

ggplot(f.m.dat) +
  aes(x = Family, fill = as.factor(Preservation_score)) +
  geom_bar(position = 'fill') +
  ylab("Proportion of total") +
  xlab("Family") +
  coord_flip() +
  labs(fill="Preservation score") +
  scale_fill_manual(values=(wes_palette("Zissou1"))) +
  theme_bw() +
  scale_y_continuous(expand = c(0,0), limits = c(0,1))

# Format into table and run Chi-Squared test
res.fam <- chisq.test(table(f.m.dat$Preservation_score, f.m.dat$Family))
res.fam
res.fam$expected # compare against what the test would have expected

vcd::mosaic(~ Preservation_score + Family,
       direction = c("v", "h"),
       data = f.m.dat,
       labeling_args = list(just_labels = "right", 
                            set_varnames = c(Preservation_score = "Preservation Score"), 
                            rot_labels = c(0,0)),
       shade = TRUE
)

################################################################################
# 3. TEMPORAL PATTERNS
################################################################################

########################################
##### PRESERVATION PER TIME PERIOD #####
########################################

##### PERIOD LEVEL #####

# Plot preservation score per time period (discreet)
ggplot(m.dat.period, aes(x = Preservation_score, fill = bin_assignment)) +
  geom_bar() +
  scale_fill_manual("Legend", values = myColours) +
  ylab("Frequency") +
  xlab("Preservation Score") +
  guides(fill=guide_legend(title="Period")) +
  theme_bw() 

# Plot proportional preservation score per time period (discreet)
ggplot(m.dat.period) +
  aes(x = Preservation_score, fill = bin_assignment) +
  geom_bar(position = 'fill') +
  scale_fill_manual("Legend", values = myColours) +
  ylab("Proportion of total") +
  guides(fill=guide_legend(title="Period")) +
  xlab("Preservation Score") +
  theme_bw() +
  scale_y_continuous(expand = c(0,0), limits = c(0,1))

# Plot preservation scores through time (continuous)
a <- as.data.frame(table(m.dat.period$Preservation_score, 
                         m.dat.period$bin_assignment))
names(a) <- c("Preservation_score", "bin", "Freq")
a <- merge(a, periods, by = 'bin')
a$mid_ma <- (a$max_ma +a$min_ma)/2
a$Preservation_score <- factor(a$Preservation_score, 
                               levels = c("5", "4", "3", "2", "1"))
(a <- ggplot(a, aes(x=mid_ma, y=Freq, fill=Preservation_score)) + 
  geom_area() +
  scale_x_reverse() +
  theme_bw() +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  guides(fill=guide_legend(title="Preservation score")) +
  ylab("Frequency") +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Zissou1"))))

# Plot proportional preservation scores through time (continuous)
a <- as.data.frame(table(m.dat.period$Preservation_score, 
                         m.dat.period$bin_assignment))
names(a) <- c("Preservation_score", "bin", "Freq")
a <- merge(a, periods, by = 'bin')
a$mid_ma <- (a$max_ma +a$min_ma)/2
a$Preservation_score <- factor(a$Preservation_score, 
                               levels = c("5", "4", "3", "2", "1"))
(a <- ggplot(a, aes(x=mid_ma, y=Freq, fill=Preservation_score)) + 
  geom_area(position = 'fill') +
  scale_x_reverse() +
  theme_bw() +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Preservation score")) +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Zissou1"))))

# Format into table and run Chi Squared test
time.res <- chisq.test(table(m.dat.period$Preservation_score, 
                             m.dat.period$bin_assignment))
time.res
time.res$expected # compare against what the test would have expected

# Mosaic plot
vcd::mosaic(~ Preservation_score + bin_assignment,
       direction = c("v", "h"),
       data = m.dat.period,
       shade = TRUE, 
       labeling_args = list(rot_labels = c(0), just_labels = "right", 
                            set_varnames = c(Preservation_score = "Preservation Score",
                                             bin_assignment = ""))
)

##### STAGE LEVEL #####

# Preservation Score
test <- as.data.frame(table(m.dat$Preservation_score, m.dat$bin_midpoint))
names(test) <- c("Preservation_score", "mid_ma", "Freq")
test$mid_ma <- as.numeric(as.character(test$mid_ma))
test$Preservation_score <- factor(test$Preservation_score, 
                                  levels = c("5", "4", "3", "2", "1"))
(a <- ggplot(test, aes(x=mid_ma, y=Freq, fill=Preservation_score)) + 
  geom_area() +
  scale_x_reverse() +
  theme_bw() +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  guides(fill=guide_legend(title="Preservation score")) +
  ylab("Frequency") +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Zissou1"))))

# Preservation Score (proportional)
a <- as.data.frame(table(m.dat.period$Preservation_score, m.dat$bin_midpoint))
names(a) <- c("Preservation_score", "mid_ma", "Freq")
a$mid_ma <- as.numeric(as.character(a$mid_ma))
a$Preservation_score <- factor(a$Preservation_score, 
                               levels = c("5", "4", "3", "2", "1"))
(a <- ggplot(a, aes(x=mid_ma, y=Freq, fill=Preservation_score)) + 
  geom_area(position = 'fill') +
  scale_x_reverse() +
  theme_bw() +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Preservation score")) +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Zissou1"))))

# Lithology
test <- as.data.frame(table(l.m.dat$Finalised_lith, l.m.dat$bin_midpoint))
names(test) <- c("Lithology", "mid_ma", "Freq")
test$mid_ma <- as.numeric(as.character(test$mid_ma))
(a <- ggplot(test, aes(x=mid_ma, y=Freq, fill=Lithology)) + 
  geom_area() +
  scale_x_reverse() +
  theme_bw() +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  guides(fill=guide_legend(title="Lithology")) +
  ylab("Frequency") +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Zissou1"))))

# Grain size
test <- as.data.frame(table(g.m.dat$Finalised_grainsize_simp, 
                            g.m.dat$bin_midpoint))
names(test) <- c("Finalised_grainsize", "mid_ma", "Freq")
test$mid_ma <- as.numeric(as.character(test$mid_ma))
(a <- ggplot(test, aes(x=mid_ma, y=Freq, fill=Finalised_grainsize)) + 
  geom_area() +
  scale_x_reverse() +
  theme_bw() +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  guides(fill=guide_legend(title="Grain size")) +
  ylab("Frequency") +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Moonrise1"))))

# Family
test <- as.data.frame(table(f.m.dat$Family, f.m.dat$bin_midpoint))
names(test) <- c("Family", "mid_ma", "Freq")
test$mid_ma <- as.numeric(as.character(test$mid_ma))
(a <- ggplot(test, aes(x=mid_ma, y=Freq, fill=Family)) + 
  geom_area() +
  scale_x_reverse() +
  theme_bw() +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  guides(fill=guide_legend(title="Family")) +
  ylab("Frequency") +
  xlab("Time (Ma)") +
  scale_fill_manual(values=c(wes_palette("Royal1"), 
                             wes_palette("Zissou1"), 
                             wes_palette("Royal2"))))

##### SERIES LEVEL #####

# Plot preservation score per time period (discreet)
ggplot(m.dat.series, aes(x = Preservation_score, fill = bin_assignment)) +
  geom_bar() +
  scale_fill_manual("Legend", values = myColours) +
  ylab("Frequency") +
  xlab("Preservation Score") +
  guides(fill=guide_legend(title="Period")) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  guides(alpha = "none")

# Plot proportional preservation score per time period (discreet)
ggplot(m.dat.series) +
  aes(x = Preservation_score, fill = bin_assignment) +
  geom_bar(position = 'fill') +
  scale_fill_manual("Legend", values = myColours) +
  ylab("Proportion of total") +
  guides(fill=guide_legend(title="Period")) +
  xlab("Preservation Score") +
  theme_bw() +
  scale_y_continuous(expand = c(0,0), limits = c(0,1))

# Plot preservation scores through time (continuous)
a <- as.data.frame(table(m.dat.series$Preservation_score, 
                         m.dat.series$bin_assignment))
names(a) <- c("Preservation_score", "bin", "Freq")
a <- merge(a, series, by = 'bin')
a$Preservation_score <- factor(a$Preservation_score, 
                               levels = c("5", "4", "3", "2", "1"))
(a <- ggplot(a, aes(x=mid_ma, y=Freq, fill=Preservation_score)) + 
  geom_area() +
  scale_x_reverse() +
  theme_bw() +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  guides(fill=guide_legend(title="Preservation score")) +
  ylab("Frequency") +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Zissou1"))))

# Plot proportional preservation scores through time (continuous)
a <- as.data.frame(table(m.dat.series$Preservation_score, 
                         m.dat.series$bin_assignment))
names(a) <- c("Preservation_score", "bin", "Freq")
a <- merge(a, series, by = 'bin')
a$Preservation_score <- factor(a$Preservation_score, 
                               levels = c("5", "4", "3", "2", "1"))
(a <- ggplot(a, aes(x=mid_ma, y=Freq, fill=Preservation_score)) + 
  geom_area(position = 'fill') +
  scale_x_reverse() +
  theme_bw() +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Preservation score")) +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Zissou1"))))

# Format into table and run Chi Squared test
res.series <- chisq.test(table(m.dat.series$Preservation_score, 
                               m.dat.series$bin_assignment))
res.series
res.series$expected # compare against what the test would have expected

# Mosaic plot
vcd::mosaic(~ Preservation_score + bin_assignment,
       direction = c("v", "h"),
       data = m.dat.series,
       shade = TRUE, 
       labeling_args = list(rot_labels = c(0), just_labels = "right", 
                            set_varnames = c(Preservation_score = "Preservation Score",
                                             bin_assignment = ""))
)

################################################################################
# 4. RUNNING MODELS
################################################################################

#################
##### SETUP #####
#################

model.data <- m.dat.rotate %>%
  filter(is.na(Finalised_grainsize) == F) %>%
  filter(is.na(Finalised_lith) == F) %>%
  filter(is.na(p_lat) == F)

# Assign specific grain size categories into either "Fine Grained" or "Coarse Grained"
model.data <- simple.grain(model.data)

# Remove remaining specimens without grain size
model.data <- model.data %>%
  filter(is.na(Finalised_grainsize) == F) 

# Bin into periods
model.data <- bin_time(model.data, bins = periods, method = 'majority')

# Create factors
order_ind <- rev(c("Permian", "Carboniferous", "Devonian", "Silurian", "Ordovician"))
model.data$bin_assignment <- as.factor(model.data$bin_assignment)
model.data$bin_assignment <- factor(model.data$bin_assignment, levels = order_ind)

model.data$Preservation_score <- factor(model.data$Preservation_score, 
                                     levels = c("1", "2", "3", "4", "5"), 
                                     ordered = TRUE)

# Remove families with low numbers (Cravenechinidae - 2 specimens) and NAs
data.set <- model.data %>%
  filter(Family != 'Cravenechinidae') %>%
  filter(!is.na(Family)) %>%
  dplyr::select(Family, Genus, Species, Rank, Museum_Number, Preservation_score, 
                Continent, Country, lat, lng, Formation, Age, Max_period, 
                Min_period, Finalised_lith, Finalised_grainsize_simp, max_ma, 
                min_ma, bin_assignment, interval_mid_ma, p_lng, p_lat)

###############################
##### LOGISTIC REGRESSION #####
###############################

# Set preservation level to explore and get dataset
pres.score <- c(1,2,3,4,5)

LRresults <- lapply(pres.score, function(f){
  new.pres.score <- f
  data.LR <- set_Pres_score(data.set, c(new.pres.score))
  colnames(data.LR)[colnames(data.LR) == "lat"] <- "Latitude"
  colnames(data.LR)[colnames(data.LR) == "p_lat"] <- "Palaeo-latitude"
  colnames(data.LR)[colnames(data.LR) == "interval_mid_ma"] <- "Age (Ma)"
  colnames(data.LR)[colnames(data.LR) == "Finalised_lith"] <- "Lithology"
  colnames(data.LR)[colnames(data.LR) == "Finalised_grainsize_simp"] <- "Grainsize"
  
  # Set full model
  full.model <- glm(formula = LR_Pres_score ~ Lithology + Grainsize +  
                      Latitude + `Age (Ma)` + `Palaeo-latitude` + Lithology*Grainsize, 
                    family = binomial(link = "logit"), 
                    data = data.LR)
  
  # Alter global actions to allow for dredge
  options(na.action = "na.fail") 
  
  # Get all models, ranked
  model.set <- MuMIn::dredge(full.model)
  
  # Revert global actions
  options(na.action = "na.omit")
  
  # Examine model set
  model.set
  
  # Get best model
  m <- MuMIn::get.models(model.set, subset = 1)
  
  ## view a summary of the model
  summary(m[[1]])
  
  # Save formula as character
  best.form <- as.character(formula(m[[1]]))
  
  # Save a table of coefficients
  (ctable <- as.data.frame(coef(summary(m[[1]]))))
  ctable$PresScore <- f
  ctable$Covariate <- rownames(ctable)
  LRresults <- list(ctable, best.form, f)
  return(LRresults)
})

model.set <- bind_rows(lapply(LRresults, function(l){data.frame(Formula = l[[2]][3], 
                                                             Score = l[[3]])}))
ctable <- bind_rows(lapply(LRresults, `[[`, 1))

# Save models
write.csv(model.set, 
          file = "Results/model.set.csv", 
          row.names = FALSE)
write.csv(ctable, 
          file = "Results/ctable.csv", 
          row.names = TRUE)

################################################################################
# 5. CORRELATION TESTS
################################################################################

#################
##### SETUP #####
#################

# Make columns for period/series, then bin time to stage level (this is so that both are in one dataset)
#m.dat.period$Period <- as.character(m.dat.period$bin_assignment)
#m.dat.period$Period.no <- as.numeric(m.dat.period$bin_assignment)
#m.dat.period <- bin_time(m.dat.period, bins = stages, method = "majority")

m.dat.series$Period <- as.character(m.dat.series$bin_assignment)
m.dat.series$Period.no <- as.numeric(m.dat.series$bin_assignment)
m.dat.period <- bin_time(m.dat.series, bins = stages, method = "majority")

# Split out species and genus level information
species.lvl <- m.dat.period %>%
  dplyr::filter(Rank == "Species") %>%
  dplyr::select(Genus, Species, Locality, bin_midpoint, Period, 
                Period.no, Preservation_score) %>%
  mutate(Combined_name = paste(Genus, Species, sep = " "))

genus.lvl <- m.dat.period %>%
  dplyr::filter(Rank == "Genus") %>%
  dplyr::select(Genus, Locality, bin_midpoint, Period, Period.no, Preservation_score) 

all.lvl <- m.dat.period %>%
  dplyr::select(Locality, bin_midpoint, Period, Period.no, Preservation_score)

##### Make tables of number of taxa per preservation score for each time frame #####

# All specimens
stage.pres.all <- as.data.frame(table(all.lvl$bin_midpoint, 
                                      all.lvl$Preservation_score))
colnames(stage.pres.all) <- c("bin_midpoint", "Preservation_score", "Freq")
period.pres.all <- as.data.frame(table(all.lvl$Period, 
                                       all.lvl$Preservation_score))
# Species
stage.pres.spec <- as.data.frame(table(species.lvl$bin_midpoint, 
                                       species.lvl$Preservation_score))
colnames(stage.pres.spec) <- c("bin_midpoint", "Preservation_score", "Freq")
period.pres.spec <- as.data.frame(table(species.lvl$Period, 
                                        species.lvl$Preservation_score))
# Genera
stage.pres.gen <- as.data.frame(table(genus.lvl$bin_midpoint, 
                                      genus.lvl$Preservation_score))
colnames(stage.pres.gen) <- c("bin_midpoint", "Preservation_score", "Freq")
period.pres.gen <- as.data.frame(table(genus.lvl$Period, 
                                       genus.lvl$Preservation_score))

# Run divDyn to get diversity/collections at genus and species level for each time frame
div.stage.spec <- binstat(species.lvl, 
                          tax="Combined_name", 
                          bin="bin_midpoint", 
                          coll = 'Locality')
div.stage.gen <- binstat(genus.lvl, 
                         tax="Genus", 
                         bin="bin_midpoint", 
                         coll = 'Locality')
div.period.spec <- binstat(species.lvl, 
                           tax= "Combined_name", 
                           bin="Period.no", 
                           coll = 'Locality')
div.period.gen <- binstat(genus.lvl, 
                          tax = "Genus", 
                          bin = "Period.no", 
                          coll = 'Locality')

# Get all collections overall
all.colls <- all.lvl %>%
  dplyr::select(bin_midpoint, Locality, Period.no, Period) %>%
  distinct()
Period.colls <- as.data.frame(table(all.colls$Period))
Stage.colls <- as.data.frame(table(all.colls$bin_midpoint))

# Create order for period
# order_ind <- c("Permian", "Carboniferous", "Devonian", "Silurian", "Ordovician")
# Create order for series
order_ind <- c("Permian", "Pennsylvanian", "Mississippian", 
               "Devonian", "Silurian", "Ordovician")

# Run function for each level
split(period.pres.gen, period.pres.gen$Var2, order = TRUE)
split(period.pres.spec, period.pres.spec$Var2, order = TRUE)
split(period.pres.all, period.pres.all$Var2, order = TRUE)
split(stage.pres.all, Var2 = stage.pres.all$Preservation_score, 
      order = FALSE, stage = TRUE)
split(stage.pres.spec, Var2 = stage.pres.spec$Preservation_score, 
      order = FALSE, stage = TRUE)
split(stage.pres.gen, Var2 = stage.pres.gen$Preservation_score, 
      order = FALSE, stage = TRUE)

temp_stages <- stages[55:92,]
temp_stages$bin_midpoint <- (temp_stages$max_ma + temp_stages$min_ma)/2 
div.stage.spec <- merge(temp_stages, div.stage.spec, 
                        by = "bin_midpoint", all.x = T)
div.stage.gen <- merge(temp_stages, div.stage.gen, 
                       by = "bin_midpoint", all.x = T)
colnames(Stage.colls)[1] <- "bin_midpoint"
col.stage.all <- merge(temp_stages, Stage.colls, 
                       by = "bin_midpoint", all.x = T)

##### Get comprehensive collections through time #####

# Load PBDB Palaeozoic invert occurrences
all.pbdb <- read.csv("Specimen_data/all_palaeozoic_binned.csv")

# Filter to correct age
all.pbdb.period <- all.pbdb %>%
  filter(max_ma < 485.4) %>%
  filter(min_ma > 251.902)

# Bin
all.pbdb.period <- bin_time(occdf = all.pbdb.period, bins = series, 
                            method = 'majority')

pbdb.stat <- binstat(all.pbdb, 
                     tax= "accepted_name", 
                     bin="bin_assignment", 
                     coll = "collection_no", 
                     noNAStart = TRUE)
pbdb.period.stat <- binstat(all.pbdb.period, 
                            tax = "accepted_name", 
                            bin="bin_midpoint", 
                            coll = "collection_no")

pbdb.stat <- pbdb.stat[7:nrow(pbdb.stat),]
pbdb.stat$bin <- pbdb.stat$bin_assignment
pbdb.stat <- merge(temp_stages, pbdb.stat, by = "bin", all.x = T)

# Load all palaeozoic echinodermata occurrences
pb.echino <- read.csv("Specimen_data/Echinodermata.csv", skip = 19)
pb.echino <- pb.echino %>%
  dplyr::filter(max_ma < 485.4000) %>%
  dplyr::filter(min_ma > 251.2000)

pb.echino.period <- bin_time(occdf = pb.echino, bins = series, method = 'majority')
pb.echino <- bin_time(occdf = pb.echino, bins = stages, method = 'majority')

pb.echino <- binstat(pb.echino, 
        tax= "accepted_name", 
        bin="bin_assignment", 
        coll = "collection_no", 
        noNAStart = TRUE)
pb.echino.period <- binstat(pb.echino.period, 
                            tax = "accepted_name", 
                            bin="bin_midpoint", 
                            coll = "collection_no")

pb.echino$bin <- pb.echino$bin_assignment
pb.echino <- merge(temp_stages, pb.echino, by = "bin", all.x = TRUE)

stage.pres.ratio <- as.data.frame(table(all.lvl$bin_midpoint))
colnames(stage.pres.ratio) <- c("bin_midpoint", "Freq")
stage.pres.ratio <- merge(temp_stages, stage.pres.ratio, by = "bin_midpoint", all.x = TRUE)
stage.pres.ratio[is.na(stage.pres.ratio)] <- 0

stage.pres.ratio$ratio_1 <- stage.pres.all.1$Freq/stage.pres.ratio$Freq
stage.pres.ratio$ratio_5 <- stage.pres.all.5$Freq/stage.pres.ratio$Freq

##################################
##### DIVERSITY CORRELATIONS #####
##################################

##### SPECIES #####

# Setup complete list of bins
div.stage.spec.comp <- merge(stage.pres.all.1, div.stage.spec, by = "bin_midpoint", all.x = TRUE)

# Stage - diversity vs.preservation score
stage.spec.div.1 <- cor.test(log10(div.stage.spec.comp$SIBs), 
                             log10(stage.pres.all.1$Freq), method = 'spearman')
stage.spec.div.2 <- cor.test(log10(div.stage.spec.comp$SIBs), 
                             log10(stage.pres.all.2$Freq), method = 'spearman')
stage.spec.div.3 <- cor.test(log10(div.stage.spec.comp$SIBs), 
                             log10(stage.pres.all.3$Freq), method = 'spearman')
stage.spec.div.4 <- cor.test(log10(div.stage.spec.comp$SIBs), 
                             log10(stage.pres.all.4$Freq), method = 'spearman')
stage.spec.div.5 <- cor.test(log10(div.stage.spec.comp$SIBs), 
                             log10(stage.pres.all.5$Freq), method = 'spearman')

stage.spec.ratio.1 <- cor.test(log10(div.stage.spec.comp$SIBs), 
                               log10(stage.pres.ratio$ratio_1), method = 'spearman')
stage.spec.ratio.5 <- cor.test(log10(div.stage.spec.comp$SIBs), 
                               log10(stage.pres.ratio$ratio_5), method = 'spearman')

stage.div.spec <- rbind(c("Stage", 1, "Diversity (species)", 
                          stage.spec.div.1$estimate, stage.spec.div.1$p.value),
                        c("Stage", 2, "Diversity (species)", 
                          stage.spec.div.2$estimate, stage.spec.div.2$p.value),
                        c("Stage", 3, "Diversity (species)", 
                          stage.spec.div.3$estimate, stage.spec.div.3$p.value),
                        c("Stage", 4, "Diversity (species)", 
                          stage.spec.div.4$estimate, stage.spec.div.4$p.value),
                        c("Stage", 5, "Diversity (species)", 
                          stage.spec.div.5$estimate, stage.spec.div.5$p.value))

# Period - diversity vs. preservation score
period.spec.div.1 <- cor.test(log10(div.period.spec$SIBs), 
                              log10(period.pres.all.1$Freq), method = 'spearman')
period.spec.div.2 <- cor.test(log10(div.period.spec$SIBs), 
                              log10(period.pres.all.2$Freq), method = 'spearman')
period.spec.div.3 <- cor.test(log10(div.period.spec$SIBs), 
                              log10(period.pres.all.3$Freq), method = 'spearman')
period.spec.div.4 <- cor.test(log10(div.period.spec$SIBs), 
                              log10(period.pres.all.4$Freq), method = 'spearman')
period.spec.div.5 <- cor.test(log10(div.period.spec$SIBs), 
                              log10(period.pres.all.5$Freq), method = 'spearman')

period.div.spec <- rbind(c("Period", 1, "Diversity (species)",
                           period.spec.div.1$estimate, period.spec.div.1$p.value),
                         c("Period", 2, "Diversity (species)",
                           period.spec.div.2$estimate, period.spec.div.2$p.value),
                         c("Period", 3, "Diversity (species)",
                           period.spec.div.3$estimate, period.spec.div.3$p.value),
                         c("Period", 4, "Diversity (species)",
                           period.spec.div.4$estimate, period.spec.div.4$p.value),
                         c("Period", 5, "Diversity (species)",
                           period.spec.div.5$estimate, period.spec.div.5$p.value))

##### GENERA #####

# Setup complete list of bins
div.stage.gen.comp <- merge(stage.pres.all.1, div.stage.gen, 
                            by = "bin_midpoint", all.x = TRUE)

# Stage - diversity vs.preservation score
stage.gen.div.1 <- cor.test(log10(div.stage.gen.comp$SIBs), 
                            log10(stage.pres.all.1$Freq), method = 'spearman')
stage.gen.div.2 <- cor.test(log10(div.stage.gen.comp$SIBs), 
                            log10(stage.pres.all.2$Freq), method = 'spearman')
stage.gen.div.3 <- cor.test(log10(div.stage.gen.comp$SIBs), 
                            log10(stage.pres.all.3$Freq), method = 'spearman')
stage.gen.div.4 <- cor.test(log10(div.stage.gen.comp$SIBs), 
                            log10(stage.pres.all.4$Freq), method = 'spearman')
stage.gen.div.5 <- cor.test(log10(div.stage.gen.comp$SIBs), 
                            log10(stage.pres.all.5$Freq), method = 'spearman')

stage.gen.ratio.1 <- cor.test(log10(div.stage.gen.comp$SIBs), 
                              log10(stage.pres.ratio$ratio_1), method = 'spearman')
stage.gen.ratio.5 <- cor.test(log10(div.stage.gen.comp$SIBs), 
                              log10(stage.pres.ratio$ratio_5), method = 'spearman')

stage.div.gen <-  rbind(c("Stage", 1, "Diversity (genera)",
                          stage.gen.div.1$estimate, stage.gen.div.1$p.value),
                        c("Stage", 2, "Diversity (genera)",
                          stage.gen.div.2$estimate, stage.gen.div.2$p.value),
                        c("Stage", 3, "Diversity (genera)",
                          stage.gen.div.3$estimate, stage.gen.div.3$p.value),
                        c("Stage", 4, "Diversity (genera)",
                          stage.gen.div.4$estimate, stage.gen.div.4$p.value),
                        c("Stage", 5, "Diversity (genera)",
                          stage.gen.div.5$estimate, stage.gen.div.5$p.value))

# Period - diversity vs. preservation score
period.gen.div.1 <- cor.test(log10(div.period.gen$SIBs), 
                             log10(period.pres.all.1$Freq), method = 'spearman')
period.gen.div.2 <- cor.test(log10(div.period.gen$SIBs), 
                             log10(period.pres.all.2$Freq), method = 'spearman')
period.gen.div.3 <- cor.test(log10(div.period.gen$SIBs), 
                             log10(period.pres.all.3$Freq), method = 'spearman')
period.gen.div.4 <- cor.test(log10(div.period.gen$SIBs), 
                             log10(period.pres.all.4$Freq), method = 'spearman')
period.gen.div.5 <- cor.test(log10(div.period.gen$SIBs), 
                             log10(period.pres.all.5$Freq), method = 'spearman')

period.div.gen <-  rbind(c("Period", 1, "Diversity (genera)", 
                           period.gen.div.1$estimate, period.gen.div.1$p.value),
                         c("Period", 2, "Diversity (genera)", 
                           period.gen.div.2$estimate, period.gen.div.2$p.value),
                         c("Period", 3, "Diversity (genera)", 
                           period.gen.div.3$estimate, period.gen.div.3$p.value),
                         c("Period", 4, "Diversity (genera)", 
                           period.gen.div.4$estimate, period.gen.div.4$p.value),
                         c("Period", 5, "Diversity (genera)", 
                           period.gen.div.5$estimate, period.gen.div.5$p.value))

####################################
##### COLLECTIONS CORRELATIONS #####
####################################

colls.stage.comp <- merge(stage.pres.all.1, Stage.colls, by = "bin_midpoint", all.x = T)

# Stage - echinoid collections vs. preservation score
stage.col.1 <- cor.test(log10(colls.stage.comp$Freq.y), 
                        log10(stage.pres.all.1$Freq), method = 'spearman')
stage.col.2 <- cor.test(log10(colls.stage.comp$Freq.y), 
                        log10(stage.pres.all.2$Freq), method = 'spearman')
stage.col.3 <- cor.test(log10(colls.stage.comp$Freq.y), 
                        log10(stage.pres.all.3$Freq), method = 'spearman')
stage.col.4 <- cor.test(log10(colls.stage.comp$Freq.y), 
                        log10(stage.pres.all.4$Freq), method = 'spearman')
stage.col.5 <- cor.test(log10(colls.stage.comp$Freq.y), 
                        log10(stage.pres.all.5$Freq), method = 'spearman')

stage.coll.ech <-  rbind(c("Stage", 1, "Echinoid collections", 
                           stage.col.1$estimate, stage.col.1$p.value),
                         c("Stage", 2, "Echinoid collections", 
                           stage.col.2$estimate, stage.col.2$p.value),
                         c("Stage", 3, "Echinoid collections", 
                           stage.col.3$estimate, stage.col.3$p.value),
                         c("Stage", 4, "Echinoid collections", 
                           stage.col.4$estimate, stage.col.4$p.value),
                         c("Stage", 5, "Echinoid collections", 
                           stage.col.5$estimate, stage.col.5$p.value))

Period.colls <- Period.colls[match(order_ind, Period.colls$Var1),]

# Period - echinoid collections vs. preservation score
period.col.1 <- cor.test(log10(Period.colls$Freq), 
                         log10(period.pres.all.1$Freq), method = 'spearman')
period.col.2 <- cor.test(log10(Period.colls$Freq), 
                         log10(period.pres.all.2$Freq), method = 'spearman')
period.col.3 <- cor.test(log10(Period.colls$Freq), 
                         log10(period.pres.all.3$Freq), method = 'spearman')
period.col.4 <- cor.test(log10(Period.colls$Freq), 
                         log10(period.pres.all.4$Freq), method = 'spearman')
period.col.5 <- cor.test(log10(Period.colls$Freq), 
                         log10(period.pres.all.5$Freq), method = 'spearman')

period.coll.ech <-  rbind(c("Period", 1, "Echinoid collections", 
                            period.col.1$estimate, period.col.1$p.value),
                          c("Period", 2, "Echinoid collections", 
                            period.col.2$estimate, period.col.2$p.value),
                          c("Period", 3, "Echinoid collections", 
                            period.col.3$estimate, period.col.3$p.value),
                          c("Period", 4, "Echinoid collections", 
                            period.col.4$estimate, period.col.4$p.value),
                          c("Period", 5, "Echinoid collections", 
                            period.col.5$estimate, period.col.5$p.value))

# Stage - global collections vs. preservation score
stage.pbdb.1 <- cor.test(log10(pbdb.stat$colls), 
                         log10(stage.pres.all.1$Freq), method = 'spearman')
stage.pbdb.2 <- cor.test(log10(pbdb.stat$colls), 
                         log10(stage.pres.all.2$Freq), method = 'spearman')
stage.pbdb.3 <- cor.test(log10(pbdb.stat$colls), 
                         log10(stage.pres.all.3$Freq), method = 'spearman')
stage.pbdb.4 <- cor.test(log10(pbdb.stat$colls), 
                         log10(stage.pres.all.4$Freq), method = 'spearman')
stage.pbdb.5 <- cor.test(log10(pbdb.stat$colls), 
                         log10(stage.pres.all.5$Freq), method = 'spearman')

stage.coll.pbdb <-  rbind(c("Stage", 1, "Global collections (pbdb)", 
                            stage.pbdb.1$estimate, stage.pbdb.1$p.value),
                          c("Stage", 2, "Global collections (pbdb)", 
                            stage.pbdb.2$estimate, stage.pbdb.2$p.value),
                          c("Stage", 3, "Global collections (pbdb)", 
                            stage.pbdb.3$estimate, stage.pbdb.3$p.value),
                          c("Stage", 4, "Global collections (pbdb)", 
                            stage.pbdb.4$estimate, stage.pbdb.4$p.value),
                          c("Stage", 5, "Global collections (pbdb)", 
                            stage.pbdb.5$estimate, stage.pbdb.5$p.value))

# Period - global collections vs. preservation score
period.pbdb.1 <- cor.test(log10(pbdb.period.stat$colls), 
                          log10(period.pres.all.1$Freq), method = 'spearman')
period.pbdb.2 <- cor.test(log10(pbdb.period.stat$colls), 
                          log10(period.pres.all.2$Freq), method = 'spearman')
period.pbdb.3 <- cor.test(log10(pbdb.period.stat$colls), 
                          log10(period.pres.all.3$Freq), method = 'spearman')
period.pbdb.4 <- cor.test(log10(pbdb.period.stat$colls), 
                          log10(period.pres.all.4$Freq), method = 'spearman')
period.pbdb.5 <- cor.test(log10(pbdb.period.stat$colls), 
                          log10(period.pres.all.5$Freq), method = 'spearman')

period.coll.pbdb <- rbind(c("Period", 1, "Global collections (pbdb)", 
                            period.pbdb.1$estimate, period.pbdb.1$p.value),
                          c("Period", 2, "Global collections (pbdb)", 
                            period.pbdb.2$estimate, period.pbdb.2$p.value),
                          c("Period", 3, "Global collections (pbdb)", 
                            period.pbdb.3$estimate, period.pbdb.3$p.value),
                          c("Period", 4, "Global collections (pbdb)", 
                            period.pbdb.4$estimate, period.pbdb.4$p.value),
                          c("Period", 5, "Global collections (pbdb)", 
                            period.pbdb.5$estimate, period.pbdb.5$p.value))

colls.stage.echino <- merge(stage.pres.all.1, pb.echino, by = "bin_midpoint", all.x = T)

# Stage - echinodermata collections vs. preservation score
stage.col.echino.1 <- cor.test(log10(colls.stage.echino$colls), 
                               log10(stage.pres.all.1$Freq), method = 'spearman')
stage.col.echino.2 <- cor.test(log10(colls.stage.echino$colls), 
                               log10(stage.pres.all.2$Freq), method = 'spearman')
stage.col.echino.3 <- cor.test(log10(colls.stage.echino$colls), 
                               log10(stage.pres.all.3$Freq), method = 'spearman')
stage.col.echino.4 <- cor.test(log10(colls.stage.echino$colls), 
                               log10(stage.pres.all.4$Freq), method = 'spearman')
stage.col.echino.5 <- cor.test(log10(colls.stage.echino$colls), 
                               log10(stage.pres.all.5$Freq), method = 'spearman')

stage.coll.echino <- rbind(c("Stage", 1, "Echinodermata collections", 
                             stage.col.echino.1$estimate, stage.col.echino.1$p.value),
                           c("Stage", 2, "Echinodermata collections", 
                             stage.col.echino.2$estimate, stage.col.echino.2$p.value),
                           c("Stage", 3, "Echinodermata collections", 
                             stage.col.echino.3$estimate, stage.col.echino.3$p.value),
                           c("Stage", 4, "Echinodermata collections", 
                             stage.col.echino.4$estimate, stage.col.echino.4$p.value),
                           c("Stage", 5, "Echinodermata collections", 
                             stage.col.echino.5$estimate, stage.col.echino.5$p.value))

# Period - echinodermata collections vs. preservation score
period.col.echino.1 <- cor.test(log10(pb.echino.period$colls), 
                                log10(period.pres.all.1$Freq), method = 'spearman')
period.col.echino.2 <- cor.test(log10(pb.echino.period$colls), 
                                log10(period.pres.all.2$Freq), method = 'spearman')
period.col.echino.3 <- cor.test(log10(pb.echino.period$colls), 
                                log10(period.pres.all.3$Freq), method = 'spearman')
period.col.echino.4 <- cor.test(log10(pb.echino.period$colls), 
                                log10(period.pres.all.4$Freq), method = 'spearman')
period.col.echino.5 <- cor.test(log10(pb.echino.period$colls), 
                                log10(period.pres.all.5$Freq), method = 'spearman')

period.coll.echino <- rbind(c("Period", 1, "Echinodermata collections", 
                              period.col.echino.1$estimate, period.col.echino.1$p.value),
                            c("Period", 2, "Echinodermata collections", 
                              period.col.echino.2$estimate, period.col.echino.2$p.value),
                            c("Period", 3, "Echinodermata collections", 
                              period.col.echino.3$estimate, period.col.echino.3$p.value),
                            c("Period", 4, "Echinodermata collections", 
                              period.col.echino.4$estimate, period.col.echino.4$p.value),
                            c("Period", 5, "Echinodermata collections", 
                              period.col.echino.5$estimate, period.col.echino.5$p.value))


(a <- ggplot(pbdb.stat, aes(x=bin_midpoint, y=colls)) + 
  geom_line() +
  scale_x_reverse() +
  theme_bw() +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  ylab("Count") +
  xlab("Time (Ma)"))

(a <- ggplot(pbdb.period.stat, aes(x=bin_midpoint, y=colls)) + 
  geom_line() +
  scale_x_reverse() +
  theme_bw() +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  ylab("Count") +
  xlab("Time (Ma)"))

(a <- ggplot(colls.stage.echino, aes(x=bin_midpoint, y=colls)) + 
  geom_line() +
  scale_x_reverse() +
  theme_bw() +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  ylab("Count") +
  xlab("Time (Ma)"))


abund.counts <- as.data.frame(t(table(m.dat$Preservation_score, m.dat$bin_midpoint)))
colnames(abund.counts) <- c("bin_midpoint", "Score", "Freq")
test1 <- merge(sili.macro.count, abund.counts, by = "bin_midpoint")
test2 <- merge(carb.macro.count, abund.counts, by = "bin_midpoint")
test <- rbind(test1, test2)

test$count <- log10(test$count + 0.0001)
test$Freq <- log10(test$Freq + 0.0001)

test2$count <- log10(test2$count + 0.0001)
test2$Freq <- log10(test2$Freq + 0.0001)

stages2 <- stages %>%
  filter(max_ma > 254.13) %>%
  filter(min_ma <485.4000)

test3 <- bin_time(occdf = stages2, bins = series, method = "mid")

test <- merge(test, test3, by = "bin_midpoint")

test3$bin_midpoint <- (test3$max_ma + test3$min_ma)/2
test1 <- merge(test1, test3, by = "bin_midpoint")
test2 <- merge(test2, test3, by = "bin_midpoint")
test4 <- rbind(test1, test2)

ggplot(test2, aes(x = count, y = Freq, colour = bin_assignment)) +
  
  geom_point() +
 # geom_smooth(method = lm) +
  facet_wrap(~ Score ) +
  theme_bw() 


ggplot(test, aes(x = count, y = Freq, colour = bin_assignment)) +
  geom_point() +
  scale_colour_manual("Legend", values = myColours) +
  #geom_smooth(method = lm) +
  facet_wrap(~ lith + Score, nrow = 2) +
  theme_bw() 


m.dat.series

a <- m.dat.series %>%
  dplyr::group_by(bin_midpoint, Finalised_lith, Rank, Preservation_score) %>%
  dplyr::summarize(test = n()) %>%
  dplyr::filter(is.na(Finalised_lith) != T)

a$Preservation_score <- as.factor(a$Preservation_score)

ggplot(a, aes(x=bin_midpoint, y=test)) + 
    geom_line(aes(colour = Rank)) +
    scale_x_reverse() +
    theme_bw() +
    coord_geo(dat = series2, 
              xlim = c(485.4, 251.902),
              rot = list(0),
              size = 4,
              pos = list("b"),
              abbrv = T) +
    ylab("Count") +
    xlab("Time (Ma)") +
  facet_wrap(~Finalised_lith + Preservation_score, nrow = 2)

#################################
##### SED TYPE CORRELATIONS #####
#################################

##### PLOTS #####

(a <- ggplot(macro.area, aes(x=bin_midpoint, y=count)) + 
  geom_line(aes(colour = lith)) +
  scale_x_reverse() +
  theme_bw() +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  ylab("Count") +
  xlab("Time (Ma)") +
  geom_point(aes(color = lith))) 

(a <- ggplot(macro.count, aes(x=bin_midpoint, y=count)) + 
  geom_line(aes(colour = lith)) +
  scale_x_reverse() +
  theme_bw() +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  ylab("Count") +
  xlab("Time (Ma)") +
  geom_point(aes(color = lith)))

##### SET UP #####

# Get appropriate specimens
all.lvl <- m.dat.period %>%
  dplyr::filter(Continent == "North America") %>%
  dplyr::select(Locality, bin_midpoint, Period, Period.no, Preservation_score)

# Find tallys of Preservation score through time
stage.pres.all <- as.data.frame(table(all.lvl$bin_midpoint, 
                                      all.lvl$Preservation_score))
colnames(stage.pres.all) <- c("bin_midpoint", "Preservation_score", "Freq")
period.pres.all <- as.data.frame(table(all.lvl$Period, 
                                       all.lvl$Preservation_score))

# Split into datasets
split(period.pres.all, period.pres.all$Var2, order = TRUE)
split(stage.pres.all, Var2 = stage.pres.all$Preservation_score, 
      order = FALSE, stage = TRUE)

##### STAGE CORRELATIONS #####

# Taphonomic grade vs. carbonate (count), stage level
stage.carb.macro.count.1 <- cor.test(log10(carb.macro.count$count), 
                                     log10(stage.pres.all.1$Freq), method = 'spearman')
stage.carb.macro.count.2 <- cor.test(log10(carb.macro.count$count), 
                                     log10(stage.pres.all.2$Freq), method = 'spearman')
stage.carb.macro.count.3 <- cor.test(log10(carb.macro.count$count), 
                                     log10(stage.pres.all.3$Freq), method = 'spearman')
stage.carb.macro.count.4 <- cor.test(log10(carb.macro.count$count), 
                                     log10(stage.pres.all.4$Freq), method = 'spearman')
stage.carb.macro.count.5 <- cor.test(log10(carb.macro.count$count), 
                                     log10(stage.pres.all.5$Freq), method = 'spearman')

stage.carb.macro.count.ratio.1 <- cor.test(log10(carb.macro.count$count), 
                                           log10(stage.pres.ratio$ratio_1), 
                                           method = 'spearman')
stage.carb.macro.count.ratio.5 <- cor.test(log10(carb.macro.count$count), 
                                           log10(stage.pres.ratio$ratio_5), 
                                           method = 'spearman')

stage.macro.carb.count <-  rbind(c("Stage", 1, "Macrostrat carbonate (count)",
                                   stage.carb.macro.count.1$estimate,
                                   stage.carb.macro.count.1$p.value),
                                 c("Stage", 2, "Macrostrat carbonate (count)",
                                   stage.carb.macro.count.2$estimate, 
                                   stage.carb.macro.count.2$p.value),
                                 c("Stage", 3, "Macrostrat carbonate (count)",
                                   stage.carb.macro.count.3$estimate, 
                                   stage.carb.macro.count.3$p.value),
                                 c("Stage", 4, "Macrostrat carbonate (count)",
                                   stage.carb.macro.count.4$estimate, 
                                   stage.carb.macro.count.4$p.value),
                                 c("Stage", 5, "Macrostrat carbonate (count)",
                                   stage.carb.macro.count.5$estimate,
                                   stage.carb.macro.count.5$p.value))

# Taphonomic grade vs. siliclastic (count), stage level
stage.sili.macro.count.1 <- cor.test(log10(sili.macro.count$count),
                                     log10(stage.pres.all.1$Freq), 
                                     method = 'spearman')
stage.sili.macro.count.2 <- cor.test(log10(sili.macro.count$count),
                                     log10(stage.pres.all.2$Freq),
                                     method = 'spearman')
stage.sili.macro.count.3 <- cor.test(log10(sili.macro.count$count),
                                     log10(stage.pres.all.3$Freq), 
                                     method = 'spearman')
stage.sili.macro.count.4 <- cor.test(log10(sili.macro.count$count),
                                     log10(stage.pres.all.4$Freq), 
                                     method = 'spearman')
stage.sili.macro.count.5 <- cor.test(log10(sili.macro.count$count),
                                     log10(stage.pres.all.5$Freq), 
                                     method = 'spearman')

stage.sili.macro.count.ratio.1 <- cor.test(log10(sili.macro.count$count), 
                                           log10(stage.pres.ratio$ratio_1), 
                                           method = 'spearman')
stage.sili.macro.count.ratio.5 <- cor.test(log10(sili.macro.count$count), 
                                           log10(stage.pres.ratio$ratio_5), 
                                           method = 'spearman')

stage.macro.sili.count <- rbind(c("Stage", 1, "Macrostrat silicilcastic (count)",
                                  stage.sili.macro.count.1$estimate, 
                                  stage.sili.macro.count.1$p.value),
                                c("Stage", 2, "Macrostrat siliciclastic (count)",
                                  stage.sili.macro.count.2$estimate, 
                                  stage.sili.macro.count.2$p.value),
                                c("Stage", 3, "Macrostrat siliciclastic (count)",
                                  stage.sili.macro.count.3$estimate, 
                                  stage.sili.macro.count.3$p.value),
                                c("Stage", 4, "Macrostrat siliciclastic (count)",
                                  stage.sili.macro.count.4$estimate, 
                                  stage.sili.macro.count.4$p.value),
                                c("Stage", 5, "Macrostrat siliciclastic (count)",
                                  stage.sili.macro.count.5$estimate, 
                                  stage.sili.macro.count.5$p.value))

# Taphonomic grade vs. carbonate (area), stage level
stage.carb.macro.area.1 <- cor.test(log10(carb.macro.area$count), 
                                    log10(stage.pres.all.1$Freq), 
                                    method = 'spearman')
stage.carb.macro.area.2 <- cor.test(log10(carb.macro.area$count), 
                                    log10(stage.pres.all.2$Freq), 
                                    method = 'spearman')
stage.carb.macro.area.3 <- cor.test(log10(carb.macro.area$count), 
                                    log10(stage.pres.all.3$Freq), 
                                    method = 'spearman')
stage.carb.macro.area.4 <- cor.test(log10(carb.macro.area$count), 
                                    log10(stage.pres.all.4$Freq), 
                                    method = 'spearman')
stage.carb.macro.area.5 <- cor.test(log10(carb.macro.area$count), 
                                    log10(stage.pres.all.5$Freq), 
                                    method = 'spearman')

stage.macro.carb.area <- rbind(c("Stage", 1, "Macrostrat carbonate (area)",
                                 stage.carb.macro.area.1$estimate, 
                                 stage.carb.macro.area.1$p.value),
                               c("Stage", 2, "Macrostrat carbonate (area)",
                                 stage.carb.macro.area.2$estimate, 
                                 stage.carb.macro.area.2$p.value),
                               c("Stage", 3, "Macrostrat carbonate (area)",
                                 stage.carb.macro.area.3$estimate, 
                                 stage.carb.macro.area.3$p.value),
                               c("Stage", 4, "Macrostrat carbonate (area)",
                                 stage.carb.macro.area.4$estimate, 
                                 stage.carb.macro.area.4$p.value),
                               c("Stage", 5, "Macrostrat carbonate (area)",
                                 stage.carb.macro.area.5$estimate, 
                                 stage.carb.macro.area.5$p.value))

# Taphonomic grade vs. siliclastic (area), stage level
stage.sili.macro.area.1 <- cor.test(log10(sili.macro.area$count), 
                                    log10(stage.pres.all.1$Freq),
                                    method = 'spearman')
stage.sili.macro.area.2 <- cor.test(log10(sili.macro.area$count), 
                                    log10(stage.pres.all.2$Freq),
                                    method = 'spearman')
stage.sili.macro.area.3 <- cor.test(log10(sili.macro.area$count), 
                                    log10(stage.pres.all.3$Freq),
                                    method = 'spearman')
stage.sili.macro.area.4 <- cor.test(log10(sili.macro.area$count), 
                                    log10(stage.pres.all.4$Freq),
                                    method = 'spearman')
stage.sili.macro.area.5 <- cor.test(log10(sili.macro.area$count), 
                                    log10(stage.pres.all.5$Freq),
                                    method = 'spearman')

stage.macro.sili.area  <-  rbind(c("Stage", 1, "Macrostrat siliciclastic (area)",
                                   stage.sili.macro.area.1$estimate,
                                   stage.sili.macro.area.1$p.value),
                                 c("Stage", 2, "Macrostrat siliciclastic (area)",
                                   stage.sili.macro.area.2$estimate,
                                   stage.sili.macro.area.2$p.value),
                                 c("Stage", 3, "Macrostrat siliciclastic (area)",
                                   stage.sili.macro.area.3$estimate,
                                   stage.sili.macro.area.3$p.value),
                                 c("Stage", 4, "Macrostrat siliciclastic (area)",
                                   stage.sili.macro.area.4$estimate,
                                   stage.sili.macro.area.4$p.value),
                                 c("Stage", 5, "Macrostrat siliciclastic (area)",
                                   stage.sili.macro.area.5$estimate,
                                   stage.sili.macro.area.5$p.value))

##### STAGE CORRELATIONS, MISSING STAGES REMOVED #####

# Remove stages based on overall occurrences

occs <- as.data.frame(table(m.dat$bin_assignment))
colnames(occs) <- c("bin", "Freq")
occs$bin <- as.numeric(as.character(occs$bin))
ref <- full_join(stages, occs, by = "bin")
ref <- ref %>%
  filter(max_ma < 485.41) %>%
  filter(min_ma > 251.901)
ref[is.na(ref) == T] <- 0
ref <- ref[ref$Freq == 0,]

# Set up carb.macro and sili.macro
stages$bin_midpoint <- (stages$max_ma + stages$min_ma)/2
carb.macro.count.removed <- full_join(stages, carb.macro.count, by = "bin_midpoint")
sili.macro.count.removed <- full_join(stages, sili.macro.count, by = "bin_midpoint")
carb.macro.area.removed <- full_join(stages, carb.macro.area, by = "bin_midpoint")
sili.macro.area.removed <- full_join(stages, sili.macro.area, by = "bin_midpoint")

# Remove rows
stage.pres.all.1.removed <- stage.pres.all.1[!stage.pres.all.1$name %in% ref$name, ]
stage.pres.all.2.removed <- stage.pres.all.2[!stage.pres.all.2$name %in% ref$name, ]
stage.pres.all.3.removed <- stage.pres.all.3[!stage.pres.all.3$name %in% ref$name, ]
stage.pres.all.4.removed <- stage.pres.all.4[!stage.pres.all.4$name %in% ref$name, ]
stage.pres.all.5.removed <- stage.pres.all.5[!stage.pres.all.5$name %in% ref$name, ]
carb.macro.count.removed <- carb.macro.count.removed[!carb.macro.count.removed$name %in% ref$name, ]
sili.macro.count.removed <- sili.macro.count.removed[!sili.macro.count.removed$name %in% ref$name, ]
carb.macro.area.removed <- carb.macro.area.removed[!carb.macro.area.removed$name %in% ref$name, ]
sili.macro.area.removed <- sili.macro.area.removed[!sili.macro.area.removed$name %in% ref$name, ]
sili.macro.area.removed <- sili.macro.area.removed[!is.na(sili.macro.area.removed$count), ]
sili.macro.count.removed <- sili.macro.count.removed[!is.na(sili.macro.count.removed$count), ]
carb.macro.area.removed <- carb.macro.area.removed[!is.na(carb.macro.area.removed$count), ]
carb.macro.count.removed <- carb.macro.count.removed[!is.na(carb.macro.count.removed$count), ]

# Taphonomic grade vs. carbonate (count), stage level
stage.carb.macro.count.removed.1 <- cor.test(log10(carb.macro.count.removed$count), 
                                     log10(stage.pres.all.1.removed$Freq), method = 'spearman')
stage.carb.macro.count.removed.2 <- cor.test(log10(carb.macro.count.removed$count), 
                                     log10(stage.pres.all.2.removed$Freq), method = 'spearman')
stage.carb.macro.count.removed.3 <- cor.test(log10(carb.macro.count.removed$count), 
                                     log10(stage.pres.all.3.removed$Freq), method = 'spearman')
stage.carb.macro.count.removed.4 <- cor.test(log10(carb.macro.count.removed$count), 
                                     log10(stage.pres.all.4.removed$Freq), method = 'spearman')
stage.carb.macro.count.removed.5 <- cor.test(log10(carb.macro.count.removed$count), 
                                     log10(stage.pres.all.5.removed$Freq), method = 'spearman')

stage.macro.carb.count.removed <-  rbind(c("Stage", 1, "Macrostrat carbonate (count)",
                                   stage.carb.macro.count.removed.1$estimate,
                                   stage.carb.macro.count.removed.1$p.value),
                                 c("Stage", 2, "Macrostrat carbonate (count)",
                                   stage.carb.macro.count.removed.2$estimate, 
                                   stage.carb.macro.count.removed.2$p.value),
                                 c("Stage", 3, "Macrostrat carbonate (count)",
                                   stage.carb.macro.count.removed.3$estimate, 
                                   stage.carb.macro.count.removed.3$p.value),
                                 c("Stage", 4, "Macrostrat carbonate (count)",
                                   stage.carb.macro.count.removed.4$estimate, 
                                   stage.carb.macro.count.removed.4$p.value),
                                 c("Stage", 5, "Macrostrat carbonate (count)",
                                   stage.carb.macro.count.removed.5$estimate,
                                   stage.carb.macro.count.removed.5$p.value))

# Taphonomic grade vs. siliclastic (count), stage level
stage.sili.macro.count.removed.1 <- cor.test(log10(sili.macro.count.removed$count),
                                     log10(stage.pres.all.1.removed$Freq), 
                                     method = 'spearman')
stage.sili.macro.count.removed.2 <- cor.test(log10(sili.macro.count.removed$count),
                                     log10(stage.pres.all.2.removed$Freq),
                                     method = 'spearman')
stage.sili.macro.count.removed.3 <- cor.test(log10(sili.macro.count.removed$count),
                                     log10(stage.pres.all.3.removed$Freq), 
                                     method = 'spearman')
stage.sili.macro.count.removed.4 <- cor.test(log10(sili.macro.count.removed$count),
                                     log10(stage.pres.all.4.removed$Freq), 
                                     method = 'spearman')
stage.sili.macro.count.removed.5 <- cor.test(log10(sili.macro.count.removed$count),
                                     log10(stage.pres.all.5.removed$Freq), 
                                     method = 'spearman')

stage.macro.sili.count.removed <- rbind(c("Stage", 1, "Macrostrat silicilcastic (count)",
                                  stage.sili.macro.count.removed.1$estimate, 
                                  stage.sili.macro.count.removed.1$p.value),
                                c("Stage", 2, "Macrostrat siliciclastic (count)",
                                  stage.sili.macro.count.removed.2$estimate, 
                                  stage.sili.macro.count.removed.2$p.value),
                                c("Stage", 3, "Macrostrat siliciclastic (count)",
                                  stage.sili.macro.count.removed.3$estimate, 
                                  stage.sili.macro.count.removed.3$p.value),
                                c("Stage", 4, "Macrostrat siliciclastic (count)",
                                  stage.sili.macro.count.removed.4$estimate, 
                                  stage.sili.macro.count.removed.4$p.value),
                                c("Stage", 5, "Macrostrat siliciclastic (count)",
                                  stage.sili.macro.count.removed.5$estimate, 
                                  stage.sili.macro.count.removed.5$p.value))

# Taphonomic grade vs. carbonate (area), stage level
stage.carb.macro.area.removed.1 <- cor.test(log10(carb.macro.area.removed$count), 
                                    log10(stage.pres.all.1.removed$Freq), 
                                    method = 'spearman')
stage.carb.macro.area.removed.2 <- cor.test(log10(carb.macro.area.removed$count), 
                                    log10(stage.pres.all.2.removed$Freq), 
                                    method = 'spearman')
stage.carb.macro.area.removed.3 <- cor.test(log10(carb.macro.area.removed$count), 
                                    log10(stage.pres.all.3.removed$Freq), 
                                    method = 'spearman')
stage.carb.macro.area.removed.4 <- cor.test(log10(carb.macro.area.removed$count), 
                                    log10(stage.pres.all.4.removed$Freq), 
                                    method = 'spearman')
stage.carb.macro.area.removed.5 <- cor.test(log10(carb.macro.area.removed$count), 
                                    log10(stage.pres.all.5.removed$Freq), 
                                    method = 'spearman')

stage.macro.carb.area.removed <- rbind(c("Stage", 1, "Macrostrat carbonate (area)",
                                 stage.carb.macro.area.removed.1$estimate, 
                                 stage.carb.macro.area.removed.1$p.value),
                               c("Stage", 2, "Macrostrat carbonate (area)",
                                 stage.carb.macro.area.removed.2$estimate, 
                                 stage.carb.macro.area.removed.2$p.value),
                               c("Stage", 3, "Macrostrat carbonate (area)",
                                 stage.carb.macro.area.removed.3$estimate, 
                                 stage.carb.macro.area.removed.3$p.value),
                               c("Stage", 4, "Macrostrat carbonate (area)",
                                 stage.carb.macro.area.removed.4$estimate, 
                                 stage.carb.macro.area.removed.4$p.value),
                               c("Stage", 5, "Macrostrat carbonate (area)",
                                 stage.carb.macro.area.removed.5$estimate, 
                                 stage.carb.macro.area.removed.5$p.value))

# Taphonomic grade vs. siliclastic (area), stage level
stage.sili.macro.area.removed.1 <- cor.test(log10(sili.macro.area.removed$count), 
                                    log10(stage.pres.all.1.removed$Freq),
                                    method = 'spearman')
stage.sili.macro.area.removed.2 <- cor.test(log10(sili.macro.area.removed$count), 
                                    log10(stage.pres.all.2.removed$Freq),
                                    method = 'spearman')
stage.sili.macro.area.removed.3 <- cor.test(log10(sili.macro.area.removed$count), 
                                    log10(stage.pres.all.3.removed$Freq),
                                    method = 'spearman')
stage.sili.macro.area.removed.4 <- cor.test(log10(sili.macro.area.removed$count), 
                                    log10(stage.pres.all.4.removed$Freq),
                                    method = 'spearman')
stage.sili.macro.area.removed.5 <- cor.test(log10(sili.macro.area.removed$count), 
                                    log10(stage.pres.all.5.removed$Freq),
                                    method = 'spearman')

stage.macro.sili.area.removed  <-  rbind(c("Stage", 1, "Macrostrat siliciclastic (area)",
                                   stage.sili.macro.area.removed.1$estimate,
                                   stage.sili.macro.area.removed.1$p.value),
                                 c("Stage", 2, "Macrostrat siliciclastic (area)",
                                   stage.sili.macro.area.removed.2$estimate,
                                   stage.sili.macro.area.removed.2$p.value),
                                 c("Stage", 3, "Macrostrat siliciclastic (area)",
                                   stage.sili.macro.area.removed.3$estimate,
                                   stage.sili.macro.area.removed.3$p.value),
                                 c("Stage", 4, "Macrostrat siliciclastic (area)",
                                   stage.sili.macro.area.removed.4$estimate,
                                   stage.sili.macro.area.removed.4$p.value),
                                 c("Stage", 5, "Macrostrat siliciclastic (area)",
                                   stage.sili.macro.area.removed.5$estimate,
                                   stage.sili.macro.area.removed.5$p.value))

# COMBINED

cor.results.removed <- as.data.frame(rbind(
                                   stage.macro.carb.count.removed,
                                   stage.macro.sili.count.removed,
                                   stage.macro.carb.area.removed,
                                   stage.macro.sili.area.removed
))

# Setup columns
colnames(cor.results.removed) <- c("Bin size", "Preservation Score", "vs.", "Rho", "p")
cor.results.removed$Rho <- signif(as.numeric(cor.results.removed$Rho), digits = 3)
cor.results.removed$p <- signif(as.numeric(cor.results.removed$p), digits = 5)

# Correct for multiple tests
cor.results.removed$BH <- p.adjust(cor.results.removed$p, method = "BH")
cor.results.removed$Signif <- ifelse(cor.results.removed$BH < 0.05, "*", "")

# Write .csv
write.csv(cor.results.removed, "Results/All_correlations_BH_corrected_removed.csv")

##### PERIOD CORRELATIONS #####

# Taphonomic grade vs. carbonate (count), period level
period.carb.macro.count.1 <- cor.test(log10(carb.macro.count.period$count), 
                                      log10(period.pres.all.1$Freq), 
                                      method = 'spearman')
period.carb.macro.count.2 <- cor.test(log10(carb.macro.count.period$count), 
                                      log10(period.pres.all.2$Freq), 
                                      method = 'spearman')
period.carb.macro.count.3 <- cor.test(log10(carb.macro.count.period$count), 
                                      log10(period.pres.all.3$Freq), 
                                      method = 'spearman')
period.carb.macro.count.4 <- cor.test(log10(carb.macro.count.period$count), 
                                      log10(period.pres.all.4$Freq), 
                                      method = 'spearman')
period.carb.macro.count.5 <- cor.test(log10(carb.macro.count.period$count), 
                                      log10(period.pres.all.5$Freq), 
                                      method = 'spearman')

period.macro.carb.count <- rbind(c("Period", 1, "Macrostrat carbonate (count)", 
                                   period.carb.macro.count.1$estimate,
                                   period.carb.macro.count.1$p.value),
                                 c("Period", 2, "Macrostrat carbonate (count)", 
                                   period.carb.macro.count.2$estimate,
                                   period.carb.macro.count.2$p.value),
                                 c("Period", 3, "Macrostrat carbonate (count)", 
                                   period.carb.macro.count.3$estimate,
                                   period.carb.macro.count.3$p.value),
                                 c("Period", 4, "Macrostrat carbonate (count)", 
                                   period.carb.macro.count.4$estimate,
                                   period.carb.macro.count.4$p.value),
                                 c("Period", 5, "Macrostrat carbonate (count)", 
                                   period.carb.macro.count.5$estimate,
                                   period.carb.macro.count.5$p.value))

# Taphonomic grade vs. siliclastic (count), period level
period.sili.macro.count.1 <- cor.test(log10(sili.macro.count.period$count),
                                      log10(period.pres.all.1$Freq), 
                                      method = 'spearman')
period.sili.macro.count.2 <- cor.test(log10(sili.macro.count.period$count),
                                      log10(period.pres.all.2$Freq), 
                                      method = 'spearman')
period.sili.macro.count.3 <- cor.test(log10(sili.macro.count.period$count),
                                      log10(period.pres.all.3$Freq), 
                                      method = 'spearman')
period.sili.macro.count.4 <- cor.test(log10(sili.macro.count.period$count),
                                      log10(period.pres.all.4$Freq), 
                                      method = 'spearman')
period.sili.macro.count.5 <- cor.test(log10(sili.macro.count.period$count),
                                      log10(period.pres.all.5$Freq), 
                                      method = 'spearman')

period.macro.sili.count <- rbind(c("Period", 1, "Macrostrat silicilcastic (count)", 
                                   period.sili.macro.count.1$estimate,
                                   period.sili.macro.count.1$p.value),
                                 c("Period", 2, "Macrostrat siliciclastic (count)", 
                                   period.sili.macro.count.2$estimate,
                                   period.sili.macro.count.2$p.value),
                                 c("Period", 3, "Macrostrat siliciclastic (count)", 
                                   period.sili.macro.count.3$estimate,
                                   period.sili.macro.count.3$p.value),
                                 c("Period", 4, "Macrostrat siliciclastic (count)", 
                                   period.sili.macro.count.4$estimate,
                                   period.sili.macro.count.4$p.value),
                                 c("Period", 5, "Macrostrat siliciclastic (count)", 
                                   period.sili.macro.count.5$estimate,
                                   period.sili.macro.count.5$p.value))

# Taphonomic grade vs. carbonate (area), period level
period.carb.macro.area.1 <- cor.test(log10(carb.macro.area.period$count), 
                                     log10(period.pres.all.1$Freq), 
                                     method = 'spearman')
period.carb.macro.area.2 <- cor.test(log10(carb.macro.area.period$count), 
                                     log10(period.pres.all.2$Freq), 
                                     method = 'spearman')
period.carb.macro.area.3 <- cor.test(log10(carb.macro.area.period$count), 
                                     log10(period.pres.all.3$Freq), 
                                     method = 'spearman')
period.carb.macro.area.4 <- cor.test(log10(carb.macro.area.period$count), 
                                     log10(period.pres.all.4$Freq), 
                                     method = 'spearman')
period.carb.macro.area.5 <- cor.test(log10(carb.macro.area.period$count), 
                                     log10(period.pres.all.5$Freq), 
                                     method = 'spearman')

period.macro.carb.area <- rbind(c("Period", 1, "Macrostrat carbonate (area)", 
                                  period.carb.macro.area.1$estimate, 
                                  period.carb.macro.area.1$p.value),
                                c("Period", 2, "Macrostrat carbonate (area)", 
                                  period.carb.macro.area.2$estimate, 
                                  period.carb.macro.area.2$p.value),
                                c("Period", 3, "Macrostrat carbonate (area)", 
                                  period.carb.macro.area.3$estimate, 
                                  period.carb.macro.area.3$p.value),
                                c("Period", 4, "Macrostrat carbonate (area)", 
                                  period.carb.macro.area.4$estimate, 
                                  period.carb.macro.area.4$p.value),
                                c("Period", 5, "Macrostrat carbonate (area)", 
                                  period.carb.macro.area.5$estimate, 
                                  period.carb.macro.area.5$p.value))

# Taphonomic grade vs. siliclastic (area), period level
period.sili.macro.area.1 <- cor.test(log10(sili.macro.area.period$count), 
                                     log10(period.pres.all.1$Freq), 
                                     method = 'spearman')
period.sili.macro.area.2 <- cor.test(log10(sili.macro.area.period$count), 
                                     log10(period.pres.all.2$Freq), 
                                     method = 'spearman')
period.sili.macro.area.3 <- cor.test(log10(sili.macro.area.period$count), 
                                     log10(period.pres.all.3$Freq), 
                                     method = 'spearman')
period.sili.macro.area.4 <- cor.test(log10(sili.macro.area.period$count), 
                                     log10(period.pres.all.4$Freq), 
                                     method = 'spearman')
period.sili.macro.area.5 <- cor.test(log10(sili.macro.area.period$count), 
                                     log10(period.pres.all.5$Freq), 
                                     method = 'spearman')

period.macro.sili.area  <- rbind(c("Period", 1, "Macrostrat siliciclastic (area)", 
                                   period.sili.macro.area.1$estimate, 
                                   period.sili.macro.area.1$p.value),
                                 c("Period", 2, "Macrostrat siliciclastic (area)", 
                                   period.sili.macro.area.2$estimate, 
                                   period.sili.macro.area.2$p.value),
                                 c("Period", 3, "Macrostrat siliciclastic (area)", 
                                   period.sili.macro.area.3$estimate, 
                                   period.sili.macro.area.3$p.value),
                                 c("Period", 4, "Macrostrat siliciclastic (area)", 
                                   period.sili.macro.area.4$estimate, 
                                   period.sili.macro.area.4$p.value),
                                 c("Period", 5, "Macrostrat siliciclastic (area)", 
                                   period.sili.macro.area.5$estimate, 
                                   period.sili.macro.area.5$p.value))

############################
##### COMPLETE RESULTS #####
############################

# Combine results
cor.results <- as.data.frame(rbind(stage.div.spec, 
                     stage.div.gen, 
                     period.div.spec, 
                     period.div.gen,
                     stage.coll.ech, 
                     period.coll.ech,
                     stage.coll.echino,
                     period.coll.echino,
                     stage.coll.pbdb, 
                     period.coll.pbdb,
                     stage.macro.carb.count,
                     stage.macro.sili.count,
                     stage.macro.carb.area,
                     stage.macro.sili.area,
                     period.macro.carb.count,
                     period.macro.sili.count,
                     period.macro.carb.area,
                     period.macro.sili.area
))

# Setup columns
colnames(cor.results) <- c("Bin size", "Preservation Score", "vs.", "Rho", "p")
cor.results$Rho <- signif(as.numeric(cor.results$Rho), digits = 3)
cor.results$p <- signif(as.numeric(cor.results$p), digits = 5)

# Correct for multiple tests
cor.results$BH <- p.adjust(cor.results$p, method = "BH")
cor.results$Signif <- ifelse(cor.results$BH < 0.05, "*", "")

# Write .csv
write.csv(cor.results, "Results/All_correlations_BH_corrected.csv")

################################################################################
# 6. MAPS
################################################################################

##################
##### MODERN #####
##################

# Set extent
e <- extent(-180, 180, -90, 90)

##### ALL SPECIMENS #####
get_grid_im(m.dat, 2,  "Echinoids", ext = e)

##### TAPH GRADES #####
get_grid_im(dplyr::filter(m.dat, Preservation_score == 1), 2,  
            "Taphonomic Grade 1 Echinoids", ext = e)
get_grid_im(dplyr::filter(m.dat, Preservation_score == 2), 2,  
            "Taphonomic Grade 2 Echinoids", ext = e)
get_grid_im(dplyr::filter(m.dat, Preservation_score == 3), 2,  
            "Taphonomic Grade 3 Echinoids", ext = e)
get_grid_im(dplyr::filter(m.dat, Preservation_score == 4), 2,  
            "Taphonomic Grade 4 Echinoids", ext = e)
get_grid_im(dplyr::filter(m.dat, Preservation_score == 5), 2,  
            "Taphonomic Grade 5 Echinoids", ext = e)
