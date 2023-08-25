################################################################################
#                                                                              #
# ECHINOID TAPHONOMY                                                           #
# J. Thompson, C.D. Dean, M. Ford, T.A.M. Ewin.                                #
# Code by C. D. Dean, started 29.09.22                                         #
#                                                                              #
################################################################################

################################################################################
# 1. SETUP
################################################################################

# Load setup file
source("0.Setup.R")

##############################
##### SETUP FOR ANALYSES #####
##############################

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
series <- read.csv("series.csv")

# Bin into series
order_ind <- c("Permian", "Pennsylvanian", "Mississippian", "Devonian", "Silurian", "Ordovician")
m.dat.series <- palaeoverse::bin_time(m.dat, bins = series, method = 'majority')
m.dat.series$bin_assignment <- as.factor(m.dat.series$bin_assignment)
m.dat.series$bin_assignment <- factor(m.dat.series$bin_assignment, levels = order_ind)

# Assign colours
myColours <- series$color
names(myColours) <- levels(m.dat.series$bin_assignment)
custom_colours <- scale_colour_manual(name = "bin_assignment", values = myColours[1:6])

################################################################################
# 2. BAR PLOTS AND CHI SQUARED
################################################################################

##########################
##### TAXONOMIC RANK #####
##########################

# Make bar plot
make.bar.plot(m.dat, m.dat$Rank, "Taxonomic Rank", colour = "Zissou1", FALSE)

# Make proportional bar plot
make.prop.bar.plot(m.dat, m.dat$Rank, "Taxonomic Rank", colour = "Zissou1", FALSE)

# Format into table and run Chi Squared test
test <- chisq.test(table(m.dat$Preservation_score, m.dat$Rank))
test
test$expected # compare against what the test would have expected

# Mosaic plot
mosaic(~ Preservation_score + Rank,
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
test <- chisq.test(table(l.m.dat$Preservation_score, l.m.dat$Finalised_lith))
test
test$expected # compare against what the test would have expected

# Mosaic plot
mosaic(~ Preservation_score + Finalised_lith,
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
make.bar.plot(g.m.dat, g.m.dat$Finalised_grainsize, "Grain Size", colour = "Darjeeling2", flip = FALSE)

# Make proportional bar plot
make.prop.bar.plot(g.m.dat, g.m.dat$Finalised_grainsize, "Grain Size", colour = "Darjeeling2", flip = FALSE)

# Format into table and run Chi Squared test
test <- chisq.test(table(g.m.dat$Preservation_score, g.m.dat$Finalised_grainsize))
test
test$expected # compare against what the test would have expected

# Mosaic plot
mosaic(~ Preservation_score + Finalised_grainsize,
       direction = c("v", "h"),
       data = g.m.dat,
       labeling_args = list(just_labels = "right", 
                            set_varnames = c(Preservation_score = "Preservation Score", 
                                             Finalised_grainsize = "Grain size")),
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

# Format into table and run Chi Squared test
test <- chisq.test(table(f.m.dat$Preservation_score, f.m.dat$Family))
test
test$expected # compare against what the test would have expected

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
a <- as.data.frame(table(m.dat.period$Preservation_score, m.dat.period$bin_assignment))
names(a) <- c("Preservation_score", "interval_name", "Freq")
a <- merge(a, df, by = 'interval_name')
a$Preservation_score <- factor(a$Preservation_score, levels = c("5", "4", "3", "2", "1"))
a <- ggplot(a, aes(x=mid_ma, y=Freq, fill=Preservation_score)) + 
  geom_area() +
  scale_x_reverse() +
  theme_bw() +
  guides(fill=guide_legend(title="Preservation score")) +
  ylab("Frequency") +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Zissou1")))
gggeo_scale(a)

# Plot proportional preservation scores through time (continuous)
a <- as.data.frame(table(m.dat.period$Preservation_score, m.dat.period$bin_assignment))
names(a) <- c("Preservation_score", "interval_name", "Freq")
a <- merge(a, df, by = 'interval_name')
a$Preservation_score <- factor(a$Preservation_score, levels = c("5", "4", "3", "2", "1"))
a <- ggplot(a, aes(x=mid_ma, y=Freq, fill=Preservation_score)) + 
  geom_area(position = 'fill') +
  scale_x_reverse() +
  theme_bw() +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Preservation score")) +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Zissou1"))) 
gggeo_scale(a)

# Format into table and run Chi Squared test
test <- chisq.test(table(m.dat.period$Preservation_score, m.dat.period$bin_assignment))
test
test$expected # compare against what the test would have expected

# Mosaic plot
mosaic(~ Preservation_score + bin_assignment,
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
test$Preservation_score <- factor(test$Preservation_score, levels = c("5", "4", "3", "2", "1"))
a <- ggplot(test, aes(x=mid_ma, y=Freq, fill=Preservation_score)) + 
  geom_area() +
  scale_x_reverse() +
  theme_bw() +
  guides(fill=guide_legend(title="Preservation score")) +
  ylab("Frequency") +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Zissou1")))
gggeo_scale(a)

# Preservation Score (proportional)
a <- as.data.frame(table(m.dat.period$Preservation_score, m.dat$bin_midpoint))
names(a) <- c("Preservation_score", "mid_ma", "Freq")
a$mid_ma <- as.numeric(as.character(a$mid_ma))
a$Preservation_score <- factor(a$Preservation_score, levels = c("5", "4", "3", "2", "1"))
a <- ggplot(a, aes(x=mid_ma, y=Freq, fill=Preservation_score)) + 
  geom_area(position = 'fill') +
  scale_x_reverse() +
  theme_bw() +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Preservation score")) +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Zissou1"))) 
gggeo_scale(a)

# Lithology
test <- as.data.frame(table(l.m.dat$Finalised_lith, l.m.dat$bin_midpoint))
names(test) <- c("Lithology", "mid_ma", "Freq")
test$mid_ma <- as.numeric(as.character(test$mid_ma))
a <- ggplot(test, aes(x=mid_ma, y=Freq, fill=Lithology)) + 
  geom_area() +
  scale_x_reverse() +
  theme_bw() +
  guides(fill=guide_legend(title="Lithology")) +
  ylab("Frequency") +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Zissou1")))
gggeo_scale(a)

# Grain size
test <- as.data.frame(table(g.m.dat$Finalised_grainsize, g.m.dat$bin_midpoint))
names(test) <- c("Finalised_grainsize", "mid_ma", "Freq")
test$mid_ma <- as.numeric(as.character(test$mid_ma))
a <- ggplot(test, aes(x=mid_ma, y=Freq, fill=Finalised_grainsize)) + 
  geom_area() +
  scale_x_reverse() +
  theme_bw() +
  guides(fill=guide_legend(title="Grain size")) +
  ylab("Frequency") +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Zissou1")))
gggeo_scale(a)

# Family
test <- as.data.frame(table(f.m.dat$Family, f.m.dat$bin_midpoint))
names(test) <- c("Family", "mid_ma", "Freq")
test$mid_ma <- as.numeric(as.character(test$mid_ma))
a <- ggplot(test, aes(x=mid_ma, y=Freq, fill=Family)) + 
  geom_area() +
  scale_x_reverse() +
  theme_bw() +
  guides(fill=guide_legend(title="Family")) +
  ylab("Frequency") +
  xlab("Time (Ma)") +
  scale_fill_manual(values=c(wes_palette("Royal1"), 
                             wes_palette("Zissou1"), 
                             wes_palette("Royal2"))) 
gggeo_scale(a)

##### SERIES LEVEL #####

# Plot preservation score per time period (discreet)
ggplot(m.dat.series, aes(x = Preservation_score, fill = bin_assignment)) +
  geom_bar() +
  scale_fill_manual("Legend", values = myColours) +
  ylab("Frequency") +
  xlab("Preservation Score") +
  guides(fill=guide_legend(title="Period")) +
  theme_bw() 

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
a <- as.data.frame(table(m.dat.series$Preservation_score, m.dat.series$bin_assignment))
names(a) <- c("Preservation_score", "bin", "Freq")
a <- merge(a, series, by = 'bin')
a$Preservation_score <- factor(a$Preservation_score, levels = c("5", "4", "3", "2", "1"))
a <- ggplot(a, aes(x=mid_ma, y=Freq, fill=Preservation_score)) + 
  geom_area() +
  scale_x_reverse() +
  theme_bw() +
  guides(fill=guide_legend(title="Preservation score")) +
  ylab("Frequency") +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Zissou1")))
gggeo_scale(a, dat = series)

# Plot proportional preservation scores through time (continuous)
a <- as.data.frame(table(m.dat.series$Preservation_score, m.dat.series$bin_assignment))
names(a) <- c("Preservation_score", "bin", "Freq")
a <- merge(a, series, by = 'bin')
a$Preservation_score <- factor(a$Preservation_score, levels = c("5", "4", "3", "2", "1"))
a <- ggplot(a, aes(x=mid_ma, y=Freq, fill=Preservation_score)) + 
  geom_area(position = 'fill') +
  scale_x_reverse() +
  theme_bw() +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Preservation score")) +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Zissou1"))) 
gggeo_scale(a, dat = series)

# Format into table and run Chi Squared test
test <- chisq.test(table(m.dat.series$Preservation_score, m.dat.series$bin_assignment))
test
test$expected # compare against what the test would have expected

# Mosaic plot
mosaic(~ Preservation_score + bin_assignment,
       direction = c("v", "h"),
       data = m.dat.series,
       shade = TRUE, 
       labeling_args = list(rot_labels = c(0), just_labels = "right", 
                            set_varnames = c(Preservation_score = "Preservation Score",
                                             bin_assignment = ""))
)

################################################################################
# 4. JACCARD SIMILARITY
################################################################################

########################################
##### PRESERVATION SCORE V. GENERA #####
########################################

# Filter to genera and select appropriate columns
jac.1 <- m.dat %>%
  dplyr::filter(Rank == "Genus" | Rank == "Species") %>%
  dplyr::select(Genus, Preservation_score) %>%
  group_by_all() %>%
  summarize(Count = n())

temp.res <- Jaccard.1(jac.1, selection = "Genus", ... = Genus)

######################################################
###### PRESERVATION SCORE V. GENERA THROUGH TIME #####
######################################################

# Setup specimens needed
jac.1 <- m.dat.period %>%
  dplyr::filter(Rank == "Genus" | Rank == "Species") %>%
  dplyr::select(Genus, Preservation_score, bin_assignment) 
jac.period <- unique(jac.1$bin_assignment)

# Run a loop
for(t in jac.period) {
  # Only select period of interest
  jac.2 <- jac.1 %>%
    dplyr::filter(bin_assignment == t) %>%
    dplyr::select(Preservation_score, Genus) %>%
    group_by_all() %>%
    summarize(Count = n())
  
  # Transform to Preservation_score by genus matrix
  jac.3 <- dcast(jac.2,
                Preservation_score ~ Genus,
                value.var = "Count", fill = 0)
  jac.3 <- as.data.frame(jac.3)
  jac.4 <- jac.3[,-1]
  rownames(jac.4) <- jac.3[,1]
  
  # Make binary (presence/absence)
  binary <- jac.4
  binary[binary > 0] <- 1
  
  # Get combinations
  combos <- combn(rownames(binary), 2)
  
  # Loop through combos
  for(a in 1:ncol(combos)){
    temp <- binary[c(combos[1, a],combos[2, a]),]
    temp <- vegdist(temp, method = "jaccard")
    if(a == 1){
      temp.res <- c(paste(combos[1, a], "v", 
                          combos[2, a], sep = ""), 
                    as.numeric(temp))
    } else{
      temp.res <- rbind(temp.res, c(paste(combos[1, a], "v", 
                                          combos[2, a], sep = ""), 
                                    as.numeric(temp)))
    }
  }
  # Save temporary results and organise
  temp.res <- as.data.frame(temp.res)
  temp.res <- cbind(temp.res, t, ncol(binary))
  colnames(temp.res) <- c("Test", "Score", "interval_name", "Genera")
  
  # If it's the first loop, make the results
  if(which(t == jac.period)[[1]] == 1){
    results <- temp.res
  }else{ # Otherwise bind to results
    results <- rbind(results, temp.res)
  }
}

# Load info for matching to intervals
df <- palaeoverse::GTS2020
df <- df[155:159,]

# Attach to results
results <- df %>%
  dplyr::select(interval_name, mid_ma) %>%
  merge(results, by = "interval_name")
results$Score <- as.numeric(results$Score)

# Organise factor levels
order_ind2 <- c("1v2", "2v3", "3v4", "4v5", "1v3", "2v4", "3v5", "1v4", "2v5", "1v5")
results$Test <- as.factor(results$Test)
results$Test <- factor(results$Test, levels = order_ind2)

# Plot the combinations for each preservation level through time
for(i in 1:5){
  test <- results %>%
    filter(grepl(i,Test) == T)
  a <- ggplot(test, aes(x=mid_ma, y=Score)) + 
    geom_line(aes(colour = Test, linetype = Test)) +
    scale_x_reverse() +
    theme_bw() +
    ylab("Jaccard Similarity") +
    xlab("Time (Ma)") +
    geom_point(aes(color = Test))
  a <- a + scale_colour_viridis_d()
  assign(paste("p", i, sep = ""), gggeo_scale(a), envir = .GlobalEnv)
}
ggarrange2(p1, p2, p3, p4, p5, ncol =3, nrow =2)

# Plot combined graph through time
a <- ggplot(results, aes(x=mid_ma, y=Score)) + 
  geom_line(aes(colour = Test, linetype = Test)) +
  scale_x_reverse() +
  theme_bw() +
  ylab("Jaccard Similarity") +
  xlab("Time (Ma)") +
  geom_point(aes(color = Test))
a <- a + scale_colour_viridis_d()
gggeo_scale(a)

########################################
##### PRESERVATION SCORE V. FAMILY #####
########################################

# Filter to genera and select appropriate columns
jac.1 <- f.m.dat %>%
  dplyr::select(Family, Preservation_score) %>%
  group_by_all() %>%
  summarize(Count = n())

temp.res <- Jaccard.1(jac.1, selection = "Family", ... = Family)

######################################################
###### PRESERVATION SCORE V. FAMILY THROUGH TIME #####
######################################################

# Setup specimens needed
jac.1 <- m.dat.period %>%
  dplyr::filter(Rank == "Family" | Rank == "Genus" | Rank == "Species") %>%
  dplyr::select(Family, Preservation_score, bin_assignment) 
jac.period <- unique(jac.1$bin_assignment)

# Run a loop
for(t in jac.period) {
  # Only select period of interest
  jac.2 <- jac.1 %>%
    dplyr::filter(bin_assignment == t) %>%
    dplyr::select(Preservation_score, Family) %>%
    group_by_all() %>%
    summarize(Count = n())
  
  # Transform to Preservation_score by Family matrix
  jac.3 <- dcast(jac.2,
                 Preservation_score ~ Family,
                 value.var = "Count", fill = 0)
  jac.3 <- as.data.frame(jac.3)
  jac.4 <- jac.3[,-1]
  rownames(jac.4) <- jac.3[,1]
  
  # Make binary (presence/absence)
  binary <- jac.4
  binary[binary > 0] <- 1
  
  # Get combinations
  combos <- combn(rownames(binary), 2)
  
  # Loop through combos
  for(a in 1:ncol(combos)){
    temp <- binary[c(combos[1, a],combos[2, a]),]
    temp <- vegdist(temp, method = "jaccard")
    if(a == 1){
      temp.res <- c(paste(combos[1, a], "v", 
                          combos[2, a], sep = ""), 
                    as.numeric(temp))
    } else{
      temp.res <- rbind(temp.res, c(paste(combos[1, a], "v", 
                                          combos[2, a], sep = ""), 
                                    as.numeric(temp)))
    }
  }
  # Save temporary results and organise
  temp.res <- as.data.frame(temp.res)
  temp.res <- cbind(temp.res, t, ncol(binary))
  colnames(temp.res) <- c("Test", "Score", "interval_name", "Family")
  
  # If it's the first loop, make the results
  if(which(t == jac.period)[[1]] == 1){
    results <- temp.res
  }else{ # Otherwise bind to results
    results <- rbind(results, temp.res)
  }
}

# Load info for matching to intervals
df <- palaeoverse::GTS2020
df <- df[155:159,]

# Attach to results
results <- df %>%
  dplyr::select(interval_name, mid_ma) %>%
  merge(results, by = "interval_name")
results$Score <- as.numeric(results$Score)

# Organise factor levels
order_ind2 <- c("1v2", "2v3", "3v4", "4v5", "1v3", "2v4", "3v5", "1v4", "2v5", "1v5")
results$Test <- as.factor(results$Test)
results$Test <- factor(results$Test, levels = order_ind2)

# Plot the combinations for each preservation level through time
for(i in 1:5){
  test <- results %>%
    filter(grepl(i,Test) == T)
  a <- ggplot(test, aes(x=mid_ma, y=Score)) + 
    geom_line(aes(colour = Test, linetype = Test)) +
    scale_x_reverse() +
    theme_bw() +
    ylab("Jaccard Similarity") +
    xlab("Time (Ma)") +
    geom_point(aes(color = Test))
  a <- a + scale_colour_viridis_d()
  assign(paste("p", i, sep = ""), gggeo_scale(a), envir = .GlobalEnv)
}
ggarrange2(p1, p2, p3, p4, p5, ncol =3, nrow =2)

# Plot combined graph through time
a <- ggplot(results, aes(x=mid_ma, y=Score)) + 
  geom_line(aes(colour = Test, linetype = Test)) +
  scale_x_reverse() +
  theme_bw() +
  ylab("Jaccard Similarity") +
  xlab("Time (Ma)") +
  geom_point(aes(color = Test))
a <- a + scale_colour_viridis_d()
gggeo_scale(a)

###############################
##### LITHOLOGY V. GENERA #####
###############################

# Filter to genera and select appropriate columns
jac.1 <- l.m.dat %>%
  dplyr::filter(Rank == "Genus" | Rank == "Species") %>% 
  dplyr::select(Genus, Finalised_lith) %>%
  group_by_all() %>%
  summarize(Count = n())

# Transform to Preservation_score by genus matrix
jac.1 <- dcast(jac.1,
              Finalised_lith ~ Genus,
              value.var = "Count", fill = 0)
jac.1 <- as.data.frame(jac.1)
jac.2 <- jac.1[,-1]
rownames(jac.2) <- jac.1[,1]

# Make binary (presence/absence)
binary <- jac.2
binary[binary > 0] <- 1

## 1 vs. 2 ##
lithology.jac <- binary[c(1,2),]
lithology.jac <- vegdist(lithology.jac, method = "jaccard")

############################################
##### LITHOLOGY V. GENERA THROUGH TIME #####
############################################

# Select relevant specimens
test <- m.dat.period %>%
  dplyr::filter(is.na(Genus) == F) %>%
  filter(is.na(Finalised_lith) == F) %>%
  dplyr::select(Genus, Finalised_lith, bin_assignment) 

# Setup loop
for(t in unique(test$bin_assignment)) {
  
  # Only select period of interest
  temp_test <- test %>%
    dplyr::filter(bin_assignment == t) %>%
    dplyr::select(Finalised_lith, Genus) %>%
    group_by_all() %>%
    summarize(Count = n())
  
  # Transform to Preservation_score by genus matrix
  temp_test <- dcast(temp_test,
                     Finalised_lith ~ Genus,
                     value.var = "Count", fill = 0)
  temp_test <- as.data.frame(temp_test)
  temp_test2 <- temp_test[,-1]
  rownames(temp_test2) <- temp_test[,1]
  
  # Make binary (presence/absence)
  binary <- temp_test2
  binary[binary > 0] <- 1
  
  if(nrow(binary) == 1){
    temp = NA
    temp.res <- cbind(temp, t, ncol(binary))
    colnames(temp.res) <- c("Score", "interval_name", "Genera")
    if(which(t == unique(test$bin_assignment))[[1]] == 1){
      results <- temp.res
    }else{
      results <- rbind(results, temp.res)
    }
    next
  }

  temp <- vegdist(binary, method = "jaccard")
  temp.res <- cbind(temp, t, ncol(binary))
  colnames(temp.res) <- c("Score", "interval_name", "Genera")
  
  if(which(t == unique(test$bin_assignment))[[1]] == 1){
    results <- temp.res
  }else{
    results <- rbind(results, temp.res)
  }
}

# Get data to plot results
df <- palaeoverse::GTS2020
df <- df[155:159,]
results <- df %>%
  dplyr::select(interval_name, mid_ma) %>%
  merge(results, by = "interval_name")
results$Score <- as.numeric(results$Score)

# Plot results
a <- ggplot(results, aes(x=mid_ma, y=Score)) + 
  geom_line() +
  scale_x_reverse() +
  theme_bw() +
  ylab("Jaccard Similarity") +
  xlab("Time (Ma)") +
  geom_point() +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.05))
gggeo_scale(a)

###############################
##### LITHOLOGY V. FAMILY #####
###############################

# Filter to genera and select appropriate columns
jac.1 <- l.m.dat %>%
  dplyr::filter(Rank == "Family" | Rank == "Genus" | Rank == "Species") %>% 
  dplyr::select(Family, Finalised_lith) %>%
  group_by_all() %>%
  summarize(Count = n())

# Transform to Preservation_score by genus matrix
jac.1 <- dcast(jac.1,
               Finalised_lith ~ Family,
               value.var = "Count", fill = 0)
jac.1 <- as.data.frame(jac.1)
jac.2 <- jac.1[,-1]
rownames(jac.2) <- jac.1[,1]

# Make binary (presence/absence)
binary <- jac.2
binary[binary > 0] <- 1

## 1 vs. 2 ##
lithology.jac <- binary[c(1,2),]
lithology.jac <- vegdist(lithology.jac, method = "jaccard")

############################################
##### LITHOLOGY V. FAMILY THROUGH TIME #####
############################################

# Select relevant specimens
test <- m.dat.period %>%
  dplyr::filter(Rank == "Family" | Rank == "Genus" | Rank == "Species") %>%
  filter(is.na(Finalised_lith) == F) %>%
  dplyr::select(Family, Finalised_lith, bin_assignment) 

# Setup loop
for(t in unique(test$bin_assignment)) {
  
  # Only select period of interest
  temp_test <- test %>%
    dplyr::filter(bin_assignment == t) %>%
    dplyr::select(Finalised_lith, Family) %>%
    group_by_all() %>%
    summarize(Count = n())
  
  # Transform to Preservation_score by genus matrix
  temp_test <- dcast(temp_test,
                     Finalised_lith ~ Family,
                     value.var = "Count", fill = 0)
  temp_test <- as.data.frame(temp_test)
  temp_test2 <- temp_test[,-1]
  rownames(temp_test2) <- temp_test[,1]
  
  # Make binary (presence/absence)
  binary <- temp_test2
  binary[binary > 0] <- 1
  
  if(nrow(binary) == 1){
    temp = NA
    temp.res <- cbind(temp, t, ncol(binary))
    colnames(temp.res) <- c("Score", "interval_name", "Family")
    if(which(t == unique(test$bin_assignment))[[1]] == 1){
      results <- temp.res
    }else{
      results <- rbind(results, temp.res)
    }
    next
  }
  
  temp <- vegdist(binary, method = "jaccard")
  temp.res <- cbind(temp, t, ncol(binary))
  colnames(temp.res) <- c("Score", "interval_name", "Family")
  
  if(which(t == unique(test$bin_assignment))[[1]] == 1){
    results <- temp.res
  }else{
    results <- rbind(results, temp.res)
  }
}

# Get data to plot results
df <- palaeoverse::GTS2020
df <- df[155:159,]
results <- df %>%
  dplyr::select(interval_name, mid_ma) %>%
  merge(results, by = "interval_name")
results$Score <- as.numeric(results$Score)

# Plot results
a <- ggplot(results, aes(x=mid_ma, y=Score)) + 
  geom_line() +
  scale_x_reverse() +
  theme_bw() +
  ylab("Jaccard Similarity") +
  xlab("Time (Ma)") +
  geom_point() +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.05))
gggeo_scale(a)

###############################
##### GRAINSIZE V. GENERA #####
###############################

# Filter to genera and select appropriate columns
test <- g.m.dat %>%
  dplyr::filter(Rank == "Genus" | Rank == "Species") %>% 
  dplyr::select(Genus, Finalised_grainsize) %>%
  group_by_all() %>%
  summarize(Count = n())

# Transform to Preservation_score by genus matrix
test <- dcast(test,
              Finalised_grainsize ~ Genus,
              value.var = "Count", fill = 0)
test <- as.data.frame(test)
test2 <- test[,-1]
rownames(test2) <- test[,1]

# Make binary (presence/absence)
binary <- test2
binary[binary > 0] <- 1

## 1 vs. 2 ##
grainsize.jacc <- binary[c(1,2),]
grainsize.jacc <- vegdist(grainsize.jacc, method = "jaccard")

############################################
##### GRAINSIZE V. GENERA THROUGH TIME #####
############################################

test <- m.dat.period %>%
  dplyr::filter(is.na(Genus) == F) %>%
  filter(is.na(Finalised_grainsize) == F) %>%
  dplyr::select(Genus, Finalised_grainsize, bin_assignment) 

test <- simple.grain(test)

for(t in unique(test$bin_assignment)) {
  
  # Only select period of interest
  temp_test <- test %>%
    dplyr::filter(bin_assignment == t) %>%
    dplyr::select(Finalised_grainsize, Genus) %>%
    group_by_all() %>%
    summarize(Count = n())
  
  # Transform to Preservation_score by genus matrix
  temp_test <- dcast(temp_test,
                     Finalised_grainsize ~ Genus,
                     value.var = "Count", fill = 0)
  temp_test <- as.data.frame(temp_test)
  temp_test2 <- temp_test[,-1]
  rownames(temp_test2) <- temp_test[,1]
  
  # Make binary (presence/absence)
  binary <- temp_test2
  binary[binary > 0] <- 1
  
  if(nrow(binary) == 1){
    temp = NA
    temp.res <- cbind(temp, t, ncol(binary))
    colnames(temp.res) <- c("Score", "interval_name", "Genera")
    
    if(which(t == unique(test$bin_assignment))[[1]] == 1){
      results <- temp.res
    }else{
      results <- rbind(results, temp.res)
    }
    next
  }
  
  temp <- vegdist(binary, method = "jaccard")
  temp.res <- cbind(temp, t, ncol(binary))
  colnames(temp.res) <- c("Score", "interval_name", "Genera")
  
  if(which(t == unique(test$bin_assignment))[[1]] == 1){
    results <- temp.res
  }else{
    results <- rbind(results, temp.res)
  }
}

# Get data to plot results
df <- palaeoverse::GTS2020
df <- df[155:159,]
results <- df %>%
  dplyr::select(interval_name, mid_ma) %>%
  merge(results, by = "interval_name")
results$Score <- as.numeric(results$Score)

# Plot results
a <- ggplot(results, aes(x=mid_ma, y=Score)) + 
  geom_line() +
  scale_x_reverse() +
  theme_bw() +
  ylab("Jaccard Similarity") +
  xlab("Time (Ma)") +
  geom_point() +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.05))
gggeo_scale(a)

###############################
##### GRAINSIZE V. FAMILY #####
###############################

# Filter to genera and select appropriate columns
test <- g.m.dat %>%
  dplyr::filter(Rank == "Family" | Rank == "Genus" | Rank == "Species") %>% 
  dplyr::select(Family, Finalised_grainsize) %>%
  group_by_all() %>%
  summarize(Count = n())

# Transform to Preservation_score by Family matrix
test <- dcast(test,
              Finalised_grainsize ~ Family,
              value.var = "Count", fill = 0)
test <- as.data.frame(test)
test2 <- test[,-1]
rownames(test2) <- test[,1]

# Make binary (presence/absence)
binary <- test2
binary[binary > 0] <- 1

## 1 vs. 2 ##
grainsize.jacc <- binary[c(1,2),]
grainsize.jacc <- vegdist(grainsize.jacc, method = "jaccard")

############################################
##### GRAINSIZE V. GENERA THROUGH TIME #####
############################################

test <- m.dat.period %>%
  dplyr::filter(Rank == "Family" | Rank == "Genus" | Rank == "Species") %>%
  filter(is.na(Finalised_grainsize) == F) %>%
  dplyr::select(Family, Finalised_grainsize, bin_assignment) 

test <- simple.grain(test)

for(t in unique(test$bin_assignment)) {
  
  # Only select period of interest
  temp_test <- test %>%
    dplyr::filter(bin_assignment == t) %>%
    dplyr::select(Finalised_grainsize, Family) %>%
    group_by_all() %>%
    summarize(Count = n())
  
  # Transform to Preservation_score by Family matrix
  temp_test <- dcast(temp_test,
                     Finalised_grainsize ~ Family,
                     value.var = "Count", fill = 0)
  temp_test <- as.data.frame(temp_test)
  temp_test2 <- temp_test[,-1]
  rownames(temp_test2) <- temp_test[,1]
  
  # Make binary (presence/absence)
  binary <- temp_test2
  binary[binary > 0] <- 1
  
  if(nrow(binary) == 1){
    temp = NA
    temp.res <- cbind(temp, t, ncol(binary))
    colnames(temp.res) <- c("Score", "interval_name", "Family")
    
    if(which(t == unique(test$bin_assignment))[[1]] == 1){
      results <- temp.res
    }else{
      results <- rbind(results, temp.res)
    }
    next
  }
  
  temp <- vegdist(binary, method = "jaccard")
  temp.res <- cbind(temp, t, ncol(binary))
  colnames(temp.res) <- c("Score", "interval_name", "Family")
  
  if(which(t == unique(test$bin_assignment))[[1]] == 1){
    results <- temp.res
  }else{
    results <- rbind(results, temp.res)
  }
}

# Get data to plot results
df <- palaeoverse::GTS2020
df <- df[155:159,]
results <- df %>%
  dplyr::select(interval_name, mid_ma) %>%
  merge(results, by = "interval_name")
results$Score <- as.numeric(results$Score)

# Plot results
a <- ggplot(results, aes(x=mid_ma, y=Score)) + 
  geom_line() +
  scale_x_reverse() +
  theme_bw() +
  ylab("Jaccard Similarity") +
  xlab("Time (Ma)") +
  geom_point() +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.05))
gggeo_scale(a)

################################################################################
# 5. RUNNING MODELS
################################################################################

###################
##### SETUP 1 #####
###################

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
  dplyr::select(Family, Genus, Species, Rank, Museum_Number, Preservation_score, Continent, 
                Country, lat, lng, Formation, Age, Max_period, Min_period, Finalised_lith, 
                Finalised_grainsize, max_ma, min_ma, bin_assignment, interval_mid_ma, p_lng, p_lat)


###################
##### SETUP 2 #####
###################

model.data <- m.dat %>%
  filter(is.na(Finalised_grainsize) == F) %>%
  filter(is.na(Finalised_lith) == F)

# Assign specific grain size categories into either "Fine Grained" or "Coarse Grained"
model.data <- simple.grain(model.data)

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
  dplyr::select(Family, Genus, Species, Rank, Museum_Number, Preservation_score, Continent, 
                Country, lat, lng, Formation, Age, Max_period, Min_period, Finalised_lith, 
                Finalised_grainsize, max_ma, min_ma, bin_assignment, interval_mid_ma)


#######################################
##### ORDINAL LOGISTIC REGRESSION #####
#######################################

# Set up modelling
samplesize <- 0.6*nrow(data.set)
index <- sample(seq_len(nrow(data.set)), size = samplesize)

#Creating training and test set 
datatrain <- data.set[index,]
table(datatrain$Finalised_lith, datatrain$Finalised_grainsize)
table(datatrain$Finalised_lith, datatrain$Preservation_score)
table(datatrain$Finalised_grainsize, datatrain$Preservation_score, datatrain$Finalised_lith)

# Get test data
datatest <- data.set[-index,]

# Set full model
full.model <- polr(Preservation_score ~ Finalised_lith + Finalised_grainsize + Family + 
             lat + p_lat + interval_mid_ma + Finalised_lith*Finalised_grainsize, 
           data = datatrain, Hess=TRUE)
summary(full.model)

# Alter global actions to allow for dredge
options(na.action = "na.fail") 

# Get all models, ranked
model.set <- dredge(full.model)

# Revert global actions
options(na.action = "na.omit")

# Get best model
m <- get.models(model.set, subset = 1)

## view a summary of the model
summary(m[[1]]) 
(ctable <- coef(summary(m[[1]])))

## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

## combined table
(ctable <- cbind(ctable, "p value" = p))
ci <- confint.default(m[[1]])

## odds ratios
exp(coef(m[[1]]))

## OR and CI
test <-exp(cbind(OR = coef(m[[1]]), ci))

# Model quality
predictPreservation_score <- predict(m[[1]],datatest)
table(datatest$Preservation_score, predictPreservation_score)
predictPreservation_score <- as.character(predictPreservation_score)
predictPreservation_score[is.na(predictPreservation_score)] <- "0"
mean(as.character(datatest$Preservation_score) != as.character(predictPreservation_score))

# Plotting 
plot_model(m[[1]], title = "Preservation Score")

plot(effects::Effect(focal.predictors = c("Finalised_grainsize", "Finalised_lith"),m[[1]]),
     main="Ordinal logistic regression; top model, showing lithology + grain size",
     axes = list(y=list(lab="Preservation Score (probability)")),
     xlab = "")
plot(effects::Effect(focal.predictors = c("lat", "Finalised_lith"),m[[1]]),
     main="Ordinal logistic regression; top model, showing latitude + lithology",
     axes = list(y=list(lab="Preservation Score (probability)")), 
     xlab = "")

###############################
##### LOGISTIC REGRESSION #####
###############################

# Set preservation level to explore and get dataset
pres.score <- 5
data.LR <- set_Pres_score(data.set, c(pres.score))
colnames(data.LR)[colnames(data.LR) == "lat"] <- "Latitude"
colnames(data.LR)[colnames(data.LR) == "p_lat"] <- "Palaeo-latitude"
colnames(data.LR)[colnames(data.LR) == "interval_mid_ma"] <- "Age (Ma)"
colnames(data.LR)[colnames(data.LR) == "Finalised_lith"] <- "Lithology"
colnames(data.LR)[colnames(data.LR) == "Finalised_grainsize"] <- "Grainsize"

# Set full model
full.model <- glm(formula = LR_Pres_score ~ Lithology + Grainsize + Family + 
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

# Save covariates to remove for failing
(to.remove <- row.names(ctable)[ctable$`Std. Error` > 10 |ctable$`Std. Error` < -10])

# Save order for changing colour
ctable$color <- 0
ctable$color[ctable$Estimate > 0] <- '#377EB8'
ctable$color[ctable$Estimate < 0] <- '#E41A1C'
ctable$color[ctable$`Pr(>|z|)` > 0.05] <- '#A9A9A9'
my.colors <- ctable$color[ctable$`Std. Error` < 10 & ctable$`Std. Error` > -10]
my.colors <- my.colors[-1]
my.colors.2 <- as.numeric(as.factor(my.colors))
my.colors.3 <- c('#377EB8','#A9A9A9','#E41A1C')

# Odds Ratio plot
(p <- plot_model(m[[1]], 
           show.values = TRUE, 
           value.offset = 0.3,
           rm.terms = to.remove,
           title = paste("Preservation Score: ", pres.score, sep = ""),
           vline.color = 'red',
           group.terms = my.colors.2,
           color = my.colors.3
           ))

p <- p + ggtitle(paste("Preservation Score: ", pres.score, sep = "")) +
  labs(subtitle = paste("Model = Preservation score ~ ", best.form[3], sep = ""))

# Save the plot
pdf(paste("Plot_Odds_ratio_pres_", pres.score, ".pdf", sep = ""), 
    width = par("din")[1], 
    height = par("din")[2])
print(p)
dev.off()

# Log-odds plot
(p <- plot_model(m[[1]], 
           show.values = TRUE, 
           value.offset = 0.3,
           rm.terms = to.remove,
           group.terms = my.colors.2,
           color = my.colors.3,
           transform = NULL
           ))

p <- p + ggtitle(paste("Preservation Score: ", pres.score, sep = "")) +
  labs(subtitle = paste("Model = Preservation score ~ ", best.form[3], sep = ""))

# Save the plot
pdf(paste("Plot_Log_odds_pres_", pres.score, ".pdf", sep = ""), 
    width = par("din")[1], 
    height = par("din")[2])
print(p)
dev.off()

# Probability plot
(p <- plot_model(m[[1]], 
                 show.values = TRUE, 
                 value.offset = 0.3,
                 rm.terms = to.remove,
                 group.terms = my.colors.2,
                 color = my.colors.3,
                 transform = "plogis"
))

p <- p + ggtitle(paste("Preservation Score: ", pres.score, sep = "")) +
  labs(subtitle = paste("Model = Preservation score ~ ", best.form[3], sep = ""))

# Save the plot
pdf(paste("Plot_Probability_pres_", pres.score, ".pdf", sep = ""), 
    width = par("din")[1], 
    height = par("din")[2])
print(p)
dev.off()

# Save models
write.csv(model.set, 
          file = paste("LR_all_models_pres_", pres.score, ".csv", sep = ""), 
          row.names = FALSE)
write.csv(ctable, 
          file = paste("LR_best_model_pres_", pres.score, ".csv", sep = ""), 
          row.names = TRUE)

################################################################################
# 7. CORRELATION TESTS
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
  dplyr::select(Genus, Species, Locality, bin_midpoint, Period, Period.no, Preservation_score) %>%
  mutate(Combined_name = paste(Genus, Species, sep = " "))

genus.lvl <- m.dat.period %>%
  dplyr::filter(Rank == "Genus") %>%
  dplyr::select(Genus, Locality, bin_midpoint, Period, Period.no, Preservation_score) 

all.lvl <- m.dat.period %>%
  dplyr::select(Locality, bin_midpoint, Period, Period.no, Preservation_score)

##### Make tables of number of taxa per preservation score for each time frame #####

# All specimens
stage.pres.all <- as.data.frame(table(all.lvl$bin_midpoint, all.lvl$Preservation_score))
colnames(stage.pres.all) <- c("bin_midpoint", "Preservation_score", "Freq")
period.pres.all <- as.data.frame(table(all.lvl$Period, all.lvl$Preservation_score))
# Species
stage.pres.spec <- as.data.frame(table(species.lvl$bin_midpoint, species.lvl$Preservation_score))
colnames(stage.pres.spec) <- c("bin_midpoint", "Preservation_score", "Freq")
period.pres.spec <- as.data.frame(table(species.lvl$Period, species.lvl$Preservation_score))
# Genera
stage.pres.gen <- as.data.frame(table(genus.lvl$bin_midpoint, genus.lvl$Preservation_score))
colnames(stage.pres.gen) <- c("bin_midpoint", "Preservation_score", "Freq")
period.pres.gen <- as.data.frame(table(genus.lvl$Period, genus.lvl$Preservation_score))

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
order_ind <- c("Permian", "Pennsylvanian", "Mississippian", "Devonian", "Silurian", "Ordovician")

# Run function for each level
split(period.pres.gen, period.pres.gen$Var2, order = TRUE)
split(period.pres.spec, period.pres.spec$Var2, order = TRUE)
split(period.pres.all, period.pres.all$Var2, order = TRUE)
split(stage.pres.all, Var2 = stage.pres.all$Preservation_score, order = FALSE, stage = TRUE)
split(stage.pres.spec, Var2 = stage.pres.spec$Preservation_score, order = FALSE, stage = TRUE)
split(stage.pres.gen, Var2 = stage.pres.gen$Preservation_score, order = FALSE, stage = TRUE)

temp_stages <- stages[55:92,]
temp_stages$bin_midpoint <- (temp_stages$max_ma + temp_stages$min_ma)/2 
div.stage.spec <- merge(temp_stages, div.stage.spec, by = "bin_midpoint", all.x = T)
div.stage.gen <- merge(temp_stages, div.stage.gen, by = "bin_midpoint", all.x = T)
colnames(Stage.colls)[1] <- "bin_midpoint"
col.stage.all <- merge(temp_stages, Stage.colls, by = "bin_midpoint", all.x = T)

##### Get comprehensive collections through time #####

# Load PBDB Palaeozoic invert occurrences
all.pbdb <- read.csv("all_palaeozoic_binned.csv")

# Filter to correct age
all.pbdb.period <- all.pbdb %>%
  filter(max_ma < 485.4) %>%
  filter(min_ma > 251.902)

# Bin
all.pbdb.period <- bin_time(all.pbdb.period, series, method = 'majority')

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

##################################
##### DIVERSITY CORRELATIONS #####
##################################

##### SPECIES #####

# Setup complete list of bins
div.stage.spec.comp <- merge(stage.pres.all.1, div.stage.spec, by = "bin_midpoint", all.x = T)

# Stage - diversity vs.preservation score
stage.spec.div.1 <- cor.test(log10(div.stage.spec.comp$SIBs), log10(stage.pres.all.1$Freq), method = 'spearman')
stage.spec.div.2 <- cor.test(log10(div.stage.spec.comp$SIBs), log10(stage.pres.all.2$Freq), method = 'spearman')
stage.spec.div.3 <- cor.test(log10(div.stage.spec.comp$SIBs), log10(stage.pres.all.3$Freq), method = 'spearman')
stage.spec.div.4 <- cor.test(log10(div.stage.spec.comp$SIBs), log10(stage.pres.all.4$Freq), method = 'spearman')
stage.spec.div.5 <- cor.test(log10(div.stage.spec.comp$SIBs), log10(stage.pres.all.5$Freq), method = 'spearman')

stage.div.spec <- rbind(c("Stage", 1, "Diversity (species)", stage.spec.div.1$estimate, stage.spec.div.1$p.value),
                        c("Stage", 2, "Diversity (species)", stage.spec.div.2$estimate, stage.spec.div.2$p.value),
                        c("Stage", 3, "Diversity (species)", stage.spec.div.3$estimate, stage.spec.div.3$p.value),
                        c("Stage", 4, "Diversity (species)", stage.spec.div.4$estimate, stage.spec.div.4$p.value),
                        c("Stage", 5, "Diversity (species)", stage.spec.div.5$estimate, stage.spec.div.5$p.value))

# Period - diversity vs. preservation score
period.spec.div.1 <- cor.test(log10(div.period.spec$SIBs), log10(period.pres.all.1$Freq), method = 'spearman')
period.spec.div.2 <- cor.test(log10(div.period.spec$SIBs), log10(period.pres.all.2$Freq), method = 'spearman')
period.spec.div.3 <- cor.test(log10(div.period.spec$SIBs), log10(period.pres.all.3$Freq), method = 'spearman')
period.spec.div.4 <- cor.test(log10(div.period.spec$SIBs), log10(period.pres.all.4$Freq), method = 'spearman')
period.spec.div.5 <- cor.test(log10(div.period.spec$SIBs), log10(period.pres.all.5$Freq), method = 'spearman')

period.div.spec <- rbind(c("Period", 1, "Diversity (species)", period.spec.div.1$estimate, period.spec.div.1$p.value),
                         c("Period", 2, "Diversity (species)", period.spec.div.2$estimate, period.spec.div.2$p.value),
                         c("Period", 3, "Diversity (species)", period.spec.div.3$estimate, period.spec.div.3$p.value),
                         c("Period", 4, "Diversity (species)", period.spec.div.4$estimate, period.spec.div.4$p.value),
                         c("Period", 5, "Diversity (species)", period.spec.div.5$estimate, period.spec.div.5$p.value))

##### GENERA #####

# Setup complete list of bins
div.stage.gen.comp <- merge(stage.pres.all.1, div.stage.gen, by = "bin_midpoint", all.x = T)

# Stage - diversity vs.preservation score
stage.gen.div.1 <- cor.test(log10(div.stage.gen.comp$SIBs), log10(stage.pres.all.1$Freq), method = 'spearman')
stage.gen.div.2 <- cor.test(log10(div.stage.gen.comp$SIBs), log10(stage.pres.all.2$Freq), method = 'spearman')
stage.gen.div.3 <- cor.test(log10(div.stage.gen.comp$SIBs), log10(stage.pres.all.3$Freq), method = 'spearman')
stage.gen.div.4 <- cor.test(log10(div.stage.gen.comp$SIBs), log10(stage.pres.all.4$Freq), method = 'spearman')
stage.gen.div.5 <- cor.test(log10(div.stage.gen.comp$SIBs), log10(stage.pres.all.5$Freq), method = 'spearman')

stage.div.gen <-  rbind(c("Stage", 1, "Diversity (genera)",stage.gen.div.1$estimate, stage.gen.div.1$p.value),
                        c("Stage", 2, "Diversity (genera)",stage.gen.div.2$estimate, stage.gen.div.2$p.value),
                        c("Stage", 3, "Diversity (genera)",stage.gen.div.3$estimate, stage.gen.div.3$p.value),
                        c("Stage", 4, "Diversity (genera)",stage.gen.div.4$estimate, stage.gen.div.4$p.value),
                        c("Stage", 5, "Diversity (genera)",stage.gen.div.5$estimate, stage.gen.div.5$p.value))

# Period - diversity vs. preservation score
period.gen.div.1 <- cor.test(log10(div.period.gen$SIBs), log10(period.pres.all.1$Freq), method = 'spearman')
period.gen.div.2 <- cor.test(log10(div.period.gen$SIBs), log10(period.pres.all.2$Freq), method = 'spearman')
period.gen.div.3 <- cor.test(log10(div.period.gen$SIBs), log10(period.pres.all.3$Freq), method = 'spearman')
period.gen.div.4 <- cor.test(log10(div.period.gen$SIBs), log10(period.pres.all.4$Freq), method = 'spearman')
period.gen.div.5 <- cor.test(log10(div.period.gen$SIBs), log10(period.pres.all.5$Freq), method = 'spearman')

period.div.gen <-  rbind(c("Period", 1, "Diversity (genera)", period.gen.div.1$estimate, period.gen.div.1$p.value),
                         c("Period", 2, "Diversity (genera)", period.gen.div.2$estimate, period.gen.div.2$p.value),
                         c("Period", 3, "Diversity (genera)", period.gen.div.3$estimate, period.gen.div.3$p.value),
                         c("Period", 4, "Diversity (genera)", period.gen.div.4$estimate, period.gen.div.4$p.value),
                         c("Period", 5, "Diversity (genera)", period.gen.div.5$estimate, period.gen.div.5$p.value))

####################################
##### COLLECTIONS CORRELATIONS #####
####################################

colls.stage.comp <- merge(stage.pres.all.1, Stage.colls, by = "bin_midpoint", all.x = T)

# Stage - echinoid collections vs. preservation score
stage.col.1 <- cor.test(log10(colls.stage.comp$Freq.y), log10(stage.pres.all.1$Freq), method = 'spearman')
stage.col.2 <- cor.test(log10(colls.stage.comp$Freq.y), log10(stage.pres.all.2$Freq), method = 'spearman')
stage.col.3 <- cor.test(log10(colls.stage.comp$Freq.y), log10(stage.pres.all.3$Freq), method = 'spearman')
stage.col.4 <- cor.test(log10(colls.stage.comp$Freq.y), log10(stage.pres.all.4$Freq), method = 'spearman')
stage.col.5 <- cor.test(log10(colls.stage.comp$Freq.y), log10(stage.pres.all.5$Freq), method = 'spearman')

stage.coll.ech <-  rbind(c("Stage", 1, "Echinoid collections", stage.col.1$estimate, stage.col.1$p.value),
                         c("Stage", 2, "Echinoid collections", stage.col.2$estimate, stage.col.2$p.value),
                         c("Stage", 3, "Echinoid collections", stage.col.3$estimate, stage.col.3$p.value),
                         c("Stage", 4, "Echinoid collections", stage.col.4$estimate, stage.col.4$p.value),
                         c("Stage", 5, "Echinoid collections", stage.col.5$estimate, stage.col.5$p.value))

Period.colls <- Period.colls[match(order_ind, Period.colls$Var1),]

# Period - echinoid collections vs. preservation score
period.col.1 <- cor.test(log10(Period.colls$Freq), log10(period.pres.all.1$Freq), method = 'spearman')
period.col.2 <- cor.test(log10(Period.colls$Freq), log10(period.pres.all.2$Freq), method = 'spearman')
period.col.3 <- cor.test(log10(Period.colls$Freq), log10(period.pres.all.3$Freq), method = 'spearman')
period.col.4 <- cor.test(log10(Period.colls$Freq), log10(period.pres.all.4$Freq), method = 'spearman')
period.col.5 <- cor.test(log10(Period.colls$Freq), log10(period.pres.all.5$Freq), method = 'spearman')

period.coll.ech <-  rbind(c("Period", 1, "Echinoid collections", period.col.1$estimate, period.col.1$p.value),
                          c("Period", 2, "Echinoid collections", period.col.2$estimate, period.col.2$p.value),
                          c("Period", 3, "Echinoid collections", period.col.3$estimate, period.col.3$p.value),
                          c("Period", 4, "Echinoid collections", period.col.4$estimate, period.col.4$p.value),
                          c("Period", 5, "Echinoid collections", period.col.5$estimate, period.col.5$p.value))

# Stage - global collections vs. preservation score
stage.pbdb.1 <- cor.test(log10(pbdb.stat$colls), log10(stage.pres.all.1$Freq), method = 'spearman')
stage.pbdb.2 <- cor.test(log10(pbdb.stat$colls), log10(stage.pres.all.2$Freq), method = 'spearman')
stage.pbdb.3 <- cor.test(log10(pbdb.stat$colls), log10(stage.pres.all.3$Freq), method = 'spearman')
stage.pbdb.4 <- cor.test(log10(pbdb.stat$colls), log10(stage.pres.all.4$Freq), method = 'spearman')
stage.pbdb.5 <- cor.test(log10(pbdb.stat$colls), log10(stage.pres.all.5$Freq), method = 'spearman')

stage.coll.pbdb <-  rbind(c("Stage", 1, "Global collections (pbdb)", stage.pbdb.1$estimate, stage.pbdb.1$p.value),
                          c("Stage", 2, "Global collections (pbdb)", stage.pbdb.2$estimate, stage.pbdb.2$p.value),
                          c("Stage", 3, "Global collections (pbdb)", stage.pbdb.3$estimate, stage.pbdb.3$p.value),
                          c("Stage", 4, "Global collections (pbdb)", stage.pbdb.4$estimate, stage.pbdb.4$p.value),
                          c("Stage", 5, "Global collections (pbdb)", stage.pbdb.5$estimate, stage.pbdb.5$p.value))

# Period - global collections vs. preservation score
period.pbdb.1 <- cor.test(log10(pbdb.period.stat$colls), log10(period.pres.all.1$Freq), method = 'spearman')
period.pbdb.2 <- cor.test(log10(pbdb.period.stat$colls), log10(period.pres.all.2$Freq), method = 'spearman')
period.pbdb.3 <- cor.test(log10(pbdb.period.stat$colls), log10(period.pres.all.3$Freq), method = 'spearman')
period.pbdb.4 <- cor.test(log10(pbdb.period.stat$colls), log10(period.pres.all.4$Freq), method = 'spearman')
period.pbdb.5 <- cor.test(log10(pbdb.period.stat$colls), log10(period.pres.all.5$Freq), method = 'spearman')

period.coll.pbdb <- rbind(c("Period", 1, "Global collections (pbdb)", period.pbdb.1$estimate, period.pbdb.1$p.value),
                          c("Period", 2, "Global collections (pbdb)", period.pbdb.2$estimate, period.pbdb.2$p.value),
                          c("Period", 3, "Global collections (pbdb)", period.pbdb.3$estimate, period.pbdb.3$p.value),
                          c("Period", 4, "Global collections (pbdb)", period.pbdb.4$estimate, period.pbdb.4$p.value),
                          c("Period", 5, "Global collections (pbdb)", period.pbdb.5$estimate, period.pbdb.5$p.value))

a <- ggplot(pbdb.stat, aes(x=bin_midpoint, y=colls)) + 
  geom_line() +
  scale_x_reverse() +
  theme_bw() +
  ylab("Count") +
  xlab("Time (Ma)") 
gggeo_scale(a)

a <- ggplot(pbdb.period.stat, aes(x=bin_midpoint, y=colls)) + 
  geom_line() +
  scale_x_reverse() +
  theme_bw() +
  ylab("Count") +
  xlab("Time (Ma)") 
gggeo_scale(a)

#################################
##### SED TYPE CORRELATIONS #####
#################################

##### SETUP: MACROSTRAT #####

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

##### PLOTS #####

a <- ggplot(macro.area, aes(x=bin_midpoint, y=count)) + 
  geom_line(aes(colour = lith)) +
  scale_x_reverse() +
  theme_bw() +
  ylab("Count") +
  xlab("Time (Ma)") +
  geom_point(aes(color = lith)) 
gggeo_scale(a)

a <- ggplot(macro.count, aes(x=bin_midpoint, y=count)) + 
  geom_line(aes(colour = lith)) +
  scale_x_reverse() +
  theme_bw() +
  ylab("Count") +
  xlab("Time (Ma)") +
  geom_point(aes(color = lith)) 
gggeo_scale(a)

##### SET UP #####

# Get appropriate specimens
all.lvl <- m.dat.period %>%
  dplyr::filter(Continent == "North America") %>%
  dplyr::select(Locality, bin_midpoint, Period, Period.no, Preservation_score)

# Find tallys of Preservation score through time
stage.pres.all <- as.data.frame(table(all.lvl$bin_midpoint, all.lvl$Preservation_score))
colnames(stage.pres.all) <- c("bin_midpoint", "Preservation_score", "Freq")
period.pres.all <- as.data.frame(table(all.lvl$Period, all.lvl$Preservation_score))

# Split into datasets
split(period.pres.all, period.pres.all$Var2, order = TRUE)
split(stage.pres.all, Var2 = stage.pres.all$Preservation_score, order = FALSE, stage = TRUE)

##### STAGE CORRELATIONS #####

# Taphonomic grade vs. carbonate (count), stage level
stage.carb.macro.count.1 <- cor.test(log10(carb.macro.count$count), log10(stage.pres.all.1$Freq), method = 'spearman')
stage.carb.macro.count.2 <- cor.test(log10(carb.macro.count$count), log10(stage.pres.all.2$Freq), method = 'spearman')
stage.carb.macro.count.3 <- cor.test(log10(carb.macro.count$count), log10(stage.pres.all.3$Freq), method = 'spearman')
stage.carb.macro.count.4 <- cor.test(log10(carb.macro.count$count), log10(stage.pres.all.4$Freq), method = 'spearman')
stage.carb.macro.count.5 <- cor.test(log10(carb.macro.count$count), log10(stage.pres.all.5$Freq), method = 'spearman')

stage.macro.carb.count <-  rbind(c("Stage", 1, "Macrostrat carbonate (count)", stage.carb.macro.count.1$estimate, stage.carb.macro.count.1$p.value),
                                 c("Stage", 2, "Macrostrat carbonate (count)", stage.carb.macro.count.2$estimate, stage.carb.macro.count.2$p.value),
                                 c("Stage", 3, "Macrostrat carbonate (count)", stage.carb.macro.count.3$estimate, stage.carb.macro.count.3$p.value),
                                 c("Stage", 4, "Macrostrat carbonate (count)", stage.carb.macro.count.4$estimate, stage.carb.macro.count.4$p.value),
                                 c("Stage", 5, "Macrostrat carbonate (count)", stage.carb.macro.count.5$estimate, stage.carb.macro.count.5$p.value))

# Taphonomic grade vs. siliclastic (count), stage level
stage.sili.macro.count.1 <- cor.test(log10(sili.macro.count$count), log10(stage.pres.all.1$Freq), method = 'spearman')
stage.sili.macro.count.2 <- cor.test(log10(sili.macro.count$count), log10(stage.pres.all.2$Freq), method = 'spearman')
stage.sili.macro.count.3 <- cor.test(log10(sili.macro.count$count), log10(stage.pres.all.3$Freq), method = 'spearman')
stage.sili.macro.count.4 <- cor.test(log10(sili.macro.count$count), log10(stage.pres.all.4$Freq), method = 'spearman')
stage.sili.macro.count.5 <- cor.test(log10(sili.macro.count$count), log10(stage.pres.all.5$Freq), method = 'spearman')

stage.macro.sili.count <-  rbind(c("Stage", 1, "Macrostrat silicilcastic (count)", stage.sili.macro.count.1$estimate, stage.sili.macro.count.1$p.value),
                           c("Stage", 2, "Macrostrat siliciclastic (count)", stage.sili.macro.count.2$estimate, stage.sili.macro.count.2$p.value),
                           c("Stage", 3, "Macrostrat siliciclastic (count)", stage.sili.macro.count.3$estimate, stage.sili.macro.count.3$p.value),
                           c("Stage", 4, "Macrostrat siliciclastic (count)", stage.sili.macro.count.4$estimate, stage.sili.macro.count.4$p.value),
                           c("Stage", 5, "Macrostrat siliciclastic (count)", stage.sili.macro.count.5$estimate, stage.sili.macro.count.5$p.value))

# Taphonomic grade vs. carbonate (area), stage level
stage.carb.macro.area.1 <- cor.test(log10(carb.macro.area$count), log10(stage.pres.all.1$Freq), method = 'spearman')
stage.carb.macro.area.2 <- cor.test(log10(carb.macro.area$count), log10(stage.pres.all.2$Freq), method = 'spearman')
stage.carb.macro.area.3 <- cor.test(log10(carb.macro.area$count), log10(stage.pres.all.3$Freq), method = 'spearman')
stage.carb.macro.area.4 <- cor.test(log10(carb.macro.area$count), log10(stage.pres.all.4$Freq), method = 'spearman')
stage.carb.macro.area.5 <- cor.test(log10(carb.macro.area$count), log10(stage.pres.all.5$Freq), method = 'spearman')

stage.macro.carb.area <-   rbind(c("Stage", 1, "Macrostrat carbonate (area)", stage.carb.macro.area.1$estimate, stage.carb.macro.area.1$p.value),
                                 c("Stage", 2, "Macrostrat carbonate (area)", stage.carb.macro.area.2$estimate, stage.carb.macro.area.2$p.value),
                                 c("Stage", 3, "Macrostrat carbonate (area)", stage.carb.macro.area.3$estimate, stage.carb.macro.area.3$p.value),
                                 c("Stage", 4, "Macrostrat carbonate (area)", stage.carb.macro.area.4$estimate, stage.carb.macro.area.4$p.value),
                                 c("Stage", 5, "Macrostrat carbonate (area)", stage.carb.macro.area.5$estimate, stage.carb.macro.area.5$p.value))

# Taphonomic grade vs. siliclastic (area), stage level
stage.sili.macro.area.1 <- cor.test(log10(sili.macro.area$count), log10(stage.pres.all.1$Freq), method = 'spearman')
stage.sili.macro.area.2 <- cor.test(log10(sili.macro.area$count), log10(stage.pres.all.2$Freq), method = 'spearman')
stage.sili.macro.area.3 <- cor.test(log10(sili.macro.area$count), log10(stage.pres.all.3$Freq), method = 'spearman')
stage.sili.macro.area.4 <- cor.test(log10(sili.macro.area$count), log10(stage.pres.all.4$Freq), method = 'spearman')
stage.sili.macro.area.5 <- cor.test(log10(sili.macro.area$count), log10(stage.pres.all.5$Freq), method = 'spearman')

stage.macro.sili.area  <-  rbind(c("Stage", 1, "Macrostrat siliciclastic (area)", stage.sili.macro.area.1$estimate, stage.sili.macro.area.1$p.value),
                                 c("Stage", 2, "Macrostrat siliciclastic (area)", stage.sili.macro.area.2$estimate, stage.sili.macro.area.2$p.value),
                                 c("Stage", 3, "Macrostrat siliciclastic (area)", stage.sili.macro.area.3$estimate, stage.sili.macro.area.3$p.value),
                                 c("Stage", 4, "Macrostrat siliciclastic (area)", stage.sili.macro.area.4$estimate, stage.sili.macro.area.4$p.value),
                                 c("Stage", 5, "Macrostrat siliciclastic (area)", stage.sili.macro.area.5$estimate, stage.sili.macro.area.5$p.value))

##### PERIOD CORRELATIONS #####

# Taphonomic grade vs. carbonate (count), period level
period.carb.macro.count.1 <- cor.test(log10(carb.macro.count.period$count), log10(period.pres.all.1$Freq), method = 'spearman')
period.carb.macro.count.2 <- cor.test(log10(carb.macro.count.period$count), log10(period.pres.all.2$Freq), method = 'spearman')
period.carb.macro.count.3 <- cor.test(log10(carb.macro.count.period$count), log10(period.pres.all.3$Freq), method = 'spearman')
period.carb.macro.count.4 <- cor.test(log10(carb.macro.count.period$count), log10(period.pres.all.4$Freq), method = 'spearman')
period.carb.macro.count.5 <- cor.test(log10(carb.macro.count.period$count), log10(period.pres.all.5$Freq), method = 'spearman')

period.macro.carb.count <- rbind(c("Period", 1, "Macrostrat carbonate (count)", period.carb.macro.count.1$estimate, period.carb.macro.count.1$p.value),
                                 c("Period", 2, "Macrostrat carbonate (count)", period.carb.macro.count.2$estimate, period.carb.macro.count.2$p.value),
                                 c("Period", 3, "Macrostrat carbonate (count)", period.carb.macro.count.3$estimate, period.carb.macro.count.3$p.value),
                                 c("Period", 4, "Macrostrat carbonate (count)", period.carb.macro.count.4$estimate, period.carb.macro.count.4$p.value),
                                 c("Period", 5, "Macrostrat carbonate (count)", period.carb.macro.count.5$estimate, period.carb.macro.count.5$p.value))

# Taphonomic grade vs. siliclastic (count), period level
period.sili.macro.count.1 <- cor.test(log10(sili.macro.count.period$count), log10(period.pres.all.1$Freq), method = 'spearman')
period.sili.macro.count.2 <- cor.test(log10(sili.macro.count.period$count), log10(period.pres.all.2$Freq), method = 'spearman')
period.sili.macro.count.3 <- cor.test(log10(sili.macro.count.period$count), log10(period.pres.all.3$Freq), method = 'spearman')
period.sili.macro.count.4 <- cor.test(log10(sili.macro.count.period$count), log10(period.pres.all.4$Freq), method = 'spearman')
period.sili.macro.count.5 <- cor.test(log10(sili.macro.count.period$count), log10(period.pres.all.5$Freq), method = 'spearman')

period.macro.sili.count <- rbind(c("Period", 1, "Macrostrat silicilcastic (count)", period.sili.macro.count.1$estimate, period.sili.macro.count.1$p.value),
                                 c("Period", 2, "Macrostrat siliciclastic (count)", period.sili.macro.count.2$estimate, period.sili.macro.count.2$p.value),
                                 c("Period", 3, "Macrostrat siliciclastic (count)", period.sili.macro.count.3$estimate, period.sili.macro.count.3$p.value),
                                 c("Period", 4, "Macrostrat siliciclastic (count)", period.sili.macro.count.4$estimate, period.sili.macro.count.4$p.value),
                                 c("Period", 5, "Macrostrat siliciclastic (count)", period.sili.macro.count.5$estimate, period.sili.macro.count.5$p.value))

# Taphonomic grade vs. carbonate (area), period level
period.carb.macro.area.1 <- cor.test(log10(carb.macro.area.period$count), log10(period.pres.all.1$Freq), method = 'spearman')
period.carb.macro.area.2 <- cor.test(log10(carb.macro.area.period$count), log10(period.pres.all.2$Freq), method = 'spearman')
period.carb.macro.area.3 <- cor.test(log10(carb.macro.area.period$count), log10(period.pres.all.3$Freq), method = 'spearman')
period.carb.macro.area.4 <- cor.test(log10(carb.macro.area.period$count), log10(period.pres.all.4$Freq), method = 'spearman')
period.carb.macro.area.5 <- cor.test(log10(carb.macro.area.period$count), log10(period.pres.all.5$Freq), method = 'spearman')

period.macro.carb.area <-   rbind(c("Period", 1, "Macrostrat carbonate (area)", period.carb.macro.area.1$estimate, period.carb.macro.area.1$p.value),
                                 c("Period", 2, "Macrostrat carbonate (area)", period.carb.macro.area.2$estimate, period.carb.macro.area.2$p.value),
                                 c("Period", 3, "Macrostrat carbonate (area)", period.carb.macro.area.3$estimate, period.carb.macro.area.3$p.value),
                                 c("Period", 4, "Macrostrat carbonate (area)", period.carb.macro.area.4$estimate, period.carb.macro.area.4$p.value),
                                 c("Period", 5, "Macrostrat carbonate (area)", period.carb.macro.area.5$estimate, period.carb.macro.area.5$p.value))

# Taphonomic grade vs. siliclastic (area), period level
period.sili.macro.area.1 <- cor.test(log10(sili.macro.area.period$count), log10(period.pres.all.1$Freq), method = 'spearman')
period.sili.macro.area.2 <- cor.test(log10(sili.macro.area.period$count), log10(period.pres.all.2$Freq), method = 'spearman')
period.sili.macro.area.3 <- cor.test(log10(sili.macro.area.period$count), log10(period.pres.all.3$Freq), method = 'spearman')
period.sili.macro.area.4 <- cor.test(log10(sili.macro.area.period$count), log10(period.pres.all.4$Freq), method = 'spearman')
period.sili.macro.area.5 <- cor.test(log10(sili.macro.area.period$count), log10(period.pres.all.5$Freq), method = 'spearman')

period.macro.sili.area  <- rbind(c("Period", 1, "Macrostrat siliciclastic (area)", period.sili.macro.area.1$estimate, period.sili.macro.area.1$p.value),
                                 c("Period", 2, "Macrostrat siliciclastic (area)", period.sili.macro.area.2$estimate, period.sili.macro.area.2$p.value),
                                 c("Period", 3, "Macrostrat siliciclastic (area)", period.sili.macro.area.3$estimate, period.sili.macro.area.3$p.value),
                                 c("Period", 4, "Macrostrat siliciclastic (area)", period.sili.macro.area.4$estimate, period.sili.macro.area.4$p.value),
                                 c("Period", 5, "Macrostrat siliciclastic (area)", period.sili.macro.area.5$estimate, period.sili.macro.area.5$p.value))

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
# 8. MAPS
################################################################################

##################
##### MODERN #####
##################

# Set extent
e <<- extent(-180, 180, -90, 90)

##### ALL SPECIMENS #####
get_grid_im(m.dat, 2,  "Echinoids", ext = e)

##### TAPH GRADES #####
get_grid_im(dplyr::filter(m.dat, Preservation_score == 1), 2,  "Taphonomic Grade 1 Echinoids", ext = e)
get_grid_im(dplyr::filter(m.dat, Preservation_score == 2), 2,  "Taphonomic Grade 2 Echinoids", ext = e)
get_grid_im(dplyr::filter(m.dat, Preservation_score == 3), 2,  "Taphonomic Grade 3 Echinoids", ext = e)
get_grid_im(dplyr::filter(m.dat, Preservation_score == 4), 2,  "Taphonomic Grade 4 Echinoids", ext = e)
get_grid_im(dplyr::filter(m.dat, Preservation_score == 5), 2,  "Taphonomic Grade 5 Echinoids", ext = e)
