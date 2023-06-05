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

################################################################################
# 2. BAR PLOTS AND CHI SQUARED
################################################################################

#####################
##### FUNCTIONS #####
#####################

# Bar plot
make.bar.plot <- function(dataset, fill, legend, colour, flip = FALSE){
  a <- ggplot(dataset) +
    aes(x = Preservation_score, fill = fill) +
    geom_bar() +
    scale_fill_manual(values=wes_palette(colour)) +
    guides(fill=guide_legend(title=legend)) +
    ylab("Frequency") +
    xlab("Preservation Score") +
    theme_bw()
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
    scale_y_continuous(expand = c(0,0), limits = c(0,1))
  if(flip == TRUE){
    a <- a + coord_flip()
  }
  return(a)
}

#####################
##### LITHOLOGY #####
#####################

# Remove data without lithological info
l.m.dat <- m.dat %>%
  filter(is.na(Finalised_lith) == F)

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

# Make dataset for grain size
g.m.dat <- m.dat %>%
  filter(is.na(Finalised_grainsize) == F) 

# Assign specific grain size categories into either "Fine Grained" or "Coarse Grained"
for (l in 1:nrow(g.m.dat)){
  if (g.m.dat$Finalised_grainsize[l]=="Fine Grained"  | g.m.dat$Finalised_grainsize[l] == "Slate" |
      g.m.dat$Finalised_grainsize[l]=="Chert" | g.m.dat$Finalised_grainsize[l] =="Mudstone/Grainstone" |
      g.m.dat$Finalised_grainsize[l]=="Shale" | g.m.dat$Finalised_grainsize[l] == "Wackestone" | g.m.dat$Finalised_grainsize[l] == "Mudstone" | 
      g.m.dat$Finalised_grainsize[l]=="Siltstone" | g.m.dat$Finalised_grainsize[l]== "Mudstone/Wackestone" |
      g.m.dat$Finalised_grainsize[l]=="Micrite" | g.m.dat$Finalised_grainsize[l]== "Mudstone/Siltstone"){
    g.m.dat$Finalised_grainsize[l]<-"Fine Grained"
  }
  else if (g.m.dat$Finalised_grainsize[l]=="Grainstone" | g.m.dat$Finalised_grainsize[l] == "Packstone"| 
           g.m.dat$Finalised_grainsize[l]=="Sandstone" | g.m.dat$Finalised_grainsize[l]=="Reefal" |
           g.m.dat$Finalised_grainsize[l] == "Coquina" |
           g.m.dat$Finalised_grainsize[l] == "Boundstone"){ 
    g.m.dat$Finalised_grainsize[l]<-"Coarse Grained"
  }
  else{
    g.m.dat$Finalised_grainsize[l] <- NA
  }
}

# Remove remaining specimens without grain size
g.m.dat <- g.m.dat %>%
  filter(is.na(Finalised_grainsize) == F) 

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

# Make plot of preservation scores against family, NA values removed
f.m.dat <- m.dat %>%
  filter(is.na(Family) == F) %>%
  filter(Family != 'Triadotiaridae') %>%
  filter(Family != 'Cravenechinidae') %>%
  filter(Family != 'Archaeocidaridae or miocidaridae')

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

#################################
##### SET UP - PERIOD LEVEL #####
#################################

# Load periods and rename columns for binning
data(periods)
colnames(periods)[1] <- "bin"
colnames(periods)[2] <- "max_ma"
colnames(periods)[3] <- "min_ma"

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

#################################
##### SET UP - SERIES LEVEL #####
#################################

# Load series data
series <- read.csv("series.csv")

# Bin into series
m.dat.series <- bin_time(m.dat, bins = series, method = 'majority')
m.dat.series$bin_assignment <- as.factor(m.dat.series$bin_assignment)
m.dat.series$bin_assignment <- factor(m.dat.series$bin_assignment, levels = order_ind)

# Assign colours
myColours <- series$color
names(myColours) <- levels(m.dat.series$bin_assignment)
custom_colours <- scale_colour_manual(name = "bin_assignment", values = myColours[1:6])

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

# Transform to Preservation_score by genus matrix
jac.1 <- dcast(jac.1,
              Preservation_score ~ Genus,
              value.var = "Count", fill = 0)
jac.1 <- as.data.frame(jac.1)
jac.2 <- jac.1[,-1]
rownames(jac.2) <- jac.1[,1]

# Make binary (presence/absence)
binary <- jac.2
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

# Final data wrangling
temp.res <- as.data.frame(temp.res)
temp.res <- cbind(temp.res, ncol(binary))
colnames(temp.res) <- c("Test", "Score", "Genera")
temp.res$Test <- factor(temp.res$Test, levels = c("1v2", "2v3", 
                                                  "3v4", "4v5", 
                                                  "1v3", "2v4", 
                                                  "3v5", "1v4", 
                                                  "2v5", "1v5"))
temp.res <- with(temp.res, temp.res[order(Test, Score, Genera),])

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

# Transform to Preservation_score by genus matrix
jac.1 <- dcast(jac.1,
               Preservation_score ~ Family,
               value.var = "Count", fill = 0)
jac.1 <- as.data.frame(jac.1)
jac.2 <- jac.1[,-1]
rownames(jac.2) <- jac.1[,1]

# Make binary (presence/absence)
binary <- jac.2
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

# Final data wrangling
temp.res <- as.data.frame(temp.res)
temp.res <- cbind(temp.res, ncol(binary))
colnames(temp.res) <- c("Test", "Score", "Family")
temp.res$Test <- factor(temp.res$Test, levels = c("1v2", "2v3", 
                                                  "3v4", "4v5", 
                                                  "1v3", "2v4", 
                                                  "3v5", "1v4", 
                                                  "2v5", "1v5"))
temp.res <- with(temp.res, temp.res[order(Test, Score, Family),])

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

for (l in 1:nrow(test)){
  if (test$Finalised_grainsize[l]=="Fine Grained"  | test$Finalised_grainsize[l] == "Slate" |
      test$Finalised_grainsize[l]=="Chert" | test$Finalised_grainsize[l] =="Mudstone/Grainstone" |
      test$Finalised_grainsize[l]=="Shale" | test$Finalised_grainsize[l] == "Wackestone" | test$Finalised_grainsize[l] == "Mudstone" | 
      test$Finalised_grainsize[l]=="Siltstone" | test$Finalised_grainsize[l]== "Mudstone/Wackestone" |
      test$Finalised_grainsize[l]=="Micrite" | test$Finalised_grainsize[l]== "Mudstone/Siltstone"){
    test$Finalised_grainsize[l]<-"Fine Grained"
  }
  else if (test$Finalised_grainsize[l]=="Grainstone" | test$Finalised_grainsize[l] == "Packstone"| 
           test$Finalised_grainsize[l]=="Sandstone" | test$Finalised_grainsize[l]=="Reefal" |
           test$Finalised_grainsize[l] == "Coquina" |
           test$Finalised_grainsize[l] == "Boundstone"){ 
    test$Finalised_grainsize[l]<-"Coarse Grained"
  }
  else{
    test$Finalised_grainsize[l] <- NA
  }
}


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

for (l in 1:nrow(test)){
  if (test$Finalised_grainsize[l]=="Fine Grained"  | test$Finalised_grainsize[l] == "Slate" |
      test$Finalised_grainsize[l]=="Chert" | test$Finalised_grainsize[l] =="Mudstone/Grainstone" |
      test$Finalised_grainsize[l]=="Shale" | test$Finalised_grainsize[l] == "Wackestone" | test$Finalised_grainsize[l] == "Mudstone" | 
      test$Finalised_grainsize[l]=="Siltstone" | test$Finalised_grainsize[l]== "Mudstone/Wackestone" |
      test$Finalised_grainsize[l]=="Micrite" | test$Finalised_grainsize[l]== "Mudstone/Siltstone"){
    test$Finalised_grainsize[l]<-"Fine Grained"
  }
  else if (test$Finalised_grainsize[l]=="Grainstone" | test$Finalised_grainsize[l] == "Packstone"| 
           test$Finalised_grainsize[l]=="Sandstone" | test$Finalised_grainsize[l]=="Reefal" |
           test$Finalised_grainsize[l] == "Coquina" |
           test$Finalised_grainsize[l] == "Boundstone"){ 
    test$Finalised_grainsize[l]<-"Coarse Grained"
  }
  else{
    test$Finalised_grainsize[l] <- NA
  }
}


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
# 5. ORDINAL LEAST REGRESSION
################################################################################

# Bin into periods
g.m.period.dat <- bin_time(g.m.dat, bins = periods, method = 'majority')

# Create factors
order_ind <- rev(c("Permian", "Carboniferous", "Devonian", "Silurian", "Ordovician"))
g.m.period.dat$bin_assignment <- as.factor(g.m.period.dat$bin_assignment)
g.m.period.dat$bin_assignment <- factor(g.m.period.dat$bin_assignment, levels = order_ind)

g.m.period.dat$Preservation_score <- factor(g.m.period.dat$Preservation_score, 
                                     levels = c("1", "2", "3", "4", "5"), 
                                     ordered = TRUE)
# Remove families with low numbers (Cravenechinidae - 2 specimens) and NAs
data <- g.m.period.dat %>%
  filter(Family != 'Cravenechinidae') %>%
  filter(!is.na(Family))

# Set up modelling
samplesize <- 0.6*nrow(data)
index <- sample(seq_len(nrow(data)), size = samplesize)

#Creating training and test set 
datatrain <- data[index,]
table(datatrain$Finalised_lith, datatrain$Finalised_grainsize)
table(datatrain$Finalised_lith, datatrain$Preservation_score)
table(datatrain$Finalised_grainsize, datatrain$Preservation_score, datatrain$Finalised_lith)

datatest <- data[-index,]

m1 <- polr(Preservation_score ~ Finalised_lith + Finalised_grainsize + Family + 
            bin_assignment + Finalised_lith*Finalised_grainsize, 
          data = datatrain, Hess=TRUE)
summary(m1)
m2 <- polr(Preservation_score ~ Finalised_lith + Finalised_grainsize + Family + 
            bin_assignment, 
          data = datatrain, Hess=TRUE)
summary(m2)
m3 <- polr(Preservation_score ~ Finalised_lith + Finalised_grainsize + Family, 
          data = datatrain, Hess=TRUE)
summary(m3)
m4 <- polr(Preservation_score ~ Finalised_lith + Finalised_grainsize + bin_assignment, 
          data = datatrain, Hess=TRUE)
summary(m4)
m5 <- polr(Preservation_score ~ Finalised_lith + bin_assignment + Family, 
          data = datatrain, Hess=TRUE)
summary(m5)
m6 <- polr(Preservation_score ~ Finalised_grainsize + Family + bin_assignment, 
          data = datatrain, Hess=TRUE)
summary(m6)
m7 <- polr(Preservation_score ~ Finalised_grainsize + bin_assignment, 
          data = datatrain, Hess=TRUE)
summary(m7)
m8 <- polr(Preservation_score ~ Finalised_grainsize + Family, 
          data = datatrain, Hess=TRUE)
summary(m8)
m9 <- polr(Preservation_score ~ Finalised_lith + bin_assignment, 
          data = datatrain, Hess=TRUE)
summary(m9)
m10 <- polr(Preservation_score ~ Finalised_lith + Family, 
          data = datatrain, Hess=TRUE)
summary(m10)
m11 <- polr(Preservation_score ~ Family + bin_assignment, 
          data = datatrain, Hess=TRUE)
summary(m11)
m12 <- polr(Preservation_score ~ bin_assignment, 
          data = datatrain, Hess=TRUE)
summary(m12)
m13 <- polr(Preservation_score ~ Family, 
          data = datatrain, Hess=TRUE)
summary(m13)
m14 <- polr(Preservation_score ~ Finalised_lith + Finalised_grainsize + Finalised_lith * Finalised_grainsize, 
          data = datatrain, Hess=TRUE)
summary(m14)
m15 <- polr(Preservation_score ~ Finalised_lith + Finalised_grainsize,
          data = datatrain, Hess=TRUE)
summary(m15)
m16 <- polr(Preservation_score ~ Finalised_lith,
          data = datatrain, Hess=TRUE)
summary(m16)
m17 <- polr(Preservation_score ~ Finalised_grainsize,
          data = datatrain, Hess=TRUE)
summary(m17)
m18 <- polr(Preservation_score ~ 1,
          data = datatrain, Hess=TRUE)
summary(m18)

cand.mod <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, 
                 m15, m16, m17, m18)

Modnames <- c()
for(i in 1:length(cand.mod)){
  temp.name <- as.character(cand.mod[[i]]$terms[[3]])
  if(length(temp.name) > 1){
    temp.name <- temp.name[-1]
    temp.name <- paste(temp.name[1], temp.name[2], sep = " + ")
  }
  Modnames <- c(Modnames, temp.name)
}
mod.tab <- AICcmodavg::aictab(cand.set = cand.mod, modnames = Modnames, sort = T)

# Best model
m <- polr(formula = Preservation_score ~ Finalised_lith + Finalised_grainsize + 
            Family + bin_assignment, data = datatrain, Hess = TRUE)

## view a summary of the model
summary(m)
(ctable <- coef(summary(m)))

## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

## combined table
(ctable <- cbind(ctable, "p value" = p))
ci <- confint.default(m)

## odds ratios
exp(coef(m))

## OR and CI
test <-exp(cbind(OR = coef(m), ci))

# Model quality
predictPreservation_score <- predict(m,datatest)
table(datatest$Preservation_score, predictPreservation_score)
predictPreservation_score <- as.character(predictPreservation_score)
predictPreservation_score[is.na(predictPreservation_score)] <- "0"
mean(as.character(datatest$Preservation_score) != as.character(predictPreservation_score))

# Plotting 
plot(Effect(focal.predictors = c("Finalised_grainsize", "Finalised_lith"),m),
     main="Ordinal logistic regression; top model, showing lithology + grain size",
     axes = list(y=list(lab="Preservation Score (probability)")),
     xlab = "")
plot(Effect(focal.predictors = c("bin_assignment", "Finalised_lith"),m),
     main="Ordinal logistic regression; top model, showing age + lithology",
     axes = list(y=list(lab="Preservation Score (probability)")), 
     xlab = "")
plot(Effect(focal.predictors = c("bin_assignment", "Finalised_grainsize"),m),
     main="Ordinal logistic regression; top model, showing age + grain size",
     axes = list(y=list(lab="Preservation Score (probability)")), 
     xlab = "")

################################################################################
# 6. CORRELATION TESTS
################################################################################

#################
##### SETUP #####
#################

# Make columns for period, then bin time to stage level (this is so that both are in one dataset)
m.dat.period$Period <- as.character(m.dat.period$bin_assignment)
m.dat.period$Period.no <- as.numeric(m.dat.period$bin_assignment)
m.dat.period <- bin_time(m.dat.period, bins = stages, method = "majority")

# Split out species and genus level information
species.lvl <- m.dat.period %>%
  dplyr::filter(Rank == "Species") %>%
  dplyr::select(Genus, Species, Locality, bin_midpoint, Period, Period.no, Preservation_score) %>%
  mutate(Combined_name = paste(Genus, Species, sep = " ")) %>%
  distinct()
genus.lvl <- m.dat.period %>%
  dplyr::filter(Rank == "Genus") %>%
  dplyr::select(Genus, Locality, bin_midpoint, Period, Period.no, Preservation_score) %>%
  distinct()

# Make tables of number of taxa per preservation score for each time frame
stage.pres.spec <- as.data.frame(table(species.lvl$bin_midpoint, species.lvl$Preservation_score))
stage.pres.gen <- as.data.frame(table(genus.lvl$bin_midpoint, genus.lvl$Preservation_score))
period.pres.spec <- as.data.frame(table(species.lvl$Period, species.lvl$Preservation_score))
period.pres.gen <- as.data.frame(table(genus.lvl$Period, genus.lvl$Preservation_score))

# Run dplyr to get diversity/collections at genus and species level for each time frame
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
# Create order
order_ind <- c("Permian", "Carboniferous", "Devonian", "Silurian", "Ordovician")

# Function to split out and reorganise data by preservation score for correlations
split <- function(data, Var2, order = FALSE){
  for(t in Var2){
    temp.data <- filter(data, Var2 == t)
    if(order == TRUE){
      temp.data <- temp.data[match(order_ind, temp.data$Var1),]
    }
    assign(paste(deparse(substitute(data)), t, sep = "."), temp.data, envir = .GlobalEnv)
  }
}

# Run function for each level
split(period.pres.gen, period.pres.gen$Var2, order = TRUE)
split(period.pres.spec, period.pres.spec$Var2, order = TRUE)
split(stage.pres.spec, stage.pres.spec$Var2, order = FALSE)
split(stage.pres.gen, stage.pres.gen$Var2, order = FALSE)

################################################
##### SPECIES LEVEL DIVERSITY CORRELATIONS #####
################################################

# Stage - diversity vs.preservation score
stage.spec.div.1 <- cor.test(log10(div.stage.spec$SIBs), log10(stage.pres.spec.1$Freq), method = 'spearman')
stage.spec.div.2 <- cor.test(log10(div.stage.spec$SIBs), log10(stage.pres.spec.2$Freq), method = 'spearman')
stage.spec.div.3 <- cor.test(log10(div.stage.spec$SIBs), log10(stage.pres.spec.3$Freq), method = 'spearman')
stage.spec.div.4 <- cor.test(log10(div.stage.spec$SIBs), log10(stage.pres.spec.4$Freq), method = 'spearman')
stage.spec.div.5 <- cor.test(log10(div.stage.spec$SIBs), log10(stage.pres.spec.5$Freq), method = 'spearman')
# Stage - collections vs. preservation score
stage.spec.col.1 <- cor.test(log10(div.stage.spec$colls), log10(stage.pres.spec.1$Freq), method = 'spearman')
stage.spec.col.2 <- cor.test(log10(div.stage.spec$colls), log10(stage.pres.spec.2$Freq), method = 'spearman')
stage.spec.col.3 <- cor.test(log10(div.stage.spec$colls), log10(stage.pres.spec.3$Freq), method = 'spearman')
stage.spec.col.4 <- cor.test(log10(div.stage.spec$colls), log10(stage.pres.spec.4$Freq), method = 'spearman')
stage.spec.col.5 <- cor.test(log10(div.stage.spec$colls), log10(stage.pres.spec.5$Freq), method = 'spearman')

stage.div.spec <- rbind(c("stage", 1, "diversity (species)", stage.spec.div.1$estimate, stage.spec.div.1$p.value),
                        c("stage", 2, "diversity (species)", stage.spec.div.2$estimate, stage.spec.div.2$p.value),
                        c("stage", 3, "diversity (species)", stage.spec.div.3$estimate, stage.spec.div.3$p.value),
                        c("stage", 4, "diversity (species)", stage.spec.div.4$estimate, stage.spec.div.4$p.value),
                        c("stage", 5, "diversity (species)", stage.spec.div.5$estimate, stage.spec.div.5$p.value))
stage.coll.spec <- rbind(c("stage", 1, "collections (species)", stage.spec.col.1$estimate, stage.spec.col.1$p.value),
                         c("stage", 2, "collections (species)", stage.spec.col.2$estimate, stage.spec.col.2$p.value),
                         c("stage", 3, "collections (species)", stage.spec.col.3$estimate, stage.spec.col.3$p.value),
                         c("stage", 4, "collections (species)", stage.spec.col.4$estimate, stage.spec.col.4$p.value),
                         c("stage", 5, "collections (species)", stage.spec.col.5$estimate, stage.spec.col.5$p.value))

# Period - diversity vs. preservation score
period.spec.div.1 <- cor.test(div.period.spec$SIBs, period.pres.spec.1$Freq, method = 'spearman')
period.spec.div.2 <- cor.test(div.period.spec$SIBs, period.pres.spec.2$Freq, method = 'spearman')
period.spec.div.3 <- cor.test(div.period.spec$SIBs, period.pres.spec.3$Freq, method = 'spearman')
period.spec.div.4 <- cor.test(div.period.spec$SIBs, period.pres.spec.4$Freq, method = 'spearman')
period.spec.div.5 <- cor.test(div.period.spec$SIBs, period.pres.spec.5$Freq, method = 'spearman')
# Period - collections vs. preservation score
period.spec.col.1 <- cor.test(div.period.spec$colls, period.pres.spec.1$Freq, method = 'spearman')
period.spec.col.2 <- cor.test(div.period.spec$colls, period.pres.spec.2$Freq, method = 'spearman')
period.spec.col.3 <- cor.test(div.period.spec$colls, period.pres.spec.3$Freq, method = 'spearman')
period.spec.col.4 <- cor.test(div.period.spec$colls, period.pres.spec.4$Freq, method = 'spearman')
period.spec.col.5 <- cor.test(div.period.spec$colls, period.pres.spec.5$Freq, method = 'spearman')

period.div.spec <- rbind(c("period", 1, "diversity (species)", period.spec.div.1$estimate, period.spec.div.1$p.value),
                         c("period", 2, "diversity (species)", period.spec.div.2$estimate, period.spec.div.2$p.value),
                         c("period", 3, "diversity (species)", period.spec.div.3$estimate, period.spec.div.3$p.value),
                         c("period", 4, "diversity (species)", period.spec.div.4$estimate, period.spec.div.4$p.value),
                         c("period", 5, "diversity (species)", period.spec.div.5$estimate, period.spec.div.5$p.value))
period.coll.spec <- rbind(c("period", 1, "collections (species)", period.spec.col.1$estimate, period.spec.col.1$p.value),
                          c("period", 2, "collections (species)", period.spec.col.2$estimate, period.spec.col.2$p.value),
                          c("period", 3, "collections (species)", period.spec.col.3$estimate, period.spec.col.3$p.value),
                          c("period", 4, "collections (species)", period.spec.col.4$estimate, period.spec.col.4$p.value),
                          c("period", 5, "collections (species)", period.spec.col.5$estimate, period.spec.col.5$p.value))


##############################################
##### GENUS LEVEL DIVERSITY CORRELATIONS #####
##############################################

# Stage - diversity vs.preservation score
stage.gen.div.1 <- cor.test(log10(div.stage.gen$SIBs), log10(stage.pres.gen.1$Freq), method = 'spearman')
stage.gen.div.2 <- cor.test(log10(div.stage.gen$SIBs), log10(stage.pres.gen.2$Freq), method = 'spearman')
stage.gen.div.3 <- cor.test(log10(div.stage.gen$SIBs), log10(stage.pres.gen.3$Freq), method = 'spearman')
stage.gen.div.4 <- cor.test(log10(div.stage.gen$SIBs), log10(stage.pres.gen.4$Freq), method = 'spearman')
stage.gen.div.5 <- cor.test(log10(div.stage.gen$SIBs), log10(stage.pres.gen.5$Freq), method = 'spearman')
# Stage - collections vs. preservation score
stage.gen.col.1 <- cor.test(log10(div.stage.gen$colls), log10(stage.pres.gen.1$Freq), method = 'spearman')
stage.gen.col.2 <- cor.test(log10(div.stage.gen$colls), log10(stage.pres.gen.2$Freq), method = 'spearman')
stage.gen.col.3 <- cor.test(log10(div.stage.gen$colls), log10(stage.pres.gen.3$Freq), method = 'spearman')
stage.gen.col.4 <- cor.test(log10(div.stage.gen$colls), log10(stage.pres.gen.4$Freq), method = 'spearman')
stage.gen.col.5 <- cor.test(log10(div.stage.gen$colls), log10(stage.pres.gen.5$Freq), method = 'spearman')

stage.div.gen <-  rbind(c("stage", 1, "diversity (genera)",stage.gen.div.1$estimate, stage.gen.div.1$p.value),
                        c("stage", 2, "diversity (genera)",stage.gen.div.2$estimate, stage.gen.div.2$p.value),
                        c("stage", 3, "diversity (genera)",stage.gen.div.3$estimate, stage.gen.div.3$p.value),
                        c("stage", 4, "diversity (genera)",stage.gen.div.4$estimate, stage.gen.div.4$p.value),
                        c("stage", 5, "diversity (genera)",stage.gen.div.5$estimate, stage.gen.div.5$p.value))
stage.coll.gen <-  rbind(c("stage", 1, "collections (genera)",stage.gen.col.1$estimate, stage.gen.col.1$p.value),
                         c("stage", 2, "collections (genera)",stage.gen.col.2$estimate, stage.gen.col.2$p.value),
                         c("stage", 3, "collections (genera)",stage.gen.col.3$estimate, stage.gen.col.3$p.value),
                         c("stage", 4, "collections (genera)",stage.gen.col.4$estimate, stage.gen.col.4$p.value),
                         c("stage", 5, "collections (genera)",stage.gen.col.5$estimate, stage.gen.col.5$p.value))

# Period - diversity vs. preservation score
period.gen.div.1 <- cor.test(log10(div.period.gen$SIBs), log10(period.pres.gen.1$Freq), method = 'spearman')
period.gen.div.2 <- cor.test(log10(div.period.gen$SIBs), log10(period.pres.gen.2$Freq), method = 'spearman')
period.gen.div.3 <- cor.test(log10(div.period.gen$SIBs), log10(period.pres.gen.3$Freq), method = 'spearman')
period.gen.div.4 <- cor.test(log10(div.period.gen$SIBs), log10(period.pres.gen.4$Freq), method = 'spearman')
period.gen.div.5 <- cor.test(log10(div.period.gen$SIBs), log10(period.pres.gen.5$Freq), method = 'spearman')
# Period - collections vs. preservation score
period.gen.col.1 <- cor.test(log10(div.period.gen$colls), log10(period.pres.gen.1$Freq), method = 'spearman')
period.gen.col.2 <- cor.test(log10(div.period.gen$colls), log10(period.pres.gen.2$Freq), method = 'spearman')
period.gen.col.3 <- cor.test(log10(div.period.gen$colls), log10(period.pres.gen.3$Freq), method = 'spearman')
period.gen.col.4 <- cor.test(log10(div.period.gen$colls), log10(period.pres.gen.4$Freq), method = 'spearman')
period.gen.col.5 <- cor.test(log10(div.period.gen$colls), log10(period.pres.gen.5$Freq), method = 'spearman')

period.div.gen <-  rbind(c("period", 1, "diversity (genera)", period.gen.div.1$estimate, period.gen.div.1$p.value),
                         c("period", 2, "diversity (genera)", period.gen.div.2$estimate, period.gen.div.2$p.value),
                         c("period", 3, "diversity (genera)", period.gen.div.3$estimate, period.gen.div.3$p.value),
                         c("period", 4, "diversity (genera)", period.gen.div.4$estimate, period.gen.div.4$p.value),
                         c("period", 5, "diversity (genera)", period.gen.div.5$estimate, period.gen.div.5$p.value))
period.coll.gen <- rbind(c("period", 1,"collections (genera)", period.gen.col.1$estimate, period.gen.col.1$p.value),
                          c("period",2, "collections (genera)",period.gen.col.2$estimate, period.gen.col.2$p.value),
                          c("period",3, "collections (genera)",period.gen.col.3$estimate, period.gen.col.3$p.value),
                          c("period",4, "collections (genera)",period.gen.col.4$estimate, period.gen.col.4$p.value),
                          c("period",5, "collections (genera)",period.gen.col.5$estimate, period.gen.col.5$p.value))
cor.results <- rbind(as.data.frame(stage.div.spec), as.data.frame(stage.coll.spec), 
      as.data.frame(period.div.spec), as.data.frame(period.coll.spec), 
      as.data.frame(stage.div.gen), as.data.frame(stage.coll.gen), 
      as.data.frame(period.div.gen), as.data.frame(period.coll.gen)
      )
colnames(cor.results) <- c("Bin size", "Preservation Score", "vs.", "Rho", "p")
cor.results$Rho <- signif(as.numeric(cor.results$Rho), digits = 3)
cor.results$p <- signif(as.numeric(cor.results$p), digits = 8)

#################################
##### SED TYPE CORRELATIONS #####
#################################

##### SETUP: PBDB #####

# Load all palaeozoic occurrences
pb.all <- read.csv("all_palaeozoic_binned.csv")

# Get distinct collections through time
All.colls <- pb.all %>%
  dplyr::select(collection_no, max_ma, min_ma, lng, lat, cc, formation, lithology1, environment, bin_assignment, bin_midpoint) %>%
  distinct()

# Make list of carbonate/siliclastic terms
carb_list <- c('"carbonate', '"limestone"', '"reef rocks"', 'bafflestone', 
               'bindstone', 'dolomite', 'floatstone', 'framestone', 'evaporite',
               'grainstone', 'lime mudstone', 'marl', 'packstone', 'rudstone',
               'wackestone')
sili_list <- c('"shale"', '"siliciclastic"', '"volcaniclastic"', 'ash', 'breccia', 
               'chert', 'claystone', 'coal', 'conglomerate', 'gravel', 'ironstone', 
               'lignite', 'mudstone', 'phosphorite', 'phyllite', 'quartzite', 
               'radiolariate', 'sandstone', 'schist', 'siderite', 'siltstone', 
               'slate', 'tuff')

# Subset collections
carb.colls <- filter(All.colls, lithology1 %in% carb_list)
sili.colls <-filter(All.colls, lithology1 %in% sili_list)

carb.colls <- carb.colls %>%
  dplyr::group_by(bin_midpoint) %>%
  dplyr::summarise(count = n(), lith = "carb") 

sili.colls <- sili.colls %>%
  dplyr::group_by(bin_midpoint) %>%
  dplyr::summarise(count = n(), lith = "sili")

All.colls <- rbind(carb.colls, sili.colls)

a <- ggplot(All.colls, aes(x=bin_midpoint, y=count)) + 
  geom_line(aes(colour = lith)) +
  scale_x_reverse() +
  theme_bw() +
  ylab("Jaccard Similarity") +
  xlab("Time (Ma)") +
  geom_point(aes(color = lith))
gggeo_scale(a)

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
gggeo_scale(a)

##### CORRELATIONS #####

# Taphonomic grade vs. carbonate (PBDB)
stage.spec.div.1 <- cor.test(log10(), log10(stage.pres.spec.1$Freq), method = 'spearman')
# Taphonomic grade vs. siliclastic (PBDB)

# Taphonomic grade vs. carbonate (Macrostrat)

# Taphonomic grade vs. siliclastic (Macrostrat)



################################################################################
# 7. MAPS
################################################################################

####################
##### FUNCTION #####
####################

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

##################
##### MODERN #####
##################

# Set extent
e <<- extent(-180, 180, -90, 90)

##### ALL SPECIMENS #####
get_grid_im(m.dat, 2,  "Museum Echinoids", ext = e)

##### TAPH GRADES #####
get_grid_im(dplyr::filter(m.dat, Preservation_score == 1), 2,  "Museum Echinoids", ext = e)
get_grid_im(dplyr::filter(m.dat, Preservation_score == 2), 2,  "Museum Echinoids", ext = e)
get_grid_im(dplyr::filter(m.dat, Preservation_score == 3), 2,  "Museum Echinoids", ext = e)
get_grid_im(dplyr::filter(m.dat, Preservation_score == 4), 2,  "Museum Echinoids", ext = e)
get_grid_im(dplyr::filter(m.dat, Preservation_score == 5), 2,  "Museum Echinoids", ext = e)

##################
##### PALAEO #####
##################
