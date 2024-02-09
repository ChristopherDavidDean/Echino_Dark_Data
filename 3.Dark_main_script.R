################################################################################
### MUSEUM `DARK DATA` ILLUMINATES THE EVOLUTIONARY HISTORY OF FOSSIL GROUPS ###
################################################################################

# Christopher D. Dean & Jeffrey R. Thompson
# 2024
# Script written by Christopher D. Dean

################################################################################
#                            FILE 3: MAIN SCRIPT                               #
################################################################################

################################################################################
# 1. SETUP
################################################################################

# Packages
library(dplyr)
library(iNEXT)
library(ggpubr)
library(rstatix)
library(DescTools)
library(kableExtra)
library(lmtest)
library(nlme)
library(MuMIn)
library(vcd)
library(gtable)
library(cowplot)

############################################
###### LOAD AND PREPARE MUSEUM DATASET #####
############################################

# Load required functions
source("1.Dark_Functions.R")

# Load setup file
source("2.Dark_Setup.R")

################################################################################
# 2. BASIC COMPARISONS
################################################################################

##############################
##### PRESERVATION SCORE #####
##############################

# Format into table and run Chi Squared test
chisq.test(table(m.dat$Preservation_score, m.dat$Museum_only))

# compare against what the test would have expected
chisq.test(table(m.dat$Preservation_score, m.dat$Museum_only))$expected

# Mosaic plot
vcd::mosaic(~ Preservation_score + Museum_only,
       direction = c("v", "h"),
       data = m.dat,
       labeling_args = list(just_labels = "right", 
                            set_varnames = c(Preservation_score = "Preservation Score", 
                                             Museum_only = "Only found in museum")),
       shade = TRUE
)

########################################
##### TABLE COMPARISONS AND GRAPHS #####
########################################

##### FUNCTIONS #####
# Make table of tallied results to plot, taking specific column.
setup.table <- function(mdata, 
                        pbdata, 
                        pubdata, 
                        alldata, 
                        column, 
                        pivot = T, 
                        useNA = F){
  if(useNA == T){
    temp <- list(as.data.frame(table(mdata[[column]], useNA = "ifany")),
                 as.data.frame(table(pbdata[[column]], useNA = "ifany")),
                 as.data.frame(table(pubdata[[column]], useNA = "ifany")),
                 as.data.frame(table(alldata[[column]], useNA = "ifany")))
  }else{
    temp <- list(as.data.frame(table(mdata[[column]])),
                 as.data.frame(table(pbdata[[column]])),
                 as.data.frame(table(pubdata[[column]])),
                 as.data.frame(table(alldata[[column]])))
  }
  temp <- temp %>%
    purrr::reduce(full_join, by = "Var1") %>%
    dplyr::rename(temp.name = Var1,
                  Museum = Freq.x,
                  PBDB = Freq.y,
                  Published = Freq.x.x,
                  All = Freq.y.y)
  temp[is.na(temp)] <- 0
  names(temp)[names(temp) == "temp.name"] <- column
  if(pivot == T){
    temp <- tidyr::pivot_longer(temp, cols = c("Museum", "PBDB", 
                                               "Published", "All"))
    temp$value <- as.numeric(temp$value)
  }
  assign(deparse(substitute(column)), temp, envir = .GlobalEnv)
}

# Plots results from setup.table, as both count and proportion
table.plots <- function(pivot, x, removeNA = T, colour){
  if(removeNA == T){
    pivot <- na.omit(pivot)
  }
  a <- ggplot(data = pivot, aes(x = name, y = value, fill = !!sym(x))) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(values=colour) +
    guides(fill=guide_legend(title="Rank")) +
    ylab("Total occurrences") +
    xlab("Dataset") +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "white"), 
          panel.grid.minor = element_line(colour = "white")) +
    guides(alpha = "none")
  b <- ggplot(data = pivot, aes(x = name, y = value, fill = !!sym(x))) +
    geom_bar(stat = 'identity', position = 'fill') +
    scale_fill_manual(values=colour) +
    guides(fill=guide_legend(title="Rank")) +
    ylab("Proportion of dataset") +
    xlab("Dataset") +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "white"), 
          panel.grid.minor = element_line(colour = "white")) +
    guides(alpha = "none")
  
  ggarrange(a, b, 
            align='hv',
            labels = c("A", "B"),
            nrow = 1, 
            ncol = 2,
            common.legend = T, 
            legend = "bottom")
}

##### FORMATIONS #####
length(unique(pb.only.dat$Formation))
length(unique(m.only.dat$Formation))
length(unique(pub.all.dat$Formation))
length(unique(all.dat$Formation))

##### TAXONOMIC RANK #####
Rank <- setup.table(m.only.dat, 
                    pb.only.dat, 
                    pub.all.dat, 
                    all.dat,
                    "Rank", 
                    pivot = T, 
                    useNA = F)
(p1 <- table.plots(Rank, "Rank", colour = c(wes_palette("Zissou1"), 
                                           wes_palette("Royal1"),
                                           wes_palette("Royal2"))))

ggsave("Dark_data_graphs/1.Rank.props.png", plot = p1, 
       device = "png", type = "cairo")

##### COUNTRIES #####

country.clean <- function(df){
  ##### Setup data for analysis #####
  df[df == "England" | df == "Wales" | df == "Scotland" | 
    df == "Northern Ireland"] <- "United Kingdom"
  df[df == "USA"] <- "United States"
  df[df == "Czechia"] <- "Czech Republic"
  df$Country <- gsub("&", "and", df$Country)
  
  df$aff_code <- countrycode(df$Museum_country, "country.name", "iso3c")
  df$samp_code <- countrycode(df$Country, "country.name", "iso3c")
  return(df)
}

m.only.dat <- country.clean(m.only.dat)
pb.only.dat <- country.clean(pb.only.dat)
pub.all.dat <- country.clean(pub.all.dat)
all.dat <- country.clean(all.dat)

cc <- setup.table(m.only.dat, 
                    pb.only.dat, 
                    pub.all.dat, 
                    all.dat,
                    "Country", 
                    pivot = F, 
                    useNA = F)
colSums(cc != 0)

##### CONTINENTS #####
Cont <- setup.table(m.only.dat, 
                    pb.only.dat, 
                    pub.all.dat, 
                    all.dat,
                    "Continent", 
                    pivot = T, 
                    useNA = F)
(p2 <- table.plots(Cont, "Continent", colour = c(wes_palette("Royal2"), 
                                                 wes_palette("Royal1"))))

ggsave("Dark_data_graphs/2.Cont.props.png", plot = p2, 
       device = "png", type = "cairo")

##### LITHOLOGY #####
lith <- setup.table(m.only.dat, 
                    pb.only.dat, 
                    pub.all.dat, 
                    all.dat,
                    "Finalised_lith", 
                    pivot = F, 
                    useNA = T)
# Cleaning data
lith <- rbind(lith[1,], lith[2,], lith[3,] + lith[4,], lith[5,])
lith <- tidyr::pivot_longer(lith, cols = c("Museum", "PBDB", "Published", "All"))
lith$value <- as.numeric(lith$value)
names(lith)[names(lith) == "Finalised_lith"] <- "Lithology"

# Plotting data
(p3 <- table.plots(lith, "Lithology", removeNA = T, colour = c(wes_palette("Zissou1"), 
                                                               wes_palette("Royal1"),
                                                               wes_palette("Royal2"))))

ggsave("Dark_data_graphs/3.Lith.props.png", plot = p3, 
       device = "png", type = "cairo")

##### GRAINSIZE #####
Grain <- setup.table(m.only.dat, 
                    pb.only.dat, 
                    pub.all.dat, 
                    all.dat,
                    "Finalised_grainsize_simp", 
                    pivot = T, 
                    useNA = T)
names(Grain)[names(Grain) == "Finalised_grainsize_simp"] <- "Grainsize"

(p4 <- table.plots(Grain, "Grainsize", removeNA = T, 
                   colour = c(wes_palette("Darjeeling2"))))

ggsave("Dark_data_graphs/4.Grain.props.png", plot = p4, 
       device = "png", type = "cairo")

################################################################################
# 4. SPATIAL OCCUPANCY THROUGH TIME
################################################################################

############################
##### COMPARATIVE MAPS #####
############################

# Make maps
a <- get_grid_im(m.only.dat, 2,  "Museum only echinoid occurrences", ext = e)
b <- get_grid_im(pb.only.dat, 2,  "PBDB echinoid occurrences", ext = e)
c <- get_grid_im(pub.all.dat, 2,  "Published echinoid occurrences", ext = e)
d <- get_grid_im(all.dat, 2,  "All echinoid occurrences data", ext = e)

p5 <- ggarrange(a, b, c, d, 
          align='hv',
          labels = c("A", "B", "C", "D"),
          nrow = 2,
          ncol = 2)

ggsave("Dark_data_graphs/5.Occ.maps.png", plot = p5, 
       device = "png", type = "cairo")

############################
##### SETUP/BASIC INFO #####
############################

# Setup data for get_grid()
m.only.dat.space <- m.only.dat %>%
  filter(!is.na(lat) & !is.na(lng)) 

# Apply get_grid()
get_grid(m.only.dat.space, res, e)
get_grid(pb.only.dat, res, e)
get_grid(pub.all.dat, res, e)
get_grid(all.dat, res, e)
pb.all <- pb.all %>%
  dplyr::rename(Locality = collection_no)
get_grid(pb.all, res, e)
pb.echino <- pb.echino %>%
  dplyr::rename(Locality = collection_no)
get_grid(pb.echino, res, e)

# Find how many unique cells each dataset occupies
totalcells <- length(unique(pb.all.spat.bin$cells)) + 17
totalecells <- length(unique(pb.echino.spat.bin$cells)) + 61
length(unique(pb.only.dat.spat.bin$cells))
length(unique(m.only.dat.space.spat.bin$cells))
length(unique(pub.all.dat.spat.bin$cells))
length(unique(all.dat.spat.bin$cells)) 

# Function to find percentage of cells covered by each dataset, in relation to PBDB collections
find.perc <- function(cell.data){
    return(paste("Total cells = ", length(unique(cell.data)), 
                 "; percentage of echino total = ", 
          ((length(unique(cell.data))/totalecells)*100), 
          "; percentage of all PBDB total = ", 
          ((length(unique(cell.data))/totalcells)*100), sep = ""))
}

length(setdiff(m.only.dat.space.spat.bin$cells, pub.all.dat.spat.bin$cells))

find.perc(pb.only.dat.spat.bin$cells)
find.perc(m.only.dat.space.spat.bin$cells)
find.perc(pub.all.dat.spat.bin$cells)
find.perc(all.dat.spat.bin$cells)

##################################
##### OCCUPANCY THROUGH TIME #####
##################################

# Loop through stages, providing occupancy estimates for each stage
results <- as.data.frame(matrix(NA, nrow = 0, ncol = 11))

for(i in unique(pb.all$bin_assignment)){
  # Get number of cells for each time bin
  m.cells <- unique(m.only.dat.space.spat.bin$cells[which(m.only.dat.space.spat.bin$bin_assignment == i)])
  pb.cells <- unique(pb.only.dat.spat.bin$cells[which(pb.only.dat.spat.bin$bin_assignment == i)])
  pub.cells <- unique(pub.all.dat.spat.bin$cells[which(pub.all.dat.spat.bin$bin_assignment == i)])
  tot.cells <- unique(all.dat.spat.bin$cells[which(all.dat.spat.bin$bin_assignment == i)])
  echino.cells <- unique(pb.echino.spat.bin$cells[which(pb.echino.spat.bin$bin_assignment == i)])
  all.cells <- unique(pb.all.spat.bin$cells[which(pb.all.spat.bin$bin_assignment == i)])
  unique.m.cells <- length(setdiff(m.only.dat.space.spat.bin$cells[which(m.only.dat.space.spat.bin$bin_assignment == i)],
                                   pub.all.dat.spat.bin$cells[which(pub.all.dat.spat.bin$bin_assignment == i)]))
  shared.cells <- length(intersect(m.only.dat.space.spat.bin$cells[which(m.only.dat.space.spat.bin$bin_assignment == i)],
                                   pub.all.dat.spat.bin$cells[which(pub.all.dat.spat.bin$bin_assignment == i)]))
  unique.pub.cells <- length(setdiff(pub.all.dat.spat.bin$cells[which(pub.all.dat.spat.bin$bin_assignment == i)], 
                                     m.only.dat.space.spat.bin$cells[which(m.only.dat.space.spat.bin$bin_assignment == i)]))
  
  # If there are museum specimens covering cells that aren't covered by all cells
  if(length(setdiff(m.cells, all.cells)) > 0){
    add.cells <- setdiff(m.cells, all.cells)
    all.cells <- c(all.cells, add.cells)
  }
  if(length(setdiff(m.cells, echino.cells)) > 0){
    add.cells <- setdiff(m.cells, echino.cells)
    echino.cells <- c(echino.cells, add.cells)
  }
  # Save results
  results <- rbind(results, 
                       c(i, 
                         stages$bin_midpoint[stages$bin == i], 
                         length(all.cells), 
                         length(echino.cells),
                         "m.dat",
                         length(m.cells), 
                         (length(m.cells)/length(all.cells))*100, 
                         (length(m.cells)/length(echino.cells))*100, 
                         unique.m.cells, 
                         shared.cells, 
                         unique.pub.cells),
                        c(i,
                         stages$bin_midpoint[stages$bin == i], 
                         length(all.cells), 
                         length(echino.cells),
                         "pb.dat",
                         length(pb.cells), 
                         (length(pb.cells)/length(all.cells))*100,
                         (length(pb.cells)/length(echino.cells))*100, 
                         unique.m.cells, 
                         shared.cells, 
                         unique.pub.cells),
                        c(i,
                         stages$bin_midpoint[stages$bin == i], 
                         length(all.cells), 
                         length(echino.cells),
                         "pub.dat",
                         length(pub.cells), 
                         (length(pub.cells)/length(all.cells))*100,
                         (length(pub.cells)/length(echino.cells))*100, 
                         unique.m.cells, 
                         shared.cells, 
                         unique.pub.cells),
                        c(i,
                          stages$bin_midpoint[stages$bin == i], 
                          length(all.cells), 
                          length(echino.cells),
                          "all.dat",
                          length(tot.cells), 
                          (length(tot.cells)/length(all.cells))*100, 
                          (length(tot.cells)/length(echino.cells))*100, 
                          unique.m.cells, 
                          shared.cells, 
                          unique.pub.cells)
                       )
}

# Cleaning and setup for plots
colnames(results) <- c("bin", "bin_midpoint", "all.cells", 
                       "echino.cells", "data", "n.cells", "occ.all", "occ.echino", 
                       "unique.m", "shared", "unique.p")

colnames(results)[5] <-"Data"
results$Data[which(results$Data == "all.dat")] <- "All data"
results$Data[which(results$Data == "m.dat")] <- "Museum"
results$Data[which(results$Data == "pb.dat")] <- "PBDB"
results$Data[which(results$Data == "pub.dat")] <- "Published"
results$Data <- as.factor(results$Data)
results$bin_midpoint <- as.numeric(results$bin_midpoint)
results$n.cells <- as.numeric(results$n.cells)
results$echino.cells <- as.numeric(results$echino.cells)
results$occ.all <- as.numeric(results$occ.all)
results$occ.echino <- as.numeric(results$occ.echino)
results$all.cells <- as.numeric(results$all.cells)
results$bin <- as.numeric(results$bin)
results$unique.m <- as.numeric(results$unique.m)
results$unique.p <- as.numeric(results$unique.p)
results$shared <- as.numeric(results$shared)

results <- results %>%
  filter(bin > 55 & bin < 93)

results %>%
  dplyr::group_by(Data) %>%
  dplyr::summarize(Mean.occ = mean(n.cells), 
                   Median.occ = median(n.cells), 
                   Mean.perc = mean(occ.echino, na.rm = T), 
                   Median.perc = median(occ.echino, na.rm = T), 
                   Mean.perc.all = mean(occ.all, na.rm = T),
                   Median.perc.all = median(occ.all, na.rm = T))

stacked.results <- results %>%
  dplyr::select(bin, bin_midpoint, unique.m, shared, unique.p) %>%
  distinct() %>%
  mutate(ratio = unique.m/unique.p) %>%
  tidyr::pivot_longer(cols = c(unique.m, shared, unique.p, ratio))

stacked.results %>%
  dplyr::group_by(name) %>%
  dplyr::summarize(Mean.m = mean(value, na.rm = T),
                   Median.m = median(value, na.rm = T))

stacked.results[stacked.results == "shared"] <- "Shared"
stacked.results[stacked.results == "unique.m"] <- "Museum only"
stacked.results[stacked.results == "unique.p"] <- "Published record only"
stacked.results[stacked.results == "ratio"] <- "Ratio"

#################
##### PLOTS #####
#################

a <- stacked.results %>%
  dplyr::filter(name == "Ratio") %>%
  ggplot(aes(x=bin_midpoint, y=value)) + 
  geom_line() +
  scale_x_reverse() +
  theme_bw() +
  ylab("Ratio of unique \nmuseum/published gridcells") +
  guides(fill=guide_legend(title="Preservation score")) +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white"), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  guides(alpha = "none")

b <- stacked.results %>%
  dplyr::filter(name != "Ratio") %>%
  ggplot(aes(x=bin_midpoint, y=value, fill=name)) + 
  geom_area() +
  scale_x_reverse() +
  theme_bw() +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Occupied Grid Cells")) +
  xlab("Time (Ma)") +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  scale_fill_manual(values=(rev(wes_palette("Chevalier1")[1:3]))) +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  guides(alpha = "none")

(p6.5 <- ggarrange(a, b, 
                 align='v',
                 labels = c("A", "B"),
                 nrow = 2,
                 ncol = 1, 
                 heights = c(0.55, 0.7),
                 legend = "bottom"))

ggsave("Dark_data_graphs/6.5.Occupancy.time.new.png", plot = p6.5, 
       device = "png", type = "cairo")

b <- ggplot(results, aes(x=bin_midpoint, y=occ.echino)) + 
  geom_line(aes(colour = Data)) +
  scale_x_reverse() +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  guides(alpha = "none") +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  ylim(0, 100) +
  ylab("Occupancy (%)") +
  xlab("Time (Ma)")  +
  scale_colour_viridis_d()

a <- ggplot(results, aes(x=bin_midpoint, y=n.cells)) + 
  geom_line(aes(colour = Data)) +
  scale_x_reverse() +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  guides(alpha = "none") +
  ylim(0, 100) +
  xlim(485.4, 251.902) +
  ylab("Number of occupied cells") +
  xlab("Time (Ma)")  +
  theme(legend.position = "none", 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_colour_viridis_d()

(p6 <- ggarrange(a, b, 
          align='v',
          labels = c("A", "B"),
          nrow = 2,
          ncol = 1,
          common.legend = T,
          legend = "bottom", 
          heights = c(0.5, 0.6)))

ggsave("Dark_data_graphs/6.Occupancy.time.png", plot = p6, 
       device = "png", type = "cairo")

#########################
###### CORRELATIONS #####
#########################

# Reshape data for tests
wide.results <- dcast(results, bin+bin_midpoint+all.cells+echino.cells~Data, 
                      value.var = c("n.cells"))

# Log data
wide.results$all.cells <- log10(wide.results$all.cells)
wide.results$Museum <- log10(wide.results$Museum)
wide.results$PBDB <- log10(wide.results$PBDB)
wide.results$Published <- log10(wide.results$Published)
wide.results$echino.cells<- log10(wide.results$echino.cells)
wide.results$`All data` <- log10(wide.results$`All data`)

# Occupancy vs. occupancy
cor1 <- cor.test(wide.results$Museum, wide.results$PBDB, method = "spearman")
cor2 <- cor.test(wide.results$Museum, wide.results$Published, method = "spearman")
cor3 <- cor.test(wide.results$Museum, wide.results$`All data`, method = "spearman")
cor4 <- cor.test(wide.results$PBDB, wide.results$Published, method = "spearman")
cor5 <- cor.test(wide.results$PBDB, wide.results$`All data`, method = "spearman")
cor6 <- cor.test(wide.results$Published, wide.results$`All data`, method = "spearman")

occ.v.occ <- rbind(c("Museum", "PBDB", cor1$estimate, cor1$p.value),
                   c("Museum", "Published", cor2$estimate, cor2$p.value),
                   c("Museum", "All", cor3$estimate, cor3$p.value),
                   c("PBDB", "Published", cor4$estimate, cor4$p.value),
                   c("PBDB", "All", cor5$estimate, cor5$p.value), 
                   c("Published", "All", cor6$estimate, cor6$p.value))

# Occupancy vs. Total invert occ
cor7 <- cor.test(wide.results$Museum, wide.results$all.cells, method = "spearman")
cor8 <- cor.test(wide.results$PBDB, wide.results$all.cells, method = "spearman")
cor9 <- cor.test(wide.results$Published, wide.results$all.cells, method = "spearman")
cor10 <- cor.test(wide.results$`All data`, wide.results$all.cells, method = "spearman")

occ.v.inv.occ <- rbind(c("Museum", "All PBDB", cor7$estimate, cor7$p.value),
                   c("PBDB", "All PBDB", cor8$estimate, cor8$p.value),
                   c("Published", "All PBDB", cor9$estimate, cor9$p.value),
                   c("All", "All PBDB", cor10$estimate, cor10$p.value))

# Occupancy vs. Total echinodermata occ
cor11 <- cor.test(wide.results$Museum, wide.results$echino.cells, method = "spearman")
cor12 <- cor.test(wide.results$PBDB, wide.results$echino.cells, method = "spearman")
cor13 <- cor.test(wide.results$Published, wide.results$echino.cells, method = "spearman")
cor14 <- cor.test(wide.results$`All data`, wide.results$echino.cells, method = "spearman")

occ.v.e.occ <- rbind(c("Museum", "Echinodermata PBDB", cor11$estimate, cor11$p.value),
                       c("PBDB", "Echinodermata PBDB", cor12$estimate, cor12$p.value),
                       c("Published", "Echinodermata PBDB", cor13$estimate, cor13$p.value),
                       c("All", "Echinodermata PBDB", cor14$estimate, cor14$p.value))

# Combine results
cor.results <- as.data.frame(rbind(occ.v.occ, 
                                   occ.v.inv.occ, 
                                   occ.v.e.occ))
colnames(cor.results) <- c("Variable 1", "vs. Variable 2", "Rho", "p")
cor.results$Rho <- signif(as.numeric(cor.results$Rho), digits = 3)
cor.results$p <- signif(as.numeric(cor.results$p), digits = 5)

# Correct for multiple tests
cor.results$BH <- p.adjust(cor.results$p, method = "BH")
cor.results$Signif <- ifelse(cor.results$BH < 0.05, "*", "")

cor.results %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  save_kable(file = "Dark_data_graphs/T1.occ.corrs.html")

################################################################################
# 5. SPATIAL RANGES THROUGH TIME
################################################################################

#################
##### SETUP #####
#################

# Function for finding ranges of organisms
make.range <- function(occdf, name, method, rank, coords = T){
  if(rank == "Genus"){
    temp.data <- occdf %>%
      dplyr::filter(Rank == "Genus" | Rank == "Species") %>%
      dplyr::filter(is.na(p_lng) == F)
    temp.data <- group_apply(temp.data, group = c("bin_assignment"), 
                             fun = tax_range_space, 
                             name = "Genus", 
                             lng = "p_lng", 
                             lat = "p_lat", 
                             method = method, 
                             coords = coords)
  }
  if(rank == "Family"){
    temp.data <- occdf %>%
      dplyr::filter(Rank == "Family" | Rank == "Genus" | Rank == "Species") %>%
      dplyr::filter(is.na(Family) == F) %>%
      dplyr::filter(is.na(p_lng) == F)
    temp.data <- group_apply(temp.data, group = c("bin_assignment"), 
                             fun = tax_range_space, 
                             name = "Family", 
                             lng = "p_lng", 
                             lat = "p_lat", 
                             method = method, 
                             coords = coords)
  }
  if(method == "con"){
    temp.data <- temp.data[temp.data$area != 0, ]
  }
  if(method == "lat"){
    temp.data <- temp.data[temp.data$range_lat != 0, ]
  }
  temp.data$Data <- name
  assign(paste(deparse(substitute(occdf)), ".range.", method, sep = ""), temp.data, 
         envir = .GlobalEnv)
}

# Function for finding latitudinal ranges
lat.range.fun <- function(pbddata, 
                          musdata, 
                          pubdata, 
                          alldata, 
                          name, rank = "Genus"){
  # Combine ranges from all datasets
  lat.ranges <- rbind(make.range(pbddata, "PBDB", "lat", rank), 
                      make.range(musdata, "Museum", "lat", rank), 
                      make.range(pubdata, "Published", "lat", rank), 
                      make.range(alldata, "All", "lat", rank))
 # lat.ranges$bin_assignment <- as.numeric(lat.ranges$bin_assignment)
  
  assign(paste("lat.ranges.", name, sep = ""), lat.ranges, envir = .GlobalEnv)
  
  lat.results <- data.frame()
  for(i in sort(unique(lat.ranges$bin_assignment))){
    # Get data for bin, and remove any taxa that can't be compared
    temp.data <- lat.ranges %>%
      dplyr::filter(bin_assignment == i) %>%
      dplyr::group_by(taxon) %>%
      dplyr::filter(n() >2)
    # If no taxa remain, skip to next bin
    if(nrow(temp.data) < 1){
      next
    }
    # Get difference between Total range and specific range for each taxon
    temp.data <- temp.data %>%
      left_join(dplyr::filter(., Data == "All"), by = "taxon") %>%
      mutate(Difference = range_lat.y- range_lat.x, 
             Perc = (range_lat.x/range_lat.y)*100, 
             Mean_shift = ((max_lat.y + min_lat.y)/2) - ((max_lat.x + min_lat.x)/2)) %>%
      dplyr::select(taxon, Data.x, Difference, Perc, Mean_shift) %>%
      dplyr::filter(Data.x != "All") %>%
      dplyr::rename(Taxon = taxon, Data = Data.x) %>%
      as.data.frame() %>%
      arrange(Taxon) %>%
      mutate(bin = i)
    # Save results
    lat.results <- rbind(lat.results, temp.data)
  }
  return(lat.results)
}

# Function for quickly summarising differences between datasets
space.results <- function(results){
  results %>%
    group_by(Data) %>%
    dplyr::summarize(Mean.diff = mean(Difference), 
                     Median.diff = median(Difference), 
                     Mean.perc = mean(Perc), 
                     Median.perc = median(Perc))
}

org.funct <- function(data, var){
  test.2 <- data %>%
    dplyr::group_by(!!sym(var), taxon) %>%
    dplyr::summarize(count= n()) %>%
    st_cast("POLYGON") %>%
    st_convex_hull() %>%
    na.omit(geometry)
  if(any(st_is_valid(test.2) == F)){
    test.2 <- test.2[-which(st_is_valid(test.2) == F),]
  }
  test.2 <- test.2 %>%
    dplyr::mutate(Ex_Area = st_area(geometry)) %>%
    as.data.frame() %>%
    dplyr::select(taxon, Ex_Area)
  return(test.2)
}

# Function that generates results of a variety of statistical tests. 
spatial.tests <- function(geo.results, focal){
  kruskal.test <- geo.results %>%
    kruskal_test(as.formula(paste(focal, " ~ Data", sep = "")))
  kruskal.effect <- geo.results %>%
    kruskal_effsize(as.formula(paste(focal, " ~ Data", sep = "")))
  dunn.test <- geo.results %>%
    dunn_test(as.formula(paste(focal, " ~ Data", sep = "")), 
              p.adjust.method = "bonferroni")
  wilcox.test <- geo.results %>%
    wilcox_test(as.formula(paste(focal, " ~ Data", sep = "")),
                p.adjust.method = "bonferroni")
  
  return(list(kruskal.test = kruskal.test, 
              kruskal.effect = kruskal.effect,
              dunn.test = dunn.test, 
              wilcox.test = wilcox.test))
}

# Calculates results for geographic area for taxa through time
geo.range.fun <- function(pbddata, 
                          musdata, 
                          pubdata, 
                          alldata, 
                          name, 
                          rank = "Genus"){
  # Combine ranges from all datasets
  geo.ranges <- rbind(make.range(pbddata, "PBDB", "con", rank, coords =  T), 
                      make.range(musdata, "Museum", "con", rank, coords =  T), 
                      make.range(pubdata, "Published", "con", rank, coords =  T),
                      make.range(alldata, "All", "con", rank, coords =  T))
  
  assign(paste("geo.ranges.", name, sep = ""), geo.ranges, envir = .GlobalEnv)
  
  geo.results <- data.frame()
  for(i in sort(unique(geo.ranges$bin_assignment))){
    # Get data for bin, and remove any taxa that can't be compared
    temp.data <- geo.ranges %>%
      dplyr::filter(bin_assignment == i) %>%
      group_by(taxon) %>%
      dplyr::filter(n_distinct(Data) > 2)
    
    # If no taxa remain, skip to next bin
    if(nrow(temp.data) < 1){
      next
    }
    temp_xy <- st_as_sf(temp.data, coords = c("p_lng", "p_lat"), crs = 4326)
    temp_xy <- st_transform(temp_xy, crs = c("+proj=eqearth"))
    polys <- temp_xy %>%
      dplyr::group_by(Data, taxon) %>%
      dplyr::summarize(count= n()) %>%
      st_cast("POLYGON") %>%
      st_convex_hull()
    
    # Setup biggest contribution
    test <- temp_xy %>%
      mutate(PBDB = case_when(Data == "Museum" ~ "PBDB",
                              Data == "Published" ~ "PBDB"), 
             Museum = case_when(Data == "PBDB" ~ "Museum",
                               Data == "Published" ~ "Museum"),
             Published = case_when(Data == "Museum" ~ "Published",
                                   Data == "PBDB" ~ "Published")
             )
  
    mus.test <- org.funct(test, "Museum")
    names(mus.test)[names(mus.test) == "Ex_Area"] <- "Museum_Ex_Area"
    PBD.test <- org.funct(test, "PBDB")
    names(PBD.test)[names(PBD.test) == "Ex_Area"] <- "PBDB_Ex_Area"
    pub.test <- org.funct(test, "Published")
    names(pub.test)[names(pub.test) == "Ex_Area"] <- "Published_Ex_Area"
    
    comb.test <- full_join(mus.test, pub.test, by = "taxon")
    comb.test <- full_join(comb.test, PBD.test, by = "taxon")
    
    print(comb.test)
    
    if(any(st_is_valid(polys) == F)){
      polys <- polys[-which(st_is_valid(polys) == F),]
      polys <- polys %>%
        group_by(taxon) %>%
        dplyr::filter(n_distinct(Data) > 2)
    }
    
    # Getting full area
    combined_polys <- polys %>%
      dplyr::group_by(taxon) %>%
      dplyr::summarise() %>%
      st_convex_hull() %>%
      dplyr::mutate(area = st_area(geometry))
      
    # Combining
    polys$area <- st_area(polys)
    temp.data <- polys %>%
      as.data.frame() %>%
      dplyr::full_join(combined_polys, by = "taxon") %>%
      mutate(Area = area.x,
             Difference = area.y - area.x, 
             Perc = (area.x/area.y)*100) %>%
      dplyr::select(taxon, Area, Data, Difference, Perc, count) %>%
      dplyr::filter(Data != "All") %>%
      dplyr::rename(Taxon = taxon) %>%
      as.data.frame() %>%
      arrange(Taxon) %>%
      mutate(bin = i)
    
    # Added area
    added.area <- combined_polys %>%
      as.data.frame() %>%
      dplyr::select(-geometry) %>%
      dplyr::full_join(comb.test, by = "taxon") %>%
      dplyr::mutate(Museum = area - Museum_Ex_Area,
             Published = area - Published_Ex_Area, 
             PBDB = area - PBDB_Ex_Area) %>%
      dplyr::select(taxon, Museum, PBDB, Published, area) %>%
      tidyr::pivot_longer(cols = c(Museum, PBDB, Published)) %>%
      dplyr::mutate(Percentage.added = (value/area)*100) %>%
      dplyr::rename(Taxon = taxon, 
                    Data = name,
                    Added.area = value)
    print(added.area)
    
    temp.data$Area <- as.numeric(temp.data$Area)
    temp.data$Difference <- as.numeric(temp.data$Difference)
    temp.data$Perc <- as.numeric(temp.data$Perc)
    temp.data <- left_join(temp.data, added.area, by = c("Taxon", "Data"))
    temp.data$Added.area <- as.numeric(temp.data$Added.area)
    temp.data$Added.area[temp.data$Added.area < 0] <- 0
    geo.results <- rbind(geo.results, temp.data)
    print(i)
  }
  return(geo.results)
}

##############################
##### LATITUDINAL RANGES #####
##############################

##### SETUP #####
# Stage
lat.results.stage <- lat.range.fun(pbddata = pb.only.dat,
              musdata = m.only.dat, 
              pubdata = pub.all.dat, 
              alldata = all.dat, 
              "stage")
# Series
lat.results.series <- lat.range.fun(pbddata = pb.only.series,
              musdata = m.only.series, 
              pubdata = pub.all.series, 
              alldata = all.series, 
              "series")

##### MAXIMUM AND MINIMUM PALAEOLATITUDES #####
max_min_lat <- pivot_longer(lat.ranges.stage, cols = c("max_lat", "min_lat"))
names(max_min_lat)[names(max_min_lat) == "name"] <- "Palaeolatitude"
max_min_lat$Palaeolatitude[which(max_min_lat$Palaeolatitude == "max_lat")] <- "Maximum"
max_min_lat$Palaeolatitude[which(max_min_lat$Palaeolatitude == "min_lat")] <- "Minimum"

(p7 <- max_min_lat %>%
  ggplot(aes(x = Data, y = value, fill = Palaeolatitude, color = Palaeolatitude)) +
  geom_boxplot(aes(alpha = 0.1), show.legend = NA) +
  geom_jitter(size=0.4, alpha=0.5) +
  scale_color_manual(name = "Palaeolatitude", values = c("#440154FF","#21908CFF")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_fill_manual(name = "Palaeolatitude", values = c("#440154FF","#21908CFF")) +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  guides(alpha = "none") +
  xlab("") +
  ylab("Palaeolatitude"))

ggsave("Dark_data_graphs/7.Lat.max.min.png", plot = p7, 
       device = "png", type = "cairo")

lat.ranges.stage %>%
  pivot_longer(-c(taxon_id, range_lat, taxon, bin_assignment, Data), 
               names_to = "lat_group") %>%
  group_by(Data, lat_group) %>%
  dplyr::summarize("avalue" = mean(value),
                   "sd" = sd(value, na.rm = F)) %>%
  dplyr::rename(value = avalue)

##### PALAEOLATITUDINAL SKEW #####

a <- lat.results.stage %>%
  ggplot(aes(x = Data, y = Mean_shift, fill = Data)) +
  guides(alpha = "none") +
  geom_boxplot(aes(alpha = 0.1), show.legend = NA) +
  scale_fill_viridis(discrete = T) + 
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  guides(alpha = "none") +
  theme(legend.position = "none") +
  ylab("Shift in range mean") 
  
b <- lat.results.stage %>%
  mutate(text = forcats::fct_reorder(Data, Mean_shift)) %>%
  ggplot(aes(x=Mean_shift, color = Data, fill=Data)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', bins = 15) +
  scale_fill_viridis(discrete = T) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  guides(alpha = "none") +
  xlab("Shift in range mean") + 
  ylab("Total") +
  facet_wrap(~ Data)

(p8 <- ggarrange(a, b, 
          align='hv',
          labels = c("A", "B"),
          nrow = 1,
          ncol = 2, 
          common.legend = T, 
          legend = "bottom"))

ggsave("Dark_data_graphs/8.Lat.shift.png", plot = p8, 
       device = "png", type = "cairo")

lat.results.stage %>%
  group_by(Data) %>%
  dplyr::summarize(Mean.shift.mean = mean(Mean_shift), 
                   Median.shift.mean = median(Mean_shift))

##### DIFFERENCES BETWEEN DATASETS #####
# Find results
space.results(lat.results.series)
space.results(lat.results.stage)

# Mean/median addition of latitudinal range from museums, as percentage of total
lat.results.stage %>%
  dplyr::filter(Data == "Museum" | Data == "Published") %>%
  dplyr::group_by(bin, Taxon) %>%
  dplyr::filter(n() >1) %>%
  dplyr::filter(Data == "Published") %>%
  mutate(Perc.new = 100 - Perc) %>%
  pull(Perc.new) %>%
  mean() # change this to get median

# Mean/median addition of latitudinal range from museums in degrees
lat.results.stage %>%
  dplyr::filter(Data == "Museum" | Data == "Published") %>%
  dplyr::group_by(bin, Taxon) %>%
  dplyr::filter(n() >1) %>%
  dplyr::filter(Data == "Published") %>%
  pull(Difference) %>%
  mean()

# Plot results
a <- lat.results.stage %>%
  ggplot(aes(x = Data, y = Difference, fill = Data)) +
  geom_boxplot() +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  guides(alpha = "none") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("Difference (°)") +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Difference between full latitudinal range and data range, stage level") +
  xlab("")
b <- lat.results.series %>%
  ggplot(aes(x = Data, y = Difference, fill = Data)) +
  geom_boxplot() +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  guides(alpha = "none") +
  ylab("Difference (°)") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Difference between full latitudinal range and data range, series level") +
  xlab("")

(p9 <- ggarrange(a, b, 
          align='hv',
          labels = c("A", "B"),
          nrow = 2,
          ncol = 1))

ggsave("Dark_data_graphs/9.Lat.diff.png", plot = p9, 
       device = "png", type = "cairo")

# Kruskal wallace test for differences between medians
spatial.tests(lat.results.stage, "Difference")
spatial.tests(lat.results.stage, "Mean_shift")
spatial.tests(lat.results.stage, "Perc")
spatial.tests(lat.ranges.stage, "range_lat")

#############################
##### GEOGRAPHIC RANGES #####
#############################

# Generate results
geo.results.stage <- geo.range.fun(pbddata = pb.only.dat,
              musdata = m.only.dat, 
              pubdata = pub.all.dat, 
              alldata = all.dat, 
              "stage")

geo.results.series <- geo.range.fun(pbddata = pb.only.series,
                                   musdata = m.only.series, 
                                   pubdata = pub.all.series, 
                                   alldata = all.series, 
                                   "series")

# Generate results
geo.results.stage.fam <- geo.range.fun(pbddata = pb.only.dat,
                                   musdata = m.only.dat,
                                   pubdata = pub.all.dat,
                                   alldata = all.dat,
                                   "stage", rank = "Family")

space.results(geo.results.stage.fam)
geo.results.stage.fam %>%
  dplyr::group_by(Data) %>%
  dplyr::summarize(Mean = mean(Percentage.added),
                   Median = median(Percentage.added))

# Find differences
space.results(geo.results.stage)
space.results(geo.results.series)

# Plot differences
a <- geo.results.stage %>%
  ggplot(aes(x = Data, y = Perc, fill = Data)) +
  geom_boxplot() +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  guides(alpha = "none") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ylab("Percentage of range overlap") +
  xlab("")

geo.results.stage$Percentage.added <- as.numeric(geo.results.stage$Percentage.added)

# Plot differences
b <- geo.results.stage %>%
    ggplot(aes(x = Data, y = Percentage.added, fill = Data)) +
    geom_boxplot() +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "white"), 
          panel.grid.minor = element_line(colour = "white")) +
    guides(alpha = "none") +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ylab(expression('Mean area added to geographic range 
                    (% of total area from convex hull)')) +
    xlab("")

(p10 <- ggarrange(a, b, 
                  align='hv',
                  labels = c("A", "B"),
                  nrow = 1))

ggsave("Dark_data_graphs/10.Range.overlap.added.area.png", plot = p10, 
       device = "png")

# Mean added area (percentage of total area)
geo.results.stage %>%
  dplyr::group_by(Data) %>%
  dplyr::summarize(Mean = mean(Percentage.added),
            Median = median(Percentage.added))

geo.results.series %>%
  dplyr::group_by(Data) %>%
  dplyr::summarize(Mean = mean(Percentage.added),
                   Median = median(Percentage.added))

# Plot of added area per data type
geo.results.stage$Percentage.added <- as.numeric(geo.results.stage$Percentage.added)
(p10.5 <- geo.results.stage %>%
  ggplot(aes(Percentage.added, fill = Data)) + 
  geom_histogram() +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  guides(alpha = "none") +
  xlab("Percentage of area added") +
  ylab("Frequency") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6))

ggsave("Dark_data_graphs/10.5.Added.area.png", plot = p10.5, 
       device = "png")

##### STATISTICAL TESTS #####
spatial.tests(geo.results.stage, "Area")
spatial.tests(geo.results.stage, "Difference")
spatial.tests(geo.results.stage, "Perc")
spatial.tests(geo.results.stage, "Added.area")

all.geo <- geo.ranges.stage %>%
  dplyr::select(-c(p_lng, p_lat)) %>%
  distinct() 

spatial.tests(all.geo, "area")

################################################################################
# 6. TEMPORAL RANGE COMPARISON
################################################################################

#################
##### SETUP #####
#################

# Function
temp_compare <- function(occdf, bin_filter){
  ranks <- c("Species", "Genus", "Family")
  comp <- data.frame()
  for(i in ranks){
    if(i == "Species"){
      test.all.dat <- all.dat %>%
        dplyr::filter(Rank == "Species") %>%
        dplyr::filter(n_bins < bin_filter)
      test.occdf <- occdf %>%
        dplyr::filter(Rank == "Species") %>%
        dplyr::filter(n_bins < bin_filter)
      a <- tax_range_time(test.all.dat, name = "Combined_name", 
                          by = "name", plot = F)
      b <- tax_range_time(test.occdf, name = "Combined_name", 
                          by = "name", plot = F)
    }
    if(i == "Genus"){
      test.all.dat <- all.dat %>%
        dplyr::filter(Rank == "Genus" | Rank == "Species") %>%
        dplyr::filter(is.na(Genus) == F) %>%
        dplyr::filter(n_bins < bin_filter)
      test.occdf <- occdf %>%
        dplyr::filter(Rank == "Genus" | Rank == "Species") %>%
        dplyr::filter(n_bins < bin_filter)
      a <- tax_range_time(test.all.dat, name = "Genus", by = "name", plot = F)
      b <- tax_range_time(test.occdf, name = "Genus", by = "name", plot = F)
    }
    if(i == "Family"){
      test.all.dat <- all.dat %>%
        dplyr::filter(Rank == "Family" | Rank == "Genus" | Rank == "Species") %>%
        dplyr::filter(is.na(Family) == F) %>%
        dplyr::filter(n_bins < bin_filter) %>%
        dplyr::filter(Family != "NO_FAMILY_SPECIFIED")
      test.occdf <- occdf %>%
        dplyr::filter(Rank == "Family" | Rank == "Genus" | Rank == "Species") %>%
        dplyr::filter(is.na(Family) == F) %>%
        dplyr::filter(n_bins < bin_filter) %>%
        dplyr::filter(Family != "NO_FAMILY_SPECIFIED")
      a <- tax_range_time(test.all.dat, name = "Family", by = "name", plot = F)
      b <- tax_range_time(test.occdf, name = "Family", by = "name", plot = F)
    }
    temp_comp <- inner_join(a, b, by = "taxon")
    temp_comp$range_diff <- temp_comp$range_myr.x - temp_comp$range_myr.y
    temp_comp$bottom_diff <- temp_comp$max_ma.x - temp_comp$max_ma.y
    temp_comp$top_diff <- temp_comp$min_ma.x - temp_comp$min_ma.y
    temp_comp$Group <- i
    temp_comp$Range_perc <- (temp_comp$range_myr.y/temp_comp$range_myr.x)*100
    comp <- rbind(comp, temp_comp)
  }
  comp$dataframe <- deparse(substitute(occdf))
  assign(paste("comp.", deparse(substitute(occdf)), sep = ""), comp, 
         envir = .GlobalEnv)
  comp_stats <- comp %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(Mean_range_diff = mean(range_diff), 
                     Mean_top_diff = mean(top_diff),
                     Mean_bot_diff = mean(bottom_diff),
                     Correct_perc = sum(range_diff == 0)/n(), 
                     Range_perc = mean(Range_perc),
                     Count = n()
    )
  comp_stats$dataframe <- deparse(substitute(occdf))
  assign(paste("comp.stats.", deparse(substitute(occdf)), sep = ""), comp_stats, 
         envir = .GlobalEnv)
}

###################
##### RESULTS #####
###################

# Run on each dataset
temp_compare(m.only.dat, bin_filter = 3)
temp_compare(pb.only.dat, bin_filter = 3)
temp_compare(pub.all.dat, bin_filter = 3)

# Combined
comp <- rbind(comp.m.only.dat, comp.pb.only.dat, comp.pub.all.dat)
(comp.stats <- rbind(comp.stats.m.only.dat, 
                     comp.stats.pb.only.dat, 
                     comp.stats.pub.all.dat))

comp %>%
  dplyr::filter(Group == "Family") %>%
  dplyr::filter(dataframe == "m.only.dat" | dataframe == "pub.all.dat") %>%
  dplyr::group_by(taxon) %>%
  dplyr::filter(n() >1) %>%
  dplyr::filter(dataframe == "pub.all.dat") %>%
  mutate(Perc.new = 100 - Range_perc) %>%
  pull(Perc.new) %>%
  mean() # change this to get median


# Tidying
comp$dataframe[comp$dataframe == "m.only.dat"] <- "Museum"
comp$dataframe[comp$dataframe == "pb.only.dat"] <- "PBDB"
comp$dataframe[comp$dataframe == "pub.all.dat"] <- "Published"
names(comp)[names(comp) == "dataframe"] <- "Data"
names(comp)[names(comp) == "range_diff"] <- "Difference"

# Plots
(p11 <- comp %>%
  ggplot(aes(x=Data, y=Difference, fill = Group)) +
  geom_boxplot() +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "white"), 
          panel.grid.minor = element_line(colour = "white")) +
    guides(alpha = "none") +
  scale_fill_viridis(discrete = T, alpha = 0.8) +
  xlab("Dataset") + 
  ylab("Difference in temporal range (Ma)"))

(p11.5 <- comp %>%
    ggplot(aes(x=Data, y=Range_perc, fill = Group)) +
    geom_boxplot() +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "white"), 
          panel.grid.minor = element_line(colour = "white")) +
    guides(alpha = "none") +
    scale_fill_viridis(discrete = T, alpha = 0.8) +
    xlab("Dataset") + 
    ylab("Percentage of temporal range overlap"))

ggsave("Dark_data_graphs/11.Temporal.diff.png", plot = p11, 
       device = "png", type = "cairo")

ggsave("Dark_data_graphs/11.5.Temporal.overlap.png", plot = p11.5, 
       device = "png", type = "cairo")

# Stat tests
comp %>%
  dplyr::filter(Group == "Species") %>%
  spatial.tests(., "Difference")
comp %>%
  dplyr::filter(Group == "Genus") %>%
  spatial.tests(., "Difference")
comp %>%
  dplyr::filter(Group == "Family") %>%
  spatial.tests(., "Difference")

################################################################################
# 7. DIVERSITY ANALYSIS
################################################################################

#####################
##### FUNCTIONS #####
#####################

# Run divDyn, either for raw results for SQS, on a chosen dataset at chosen taxonomic rank.
run.div <- function(dataset, Rank = "Genus", type = "raw"){
  if(Rank == "Genus"){
    temp.data <- dataset %>%
      dplyr::filter(Rank == "Genus" | Rank == "Species")
  }
  if(Rank == "Species"){
    temp.data <- dataset %>%
      dplyr::filter(Rank == "Species")
  }
  if(type == "raw"){
    if(Rank == "Genus"){
      data.div <- divDyn::divDyn(x = temp.data, bin = "bin_assignment", 
                                 tax = "Genus", noNAStart = T)
    }else{
      data.div <- divDyn::divDyn(x = temp.data, bin = "bin_assignment", 
                                 tax = "Combined_name", noNAStart = T)
    }
  }
  if(type == "SQS"){
    quorum_levels <- round(seq(from = 0.3, to = 0.6, by = 0.1), 1)
    data.div <- list()
    for(i in 1:length(quorum_levels)){
      if(Rank == "Genus"){
        data.div.temp <- divDyn::subsample(x = temp.data, 
                                           bin = "bin_assignment", 
                                           tax = "Genus", 
                                           type = "sqs", 
                                           q = quorum_levels[i],
                                           iter = 500, 
                                           noNAStart = T)
        data.div.temp$quorum_level <- quorum_levels[i]
      }else{
        data.div.temp <- divDyn::subsample(x = temp.data, 
                                           bin = "bin_assignment", 
                                           tax = "Combined_name", 
                                           type = "sqs", 
                                           q = quorum_levels[i],
                                           iter = 500, 
                                           noNAStart = T)
        data.div.temp$quorum_level <- quorum_levels[i]
      }
      data.div[[i]] <- data.div.temp
    }
  }
  return(data.div)
}

run.col <- function(dataset, Rank = "Genus"){
  if(Rank == "Genus"){
    temp.data <- dataset %>%
      dplyr::filter(Rank == "Genus" | Rank == "Species")
    data.col <- divDyn::binstat(x = temp.data, bin = "bin_assignment",
                                coll = "Locality", 
                               tax = "Genus", noNAStart = T)
  }
  if(Rank == "Species"){
    temp.data <- dataset %>%
      dplyr::filter(Rank == "Species")
    data.col <- divDyn::binstat(x = temp.data, bin = "bin_assignment", 
                                coll = "Locality", 
                                tax = "Combined_name", noNAStart = T)
  }
  return(data.col)
}

# Combines outputs from divDyn to provide a dataframe of diversity values. 
setup.div <- function(pbddata,
                      musdata,
                      pubdata,
                      alldata){
  # Only select relevant info (SIB diversity)
  m.only.div <- musdata %>%
    dplyr::select(bin_assignment, divSIB)
  pb.only.div <- pbddata %>%
    dplyr::select(bin_assignment, divSIB)
  pub.all.div <- pubdata %>%
    dplyr::select(bin_assignment, divSIB)
  all.div <- alldata %>%
    dplyr::select(bin_assignment, divSIB)
  
  # Set up dataset for merging, rename bin column
  div.tot <- stages
  names(div.tot)[names(div.tot) == "bin"]  <- "bin_assignment"
  
  # Merge data
  div.tot <- left_join(div.tot, m.only.div, by = "bin_assignment")
  names(div.tot)[names(div.tot) == "divSIB"] <- "Museum_SIB"
  div.tot <- left_join(div.tot, pb.only.div, by = "bin_assignment")
  names(div.tot)[names(div.tot) == "divSIB"] <- "PBDB_SIB"
  div.tot <- left_join(div.tot, pub.all.div, by = "bin_assignment")
  names(div.tot)[names(div.tot) == "divSIB"] <- "Pub_SIB"
  div.tot <- left_join(div.tot, all.div, by = "bin_assignment")
  names(div.tot)[names(div.tot) == "divSIB"] <- "All_SIB"
  
  # Remove rows with NAs, then remove unneccessary columns
  div.tot <- div.tot[rowSums(is.na(div.tot[,c("Museum_SIB", "PBDB_SIB", 
                                              "Pub_SIB",  "All_SIB")]))!=4,]
  div.tot <- div.tot %>%
    dplyr::mutate(bin_midpoint = (max_ma + min_ma)/2) %>%
    dplyr::select(-c(max_ma, min_ma, abbr, color, name))
  return(div.tot)
}

# Function for plotting the output of divDyn.
divDyn.plot <- function(estD_output, leg = "none"){
  ## The output is a 'list' object so we'll need to convert it into a dataframe and clean it up before plotting
  estD_plotting <- bind_rows(estD_output) # binds rows of a list
  
  ## Ensure that the quorum level column is being treated as a 'factor' to avoid errors while plotting:
  estD_plotting$quorum_level <- as.factor(estD_plotting$quorum_level)
  
  ## Create a colour gradient for as many colours as you have quorum levels:
  teal_gradient <- scales::seq_gradient_pal("turquoise", "darkslategrey", "Lab")(seq(0, 1, length.out = 4))
  
  stages.red <- stages %>%
    dplyr::select(bin, bin_midpoint) %>%
    dplyr::rename(bin_assignment = bin)
  
  estD_plotting <- left_join(estD_plotting, stages.red, by = "bin_assignment")
  
  divDyn_plot <- ggplot(estD_plotting, aes(x = bin_midpoint, 
                                           y = divSIB, 
                                           colour = quorum_level)) + 
    ## Each quorum level is called individually to be plotted:
    ## Set our line and point sizes (and shapes):
    geom_line(linewidth = 1) +
    scale_shape_manual(values=c(15, 16, 17)) +
    ## Add our colours, theme, and axes labels:
    scale_colour_manual(values = teal_gradient) +
    scale_x_reverse() +
    ylim(0, 6) +
    xlim(477.7000, 251.9020) +
    geom_point(size = 3) +
    labs(x = "Time (Ma)", y = "Coverage rarified richness") +
    theme_bw() +
    labs(color='Quorum level') +
    theme(panel.grid.major = element_line(colour = "white"), 
          panel.grid.minor = element_line(colour = "white")) +
    theme(legend.position = "right") + 
    theme(plot.margin = margin(0.1,0.2,0, 1, "cm"), 
          axis.title.y = element_text(hjust= 0.8, size = 9),
          axis.title.x = element_text(size = 9))
  if(leg == "none"){
    divDyn_plot <- divDyn_plot + theme(legend.position = "none", 
                                       axis.text.x = element_blank(), 
                                       axis.ticks.x = element_blank(), 
                                       axis.title.x = element_blank())
  }
  return(gggeo_scale(divDyn_plot))
}


# Function to generate frequency data for use in iNEXT
freq.data <- function(dataset, Rank = "Genus", bin.dataset, filter_data = T){
  if(Rank == "Genus"){
    temp.data <- dataset %>%
      dplyr::filter(Rank == "Genus" | Rank == "Species")
    names(temp.data)[names(temp.data) == "Genus"] <- "accepted_name"
  }
  if(Rank == "Species"){
    temp.data <- dataset %>%
      dplyr::filter(Rank == "Species")
    names(temp.data)[names(temp.data) == "Combined_name"] <- "accepted_name"
  }
  taxa_data <- subset(temp.data, select=c(accepted_name, Locality, bin_assignment))
  freq_data <- lapply(bin.dataset$bin, function(f) {
    tmp <- taxa_data %>% dplyr::filter(bin_assignment == f) %>% 
      dplyr::count(., accepted_name) %>% 
      dplyr::arrange(desc(n)) %>% 
      tibble::add_row(accepted_name = "n", n = sum(.$n), .before = 1) %>%
      dplyr::select(n)
    freq_raw <- as.numeric(tmp$n)
    freq_raw
  })
  if(deparse(substitute(bin.dataset)) == "series"){
    names(freq_data) <- bin.dataset$bin # give each list element its correct interval name
  }else{
    names(freq_data) <- bin.dataset$name # give each list element its correct interval name
  }
 if(filter_data == T){
   # Set to 2 to avoid issue of singletons
   freq_data <- freq_data[lapply(freq_data,length)>2] 
   # Remove any time periods with only 1 occurrence of each taxa
   #freq_data <- freq_data[names(freq_data) %in% names(Filter(function(vector) all(vector[-1] == 1), freq_data)) == F] 
   # Keep freq_data with above 15 total occurrences (otherwise fails)
   freq_data <- freq_data[lapply(freq_data,function(l)l[[1]]>15) == T]
   freq_data <- freq_data[lapply(freq_data,length)>3] 
 }
  return(freq_data)
}

# Function to run iNEXT, esting estimateD. Uses set quora of 0.3-0.6
run.iNEXT <- function(freq_data, bin.dataset){
  ## And create a new list to store output
  estD_output <- list() 
  ## Create a vector of quorum levels that we want to compute
  ## 0.4 is considered the 'standard', but the fashion now is to plot multiple quorum levels
  quorum_levels <- round(seq(from = 0.3, to = 0.6, by = 0.1), 1)
  
  for(t in 1:1000){
    tryCatch({
      for(i in 1:length(quorum_levels)) {
        # estimateD function in iNEXT, running over each quorum level
        estD_tmp <- iNEXT::estimateD(freq_data, q = c(0,1,2), datatype = "incidence_freq", 
                                     base = "coverage", level = quorum_levels[i])
        # filter to the diversity estimates (order = 1):
        estD_tmp <- dplyr::filter(estD_tmp, Order.q == 0)
        # organise the output:
        estD_tmp$quorum_level <- quorum_levels[i]
        names(estD_tmp)[names(estD_tmp) == "Assemblage"] <- "name"
        
        estD_tmp <- left_join(estD_tmp, bin.dataset, by = "name")
        # add the output to the newly created list
        estD_output[[i]] <- estD_tmp
      }
      print(estD_output)
      if(length(estD_output) > 0){
        print("success")
        break
      }
    }, error = function(e){print(paste("Failed run ", t, sep = "")) 
      print(conditionMessage(e))})
  }
  return(estD_output)
}


#########################
##### RAW DIVERSITY #####
#########################

##### DIVERSITY ANALYSIS #####
# Run function to get diversity at chosen rank
m.only.div <- run.div(m.only.dat, "Genus", type = "raw")
pb.only.div <- run.div(pb.only.dat, "Genus", type = "raw")
pub.all.div <- run.div(pub.all.dat, "Genus", type = "raw")
all.div <- run.div(all.dat, "Genus", type = "raw") 

comb.div <- setup.div(pbddata = pb.only.div,
                  musdata = m.only.div, 
                  pubdata = pub.all.div, 
                  alldata = all.div)

# Pivot to long to make plot
div.tot.long <- tidyr::pivot_longer(comb.div, cols = c("Museum_SIB", 
                                                      "PBDB_SIB", 
                                                      "Pub_SIB", 
                                                      "All_SIB"))

names(div.tot.long)[names(div.tot.long) == "name"] <- "Data"
div.tot.long$Data[div.tot.long$Data == "Museum_SIB"] <- "Museum"
div.tot.long$Data[div.tot.long$Data == "PBDB_SIB"] <- "PBDB"
div.tot.long$Data[div.tot.long$Data == "Pub_SIB"] <- "Published"
div.tot.long$Data[div.tot.long$Data == "All_SIB"] <- "All"

# Plot
(p12 <- ggplot(div.tot.long, aes(x=bin_midpoint, y=value)) + 
  geom_line(aes(colour = Data)) +
  scale_x_reverse() +
  theme_bw() +
  ylab("Raw Diversity") +
  xlab("Time (Ma)") +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  scale_colour_viridis_d() +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            ylim = c(0, 20),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  theme(legend.position = "bottom")+
  geom_point(aes(color = Data)) )

ggsave("Dark_data_graphs/12.Raw.div.png", plot = p12, 
       device = "png", type = "cairo")

# Work out collections 
comb.div <- all.dat %>%
  dplyr::group_by(bin_assignment) %>%
  dplyr::summarize(Echinoid_collections = n_distinct(Locality)) %>%
  dplyr::full_join(comb.div, by = "bin_assignment")
comb.div <- pb.echino %>%
  dplyr::group_by(bin_assignment) %>%
  dplyr::summarize(PB_Echino_collections = n_distinct(collection_name)) %>%
  dplyr::full_join(comb.div, by = "bin_assignment")
comb.div <- pb.all %>%
  dplyr::group_by(bin_assignment) %>%
  dplyr::summarize(PBDB_collections = n_distinct(collection_name)) %>%
  dplyr::full_join(comb.div, by = "bin_assignment") %>%
  dplyr::filter(!is.na(bin_midpoint)) %>%
  replace(is.na(.), 0)

##### CORRELATIONS #####
cor1 <- cor.test(log10(comb.div$Museum_SIB), log10(comb.div$PBDB_SIB), method = "spearman")
cor2 <- cor.test(log10(comb.div$Museum_SIB), log10(comb.div$Pub_SIB), method = "spearman")
cor3 <- cor.test(log10(comb.div$Museum_SIB), log10(comb.div$All_SIB), method = "spearman")
cor4 <- cor.test(log10(comb.div$Pub_SIB), log10(comb.div$PBDB_SIB), method = "spearman")
cor5 <- cor.test(log10(comb.div$All_SIB), log10(comb.div$PBDB_SIB), method = "spearman")
cor6 <- cor.test(log10(comb.div$All_SIB), log10(comb.div$Pub_SIB), method = "spearman")

div.corr.1 <-rbind(c("Museum", "PBDB", cor1$estimate, cor1$p.value),
                   c("Museum", "Published", cor2$estimate, cor2$p.value),
                   c("Museum", "All", cor3$estimate, cor3$p.value),
                   c("PBDB", "Published", cor4$estimate, cor4$p.value),
                   c("PBDB", "All", cor5$estimate, cor5$p.value), 
                   c("Published", "All", cor6$estimate, cor6$p.value))

cor7 <- cor.test(log10(comb.div$Museum_SIB), log10(comb.div$Echinoid_collections), method = "spearman")
cor8 <- cor.test(log10(comb.div$Museum_SIB), log10(comb.div$PB_Echino_collections), method = "spearman")
cor9 <- cor.test(log10(comb.div$Museum_SIB), log10(comb.div$PBDB_collections), method = "spearman")

div.corr.2 <-rbind(c("Museum", "Dataset collections", cor7$estimate, cor7$p.value),
                   c("Museum", "PBDB Echinodermata collections", cor8$estimate, cor8$p.value),
                   c("Museum", "PBDB invertebrate collections", cor9$estimate, cor9$p.value))

cor10 <- cor.test(log10(comb.div$PBDB_SIB), log10(comb.div$Echinoid_collections), method = "spearman")
cor11 <- cor.test(log10(comb.div$PBDB_SIB), log10(comb.div$PB_Echino_collections), method = "spearman")
cor12 <- cor.test(log10(comb.div$PBDB_SIB), log10(comb.div$PBDB_collections), method = "spearman")

div.corr.3 <-rbind(c("PBDB", "Dataset collections", cor10$estimate, cor10$p.value),
                   c("PBDB", "PBDB Echinodermata collections", cor11$estimate, cor11$p.value),
                   c("PBDB", "PBDB invertebrate collections", cor12$estimate, cor12$p.value))

cor13 <- cor.test(log10(comb.div$Pub_SIB), log10(comb.div$Echinoid_collections), method = "spearman")
cor14 <- cor.test(log10(comb.div$Pub_SIB), log10(comb.div$PB_Echino_collections), method = "spearman")
cor15 <- cor.test(log10(comb.div$Pub_SIB), log10(comb.div$PBDB_collections), method = "spearman")

div.corr.4 <-rbind(c("Published", "Dataset collections", cor13$estimate, cor13$p.value),
                   c("Published", "PBDB Echinodermata collections", cor14$estimate, cor14$p.value),
                   c("Published", "PBDB invertebrate collections", cor15$estimate, cor15$p.value))

cor16 <- cor.test(log10(comb.div$All_SIB), log10(comb.div$Echinoid_collections), method = "spearman")
cor17 <- cor.test(log10(comb.div$All_SIB), log10(comb.div$PB_Echino_collections), method = "spearman")
cor18 <- cor.test(log10(comb.div$All_SIB), log10(comb.div$PBDB_collections), method = "spearman")

div.corr.5 <-rbind(c("All", "Dataset collections", cor16$estimate, cor16$p.value),
                   c("All", "PBDB Echinodermata collections", cor17$estimate, cor17$p.value),
                   c("All", "PBDB invertebrate collections", cor18$estimate, cor18$p.value))

# Combine results
cor.results <- as.data.frame(rbind(div.corr.1, div.corr.2, div.corr.3, 
                                   div.corr.4, div.corr.5))
colnames(cor.results) <- c("Variable 1", "vs. Variable 2", "Rho", "p")
cor.results$Rho <- signif(as.numeric(cor.results$Rho), digits = 3)
cor.results$p <- signif(as.numeric(cor.results$p), digits = 5)

# Correct for multiple tests
cor.results$BH <- p.adjust(cor.results$p, method = "BH")
cor.results$Signif <- ifelse(cor.results$BH < 0.05, "*", "")

cor.results %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  save_kable(file = "Dark_data_graphs/T2.raw.div.corrs.html")

##################################
##### SQS DIVERSITY - DIVDYN #####
##################################

##### DIVERSITY ANALYSIS #####
# Generate results
divdyn.m <- run.div(m.only.dat, Rank = "Genus", type = "SQS") 
divdyn.pb <- run.div(pb.only.dat, Rank = "Genus", type = "SQS") 
divdyn.p <- run.div(pub.all.dat, Rank = "Genus", type = "SQS") 
divdyn.a <- run.div(all.dat, Rank = "Genus", type = "SQS") 

# Setup for plotting
a <- divDyn.plot(divdyn.m, "none")
b <- divDyn.plot(divdyn.pb,"none")
c <- divDyn.plot(divdyn.p, "none")
d <- divDyn.plot(divdyn.a, "A")

# Plot together
(p13 <- ggarrange(a, b, c, d, 
          align='hv',
          labels = c("A", "B", "C", "D"),
          nrow = 4))

ggsave("Dark_data_graphs/13.SQS.div.divDyn.png", plot = p13, 
       device = "png", type = "cairo")

##### CORRELATIONS #####

#################################
##### SQS DIVERSITY - INEXT #####
#################################

##### DIVERSITY #####
# Get frequency data
freq.p <- freq.data(pub.all.dat, Rank = "Genus", stages)
freq.m <- freq.data(m.only.dat, Rank = "Genus", stages)
freq.a <- freq.data(all.dat, Rank = "Genus", stages)
freq.pb <- freq.data(pb.only.dat, Rank = "Genus", stages)

# Run iNEXT
iNEXT.p <- run.iNEXT(freq.p, stages)
iNEXT.m <- run.iNEXT(freq.m, stages)
iNEXT.a <- run.iNEXT(freq.a, stages)
iNEXT.pb <- run.iNEXT(freq.pb, stages)

## The output is a 'list' object so we'll need to convert it into a dataframe and clean it up before plotting
test.list <- list(All = iNEXT.a,
                  Published = iNEXT.p,
                  Museum = iNEXT.m,
                  PBDB = iNEXT.pb)
estD_plotting <- bind_rows(lapply(seq_along(test.list), function(x){
  name <- names(test.list[x])
  test <- bind_rows(test.list[x])
  test$Data <- name
  return(test)
}))

## Ensure that the quorum level column is being treated as a 'factor' to avoid errors while plotting:
estD_plotting$quorum_level <- as.factor(estD_plotting$quorum_level)

## Create a colour gradient for as many colours as you have quorum levels:
teal_gradient <- scales::seq_gradient_pal("turquoise", "darkslategrey", "Lab")(seq(0, 1, length.out = 4))

(p14 <- ggplot(estD_plotting, aes(x = bin_midpoint, y = qD, ymin = qD.LCL, ymax = qD.UCL, colour = quorum_level)) + 
  ## Each quorum level is called individually to be plotted:
  geom_ribbon(data = subset(estD_plotting, quorum_level == 0.3), 
              aes(x = bin_midpoint, ymin = qD.LCL, ymax = qD.UCL), 
              inherit.aes = FALSE, fill = teal_gradient[1], alpha = 0.2) +
  geom_ribbon(data = subset(estD_plotting, quorum_level == 0.4), 
              aes(x = bin_midpoint, ymin = qD.LCL, ymax = qD.UCL), 
              inherit.aes = FALSE, fill = teal_gradient[2], alpha = 0.2) +
  geom_ribbon(data = subset(estD_plotting, quorum_level == 0.5), 
              aes(x = bin_midpoint, ymin = qD.LCL, ymax = qD.UCL), 
              inherit.aes = FALSE, fill = teal_gradient[3], alpha = 0.2) +
  geom_ribbon(data = subset(estD_plotting, quorum_level == 0.6), 
              aes(x = bin_midpoint, ymin = qD.LCL, ymax = qD.UCL), 
              inherit.aes = FALSE, fill = teal_gradient[4], alpha = 0.2) +
  ## Set our line and point sizes (and shapes):
  geom_line(linewidth = 1) +
  ylim(0, 6) +
  geom_point(aes(pch = Method), size = 3) +
  scale_shape_manual(values=c(15, 16, 17)) +
  ## Add our colours, theme, and axes labels:
  scale_colour_manual(values = teal_gradient) +
  labs(color='Quorum level') +
  scale_x_reverse() +
  xlim(477.7000, 251.9020) +
  labs(x = "Time (Ma)", y = "Coverage rarified richness") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  theme(legend.position = "bottom") + 
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  facet_wrap(~Data, nrow = 4))

ggsave("Dark_data_graphs/14.SQS.div.iNEXT.png", plot = p14, 
       device = "png", type = "cairo")

###############################
##### SAMPLE COMPLETENESS #####
###############################

# Setup
#series$bin <- seq(1, 6, 1)
series$bin_midpoint <- series$mid_ma

# Get frequency data
freq.pb.s <- freq.data(pb.only.series, "Genus", series, filter_data = T)
freq.m.s <- freq.data(m.only.series, "Genus", series, filter_data = F)
freq.p.s <- freq.data(pub.all.series, "Genus", series, filter_data = F)
freq.a.s <- freq.data(all.series, "Genus", series, filter_data = F)

# Run iNEXT
estD_m <- iNEXT::iNEXT(freq.m.s, q = 1, datatype = "incidence_freq")
estD_pb <- iNEXT::iNEXT(freq.pb.s, q = 1, datatype = "incidence_freq")
estD_p <- iNEXT::iNEXT(freq.p.s , q = 1, datatype = "incidence_freq")
estD_a <- iNEXT::iNEXT(freq.a.s , q = 1, datatype = "incidence_freq")

a <- ggiNEXT(estD_m, type=1, color.var = "Assemblage")
a <- a + scale_fill_manual(values=c(wes_palette("Zissou1"), 
                                   wes_palette("Royal1"))) + 
  scale_color_manual(values=c(wes_palette("Zissou1"), 
                             wes_palette("Royal1")))  +
  theme_bw() +
  theme(legend.text = element_text(size = 10), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white"))
b <- ggiNEXT(estD_pb, type=1, color.var = "Assemblage")
b <- b + scale_fill_manual(values=c("#3B9AB2", "#78B7C5", "#EBCC2A", "#F21A00")) + 
  scale_color_manual(values=c("#3B9AB2", "#78B7C5", "#EBCC2A", "#F21A00")) +
  theme_bw() +
  theme(legend.text = element_text(size = 10), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white"))

c <- ggiNEXT(estD_p, type=1, color.var = "Assemblage")
c <- c + scale_fill_manual(values=c(wes_palette("Zissou1"), 
                               wes_palette("Royal1"))) + 
  scale_color_manual(values=c(wes_palette("Zissou1"), 
                              wes_palette("Royal1")))  +
  theme_bw() +
  theme(legend.text = element_text(size = 10), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white"))
d <- ggiNEXT(estD_a, type=1, color.var = "Assemblage")
d <- d + scale_fill_manual(values=c(wes_palette("Zissou1"), 
                               wes_palette("Royal1"))) + 
  scale_color_manual(values=c(wes_palette("Zissou1"), 
                              wes_palette("Royal1")))  +
  theme_bw() +
  theme(legend.text = element_text(size = 10), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white"))

# Combine plots
(p15 <- ggarrange(a, b, c, d, 
          align='v',
          labels = c("A", "B", "C", "D"),
          nrow = 2, 
          ncol = 2, 
          common.legend = T,
          legend = "bottom"))

ggsave("Dark_data_graphs/15.Sample.completeness.png", plot = p15, 
       device = "png", type = "cairo")

################################################################################
# 8. SEDIMENT AFFINITY
################################################################################

#####################
##### FUNCTIONS #####
#####################

###############
##### SRA #####
###############

make.fam.lith <- function(data){
  fam.lith.data <- data %>%
    dplyr::filter(Rank == "Family"| Rank == "Genus" | Rank == "Species") %>%
    dplyr::filter(Family != "NO_FAMILY_SPECIFIED" & 
                  Family != "Archaeocidaridae or miocidaridae") %>%
    dplyr::filter(!is.na(Family)) %>%
    dplyr::filter(!is.na(Finalised_lith)) %>%
    dplyr::filter(Finalised_lith != "Mixed") %>%
    dplyr::filter(Finalised_lith != "")
  return(fam.lith.data)
}

setup.SRA <- function(fam.lith.data, bin.dataset, data.source){
  families <- c("Archaeocidaridae", "Lepidesthidae", "Lepidocentridae", 
                "Palaechinidae", "Proterocidaridae")
  fam_list <- list()
  
  for(f in families){
    print(f)
    fam_df <- fam.lith.data %>%
      dplyr::filter(Family == f)
    num.rows <- nrow(fam_df)
    if(num.rows < 2){
      next
    }
    fam_df <- fam_df %>%
      dplyr::group_by(bin_assignment) %>%
      {as.data.frame(table(.$bin_assignment, .$Finalised_lith))} %>%
      filter(Var2 != "") %>%
      tidyr::pivot_wider(names_from = Var2, 
                         values_from = Freq) %>%
      dplyr::select(Var1, Carbonate, Siliciclastic) %>%
      dplyr::rename(bin = Var1) %>%
      dplyr::mutate(Total = Carbonate + Siliciclastic, 
                    Ac = Carbonate/Total, 
                    As = Siliciclastic/Total)
    
    bin.dataset2 <- bin.dataset %>%
      dplyr::filter(bin %in% fam_df$bin)
    
    fam_df_time <- dplyr::full_join(bin.dataset2, fam_df)
    
    bin.list <- list()
    for(b in fam_df_time$bin){
      print(b)
      bin.list[[match(b, bin.dataset2$bin)]] <- test.func(b, 
                                                          fam_df_time, 
                                                          fam.lith.data = fam.lith.data, 
                                                          f)
    }
                      
    sim.aff <- bind_rows(bin.list)
    names(sim.aff)[names(sim.aff) == "Carbonate"] <- "sim.Ac"
    names(sim.aff)[names(sim.aff) == "Siliciclastic"] <- "sim.As"
    sim.aff <- sim.aff %>% dplyr::mutate_all(~replace(., is.nan(.), 0))
    sim.aff[is.na(sim.aff)] <- 0
    fam_df_time <- dplyr::full_join(fam_df_time, sim.aff, by = "bin")
    fam_df_time <- fam_df_time %>%
      dplyr::mutate(Carb_SRA = (Ac - sim.Ac)/SDrndC, 
                    Family = f)
    fam_df_time$Data <- data.source
    fam_list[[match(f, families)]] <- fam_df_time
    names(fam_list)[[match(f, families)]] <- f
  }
  return(fam_list)
}

test.func <- function(x, fam_df_time, fam.lith.data, f) {
  bin_total <- fam_df_time %>%
    dplyr::filter(bin == x)
  total.occs <- as.numeric(bin_total[names(bin_total) == "Total"])
  combined <- data.frame(Finalised_lith = c("Carbonate", "Siliciclastic"))
  for(i in 1:1000){
    tmp <- fam.lith.data %>% 
      dplyr::filter(bin_assignment == x) %>%
      dplyr::slice_sample(., n = total.occs, replace = TRUE) %>%
      dplyr::count(., Finalised_lith) %>%
      dplyr::mutate(Affinity = n/total.occs) %>%
      dplyr::select(-n)
    colnames(tmp)[2] <- i
    combined <- dplyr::full_join(combined, tmp, by = "Finalised_lith")
  }
  rownames(combined) <- combined[,1]
  combined[,1] <- NULL
  combined[is.na(combined)] <- 0
  comb.sd <- data.frame(bin = x, 
                        SDrndC = sd(combined[1,], na.rm = T), 
                        SDrndS = sd(combined[2,], na.rm = T))
  combined <- t(as.data.frame(rowMeans(combined, na.rm = T)))
  combined <- as.data.frame(combined)
  combined$bin <- x
  combined <- dplyr::full_join(combined, comb.sd, by = "bin")
  
  return(combined)
}

plot.SRA <- function(fam_list, bin.choice, fill_var, 
                     group_var, colour){
  test <- test %>%
    dplyr::select(Carb_SRA, Family, bin, Data) %>%
    dplyr::filter(!is.na(Carb_SRA)) %>%
    dplyr::filter(Carb_SRA != -Inf) %>%
    dplyr::filter(Carb_SRA != Inf)
  test <- left_join(test, bin.choice, by = "bin")
  a <- ggplot(test,   aes_string(x = "mid_ma", y= "Carb_SRA", color = fill_var)) + 
    ## Each quorum level is called individually to be plotted:
    ## Set our line and point sizes (and shapes):
    geom_line() +
    geom_point(size = 2) +
    ## Add our colours, theme, and axes labels:
    scale_colour_manual(values = wes_palette(colour)) +
    scale_x_reverse() +
    labs(x = "Time (Ma)", y = "Standardized Relative Affinity") +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "white"), 
          panel.grid.minor = element_line(colour = "white")) +
    geom_hline(yintercept = 0, 
               linetype = "dashed", 
               color = "light grey") +
    theme(plot.margin = margin(0.1,0.2,0, 1, "cm")) +
    coord_geo(dat = series2, 
              xlim = c(485.4, 251.902),
              ylim = c(-10, 6),
              rot = list(0),
              size = 4,
              pos = list("b"),
              abbrv = T) +
    facet_wrap(group_var)
  a <- shift_legend(a)
  a <- as_ggplot(a)
  return(a)
}


shift_legend <- function(p){
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}

# Series
lith.m <- make.fam.lith(m.only.series)
lith.pb <- make.fam.lith(pb.only.series)
lith.p <- make.fam.lith(pub.all.series)
lith.a <- make.fam.lith(all.series)

temp.museum <- setup.SRA(fam.lith.data = lith.m, bin.dataset = series, 
                         data.source = "Museum")
temp.pbdb <- setup.SRA(fam.lith.data = lith.pb, bin.dataset = series, 
                       data.source = "PBDB")
temp.published <- setup.SRA(fam.lith.data = lith.p, bin.dataset = series,
                            data.source = "Published")
temp.all <- setup.SRA(fam.lith.data = lith.a, bin.dataset = series, 
                      data.source = "All")

combined.SRA <- list(Museum = temp.museum, PBDB = temp.pbdb, Published = temp.published, 
             All = temp.all)

combined.SRA <- bind_rows(lapply(combined.SRA, function(f){f <- bind_rows(f)}))

(p16 <- plot.SRA(combined.SRA, series, 
                 fill_var = "Data", 
                 group_var = "Family", 
                 colour = "AsteroidCity1"))
(p16.5 <- plot.SRA(combined.SRA, series, 
                   fill_var = "Family", 
                   group_var = "Data", 
                   colour = "Darjeeling1"))

ggsave("Dark_data_graphs/16.SRA.png", plot = p16, 
       device = "png", type = "cairo")
ggsave("Dark_data_graphs/16.5.SRA.png", plot = p16.5, 
       device = "png", type = "cairo")

################################################################################
# 9. COLONIALISM
################################################################################

# For this section, code was reworked from that provided for the following publication:

# Raja, N.B., Dunne, E.M., Matiwane, A., Khan, T.M., Nätscher, P.S., Ghilardi, A.M. 
# and Chattopadhyay, D., 2022. Colonial history and global economics distort our 
# understanding of deep-time biodiversity. Nature ecology & evolution, 6(2), pp.145-154.

# We thank the authors for providing this code and for enabling this research. 

#####################
##### FUNCTIONS #####
#####################

# Function to run and plot network analysis
setup.col.sec <- function(data){
  ##### Setup data for analysis #####
  df <-  data %>%
    dplyr::select(Locality, Museum_country, Country, lat, lng) %>%
    distinct() %>%
    dplyr::filter(!is.na(Museum_country)) %>%
    dplyr::filter(Museum_country != "") %>%
    dplyr::filter(!is.na(Country))  %>%
    dplyr::rename(samp_country = Country, 
                  aff_country = Museum_country) %>%
    dplyr::mutate(collection_no = 1:nrow(.)) %>%
    dplyr::select(-c(Locality, lat, lng))
  
  df[df == "England" | df == "Wales" | df == "Scotland" | df == "Northern Ireland"] <- "United Kingdom"
  df[df == "USA"] <- "United States"
  df[df == "Czechia"] <- "Czech Republic"
  
  df$aff_code <- countrycode(df$aff_country, "country.name", "iso3c")
  df$samp_code <- countrycode(df$samp_country, "country.name", "iso3c")
  return(df)
}

network.plots <- function(df){
  # each country gets an id
  countries <- unique(c(df$aff_country, df$samp_country))
  countries <- gsub("&", "and", countries)
  countries <- data.frame(id=1:length(countries),
                          label=countries)
  countries$region <- countrycode::countryname(countries$label, 
                                               destination = "continent") #may change
  countries <- countries[order(countries$region),]
  
  regs <- unique(countries$region)
  cols <- palette.colors(length(regs))
  countries$cols <- plyr::mapvalues(countries$region, from=regs, to=cols)
  
  # By region
  df$samp_region <- countrycode::countryname(df$samp_country, destination="region23")
  df$aff_region <- countrycode::countryname(df$aff_country, destination="region23")
  
  df$samp_continent <- countrycode::countryname(df$samp_country, destination="continent")
  df$aff_continent <- countrycode::countryname(df$aff_country, destination="continent")
  
  # each region gets an id
  regions <- unique(c(df$samp_region, df$aff_region))
  regions <- data.frame(id=1:length(regions),
                        label=regions)
  
  ##### Create edge list #####
  edges <- df %>%  
    dplyr::group_by(samp_country, aff_country) %>% 
    dplyr::summarise(weight = n()) %>% 
    dplyr::ungroup()
  
  #add node id
  edges <- edges %>% 
    dplyr::left_join(countries, by = c("aff_country" = "label")) %>% 
    dplyr::rename(from = id)
  
  edges <- edges %>% 
    dplyr::left_join(countries %>%  dplyr::select(id, label), 
                     by = c("samp_country" = "label")) %>% 
    dplyr::rename(to = id)
  
  edges <- na.omit(edges)
  
  ###### Plot #####
  edges <- edges[edges$samp_country != "",] 
  edges2 <- edges[edges$samp_country != edges$aff_country,]
  
  df2 <- countries %>%  dplyr::filter(!is.na(region))
  df2 <- df2 %>% dplyr::left_join(edges %>%  dplyr::select(aff_country, weight), 
                                  by=c("label"="aff_country"))
  
  df2$weight[is.na(df2$weight)] <- 0
  #df2 <- df2[df2$label %in% unique(edges2$samp_country, edges2$aff_country),]
  
  df2 <- df2 %>%  dplyr::group_by(label, region) %>% 
    dplyr::summarise(weight=sum(weight))
  
  df2 <- df2[order(df2$region, df2$weight, decreasing = T),]
  regs <- unique(df2$region)
  
  df2$y <- 0
  df2$x <- 0
  
  seq_max <- c(20,22,21, 22, 18)
  x_mid <- c(5,4,2,3,1)
  f <- c(3,5,8,5,10)
  
  for(i in 1:length(regs)){
    temp <- df2[df2$region == regs[i],]
    df2[df2$region == regs[i],]$y <- seq(seq_max[i],1, length.out = nrow(temp))
    set.seed(42)
    df2[df2$region == regs[i],]$x <- jitter(rep(x_mid[i], nrow(temp)), factor=f[i])
  }
  
  df2[df2$label=="United States",]$y <-  24
  df2[df2$label=="China",]$y <- 22
  df2[df2$label=="Germany",]$y <- 22.8
  
  edges2 <- edges2 %>% 
    dplyr::left_join(df2 %>%  dplyr::select(label, y,x), 
                     by=c("samp_country" = "label")) %>% 
    dplyr::rename(x2 = x, y2=y) %>% 
    dplyr::left_join(df2 %>%  dplyr::select(label, y,x), 
                     by=c("aff_country" = "label")) %>% 
    dplyr::rename(x1=x, y1=y)
  
  p.breaks <- df2 %>% dplyr::group_by(region) %>% 
    dplyr::summarise(x=median(x)) %>% 
    arrange(x)
  
  edges2$region[edges2$weight < 26] <- "Other"
  
  edges2$region <- factor(edges2$region,levels=c(regs, "Other"))
  df2$weight_alpha <- df2$weight
  df2$weight_alpha[df2$label %in% df2$label[order(df2$weight, decreasing = T)][1:25]] <- 9999
  
  curvature = 0.3
  alpha = 0.8
  breaks <- c(50,100,500,1000,
              2000, 5000, 10000)
  lab_size <- 3
  
  palv <- viridis::viridis(5)
  names(palv) <- c("Europe", "Americas", "Oceania", "Asia", "Africa")
  palv["Europe"] <- "#75428f"
  palv["Americas"] <- "#365091"
  palv["Africa"] <- "#fdaa25"
  palv["Asia"] <- "#027218"
  
  p1 <- ggplot() + geom_curve(data=na.omit(edges2[edges2$x2 > edges2$x1,]), 
                              aes(x=x1, y=y1, xend=x2, yend=y2, col=region
                              ), 
                              curvature = curvature, 
                              alpha=alpha,
                              size=0.5) +
    geom_curve(data=na.omit(edges2[edges2$x1 > edges2$x2,]), 
               aes(x=x1, y=y1, xend=x2, yend=y2, col=region), 
               curvature = -curvature, alpha=alpha) +
    geom_point(data=df2, 
               aes(x=x, y=y, 
                   size=weight, col=region,
                   alpha=weight_alpha)) +
    geom_text(data=df2[df2$x < 3,], aes(x=x, y=y, 
                                        label=ifelse(label %in% df2$label[order(df2$weight, decreasing = T)][1:25], label, NA)),
              hjust=1, vjust=0, size=lab_size, fontface="bold", nudge_x = -0.1) + 
    geom_text(data=df2[df2$x > 3,], aes(x=x, y=y, 
                                        label=ifelse(label %in% df2$label[order(df2$weight, decreasing = T)][1:25], label, NA)),
              hjust=0, vjust=0, size=lab_size, nudge_y = 0.2, nudge_x = 0.1, fontface="bold") +
    scale_x_continuous(breaks=p.breaks$x, limits=c(0,6)) +
    scale_size_continuous(breaks=breaks, 
                          range=c(1,15)) +
    scale_color_manual(values=c("Other"="#e5e5e560", 
                                palv)) +
    scale_alpha_continuous(breaks=c(0, 50,100,200,500, 1000,2000, 9999), 
                           range=c(0.3,1)) +
    guides(alpha="none") +
    labs(col="Region", size="Number of \noutgoing nodes") +
    theme_void() 
  return(p1)
}

# Function for setting up network statistics
reshape.df <- function(x){
  x <- sort(x, decreasing = TRUE)[1:15]
  names <- names(x)
  x <- data.frame(value=x)
  if(rownames(x)[1] == "1"){
    x <- na.omit(x)
    rownames(x) <-  names[complete.cases(names)]
  }
  x$code <- rownames(x)
  x$country <- countrycode(x$code, "iso3c", "country.name")
  return(x)
}

# Function to plot foreign vs. domestic contributions to fossil collections
research.prep <- function(df, data){
  colls_n <- data.frame(table(df$aff_code), stringsAsFactors = FALSE)
  colnames(colls_n) <- c("code", "freq")
  total_colls <- length(unique(df$collection_no))
  colls_n$freq <- colls_n$freq/total_colls
  colls_n <- colls_n[order(colls_n$freq, decreasing = TRUE),]
  topcountries <- colls_n[1:15,]
  
  topcountries <- na.omit(topcountries)
  topcountries$country <- countrycode(topcountries$code, 
                                      origin = "iso3c", 
                                      destination = "country.name")
  
  dat2 <- df[df$aff_code != df$samp_code,]
  dat2 <- unique(dat2[,c("collection_no", "aff_code")])
  
  colls_n2 <- data.frame(table(dat2$aff_code), stringsAsFactors = FALSE)
  colnames(colls_n2) <- c("code", "foreign")
  colls_n2$foreign <- colls_n2$foreign/total_colls
  
  colls_foreign <-colls_n2[colls_n2$code %in% topcountries$code,]
  
  topcountries <- merge(topcountries, colls_n2, all = TRUE)
  topcountries[is.na(topcountries)] <- 0
  topcountries$local <- topcountries$freq - topcountries$foreign
  
  if(nrow(topcountries) > 15){
    topcountries <- topcountries[order(topcountries$freq, decreasing = TRUE),]
    topcountries <- topcountries[1:15,]
  }
  
  topcountries <- topcountries %>%  dplyr::select(local, foreign, country) %>% 
    pivot_longer(cols=c("local", "foreign"), names_to = c("type"), values_to="freq")
  
  topcountries$type <- factor(topcountries$type, levels=c("local", "foreign"))
  topcountries$Data <- data
  
  return(topcountries)
}

continent.stats <- function(df, data){
  codes <- c("samp_code", "aff_code")
  iso2c_to_continents <- c(CAN = "North America" , USA = "North America")
  
  
  a <- bind_rows(lapply(codes, function(x){
    dat2 <- df[c(x, "collection_no")]
    dat2 <- unique(dat2[,c("collection_no", x)])
    dat2$Continent <- countrycode(dat2[,x], 
                                  origin = "iso3c", 
                                  destination = "continent", 
                                  custom_match = iso2c_to_continents)
    
    dat2 <- data.frame(table(dat2$Continent), stringsAsFactors = FALSE)
    
    colnames(dat2) <- c("Continent", "Total")
    
    dat2$Frequency <- dat2$Total/sum(dat2$Total)
    dat2$Origin <- ifelse(x == "samp_code", "Country", "Museum")
    dat2 <- dat2 %>%
      mutate(Continent = str_replace(Continent, "Americas", "South America"))
    return(dat2)
  }))
  a$Data <- data
  return(a)
}

top.countries <- function(df, code, sum = NA){
  colls_n <- data.frame(table(df[code]), stringsAsFactors = FALSE)
  colnames(colls_n) <- c(code, "freq")
  total_colls <- length(unique(df$collection_no))
  colls_n$freq <- colls_n$freq/total_colls
  colls_n <- colls_n[order(colls_n$freq, decreasing = TRUE),]
  topcountries <- colls_n[1:15,]
  
  topcountries <- na.omit(topcountries)
  if(is.na(sum) == T){
    return(topcountries)
  }else{
    perc <- sum(topcountries$freq[1:sum])*100
    return(list(topcountries[1:sum,], perc))
  }
}

#################
##### SETUP #####
#################

m.df <- setup.col.sec(m.only.dat)
p.df <- setup.col.sec(pub.all.dat)
a.df <- setup.col.sec(all.dat)

#############################
##### TOP CONTRIBUTIONS #####
#############################

top.countries(m.df, code = "aff_code", 5)
top.countries(m.df, code = "samp_code", 5)
top.countries(p.df, code = "aff_code", 5)
top.countries(p.df, code = "samp_code", 5)
top.countries(a.df, code = "aff_code", 5)
top.countries(a.df, code = "samp_code", 5)

All.cont <- bind_rows(list(continent.stats(a.df, "All"),
               continent.stats(p.df, "Published"),
               continent.stats(m.df, "Museum")))

ggplot(All.cont, aes(fill = Origin, y = Frequency, x = Continent)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  coord_flip() +
  facet_wrap(~Data, nrow = 3)

All.cont %>%
  dplyr::group_by(Data, Origin) %>%
  dplyr::filter(Continent == "North America" | Continent == "Europe") %>%
  dplyr::summarise(Total_freq = sum(Frequency))

Research <- bind_rows(list(research.prep(a.df, data = "All"), 
                    research.prep(p.df, data = "Published"),
                    research.prep(m.df, data = "Museum")))

reorder(geo, values, sum, na.rm=T)

(p17 <- ggplot(Research, aes(x=reorder(country, freq, sum), y=freq*100, fill=type)) +
  geom_bar(stat="identity") +
  labs(x="", y=" Percentage contribution to fossil collections", fill="Fieldwork") +
  scale_fill_manual(values=wes_palette("Moonrise2")[c(3,1)], 
                    labels=c("Domestic research", "Foreign research")) +
  theme_bw() +
    theme(legend.position = "bottom", 
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  coord_flip() +
  ylim(0, 60) +
  guides(fill=guide_legend(ncol=2)) +
  facet_wrap(~Data, nrow = 3))

ggsave("Dark_data_graphs/17.Contribution.png", plot = p17, 
       device = "png", type = "cairo")

######################################
##### NETWORK ANALYSIS AND PLOTS #####
######################################

a <- network.plots(m.df)
b <- network.plots(p.df)
c <- network.plots(a.df)

(p18 <- ggarrange(a, b, c, 
          align='hv',
          labels = c("A", "B", "C"),
          nrow = 1, 
          ncol = 3,
          common.legend = T, 
          legend = "bottom"))

ggsave("Dark_data_graphs/18.Extraction.png", plot = p18, 
       device = "png", type = "cairo")

pdf("Dark_data_graphs/18.Extraction.pdf", 
    width = par("din")[1], 
    height = par("din")[2])
print(p18)
dev.off()

########################################
##### NETWORK STATISTICS AND PLOTS #####
########################################

network.stats <- function(df){
  
  pal <- c("#f0ffe9", "#ffe599", "#bbe487", "#4e9755", "#173109")
  
  # Setup
  gr <- df[,c("aff_code", "samp_code")]
  gr <- gr[gr$aff_code != gr$samp_code,] # remove self nodes
  gr <- igraph::graph_from_data_frame(gr, directed = T)
  
  # The degree of a node in a network is the number of connections it has to other 
  # nodes and the degree distribution is the probability distribution of these degrees 
  # over the whole network. 
  gr$degree <- igraph::degree(gr) 
  
  # Detecting the amount of influence a node has over the flow of information in a graph. 
  gr$betweenness <- igraph::betweenness(gr) 
  
  # The more central a node is, the closer it is to all other nodes. 
  gr$closeness <- igraph::closeness(gr) 
  
  # Run function
  degree_top15 <- reshape.df(gr$degree)
  betweenness_top15 <- reshape.df(gr$betweenness)
  closeness_top15 <- reshape.df(gr$closeness)
  
  # Plot Degree
  p3 <- ggplot(degree_top15, aes(x=reorder(country, value), y=value)) +
    geom_bar(stat="identity", fill=pal[4]) + coord_flip(expand=FALSE) +
    labs(y="Degree", x="") +
    ggthemes::theme_hc()
  
  # Plot betweenness
  p4 <- ggplot(betweenness_top15, aes(x=reorder(country, value), y=value)) +
    geom_bar(stat="identity", fill=pal[4]) + coord_flip(expand=FALSE) +
    labs(y="Betweenness", x="") +
    ggthemes::theme_hc()
  
  # Plot closeness
  p5 <- ggplot(closeness_top15, aes(x=reorder(country, value), y=value)) +
    geom_bar(stat="identity", fill=pal[4]) + coord_flip(expand=FALSE) +
    labs(y="Closeness", x="") +
    ggthemes::theme_hc()
  
  # Combine plots
  p6 <- ggarrange(p3, p4, p5, 
                  align='hv',
                  labels = c("A", "B", "C"),
                  nrow = 1, 
                  ncol = 3)
  return(p6)
}

(p19 <- network.stats(m.df))
(p20 <- network.stats(p.df))
(p21 <- network.stats(a.df))

ggsave("Dark_data_graphs/19.Museum.network.stats.png", plot = p19, 
       device = "png", type = "cairo")
ggsave("Dark_data_graphs/20.Published.network.stats.png", plot = p20, 
       device = "png", type = "cairo")
ggsave("Dark_data_graphs/21.All.data.network.stats.png", plot = p21, 
       device = "png", type = "cairo")


################################################################################
# 10. MODELLING COMPARISONS
################################################################################

#####################
##### FUNCTIONS #####
#####################

run.models <- function(dataset, bin.dataset, div = "raw", na = T){
 
  if(na == T){
    # Reduce to N. America
    dataset <- dataset %>%
      dplyr::filter(Continent == "North America")
  } 
  
  ##### DIVERSITY #####
  # Get diversity estimates
  if(div == "raw"){
    model.div <- run.div(dataset, "Genus", type = "raw")
    model.div <- model.div %>%
      dplyr::select(bin_assignment, divSIB) %>%
      dplyr::rename(bin = bin_assignment)
  } 
  if(div == "SQS"){
    model.div <- run.div(dataset, "Genus", type = "SQS")
    # Keep relevant data only
    model.div <- model.div[[4]] %>%
      dplyr::select(bin_assignment, divCSIB) %>%
      dplyr::rename(bin = bin_assignment)
  }
  comb.mod <- full_join(comb.mod, model.div, by = "bin")
  if(div == "raw"){
    comb.mod$divSIB <- log10(comb.mod$divSIB+1) # Add 1 to stop infinity appearing
  }else{
    comb.mod$divSIB <- log10(comb.mod$divCSIB)
  }

  comb.mod <- comb.mod %>%
    dplyr::filter(bin_midpoint > 253 & bin_midpoint < 481.55000)
  
  div.model <- comb.mod %>%
    dplyr::filter(!is.na(divSIB)) %>%
    dplyr::filter(!is.na(colls)) %>%
    dplyr::filter(!is.na(hardie)) 
  
  if(na == T){
    mod.form <- as.formula("divSIB ~ mean_sl + carb + sili + u + colls + hardie + wilkinson")
  }else{
    mod.form <- as.formula("divSIB ~ mean_sl + colls + u + hardie + wilkinson")
  }

  # Full model
  full.model <- gls(mod.form, data = div.model, method = "ML", 
                    correlation = corARMA(p=1))

  dredged <- dredge(full.model, rank = "AICc")
  top_model <- get.models(dredged, subset = 1)[[1]]
  
  # model results
  model.results <- list()
  model.results[[1]] <- summary(top_model)
  if(dredged$delta[2] < 2){
    model.results[[2]] <- summary(model.avg(dredged, subset = delta <= 2))
  }else{
    model.results[[2]] <- summary(model.avg(dredged, subset = delta <= dredged$delta[2]))
  }
  model.results[[3]] <- dredged
  model.results[[4]] <- r.squaredLR(top_model)
  names(model.results) <- c("Top.model", "Mod.Av", "All.mods", "R.squared")
  return(model.results)
}

ca.binning <- function(data){
  test <- data %>%
    dplyr::filter(max_ma < max(stages$max_ma)) %>%
    dplyr::filter(min_ma > min(stages$min_ma))
  test <- bin_time(test, stages, method = 'mid')
  test <- test %>%
    dplyr::group_by(bin_midpoint) %>%
    dplyr::summarize(mean = mean(MgCa))
  return(test)
}

model.results <- function(model.list){
  # Getting top model data
  a <- bind_rows(lapply(seq_along(model.list), function(x){
    name <- names(model.list[x])
    test <- as.data.frame(model.list[[x]][["All.mods"]][1,])
    test$Data <- name
    return(test)
  }))
  # Getting top model covariates/signif
  b <- bind_rows(lapply(seq_along(model.list), function(x){
    name <- names(model.list[x])
    test <- as.data.frame(model.list[[x]][["Top.model"]]["tTable"])
    colnames(test) <- c("Estimate", "Std. Error", 
                        "T value", "p value")
    test <- cbind(data.frame(Data = name, Covariate = rownames(test)), test)
    rownames(test) <- NULL
    test$Signif <- ifelse(test$`p value`< 0.05, "*", "")
    return(test)
  }))
  # Getting model average data
  c <- bind_rows(lapply(seq_along(model.list), function(x){
    name <- names(model.list[x])
    test <- as.data.frame(model.list[[x]][["Mod.Av"]]["coefmat.full"])
    colnames(test) <- c("Estimate", "Std. Error", 
                        "Adjusted SE", "Z value", "p value")
    test <- cbind(data.frame(Data = name, Covariate = rownames(test)), test)
    rownames(test) <- NULL
    test$Signif <- ifelse(test$`p value`< 0.05, "*", "")
    return(test)
  }))
  d <- bind_rows(lapply(seq_along(model.list), function(x){
    name <- names(model.list[x])
    test <- as.data.frame(model.list[[x]][["R.squared"]][1])
    test <- cbind(data.frame(Data = name), test)
    colnames(test) <- c("Data", "R.squared")
    return(test)
  }))
  return(list(a, b, c, d))
}

###########################
##### COVARIATE SETUP #####
###########################

# Get echinoderm collections
NA.echino <- pb.echino %>%
  dplyr::filter(cc == "US" | cc == "CA" | cc == "MX")
echino.NA.col <- binstat(NA.echino, tax = 'genus', coll = "Locality", 
                         bin = "bin_assignment", noNAStart = T)
names(echino.NA.col)[names(echino.NA.col) == "bin_assignment"] <- "bin"

# Sea level
sea.lvl <- read.csv("vanderMeer_2022.csv")
sea.lvl$max_ma <- sea.lvl$Ma+0.0001
names(sea.lvl)[names(sea.lvl) == "Ma"] <- "min_ma"
sea.lvl <- sea.lvl %>%
  dplyr::filter(max_ma < 538) %>%
  dplyr::filter(min_ma > 0)
sea.lvl <- bin_time(sea.lvl, stages, method = "mid")
sea.lvl <- sea.lvl %>%
  dplyr::group_by(bin_midpoint) %>%
  dplyr::summarize(mean_sl = mean(TGE_SL_isocorr_m))

# Macrostrat
macro <- as.data.frame(pivot_wider(macro.area, names_from = lith, 
                                   values_from = count))

# MgCa ratios
ca1 <- read.csv("Horita_data_extracted.csv")
ca2 <- read.csv("Hardie_data_extracted.csv")
ca3 <- read.csv("Wilkinson+Algeo_data_extracted.csv")

ca1 <- ca.binning(ca1)
ca2 <- ca.binning(ca2)
ca3 <- ca.binning(ca3)

MgCa <- full_join(ca1, ca2, by = 'bin_midpoint')
MgCa <- full_join(MgCa, ca3, by = 'bin_midpoint')
colnames(MgCa) <- c("bin_midpoint", "horita", "hardie", "wilkinson")

# Setup combined data
comb.mod <- stages %>%
  dplyr::select(name, bin, bin_midpoint)

# Combine data
comb.mod <- full_join(comb.mod, sea.lvl, by = "bin_midpoint")
comb.mod <- full_join(comb.mod, macro, by = "bin_midpoint")
comb.mod <- full_join(comb.mod, echino.NA.col, by = "bin")
comb.mod <- full_join(comb.mod, MgCa, by = 'bin_midpoint')

comb.mod <- comb.mod %>%
  dplyr::filter(bin_midpoint > 253 & bin_midpoint < 481.55000)

# Standardise data
comb.mod$mean_sl <- log10(comb.mod$mean_sl)
comb.mod$carb <- log10(comb.mod$carb)
comb.mod$sili <- log10(comb.mod$sili)
comb.mod$colls <- log10(comb.mod$colls)
comb.mod$hardie <- log10(comb.mod$hardie)
comb.mod$horita <- log10(comb.mod$horita)
comb.mod$wilkinson <- log10(comb.mod$wilkinson)

########################################
##### MODELLING AND SAVING RESULTS #####
########################################

# Run models
a1 <- run.models(pub.all.dat, stages, div = "raw")
a2 <- run.models(pub.all.dat, stages, div = "raw", na = F)

b1 <- run.models(all.dat, stages, div = "raw")
b2 <- run.models(all.dat, stages, div = "raw", na = F)

c1 <- run.models(m.only.dat, stages, div = "raw")
c2 <- run.models(m.only.dat, stages, div = "raw", na = F)

d1 <- run.models(pb.only.dat, stages, div = "raw")
d2 <- run.models(pb.only.dat, stages, div = "raw", na = F)

NA.model.results <- list("Published" = a1, 
                         "All" = b1, 
                         "Museum" = c1, 
                         "PBDB" = d1)
glob.model.results <- list("Published" = a2,
                             "All" = b2, 
                           "Museum" = c2,
                           "PBDB" = d2)

# Results Viewed
model.results(NA.model.results)
model.results(glob.model.results)
