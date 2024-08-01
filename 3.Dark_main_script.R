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
# 0. INDEX
################################################################################

# 1. SETUP
# 2. BASIC COMPARISONS
# 3. SPATIAL OCCUPANCY THROUGH TIME
# 4. SPATIAL RANGES THROUGH TIME
# 5. TEMPORAL RANGE COMPARISON
# 6. DIVERSITY ANALYSIS
# 7. SEDIMENT AFFINITY
# 8. SCIENTIFIC COLONIALISM
# 9. MODELLING COMPARISONS

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
m.dat$Museum_only[m.dat$Museum_only == 0] <- "Published"
m.dat$Museum_only[m.dat$Museum_only == 1] <- "Dark Data"

chisq.test(table(m.dat$Preservation_score, m.dat$Museum_only))

# compare against what the test would have expected
chisq.test(table(m.dat$Preservation_score, m.dat$Museum_only))$expected

# Mosaic plot
vcd::mosaic(~ Preservation_score + Museum_only,
       direction = c("v", "h"),
       data = m.dat,
       labeling_args = list(set_varnames = c(Preservation_score = "Preservation Score", 
                                             Museum_only = "")),
       shade = TRUE
)

########################################
##### TABLE COMPARISONS AND GRAPHS #####
########################################

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
                    "Rank_simp", 
                    pivot = T, 
                    useNA = F)
names(Rank)[names(Rank) == "Rank_simp"] <- "Rank"

(p1 <- table.plots(Rank, "Rank", colour = c(wes_palette("Zissou1"), 
                                           wes_palette("Royal1"),
                                           wes_palette("Royal2")), 
                   labs = c("A", "B")))

ggsave("Dark_data_graphs/1.Rank.props.png", plot = p1, 
       device = "png", type = "cairo")

##### COUNTRIES #####
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

# Add continents
cc$Continent <- countrycode(sourcevar = cc$Country, origin = "country.name", 
                            destination = "continent")

# Reorder columns
cc <- cc[,c(6,1,2,3,4,5)]

# Sum rows
colSums(cc != 0)

# Write table
write.csv(cc, file = "Dark_data_graphs/Final_Figs/Supp_Figs/SI_Table_1.csv")

##### CONTINENTS #####
Cont <- setup.table(m.only.dat, 
                    pb.only.dat, 
                    pub.all.dat, 
                    all.dat,
                    "Continent", 
                    pivot = T, 
                    useNA = F)
(p2 <- table.plots(Cont, "Continent", colour = c(wes_palette("Royal2"), 
                                                 wes_palette("Royal1")), 
                   labs = c("C", "D")))

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
lith <- tidyr::pivot_longer(lith, cols = c("Dark Data", "PBDB", "Published", "All"))
lith$value <- as.numeric(lith$value)
names(lith)[names(lith) == "Finalised_lith"] <- "Lithology"

# Plotting data
(p3 <- table.plots(lith, "Lithology", removeNA = T, colour = c(wes_palette("Zissou1"), 
                                                               wes_palette("Royal1"),
                                                               wes_palette("Royal2")), 
                   labs = c("E", "F")))

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
                   colour = c(wes_palette("Darjeeling2")), 
                   labs = c("G", "H")))

ggsave("Dark_data_graphs/4.Grain.props.png", plot = p4, 
       device = "png", type = "cairo")

##### COMBINED #####

(f1 <- ggarrange(p1, p2, p3, p4, 
                align='hv',
                nrow = 2, 
                ncol = 2))

ggsave("Dark_data_graphs/Final_Figs/1.Combined.plot.png", plot = f1, 
       device = "png", type = "cairo")

################################################################################
# 3. SPATIAL OCCUPANCY THROUGH TIME
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
pb.cells <- length(unique(pb.only.dat.spat.bin$cells))
m.cells <- length(unique(m.only.dat.space.spat.bin$cells))
pub.cells <- length(unique(pub.all.dat.spat.bin$cells))
all.cells <- length(unique(all.dat.spat.bin$cells)) 

# Cell increase through addition of DD
(all.cells - pub.cells)/pub.cells

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
results$Data[which(results$Data == "m.dat")] <- "Dark Data"
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

test <- results %>%
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
stacked.results[stacked.results == "unique.m"] <- "Dark Data"
stacked.results[stacked.results == "unique.p"] <- "Published record"
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
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  ylab("Ratio of unique \nDark Data/published gridcells") +
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


b <- stacked.results %>%
  dplyr::filter(name == "Ratio") %>%
  ggplot(aes(x=bin_midpoint, y=value)) + 
  geom_line() +
  scale_x_reverse() +
  theme_bw() +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = T) +
  xlab("Time (Ma)") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  ylab("Ratio of unique \nDark Data/published gridcells") +
  guides(fill=guide_legend(title="Preservation score")) +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  guides(alpha = "none")

a <- stacked.results %>%
  dplyr::filter(name != "Ratio") %>%
  ggplot(aes(x=bin_midpoint, y=value, fill=name)) + 
  geom_area() +
  scale_x_reverse() +
  theme_bw() +
  ylab("Number of occupied 1° grid cells") +
  guides(fill=guide_legend(title="Occupied Grid Cells")) +
  scale_fill_manual(values=(rev(wes_palette("Chevalier1")[1:3]))) +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white"), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  guides(alpha = "none")

(p6.5 <- ggarrange(a, b, 
                   align='v',
                   labels = c("A", "B"),
                   nrow = 2,
                   ncol = 1, 
                   heights = c(0.55, 0.7),
                   legend = "bottom", 
                   common.legend = T))

ggsave("Dark_data_graphs/6.5.Occupancy.time.new.png", plot = p6.5, 
       device = "png", type = "cairo")

ggsave("Dark_data_graphs/Final_Figs/Fig_3.png", plot = p6.5, 
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
wide.results$`Dark Data` <- log10(wide.results$`Dark Data`)
wide.results$PBDB <- log10(wide.results$PBDB)
wide.results$Published <- log10(wide.results$Published)
wide.results$echino.cells<- log10(wide.results$echino.cells)
wide.results$`All data` <- log10(wide.results$`All data`)

# Occupancy vs. occupancy
cor1 <- cor.test(wide.results$`Dark Data`, wide.results$PBDB, method = "spearman")
cor2 <- cor.test(wide.results$`Dark Data`, wide.results$Published, method = "spearman")
cor3 <- cor.test(wide.results$`Dark Data`, wide.results$`All data`, method = "spearman")
cor4 <- cor.test(wide.results$PBDB, wide.results$Published, method = "spearman")
cor5 <- cor.test(wide.results$PBDB, wide.results$`All data`, method = "spearman")
cor6 <- cor.test(wide.results$Published, wide.results$`All data`, method = "spearman")

occ.v.occ <- rbind(c("Dark Data", "PBDB", cor1$estimate, cor1$p.value),
                   c("Dark Data", "Published", cor2$estimate, cor2$p.value),
                   c("Dark Data", "All", cor3$estimate, cor3$p.value),
                   c("PBDB", "Published", cor4$estimate, cor4$p.value),
                   c("PBDB", "All", cor5$estimate, cor5$p.value), 
                   c("Published", "All", cor6$estimate, cor6$p.value))

# Occupancy vs. Total invert occ
cor7 <- cor.test(wide.results$`Dark Data`, wide.results$all.cells, method = "spearman")
cor8 <- cor.test(wide.results$PBDB, wide.results$all.cells, method = "spearman")
cor9 <- cor.test(wide.results$Published, wide.results$all.cells, method = "spearman")
cor10 <- cor.test(wide.results$`All data`, wide.results$all.cells, method = "spearman")

occ.v.inv.occ <- rbind(c("Dark Data", "All PBDB", cor7$estimate, cor7$p.value),
                   c("PBDB", "All PBDB", cor8$estimate, cor8$p.value),
                   c("Published", "All PBDB", cor9$estimate, cor9$p.value),
                   c("All", "All PBDB", cor10$estimate, cor10$p.value))

# Occupancy vs. Total echinodermata occ
cor11 <- cor.test(wide.results$`Dark Data`, wide.results$echino.cells, method = "spearman")
cor12 <- cor.test(wide.results$PBDB, wide.results$echino.cells, method = "spearman")
cor13 <- cor.test(wide.results$Published, wide.results$echino.cells, method = "spearman")
cor14 <- cor.test(wide.results$`All data`, wide.results$echino.cells, method = "spearman")

occ.v.e.occ <- rbind(c("Dark Data", "Echinodermata PBDB", cor11$estimate, cor11$p.value),
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
# 4. SPATIAL RANGES THROUGH TIME
################################################################################

##############################
##### LATITUDINAL RANGES #####
##############################

##### SETUP #####
# Stage, genera
lat.results.stage <- lat.range.fun(pbddata = pb.only.dat,
              musdata = m.only.dat, 
              pubdata = pub.all.dat, 
              alldata = all.dat, 
              "stage")
# Series, genera
lat.results.series <- lat.range.fun(pbddata = pb.only.series,
              musdata = m.only.series, 
              pubdata = pub.all.series, 
              alldata = all.series, 
              "series")
# Stage, family
lat.results.stage.fam <- lat.range.fun(pbddata = pb.only.dat,
                                   musdata = m.only.dat, 
                                   pubdata = pub.all.dat, 
                                   alldata = all.dat, 
                                   "stage", 
                                   rank = "Family")
lat.results.stage.fam <- lat.results.stage.fam %>%
  dplyr::filter(Taxon != "NO_FAMILY_SPECIFIED")

##### MAXIMUM AND MINIMUM PALAEOLATITUDES #####

### Select which results to test ###
# Stage, Genera
lat.ranges.stage.Genus$Data[lat.ranges.stage.Genus$Data == "Museum"] <- "Dark Data"
max_min_lat <- pivot_longer(lat.ranges.stage.Genus, cols = c("max_lat", "min_lat"))
# Stage, Family
lat.ranges.stage.Family$Data[lat.ranges.stage.Family$Data == "Museum"] <- "Dark Data"
max_min_lat <- pivot_longer(lat.ranges.stage.Family, cols = c("max_lat", "min_lat"))

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


lat.ranges.stage.Genus %>%
  pivot_longer(-c(taxon_id, range_lat, taxon, bin_assignment, Data), 
               names_to = "lat_group") %>%
  group_by(Data, lat_group) %>%
  dplyr::summarize("avalue" = mean(value),
                   "sd" = sd(value, na.rm = F)) %>%
  dplyr::rename(value = avalue)

##### PALAEOLATITUDINAL SKEW #####

lat.results.stage$Rank <- "Genus"
lat.results.stage.fam$Rank <- "Family"
combined <- rbind(lat.results.stage, lat.results.stage.fam)

# Save
write.csv(combined, file = "Dark_data_graphs/Final_Figs/Supp_Figs/SI_Table_4.csv")

a <- combined %>%
  ggplot(aes(x = Rank, y = Mean_shift, fill = Data)) +
  guides(alpha = "none") +
  geom_boxplot(aes(alpha = 0.1), show.legend = NA) +
  scale_fill_viridis(discrete = T) + 
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  guides(alpha = "none") +
  ylab("Shift in range mean (°)") 
  
b <- combined %>%
  mutate(text = forcats::fct_reorder(Data, Mean_shift)) %>%
  ggplot(aes(x=Mean_shift, color = Data, fill=Data)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', bins = 15) +
  scale_fill_viridis(discrete = T) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  guides(alpha = "none") +
  xlab("Shift in range mean (°)") + 
  ylab("Total") +
  facet_wrap(~Rank + Data)

(p8 <- ggarrange(a, b, 
          align='hv',
          labels = c("A", "B"),
          nrow = 1,
          ncol = 2, 
          common.legend = T, 
          legend = "bottom"))

ggsave("Dark_data_graphs/8.Lat.shift.png", plot = p8, 
       device = "png", type = "cairo")

ggsave("Dark_data_graphs/Final_Figs/X.Lat.shift.png", plot = p8, 
       device = "png", type = "cairo")

test <- combined %>%
  group_by(Data, Rank) %>%
  dplyr::summarize(Mean.shift.mean = mean(Mean_shift), 
                   Median.shift.mean = median(Mean_shift))

##### DIFFERENCES BETWEEN DATASETS #####
# Find results
space.results(lat.results.series)
space.results.gen <- space.results(lat.results.stage)
space.results.fam <- space.results(lat.results.stage.fam)
space.results.gen$Rank <- "Genus"
space.results.fam$Rank <- "Family"
comb.space.results <- rbind(space.results.gen, space.results.fam)
write.csv(comb.space.results, file = "Dark_data_graphs/Final_Figs/Supp_Figs/SI_Table_5.csv")

# Mean/median addition of latitudinal range from museums, as percentage of total
lat.results.stage %>%
  dplyr::filter(Data == "Dark Data" | Data == "Published") %>%
  dplyr::group_by(bin, Taxon) %>%
  dplyr::filter(n() >1) %>%
  dplyr::filter(Data == "Published") %>%
  mutate(Perc.new = 100 - Perc) %>%
  pull(Perc.new) %>%
  mean() # change this to get median

# Mean/median addition of latitudinal range from museums, as percentage of total
lat.results.stage.fam %>%
  dplyr::filter(Data == "Dark Data" | Data == "Published") %>%
  dplyr::group_by(bin, Taxon) %>%
  dplyr::filter(n() >1) %>%
  dplyr::filter(Data == "Published") %>%
  mutate(Perc.new = 100 - Perc) %>%
  pull(Perc.new) %>%
  mean() # change this to get median

# Mean/median addition of latitudinal range from museums in degrees
lat.results.stage %>%
  dplyr::filter(Data == "Dark Data" | Data == "Published") %>%
  dplyr::group_by(bin, Taxon) %>%
  dplyr::filter(n() >1) %>%
  dplyr::filter(Data == "Published") %>%
  pull(Difference) %>%
  mean()

# Mean/median addition of latitudinal range from museums in degrees
lat.results.stage.fam %>%
  dplyr::filter(Data == "Dark Data" | Data == "Published") %>%
  dplyr::group_by(bin, Taxon) %>%
  dplyr::filter(n() >1) %>%
  dplyr::filter(Data == "Published") %>%
  pull(Difference) %>%
  mean()

# Plot results
(p9 <- combined %>%
  ggplot(aes(x = Rank, y = Difference, fill = Data)) +
  geom_boxplot() +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
    theme(legend.position = "none") +
  guides(alpha = "none") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  scale_colour_viridis_d() +
    geom_jitter(shape=16, position=position_jitterdodge(jitter.width = 0.2), aes(color = Data), alpha = 0.8) +
  ylab("Difference between full latitudinal range and data (°)") +
  xlab(""))

(p9 <- ggarrange(a, b, 
          align='hv',
          labels = c("A", "B"),
          nrow = 2,
          ncol = 1))

ggsave("Dark_data_graphs/9.Lat.diff.png", plot = p9, 
       device = "png", type = "cairo")

ggsave("Dark_data_graphs/Final_Figs/X.Lat.diff.png", plot = p9, 
       device = "png", type = "cairo")

(p8 <- ggarrange(p9, ggarrange(a, b, ncol = 2, 
                 labels = c("C", "D"), common.legend = T, 
                 legend = "bottom"), labels = "B",
                 nrow = 2
                 ))

ggsave("Dark_data_graphs/Final_Figs/X.Lat.png", plot = p8, 
       device = "png", type = "cairo")

# Kruskal wallace test for differences between medians
spatial.tests(lat.results.stage, "Difference")
spatial.tests(lat.results.stage, "Mean_shift")
spatial.tests(lat.results.stage, "Perc")
spatial.tests(lat.ranges.stage.Genus, "range_lat")
spatial.tests(lat.ranges.stage.Family, "range_lat")

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

geo.results.stage.fam <- geo.range.fun(pbddata = pb.only.dat,
                                   musdata = m.only.dat,
                                   pubdata = pub.all.dat,
                                   alldata = all.dat,
                                   "stage", rank = "Family")
geo.results.stage.fam <- geo.results.stage.fam %>%
  dplyr::filter(Taxon != "NO_FAMILY_SPECIFIED")

# Find differences
space.results(geo.results.stage)
space.results(geo.results.series)
space.results(geo.results.stage.fam)

geo.results.stage$Rank <- "Genus"
geo.results.stage.fam$Rank <- "Family"

geo.results.combined <- rbind(geo.results.stage, geo.results.stage.fam)

# Save
write.csv(geo.results.combined, file = "Dark_data_graphs/Final_Figs/Supp_Figs/SI_Table_2.csv")

# Plot differences at generic level
a <- geo.results.combined %>%
  ggplot(aes(x = Rank, y = Perc, fill = Data)) +
  geom_boxplot(alpha = 0.4) +
  geom_jitter(shape=16, position=position_jitterdodge(jitter.width = 0.2), aes(color = Data), alpha = 0.8) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  guides(alpha = "none") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  scale_colour_viridis_d() +
  #geom_jitter(aes(color=Data), size=0.4, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ylab("Percentage of range overlap") +
  xlab("")

geo.results.combined$Percentage.added <- as.numeric(geo.results.combined$Percentage.added)

b <- geo.results.combined %>%
    ggplot(aes(x = Rank, y = Percentage.added, fill = Data)) +
    geom_boxplot(alpha = 0.4) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "white"), 
          panel.grid.minor = element_line(colour = "white")) +
    guides(alpha = "none") +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  scale_colour_viridis_d() +
  geom_jitter(shape=16, position=position_jitterdodge(jitter.width = 0.2), aes(color = Data), alpha = 0.8) +
    theme(
      plot.title = element_text(size=11)
    ) +
    ylab(expression('Area added to geographic range (% of total area from convex hull)')) +
    xlab("")

(p10 <- ggarrange(a, b, 
                  align='hv',
                  labels = c("A", "B"),
                  legend = "bottom",
                  common.legend = T,
                  nrow = 1))

ggsave("Dark_data_graphs/10.Range.overlap.added.area.genera.png", plot = p10, 
       device = "png")

ggsave("Dark_data_graphs/Final_Figs/x.Geographic.range.png", plot = p10, 
       device = "png")

# Mean added area (percentage of total area)
a <- geo.results.stage %>%
  dplyr::group_by(Data) %>%
  dplyr::summarize(`Mean added area` = mean(Percentage.added),
            `Median added area` = median(Percentage.added)) %>%
  dplyr::mutate(Rank = "Genus")

geo.results.series %>%
  dplyr::group_by(Data) %>%
  dplyr::summarize(Mean = mean(Percentage.added),
                   Median = median(Percentage.added))

b <- geo.results.stage.fam %>%
  dplyr::group_by(Data) %>%
  dplyr::summarize(`Mean added area` = mean(Percentage.added),
                   `Median added area` = median(Percentage.added)) %>%
  dplyr::mutate(Rank = "Family")

SI_Table_3 <- rbind(a, b)
write.csv(SI_Table_3, file = "Dark_data_graphs/Final_Figs/Supp_Figs/SI_Table_3.csv")

# Plot of added area per data type
geo.results.stage$Percentage.added <- as.numeric(geo.results.stage$Percentage.added)
(p10.51 <- geo.results.stage %>%
  ggplot(aes(Percentage.added, fill = Data)) + 
  geom_histogram() +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")) +
  guides(alpha = "none") +
  xlab("Percentage of area added") +
  ylab("Frequency") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6))

ggsave("Dark_data_graphs/10.5.Added.area.png", plot = p10.51, 
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

spatial.tests(geo.results.stage.fam, "Area")
spatial.tests(geo.results.stage.fam, "Difference")
spatial.tests(geo.results.stage.fam, "Perc")
spatial.tests(geo.results.stage.fam, "Added.area")

################################################################################
# 5. TEMPORAL RANGE COMPARISON
################################################################################

################################
##### TOTAL TEMPORAL RANGE #####
################################

# Run on each dataset
temp_compare(m.only.dat, bin_filter = 3)
temp_compare(pb.only.dat, bin_filter = 3)
temp_compare(pub.all.dat, bin_filter = 3)

# Combined
comp <- rbind(comp.m.only.dat, comp.pb.only.dat, comp.pub.all.dat)
(comp.stats <- rbind(comp.stats.m.only.dat, 
                     comp.stats.pb.only.dat, 
                     comp.stats.pub.all.dat))
comp.stats$dataframe[comp.stats$dataframe == "m.only.dat"] <- "Dark Data"
comp.stats$dataframe[comp.stats$dataframe == "pb.only.dat"] <- "PBDB"
comp.stats$dataframe[comp.stats$dataframe == "pub.all.dat"] <- "Published"
names(comp.stats)[names(comp.stats) == "dataframe"] <- "Data"
names(comp.stats)[names(comp.stats) == "Mean_range_diff"] <- "Mean Range Difference"
names(comp.stats)[names(comp.stats) == "Mean_top_diff"] <- "Mean Top Difference"
names(comp.stats)[names(comp.stats) == "Mean_bot_diff"] <- "Mean Bottom Difference"
names(comp.stats)[names(comp.stats) == "Range_perc"] <- "Percentage of total range"
names(comp.stats)[names(comp.stats) == "Correct_perc"] <- "Percentage of taxa with correct range"

write.csv(comp.stats, file = "Dark_data_graphs/Final_Figs/Supp_Figs/SI_Table_7.csv")

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
comp$dataframe[comp$dataframe == "m.only.dat"] <- "Dark Data"
comp$dataframe[comp$dataframe == "pb.only.dat"] <- "PBDB"
comp$dataframe[comp$dataframe == "pub.all.dat"] <- "Published"
names(comp)[names(comp) == "dataframe"] <- "Data"
names(comp)[names(comp) == "range_diff"] <- "Difference"

write.csv(comp, file = "Dark_data_graphs/Final_Figs/Supp_Figs/SI_Table_6.csv")

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
    ggplot(aes(x=Group, y=Range_perc, fill = Data)) +
    geom_boxplot() +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "white"), 
          panel.grid.minor = element_line(colour = "white"), 
          axis.title.x=element_blank()) +
    guides(alpha = "none") +
    scale_fill_viridis(discrete = T, alpha = 0.8) +
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

#######################################
##### TEMPORAL RANGE DISTRIBUTION #####
#######################################

# Get position within temporal range of each Genera
m.pos.gen <- perc.range(occdf = m.only.dat, label = "Museum", rank = "Genus")
pub.pos.gen <- perc.range(occdf = pub.all.dat, label = "Published", rank = "Genus")
pb.pos.gen <- perc.range(occdf = pb.only.dat, label = "PBDB", rank = "Genus")

# Get position within temporal range of each Genera
m.pos.sp <- perc.range(occdf = m.only.dat, label = "Museum", rank = "Species")
pub.pos.sp <- perc.range(occdf = pub.all.dat, label = "Published", rank = "Species")
pb.pos.sp <- perc.range(occdf = pb.only.dat, label = "PBDB", rank = "Species")

# Get position within temporal range of each Genera
m.pos.fam <- perc.range(occdf = m.only.dat, label = "Museum", rank = "Family")
pub.pos.fam <- perc.range(occdf = pub.all.dat, label = "Published", rank = "Family")
pb.pos.fam <- perc.range(occdf = pb.only.dat, label = "PBDB", rank = "Family")

# Combine dataset
pos <- rbind(m.pos.gen, pub.pos.gen, pb.pos.gen, m.pos.fam, pub.pos.fam, pb.pos.fam, 
             m.pos.sp, pub.pos.sp, pb.pos.sp)

pos$Rank <- as.factor(pos$Rank)
pos$Data[pos$Data == "Museum"] <- "Dark Data"

pos <- pos %>%
  dplyr::filter(!(Rank == "Species" & Data == "PBDB"))

# Plot
(p11.7 <- ggplot(pos, aes(x = Rank, y = Position, group = interaction(Rank, Data))) +
  geom_violin(width = 1, position = position_dodge(1), aes(fill = Data), scale = "area") +
  scale_fill_viridis(discrete = T, alpha = 0.8) + 
  scale_colour_viridis_d() +
  ylab("Position of occurrence within range (%)") +
  geom_jitter(shape=16, position=position_jitterdodge(dodge.width = 1, 
                                                      jitter.width = 0.25), aes(color = Data), alpha = 0.5) +
  geom_boxplot(width = 0.3, position = position_dodge(1), alpha = 0.8) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white")))

# Save
ggsave("Dark_data_graphs/11.7.Temporal.range.position.png", plot = p11.7, 
       device = "png", type = "cairo")

(p11 <- ggarrange(p11.5, p11.7, 
          align = 'hv',
          labels = c("A", "B"),
          legend = "bottom",
          nrow = 2, 
          common.legend = T))

ggsave("Dark_data_graphs/Final_Figs/X.Temporal.range.png", plot = p11, 
       device = "png", type = "cairo")

################################################################################
# 6. DIVERSITY ANALYSIS
################################################################################

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
div.tot.long$Data[div.tot.long$Data == "Museum_SIB"] <- "Dark Data"
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
cor.results$`Variable 1`[cor.results$`Variable 1` == "Museum"] <- "Dark Data"

cor.results %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  save_kable(file = "Dark_data_graphs/Final_Figs/Supp_Figs/SI_Table_8.html")

write.csv(cor.results, file = "Dark_data_graphs/Final_Figs/Supp_Figs/SI_Table_8.csv")

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
estD_plotting$Data[estD_plotting$Data == "Museum"] <- "Dark Data"
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
  ylim(0, 7) +
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
  facet_wrap(~Data, nrow = 4) +
  guides(shape = "none"))

ggsave("Dark_data_graphs/14.SQS.div.iNEXT.png", plot = p14, 
       device = "png", type = "cairo")
ggsave("Dark_data_graphs/Final_Figs/X.SQS.png", plot = p14, 
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
# 7. SEDIMENT AFFINITY
################################################################################

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

combined.SRA$Data[combined.SRA$Data == "Museum"] <- "Dark Data"

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

ggsave("Dark_data_graphs/Final_Figs/X.SRA.png", plot = p16, 
       device = "png", type = "cairo")

ggsave("Dark_data_graphs/16.5.SRA.png", plot = p16.5, 
       device = "png", type = "cairo")


(p.comb <- ggarrange(p14, p16,
          labels = c("A", "B"), 
          widths = c(0.7, 1),
          nrow = 1, 
          ncol = 2))

ggsave("Dark_data_graphs/Final_Figs/SQS.SRA.png", plot = p.comb, 
       device = "png", type = "cairo")

################################################################################
# 8. SCIENTIFIC COLONIALISM
################################################################################

# For this section, code was reworked from that provided for the following publication:

# Raja, N.B., Dunne, E.M., Matiwane, A., Khan, T.M., Nätscher, P.S., Ghilardi, A.M. 
# and Chattopadhyay, D., 2022. Colonial history and global economics distort our 
# understanding of deep-time biodiversity. Nature ecology & evolution, 6(2), pp.145-154.

# We thank the authors for providing this code and for enabling this research. 

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

All.cont %>%
  dplyr::group_by(Data, Origin) %>%
  dplyr::filter(Continent == "North America" | Continent == "Europe") %>%
  dplyr::summarise(Total_freq = sum(Frequency))

Research <- bind_rows(list(research.prep(a.df, data = "All"), 
                    research.prep(p.df, data = "Published"),
                    research.prep(m.df, data = "Museum")))

Research$Data[Research$Data == "Museum"] <- "Dark Data"

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

ggsave("Dark_data_graphs/Final_Figs/2.Contribution.png", plot = p17, 
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

(p19 <- network.stats(m.df))
(p20 <- network.stats(p.df))
(p21 <- network.stats(a.df))

ggsave("Dark_data_graphs/19.DarkData.network.stats.png", plot = p19, 
       device = "png", type = "cairo")
ggsave("Dark_data_graphs/20.Published.network.stats.png", plot = p20, 
       device = "png", type = "cairo")
ggsave("Dark_data_graphs/21.All.data.network.stats.png", plot = p21, 
       device = "png", type = "cairo")

################################################################################
# 9. MODELLING COMPARISONS
################################################################################

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
sea.lvl <- read.csv("Additional_data/vanderMeer_2022.csv")
sea.lvl$max_ma <- sea.lvl$Ma+0.0001
names(sea.lvl)[names(sea.lvl) == "Ma"] <- "min_ma"
sea.lvl <- sea.lvl %>%
  dplyr::filter(max_ma < 538) %>%
  dplyr::filter(min_ma > 0)
sea.lvl <- bin_time(occdf = sea.lvl, bins = stages, method = "mid")
sea.lvl <- sea.lvl %>%
  dplyr::group_by(bin_midpoint) %>%
  dplyr::summarize(mean_sl = mean(TGE_SL_isocorr_m))

# Macrostrat
macro <- as.data.frame(pivot_wider(macro.area, names_from = lith, 
                                   values_from = count))

# MgCa ratios
ca1 <- read.csv("Additional_data/Horita_data_extracted.csv")
ca2 <- read.csv("Additional_data/Hardie_data_extracted.csv")
ca3 <- read.csv("Additional_data/Wilkinson+Algeo_data_extracted.csv")

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
                         "Dark Data" = c1, 
                         "PBDB" = d1)
glob.model.results <- list("Published" = a2,
                             "All" = b2, 
                           "Dark Data" = c2,
                           "PBDB" = d2)
# Results Viewed
model.results.1 <- model.results(NA.model.results)
model.results.2 <-model.results(glob.model.results)

pub.mod <- paste(names(NA.model.results$Published$Top.model$coefficients), collapse = " + ")
all.mod <- paste(names(NA.model.results$All$Top.model$coefficients), collapse = " + ")
DD.mod <- paste(names(NA.model.results$`Dark Data`$Top.model$coefficients), collapse = " + ")
PBDB.mod <- paste(names(NA.model.results$PBDB$Top.model$coefficients), collapse = " + ")

NA.table <- model.results.1[[1]][,c(14,9,10,11,13)]
NA.table$`Top model` <- c(pub.mod, all.mod, DD.mod, PBDB.mod)
NA.table$Model <- "North America only"
NA.table <- NA.table[,c(7, 1, 6, 2, 3, 4)]
NA.table <- NA.table[order(NA.table$Data), ]

pub.mod <- paste(names(glob.model.results$Published$Top.model$coefficients), collapse = " + ")
all.mod <- paste(names(glob.model.results$All$Top.model$coefficients), collapse = " + ")
DD.mod <- paste(names(glob.model.results$`Dark Data`$Top.model$coefficients), collapse = " + ")
PBDB.mod <- paste(names(glob.model.results$PBDB$Top.model$coefficients), collapse = " + ")

glob.table <- model.results.2[[1]][,c(12,7,8,9)]
glob.table$`Top model` <- c(pub.mod, all.mod, DD.mod, PBDB.mod)
glob.table$Model <- "Global"
glob.table <- glob.table[,c(6, 1, 5, 2, 3, 4)]
glob.table <- glob.table[order(glob.table$Data), ]

Table1<- rbind(NA.table, glob.table)

Table1 %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  save_kable(file = "Dark_data_graphs/Final_Figs/Table_1.html")

write.csv(Table1, file = "Dark_data_graphs/Final_Figs/Table_1.csv")

model.results.1[[2]]$Model <- "North America only"
model.results.2[[2]]$Model <- "Global"

Table2 <- rbind(model.results.1[[2]], model.results.2[[2]])
Table2 <- Table2[,c(8, 1, 2, 3, 4, 5, 6, 7)]

Table2 %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  save_kable(file = "Dark_data_graphs/Final_Figs/Table_2.html")

write.csv(Table2, file = "Dark_data_graphs/Final_Figs/Table_2.csv")
