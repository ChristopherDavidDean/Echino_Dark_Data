################################################################################
######## TAPHONOMIC CONTROLS ON A MULTI-ELEMENT SKELETAL FOSSIL RECORD #########
################################################################################

# Jeffrey R. Thompson, Christopher D. Dean, Madeline Ford, Timothy A. M. Ewin
# 2023
# Script written by Christopher D. Dean

################################################################################
#                            FILE 4: MAIN FIGURES                              #
################################################################################

################################################################################
# 1. SETUP
################################################################################

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Load setup file
source("2.Setup.R")

################################################################################
# 2. FIGURE 1 - INTRO FIGURE
################################################################################

# Maps
# Set extent
e <- extent(-180, 180, -90, 90)
res <- 2

rasterise.dat <- function(data, res, ext){
  xy <- cbind(as.double(data$lng), as.double(data$lat))
  #xy <- unique(xy)
  r <- raster::raster(ext = ext, res = res)
  r <- raster::rasterize(xy, r, fun = 'count')
  return(r)
}

stacked <- stack(rasterise.dat(dplyr::filter(m.dat, Preservation_score == 1), 2, 
                               ext = e), 
      rasterise.dat(dplyr::filter(m.dat, Preservation_score == 5), 2, ext = e))
raster.names <- c("TG1", "TG5")

# find map to use as backdrop
countries <- maps::map("world", plot=FALSE, fill = TRUE) 
# Turn map into spatialpolygons
countries <- maptools::map2SpatialPolygons(countries, 
                                           IDs = countries$names, 
                                           proj4string = CRS("+proj=longlat")) 
# Customize the colorkey
my.theme <- rasterVis::rasterTheme(region=viridis(8))
a <- rasterVis::levelplot(stacked, layout = c(1, 2), margin=T, 
                          par.settings=my.theme,
                          names.attr = raster.names) + 
  # Plots background colour
  latticeExtra::layer(sp.polygons(countries, col = 0, fill = "light grey"), 
                      under = T)

b <- ggplot(f.m.dat) +
  aes(x = Preservation_score, fill = Family) +
  geom_bar() +
  ylab("Frequency") +
  xlab("Taphonomic Grade (TG)") +
  scale_fill_manual(values=c(wes_palette("Royal1"), 
                             wes_palette("Zissou1"), 
                             wes_palette("Royal2"))) +
  theme_bw() + 
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white"), 
        legend.position = "none") +
  guides(alpha = "none")


# Lithology
test <- as.data.frame(table(l.m.dat$Finalised_lith, l.m.dat$bin_midpoint))
names(test) <- c("Lithology", "mid_ma", "Freq")
test$mid_ma <- as.numeric(as.character(test$mid_ma))
c <- ggplot(test, aes(x=mid_ma, y=Freq, fill=Lithology)) + 
  geom_area() +
  scale_x_reverse() +
  theme_bw() +
  guides(fill=guide_legend(title="Lithology")) +
  ylab("Frequency") +
  xlab("Time (Ma)") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white"), 
        legend.position = "bottom") +
  guides(alpha = "none") +
  scale_fill_manual(values=(wes_palette("Zissou1"))) +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = list(T))

# Grain size
test <- as.data.frame(table(g.m.dat$Finalised_grainsize_simp, 
                            g.m.dat$bin_midpoint))
names(test) <- c("Finalised_grainsize", "mid_ma", "Freq")
test$mid_ma <- as.numeric(as.character(test$mid_ma))
d <- ggplot(test, aes(x=mid_ma, y=Freq, fill=Finalised_grainsize)) + 
  geom_area() +
  scale_x_reverse() +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white"), 
        legend.position = "bottom") +
  guides(alpha = "none") +
  guides(fill=guide_legend(title="Grain size")) +
  ylab("Frequency") +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Darjeeling2"))) +
  coord_geo(dat = series2, 
            xlim = c(485.4, 251.902),
            rot = list(0),
            size = 4,
            pos = list("b"),
            abbrv = list(T))

p1 <- ggarrange(a, b, c, d, 
                align='hv',
                labels = c("B", "C", "D", "E"),
                nrow = 2, 
                ncol = 2)

################################################################################
# 3. FIGURE 2 - GRAIN SIZE/LITHOLOGY
################################################################################

# Make bar plot
a <- make.bar.plot(l.m.dat, l.m.dat$Finalised_lith, "Lithology", 
                   colour = "Zissou1", FALSE)

# Make bar plot
b <- make.bar.plot(g.m.dat, g.m.dat$Finalised_grainsize_simp, "Grain Size", 
                   colour = "Darjeeling2", flip = FALSE)

ggarrange(a, b, 
          align='hv',
          labels = c("A", "B"),
          nrow = 1, 
          ncol = 2)

# Mosaic plot
vcd::mosaic(~ Preservation_score + Finalised_lith,
       direction = c("v", "h"),
       data = l.m.dat,
       labeling_args = list(
         set_varnames = c(Preservation_score = "Preservation Score", 
                          Finalised_lith = "Lithology")),
       shade = TRUE
)

# Mosaic plot
vcd::mosaic(~ Preservation_score + Finalised_grainsize_simp,
       direction = c("v", "h"),
       data = g.m.dat,
       labeling_args = list(just_labels = "right", 
                            set_varnames = c(Preservation_score = "Preservation Score", 
                                             Finalised_grainsize_simp = "Grain size")),
       shade = TRUE
)


################################################################################
# 4. FIGURE 3 - TAXONOMIC RANK AND PRESERVATION SCORE THROUGH TIME
################################################################################

a <- make.bar.plot(m.dat, m.dat$Rank, "Taxonomic Rank", colour = "Zissou1", FALSE)

b <- ggplot(m.dat.series, aes(x = Preservation_score, fill = bin_assignment)) +
  geom_bar() +
  scale_fill_manual("Legend", values = myColours) +
  ylab("Frequency") +
  xlab("Preservation Score") +
  guides(fill=guide_legend(title="Period")) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white"), 
        legend.position = "bottom") +
  guides(alpha = "none")

ggarrange(a, b, 
          align='hv',
          labels = c("A", "B"),
          nrow = 1, 
          ncol = 2)

# Mosaic plot
vcd::mosaic(~ Preservation_score + Rank,
       direction = c("v", "h"),
       data = m.dat,
       labeling_args = list(
         set_varnames = c(Preservation_score = "Preservation Score", 
                          Rank = "Taxonomic Rank")),
       shade = TRUE
)

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
# 5. FIGURE 4 - LOGISTIC REGRESSION
################################################################################

library(ggforestplot)

ctable <- read.csv("Results/ctable.csv")

ctable$PresScore <- as.factor(ctable$PresScore)
ctable <- ctable %>%
  filter(!(`Std..Error` > 10))

ctable[ctable == "`Age (Ma)`"] <- "Age (Ma)"
ctable[ctable == "`Palaeo-latitude`"] <- "Palaeo-latitude"
ctable[ctable == "(Intercept)"] <- "Intercept"
ctable[ctable == "GrainsizeFine Grained"] <- "Grainsize: Fine Grained"
ctable[ctable == "LithologySiliciclastic"] <- "Lithology: Siliciclastic"
ctable[ctable == "GrainsizeFine Grained:LithologySiliciclastic"] <- 
  "Grainsize x Lithology: Fine Grained x Siliciclastic"
 
ctable <- ctable[order(ctable$Covariate),]
ctable <- ctable[order(ctable$PresScore),]

forestplot(
  df = ctable,
  name = Covariate,
  estimate = Estimate,
  se = `Std..Error`,
  title = "",
  colour = PresScore,
  pvalue = `Pr...z..`,
  psignif = 0.05
) + ggplot2::facet_wrap(
  facets = ~PresScore,
  nrow = 5, 
  ncol = 1,
  scales = "free_y"
) + theme(strip.text.x = element_blank(), 
          panel.border = element_rect(colour = "black", 
                                      fill=NA, linewidth=0.5)) + 
  scale_colour_manual(values = wesanderson::wes_palette("Zissou1", 
                                                        type = 'discrete'))

