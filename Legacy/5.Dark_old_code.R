################################################################################
### MUSEUM `DARK DATA` ILLUMINATES THE EVOLUTIONARY HISTORY OF FOSSIL GROUPS ###
################################################################################

# Christopher D. Dean & Jeffrey R. Thompson
# 2024
# Script written by Christopher D. Dean

################################################################################
#                              FILE 5: SOLD CODE                               #
################################################################################


#============================= PRESENCE.ABSENCE.RASTER ===================================

# Function from https://amywhiteheadresearch.wordpress.com/2013/05/27/creating-a-presence-absence-raster-from-point-data/
# Credit to Amy Whitehead

presence.absence.raster <- function (mask.raster,species.data,raster.label="") {
  require(raster)
  
  # set the background cells in the raster to 0
  mask.raster[!is.na(mask.raster)] <- 0
  
  #set the cells that contain points to 1
  speciesRaster <- rasterize(species.data,mask.raster,field=1)
  speciesRaster <- merge(speciesRaster,mask.raster)
  
  #label the raster
  names(speciesRaster) <- raster.label
  return(speciesRaster)
}

#================================== OCCUPANCY_EST ========================================

occupancy_est <- function(occdf, totalfoss, res, e, bins){
  
  total_occupancy <- data.frame(stringsAsFactors = FALSE)
  all_total_cells <- c()
  
  for (t in rev(unique(occdf$bin))){ # For each bin
    temp_data <- totalfoss %>%
      dplyr::filter(bin_no == t) 
    grid_data <- get_grid(temp_data, res, e) # Make raster and grid up cells
    total_cells <- length(unique(grid_data$cells))
    all_total_cells <- c(all_total_cells, total_cells)
    
    temp_stack <- stack()
    
    if (level == "family"){
      family_occupancy <- data.frame(stringsAsFactors = FALSE)
      for (f in names){
        temp_occ_data <- grid_data %>%
          dplyr::filter(family_name == f)
        temp_vec <- c(t, f, length(unique(temp_occ_data$cells)), length(unique(temp_occ_data$cells))/total_cells*100)
        family_occupancy <- rbind(family_occupancy, temp_vec, stringsAsFactors = FALSE)
        
        r <- raster(res = res, val = 1, ext = e)
        xy <- cbind(as.double(temp_occ_data$lng), as.double(temp_occ_data$lat))  
        occ_raster <- presence.absence.raster(mask.raster = r, species.data = xy, raster.label = f)
        temp_stack <- stack(temp_stack, occ_raster)
      }
      colnames(family_occupancy) <- c("Bin", "Total_cells", "Occupied_cells", "Perc_Occupied_cells")
      total_family_occupancy <- rbind(total_family_occupancy, family_occupancy, stringsAsFactors = FALSE)
    }
    stack_name <- paste("Bin", t, "Stack", sep = "_")
    assign(stack_name, temp_stack, envir = .GlobalEnv)
  }
  if (level == "family"){
    reshaped <- total_family_occupancy[,-3]
    perc_fam_occ <- reshape(reshaped, idvar = "Bin", timevar = "Family", direction = "wide")
    
    reshaped <- total_family_occupancy[,-4]
    total_fam_occ <- reshape(reshaped, idvar = "Bin", timevar = "Family", direction = "wide")
    
    f.names <- colnames(perc_fam_occ)[2:length(perc_fam_occ)]
    f.names <- c("Bin", gsub(".*\\.","",f.names))
    colnames(perc_fam_occ) <- f.names
    colnames(total_fam_occ) <- f.names
    perc_fam_occ <<- perc_fam_occ
    total_fam_occ <<- total_fam_occ
    
    total_family_occupancy$Bin <- as.numeric(total_family_occupancy$Bin)
    total_family_occupancy$Perc_Occupied_cells <- as.numeric(total_family_occupancy$Perc_Occupied_cells)
    total_family_occupancy <<- total_family_occupancy
  }
  all_total_cells <<- all_total_cells
}


#================================== OCCUPANCY_EST ========================================

occupancy_est <- function(fossils, res, bins, formations, level = "genus"){
  e <- get_extent(fossils) # Get extent
  binned_occs <- FormBin_M1(fossils, formations, bins) # Bin occurrences
  if(level == "family"){
    total_family_occupancy <- data.frame(stringsAsFactors = FALSE)
  }
  if(level == "genus"){
    total_genus_occupancy <- data.frame(stringsAsFactors = FALSE)
  }
  
  all_total_cells <- c()
  
  for (t in rev(unique(binned_occs$bin_no))){ # For each bin
    temp_data <- binned_occs %>%
      dplyr::filter(bin_no == t) 
    if(level == "family"){
      names <- unique(temp_data$family_name)
    }
    if (level == "genus"){
      names <- unique(temp_data$occurrence.genus_name)
    }
    
    grid_data <- get_grid(temp_data, res, e) # Make raster and grid up cells
    total_cells <- length(unique(grid_data$cells))
    all_total_cells <- c(all_total_cells, total_cells)
    
    temp_stack <- stack()
    
    if (level == "family"){
      family_occupancy <- data.frame(stringsAsFactors = FALSE)
      for (f in names){
        temp_occ_data <- grid_data %>%
          dplyr::filter(family_name == f)
        temp_vec <- c(t, f, length(unique(temp_occ_data$cells)), length(unique(temp_occ_data$cells))/total_cells*100)
        family_occupancy <- rbind(family_occupancy, temp_vec, stringsAsFactors = FALSE)
        
        r <- raster(res = res, val = 1, ext = e)
        xy <- cbind(as.double(temp_occ_data$lng), as.double(temp_occ_data$lat))  
        occ_raster <- presence.absence.raster(mask.raster = r, species.data = xy, raster.label = f)
        temp_stack <- stack(temp_stack, occ_raster)
      }
      colnames(family_occupancy) <- c("Bin", "Family", "Occupied_cells", "Perc_Occupied_cells")
      total_family_occupancy <- rbind(total_family_occupancy, family_occupancy, stringsAsFactors = FALSE)
    }
    if (level == "genus"){
      genus_occupancy <- data.frame(stringsAsFactors = FALSE)
      for (f in names){
        temp_occ_data <- grid_data %>%
          dplyr::filter(occurrence.genus_name == f)
        temp_vec <- c(t, as.character(temp_occ_data$family_name[[1]]), f, length(unique(temp_occ_data$cells)), length(unique(temp_occ_data$cells))/total_cells*100)
        genus_occupancy <- rbind(genus_occupancy, temp_vec, stringsAsFactors = FALSE)
        
        r <- raster(res = res, val = 1, ext = e)
        xy <- cbind(as.double(temp_occ_data$lng), as.double(temp_occ_data$lat))  
        occ_raster <- presence.absence.raster(mask.raster = r, species.data = xy, raster.label = f)
        temp_stack <- stack(temp_stack, occ_raster)
      }
      colnames(genus_occupancy) <- c("Bin", "Family", "Genus", "Occupied_cells", "Perc_Occupied_cells")
      total_genus_occupancy <- rbind(total_genus_occupancy, genus_occupancy, stringsAsFactors = FALSE)
    }
    stack_name <- paste("Bin", t, "Stack", sep = "_")
    assign(stack_name, temp_stack, envir = .GlobalEnv)
  }
  if (level == "family"){
    reshaped <- total_family_occupancy[,-3]
    perc_fam_occ <- reshape(reshaped, idvar = "Bin", timevar = "Family", direction = "wide")
    
    reshaped <- total_family_occupancy[,-4]
    total_fam_occ <- reshape(reshaped, idvar = "Bin", timevar = "Family", direction = "wide")
    
    f.names <- colnames(perc_fam_occ)[2:length(perc_fam_occ)]
    f.names <- c("Bin", gsub(".*\\.","",f.names))
    colnames(perc_fam_occ) <- f.names
    colnames(total_fam_occ) <- f.names
    perc_fam_occ <<- perc_fam_occ
    total_fam_occ <<- total_fam_occ
    
    
    total_family_occupancy$Bin <- as.numeric(total_family_occupancy$Bin)
    total_family_occupancy$Perc_Occupied_cells <- as.numeric(total_family_occupancy$Perc_Occupied_cells)
    total_family_occupancy <<- total_family_occupancy
    
  }
  if (level == "genus"){
    names <- unique(total_genus_occupancy$Family)
    family_grouped_occupancy <- list()
    for (f in names){
      temp_fam_split <- total_genus_occupancy %>%
        dplyr::filter(Family == f) %>%
        dplyr::select(Bin, Genus, Perc_Occupied_cells) 
      temp_fam_split <- reshape(temp_fam_split, idvar = "Bin", timevar = "Genus", direction = "wide")
      temp_fam_split$Bin <- as.numeric(temp_fam_split$Bin)
      if (nrow(temp_fam_split) < nrow(bins)){
        temp_fam_split <- temp_fam_split %>%
          tidyr::complete(Bin = seq(1, nrow(bins), 1)) %>%
          purrr::map_df(rev)
      }
      g.names <- colnames(temp_fam_split)[2:length(temp_fam_split)]
      g.names <- c("Bin", gsub(".*\\.","",g.names))
      colnames(temp_fam_split) <- g.names
      family_grouped_occupancy[[f]] <- temp_fam_split
    }
    names(family_grouped_occupancy) <- names
    family_grouped_occupancy <<- family_grouped_occupancy
    
    total_genus_occupancy$Bin <- as.numeric(total_genus_occupancy$Bin)
    total_genus_occupancy$Perc_Occupied_cells <- as.numeric(total_genus_occupancy$Perc_Occupied_cells)
    total_genus_occupancy <<- total_genus_occupancy
  }
  all_total_cells <<- all_total_cells
}

#=================================== OCC_RASTER ==========================================

# Set up background info

occ_raster <- function(tax_name){
  regexp <- "[[:digit:]]+"
  stack_names <- rev(sort(grep("Stack",names(.GlobalEnv),value=TRUE)))
  tax_stack <- stack()
  for (t in stack_names){
    temp_stack <- get(t)
    if (tax_name %in% names(temp_stack) == TRUE){
      named_layer <- raster::subset(temp_stack, grep(tax_name, names(temp_stack), value = T))
      values(named_layer)[values(named_layer) == 0] = NA
      bin <- str_extract(t, regexp)
      names(named_layer) <- paste(tax_name, "_Bin_", bin, sep = "")
      tax_stack <- stack(tax_stack, named_layer)
    }
  }
  tax_stack <<- tax_stack
}


occ_raster_plotter <- function(tax_stack){
  countries <- maps::map("world", plot=FALSE, fill = TRUE) # find map to use as backdrop
  states <- maps::map("state", plot = FALSE, fill = TRUE)
  countries <<- maptools::map2SpatialPolygons(countries, IDs = countries$names, proj4string = CRS("+proj=longlat")) # Turn map into spatialpolygons
  states <<- maptools::map2SpatialPolygons(states, IDs = states$names, proj4string = CRS("+proj=longlat")) # Turn map into spatialpolygons
  mapTheme <- rasterVis::rasterTheme(region=brewer.pal(8,"Greens")) # MGVF
  
  title_name <- gsub("_.*","",names(tax_stack)[1])
  subbing <- paste(".*", title_name, "_", sep = "")
  names(tax_stack) <- sub(subbing, '', names(tax_stack))
  rasterVis::levelplot(tax_stack, margin=F, par.settings=mapTheme,  main = title_name) + #create levelplot for raster
    latticeExtra::layer(sp.polygons(states, col = "white", fill = NA), under = T)  +                   
    latticeExtra::layer(sp.polygons(countries, col = 0, fill = "light grey"), under = T)
}


#=================================== GET_EXTENT ==========================================

# Setup raster for resolution and extent. Note: these values should be the same ones used for the file 1.Setup_occupancy_DR.
get_extent <- function(data){
  maxLat <- round_any((max(data$lat) + 1), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
  minLat <- round_any((min(data$lat) - 1), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
  maxLng <- round_any((max(data$lng) + 1), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
  minLng <- round_any((min(data$lng) - 1), 0.5) #get value for defining extent, and increase by x for visualisation purposes
  e <<- extent(minLng, maxLng, minLat, maxLat) # build extent object
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
