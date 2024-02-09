################################################################################
######## TAPHONOMIC CONTROLS ON A MULTI-ELEMENT SKELETAL FOSSIL RECORD #########
################################################################################

# Jeffrey R. Thompson, Christopher D. Dean, Madeline Ford, Timothy A. M. Ewin
# 2023
# Script written by Christopher D. Dean

################################################################################
#                              FILE 0: OLD CODE                                #
################################################################################

# Function for sorting jaccard similarity combined
Jaccard.1 <- function(data.set, dependent, selection, ...){
  # Transform to Preservation_score by genus matrix
  jac.1 <- dcast(data.set,
                 dependent ~ ...,
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
  colnames(temp.res) <- c("Test", "Score", selection)
  temp.res$Test <- factor(temp.res$Test, levels = c("1v2", "2v3", 
                                                    "3v4", "4v5", 
                                                    "1v3", "2v4", 
                                                    "3v5", "1v4", 
                                                    "2v5", "1v5"))
  temp.res <- with(temp.res, temp.res[order(Test, Score),])
  return(temp.res)
}

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

pres.score <- 1
data.LR <- set_Pres_score(data.set, c(pres.score))
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
