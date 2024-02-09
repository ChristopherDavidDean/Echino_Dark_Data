# The following lines of code perform all analyses from Thompson and Bottjer (2018), 
# "Quantitative analysis of substrate preference in Carboniferous stem group echinoids"
# as long as you have the files "Supplemental_Table_2" and "Supplemental_Table_1", which are in the supplement of the paper.

# The directories for import and export are based on my computer (JRT), and should 
# be edited to whatever directories on your machine you wish to use. Have fun!

#The following three lines install packages necessary for data manipulation and exporting data
install.packages("dplyr")
install.packages("geoscale")
install.packages("xlsx")
install.packages("vegan")
Collections<-read.csv("/Users/jeffreythompson/Research/Substrate_Affinity/Submission_1/Supplemental_Table_2.csv") #Imports Carboniferous PaleobioDB collections

#Makes tables of all siliciclastic and carbonate collections from "Supplemental_Table_2" file downloaded form the PaleoBioDB
Siliciclastic_Collections<-subset(Collections, Collections[,19]=="conglomerate"| Collections[,19]=="claystone"| Collections[,19]=="sandstone"| Collections[,19]=="shale"| Collections[,19]=="mudstone"| Collections[,19]=="siltstone"| Collections[,18]==" silty sandstone")
Num_Silic_Collections<-nrow(Siliciclastic_Collections)

Carbonate_Collections<-subset(Collections, Collections[,19]=="chert" | Collections[,19]=="floatstone"| Collections[,19]=="framestone"| Collections[,19]=="rudstone"| Collections[,19]=="grainstone"| Collections[,19]=="reef rocks"| Collections[,19]=="bindstone"| Collections[,19]=="bafflestone" | Collections[,19]=="dolomite"| Collections[,19]=="lime mudstone"| Collections[,19]=="limestone"| Collections[,19]=="carbonate"| Collections[,19]=="packstone"| Collections[,19]=="wackestone")
Num_Carb_Collections<-nrow(Carbonate_Collections)

#This imports the raw database of echinoid occurrences from the Supplemental_Table_1 file
setwd("/Users/jeffreythompson/Documents/Research/Palaeozoic_Echinoid_Diversity")

Total<-read.delim("Carboniferous_Occurrences.txt")

# New is the primary working database for all analyses performed in the manuscript.
# The next few lines of code remove duplicate occurrences such that only unique occurrences remain for analyases
New<-distinct(Total, Species, Genus, Family, Locality, Formation, Lithology, Grain.Size, .keep_all=TRUE)
Try <-New[!duplicated(New[c(3, 6, 7, 4, 11, 10)]),]
Replicates_1<-Try[Try[,3]=="?" & duplicated(Try[,c(6, 7, 4, 10)]),]
Try<-anti_join(Try, Replicates_1, by=c("Genus", "Locality", "Formation", "Family", "Lithology", "Grain.Size"))
Replicates_2<-Try[Try[,11]=="?" & duplicated(Try[,c(3, 6, 7, 4, 10)]),] 
New<-anti_join(Try, Replicates_2, by=c("Genus", "Locality", "Formation", "Family", "Lithology", "Grain.Size")) 
New<-as.matrix(New)


#This loop changes the more specific grain size categories in New into either "Fine_Grained" or "Coarse_Grained"
for (l in 1:nrow(New)){
  if (New[l,11]=="Fine Grained" | New[l,11]=="Fine Grained Sandstone" | New[l,11] =="very fine grained sandstone" | New[l,11] =="Silty mudstone" | New[l,11]=="Shale" | New[l,11]=="Wackestone" | New[l,11]=="mudstone" | New[l,11]=="Siltstone" | New[l,11]== "mudstone/wackestone"){
        New[l,11]<-"Fine_Grained"
  }
        else if (New[l,11]=="Grainstone" | New[l,11] == "Packstone"| New[l,11] == "Coarse Grained Sandstone"){ 
          New[l,11]<-"Coarse_Grained"
        }
}

#Creates vectors for names of stages and Subperiods
TimeBins<-c("Tournaisian", "Visean", "Serpukhovian", "Bashkirian", "Moscovian", "Kasimovian", "Gzhelian")
SubPeriods<-c("Mississippian", "Pennsylvanian")
Taxon_Names<-c("Palaechinidae", "Archaeocidaridae", "Proterocidaridae", "Lepidesthidae", "Lepidocentridae")

#Computes Number of Carbonate Collections from collection matrix in Each Stage
Carb_Stages<-c(length(TimeBins), 0)
for (l in 1:length(TimeBins)){
  Carb_Stages[l]<-nrow(subset(Carbonate_Collections, Carbonate_Collections[,10]==TimeBins[l]))
}

#Computes Number of Carbonate Collections from collection matrix in Each Subperiod
Carb_SubPeriods<-c(length(SubPeriods), 0)
Carb_SubPeriods[1]<-sum(Carb_Stages[1:3])
Carb_SubPeriods[2]<-sum(Carb_Stages[4:7])
Num_Carb_Tot<-sum(Carb_SubPeriods)

#Computes Number of Siliciclastic Collections from collection matrix in Each Stage
Silic_Stages<-c(length(TimeBins), 0)
for (l in 1:length(TimeBins)){
  Silic_Stages[l]<-nrow(subset(Siliciclastic_Collections, Siliciclastic_Collections[,10]==TimeBins[l]))
}

#Computes Number of Siliciclastic Collections from collection matrix in Each Subperiod
Silic_SubPeriods<-c(length(SubPeriods), 0)
Silic_SubPeriods[1]<-sum(Silic_Stages[1:3])
Silic_SubPeriods[2]<-sum(Silic_Stages[4:7])

#Computes Number of Siliciclastic Collections in entire Carboniferous
Num_Silic_Tot<-sum(Silic_SubPeriods)
Proportion_Silic_Carb_Tot<-Num_Silic_Tot/(Num_Carb_Tot+Num_Silic_Tot)
Proportion_Carb_Silic_Tot<-Num_Carb_Tot/(Num_Silic_Tot+Num_Carb_Tot)


#Calculates Stage Level Proportions for ratio of carbonate collections to siliciclastic collections and vice versa
Proportions_Stage_Carb_Silic<-c(length(TimeBins), 0)
Proportions_Stage_Silic_Carb<-c(length(TimeBins), 0)
for (l in 1:length(TimeBins)){
  Proportions_Stage_Carb_Silic[l]<-Carb_Stages[l]/(Carb_Stages[l]+Silic_Stages[l])
  Proportions_Stage_Silic_Carb[l]<-Silic_Stages[l]/(Silic_Stages[l]+Carb_Stages[l])
}

#Plot Stage Level Number of Carb and Silic Collections
geoscalePlot(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Carb_Stages, units=c( "Age", "Epoch"), type="b", lwd=2, data.lim=c(0, 700), label="Number of Collections", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
lines(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Silic_Stages, type="b", col=b, lwd=2, pch=2)
legend(325, 700, c("Carbonates", "Siliciclastics", "Proportion of Carbonates"), cex=.6, pch=c(1:2), col=c(1:2,"azure4"), text.width=15)
par(new=T)
geoscalePlot(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Proportions_Stage_Carb_Silic, units=c( "Age", "Epoch"), type="b", lwd=2, data.lim=c(0, 1), label="Number of Collections", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9, no.axis=TRUE, col="azure4")
axis(4, labels=TRUE, ylab="Proportion of Occurrences")


#Calculates Subperiod Level Proportions for ratio of carbonates to siliciclastics and vice versa
Proportions_Subperiod_Carb_Silic<-c(length(SubPeriods), 0)
Proportions_Subperiod_Silic_Carb<-c(length(SubPeriods), 0)
for (l in 1:length(SubPeriods)){
  Proportions_Subperiod_Carb_Silic[l]<-Carb_SubPeriods[l]/(Carb_SubPeriods[l]+Silic_SubPeriods[l])
  Proportions_Subperiod_Silic_Carb[l]<-Silic_SubPeriods[l]/(Silic_SubPeriods[l]+Carb_SubPeriods[l])
}

#Plot Subperiod Level Number of Carb and Silic Collections
geoscalePlot(c(341.05, 311.05), Carb_SubPeriods, units=c( "Epoch"), type="b", lwd=2, data.lim=c(0, 1600), label="Number of Collections", erotate=0, direction="horizontal", abbrev=c("Subperiod"), age.lim=c(357.1,301.1), cex.ts=.9)
  lines(c(341.05, 311.05), Silic_SubPeriods, type="b", col=b, lwd=2, pch=2)
legend(325, 500, c("Carbonates", "Siliciclastics","Proportion of Carbonates"), cex=.6, pch=c(1:2), col=c(1:2,"azure4"), text.width=15)
par(new=T)
geoscalePlot(c(341.05, 311.05), Proportions_Subperiod_Carb_Silic, units=c("Epoch"), type="b", lwd=2, data.lim=c(0, 1), label="Number of Collections", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9, no.axis=TRUE, col="azure4")
axis(4, labels=TRUE, ylab="Proportion of Occurrences")

#Plot Number of Carb and Silic Collections For Entire Carboniferous
geoscalePlot(c(328.9), Num_Carb_Tot, units=c( "Period"), type="b", lwd=2, data.lim=c(0, 3000), label="Number of Collections", erotate=0, direction="horizontal", abbrev=c("Subperiod"), age.lim=c(357.1,301.1), cex.ts=.9)
lines(c(328.9), Num_Silic_Tot, type="b", col=b, lwd=2, pch=2)
legend(325, 3000, c("Carbonates", "Siliciclastics","Proportion of Carbonates"), cex=.6, pch=c(1:2), col=c(1:2,"azure4"), text.width=15)
par(new=T)
geoscalePlot(c(328.9), Proportion_Carb_Silic_Tot, units=c( "Period"), type="b", lwd=2, data.lim=c(0, 1), label="Number of Collections", erotate=0, direction="horizontal", abbrev=c("Subperiod"), age.lim=c(357.1,301.1), cex.ts=.9, no.axis=TRUE, col="azure4")
axis(4, labels=TRUE, ylab="Proportion of Occurrences")

#Stage Level Occurrences
#Calculates Total Number of Taxa in Each Stage
Total_Matrix_Stage<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
rownames(Total_Matrix_Stage)<-Taxon_Names
colnames(Total_Matrix_Stage)<-TimeBins
for (k in 1:length(Taxon_Names)){
  for (p in 1:length(TimeBins)){
    Total_Taxon_Stage<-subset(New, New[,4]==Taxon_Names[k] & New[,8]==TimeBins[p] )
    Total_Matrix_Stage[k,p]<-nrow(Total_Taxon_Stage)
  }
}



#This imports the raw database of echinoid occurrences from the Supplemental_Table_1 file
setwd("/Users/jeffreythompson/Documents/Research/Palaeozoic_Echinoid_Diversity")

Total<-read.delim("Carboniferous_Occurrences.txt")





# New is the primary working database for all analyses performed in the manuscript.
# The next few lines of code remove duplicate occurrences such that only unique occurrences remain for analyases
New<-distinct(Total, Species, Genus, Family, Locality, Formation, .keep_all=TRUE)
#Try<-New[!duplicated(New[c(3, 6, 7, 4, 11, 10)]),]
#Replicates_1<-Try[Try[,3]=="?" & duplicated(Try[,c(6, 7, 4, 10)]),]
#Try<-anti_join(Try, Replicates_1, by=c("Genus", "Locality", "Formation", "Family", "Lithology", "Grain.Size"))
#Replicates_2<-Try[Try[,11]=="?" & duplicated(Try[,c(3, 6, 7, 4, 10)]),] 
#New<-anti_join(Try, Replicates_2, by=c("Genus", "Locality", "Formation", "Family", "Lithology", "Grain.Size")) 
#New<-as.matrix(New)



New<-data.frame(New)



#Make presence abundance lists for multiple iterations to account for differences in stratigraphic uncertainty
#This doesn't currently account for range through, it should
Diversity_Iterations<-matrix(nrow=100, ncol=length(TimeBins), 0)

for (r in 1:100){

New_iterate<-New
for (l in 1:nrow(New)){
  if (New_iterate[l,]$Age=="Bashkirian-Moscovian"){
    New_iterate[l,]$Age<-sample(c("Bashkirian", "Moscovian"),1)
  }
  if (New_iterate[l,]$Age=="Moscovian or Kasimovian"){
    New_iterate[l,]$Age<-sample(c("Kasimovian", "Moscovian"),1)
  }
  if (New_iterate[l,]$Age=="Visean or Serpukhovian"){
    New_iterate[l,]$Age<-sample(c("Serpukhovian", "Visean"),1)
  }
  if (New_iterate[l,]$Age=="Tournaisian or Visean"){
    New_iterate[l,]$Age<-sample(c("Tournaisian", "Visean"),1)
  }
  if (New_iterate[l,]$Age=="Tournaisian to Visean"){
    New_iterate[l,]$Age<-sample(c("Tournaisian", "Visean"),1)
  }
  if (New_iterate[l,]$Age=="?Tournaisian-Visean"){
    New_iterate[l,]$Age<-sample(c("Tournaisian", "Visean"),1)
  }
  if (New_iterate[l,]$Age=="Serpukhovian or Bashkirian"){
    New_iterate[l,]$Age<-sample(c("Serpukhovian", "Bashkirian"),1)
  }
  if (New_iterate[l,]$Age=="Serpukhovian, Bashkirian or Moscovian"){
    New_iterate[l,]$Age<-sample(c("Serpukhovian", "Bashkirian", "Moscovian"),1)
  }
  if (New_iterate[l,]$Age=="Kasimovian or Gzhelian"){
    New_iterate[l,]$Age<-sample(c("Kasimovian", "Gzhelian"),1)
  }
  if (New_iterate[l,]$Age=="Pennsylvanian"){
    New_iterate[l,]$Age<-sample(c("Kasimovian", "Bashkirian", "Moscovian", "Gzhelian"),1)
  }
  if (New_iterate[l,]$Age=="Mississippian"){
    New_iterate[l,]$Age<-sample(c("Tournaisian", "Visean", "Serpukhovian"),1)
  }
  if (New_iterate[l,]$Age=="Carboniferous"){
    New_iterate[l,]$Age<-sample(c("Tournaisian", "Visean", "Serpukhovian", "Kasimovian", "Bashkirian", "Moscovian", "Gzhelian"),1)
  }
  if (New_iterate[l,]$Age=="Mississippian-Pennsylvanian"){
    New_iterate[l,]$Age<-sample(c("Tournaisian", "Visean", "Serpukhovian", "Kasimovian", "Bashkirian", "Moscovian", "Gzhelian"),1)
  }
  if (New_iterate[l,]$Age=="Mississippian or Pennsylvanian"){
    New_iterate[l,]$Age<-sample(c("Tournaisian", "Visean", "Serpukhovian", "Kasimovian", "Bashkirian", "Moscovian", "Gzhelian"),1)
  }
}
Genus_Names<-unique(New_iterate$Genus)
Total_Matrix<-matrix(nrow=length(Genus_Names), ncol=length(TimeBins), 0)
Total_Matrix_PA<-matrix(nrow=length(Genus_Names), ncol=length(TimeBins), 0)
rownames(Total_Matrix)<-Genus_Names
colnames(Total_Matrix)<-TimeBins
rownames(Total_Matrix_PA)<-Genus_Names
colnames(Total_Matrix_PA)<-TimeBins

for (k in 1:length(Genus_Names)){
  for (p in 1:length(TimeBins)){
    Total_Taxon_Stage<-subset(New_iterate, New_iterate[,3]==Genus_Names[k] & New_iterate[,8]==TimeBins[p] )
    Total_Matrix[k,p]<-nrow(Total_Taxon_Stage)
    if (Total_Matrix[k,p]>=1){
      Total_Matrix_PA[k,p]<-1
    }
  }
}

#Remove question mark genera
Total_Matrix_PA<-Total_Matrix_PA[c(-16),]
Total_Matrix<-Total_Matrix[c(-16),]

#Add range through modification
for (l in 1:nrow(Total_Matrix_PA)){
  for (d in 1:ncol(Total_Matrix_PA)){
    if (Total_Matrix_PA[l,d-1]==1 && Total_Matrix_PA[l,d]==0 && Total_Matrix_PA[l,d:ncol(Total_Matrix_PA)]==1){
      Total_Matrix_PA[l,d]<-1
    }
  }
}


#Make sums
Sums<-matrix(ncol=length(TimeBins), 0)
for (i in 1:ncol(Total_Matrix_PA)){
  Sums[i]<-sum(Total_Matrix_PA[,i])
}

Diversity_Iterations[r,]<-Sums
}


Medians<--matrix(ncol=length(TimeBins), 0)
ninety_five_high<--matrix(ncol=length(TimeBins), 0)
ninety_five_low<--matrix(ncol=length(TimeBins), 0)
for (a in 1:length(TimeBins)){
  Medians[a]<-median(Diversity_Iterations[,a])
  ninety_five_high[a]<-quantile(Diversity_Iterations[,a], .975)
  ninety_five_low[a]<-quantile(Diversity_Iterations[,a], .025)
}


Carb_Stage_Midpoints<-c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3)

geoscalePlot(Carb_Stage_Midpoints, Medians, units=c( "Age", "Epoch"), type="b", lwd=2, data.lim=c(0, 20), label="Genus Richness", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)


polygon(c(Carb_Stage_Midpoints,rev(Carb_Stage_Midpoints)), c(ninety_five_high,rev(ninety_five_low)), col=rgb(1, 0, 0,0.5), border = FALSE)






#Make presence abundance list for single iteration without including stratigraphically uncertain occurrences
Genus_Names<-unique(New$Genus)
Total_Matrix<-matrix(nrow=length(Genus_Names), ncol=length(TimeBins), 0)
Total_Matrix_PA<-matrix(nrow=length(Genus_Names), ncol=length(TimeBins), 0)
rownames(Total_Matrix)<-Genus_Names
colnames(Total_Matrix)<-TimeBins
rownames(Total_Matrix_PA)<-Genus_Names
colnames(Total_Matrix_PA)<-TimeBins

for (k in 1:length(Genus_Names)){
  for (p in 1:length(TimeBins)){
    Total_Taxon_Stage<-subset(New, New[,3]==Genus_Names[k] & New[,8]==TimeBins[p] )
    Total_Matrix[k,p]<-nrow(Total_Taxon_Stage)
    if (Total_Matrix[k,p]>=1){
      Total_Matrix_PA[k,p]<-1
    }
  }
}

#Remove question mark genera
Total_Matrix_PA<-Total_Matrix_PA[c(-16),]
Total_Matrix<-Total_Matrix[c(-16),]

#Add range through modification. Think this works, but need to check

Total_Matrix_PA_Range_Through<-Total_Matrix_PA
for (l in 1:nrow(Total_Matrix_PA)){
  for (d in 2:(ncol(Total_Matrix_PA)-1)){
    if (Total_Matrix_PA[l,d]==0 && Total_Matrix_PA[l,d-1]==1 && Total_Matrix_PA[l,(d+1):ncol(Total_Matrix_PA)]==1){
      Total_Matrix_PA_Range_Through[l,d]<-1
      #print(Total_Matrix_PA[l,d])
    }
  }
}

#Take sums of each time bin
Sums<-matrix(ncol=length(TimeBins), 0)
for (i in 1:ncol(Total_Matrix_PA)){
  Sums[i]<-sum(Total_Matrix_PA[,i])
}
 

geoscalePlot(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Sums, units=c( "Age", "Epoch"), type="b", lwd=2, data.lim=c(0, 20), label="Genus Richness", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)




for (x in 2:nrow(Carbonate_Matrix_Plot)){
  lines(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Carbonate_Matrix_Plot[x,], type="b", pch=x, col=x, lwd=2)
}
legend(320, 70, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)













#Code for paleoecology analyses. NEED TO SORT OUT HOW TO DEAL WITH SP.'s HERE!
Localities<-unique(New$Locality)

Formations<-unique(New$Formation)

Subperiods<-unique(New$Subperiod)

Ages<-unique(New$Age)

#Select the taxonomic level for analyses
Taxon_Names<-unique(New$Family)
Taxon_Column<-4

#Taxon_Names<-unique(New$Genus)
#Taxon_Column<-3

#Select the sites level for analyses (localities or ages)
#Sites<-Localities
#Sites_Column<-6

Sites<-Ages
Sites_Column<-8

#Select the cutoff threshold for number of taxa to include
Occurrence_Threshold<-3

#Make Empty Sites by Taxon matrices using abundance and presence absence data
Total_Matrix_Localities<-matrix(nrow=length(Taxon_Names), ncol=length(Sites), 0)
Total_Matrix_Localities_PA<-matrix(nrow=length(Taxon_Names), ncol=length(Sites), 0)
rownames(Total_Matrix_Localities)<-Taxon_Names
colnames(Total_Matrix_Localities)<-Sites
rownames(Total_Matrix_Localities_PA)<-Taxon_Names
colnames(Total_Matrix_Localities_PA)<-Sites
Matrix_ages<-matrix(nrow=1, ncol=length(Sites), 0)
Matrix_subperiods<-matrix(nrow=1, ncol=length(Sites), 0)


#Fill abundance and P/A matrices based off of file and corresponding file for age and subperiod classification of groups
for (k in 1:length(Taxon_Names)){
  for (p in 1:length(Sites)){
    Total_Taxon_Localities<-subset(New, New[,Taxon_Column]==Taxon_Names[k] & New[,Sites_Column]==Sites[p] )
    Total_Matrix_Localities[k,p]<-nrow(Total_Taxon_Localities)
    if (Total_Matrix_Localities[k,p]>=1){
      Total_Matrix_Localities_PA[k,p]<-1
      Matrix_ages[p]<-Total_Taxon_Localities$Age
      Matrix_subperiods[p]<-Total_Taxon_Localities$Subperiod
    }
  }
}


#Remove question mark genera and localities with just question marks or localities named ?
Question_Mark_Taxa<-rownames(Total_Matrix_Localities)=="?"
Question_Mark_Taxa_Index<-which(Question_Mark_Taxa==TRUE)
Total_Matrix_Localities<-Total_Matrix_Localities[-c(Question_Mark_Taxa_Index),]

#Remove sites named ? and corresponding sites in group matrices
Question_Mark_Sites<-colnames(Total_Matrix_Localities)=="?"
Question_Mark_Sites_Index<-which(Question_Mark_Sites==TRUE)
Total_Matrix_Localities<-Total_Matrix_Localities[,-c(Question_Mark_Sites_Index)]
Matrix_ages<-Matrix_ages[-c(Question_Mark_Sites_Index)]
Matrix_subperiods<-Matrix_subperiods[-c(Question_Mark_Sites_Index)]

#The below is for use if I don't use a threshold
#Total_Matrix_Localities_index<-colSums(Total_Matrix_Localities)>0
#Total_Matrix_Localities<-Total_Matrix_Localities[,colSums(Total_Matrix_Localities)>0]
#Localities_above_zero<-which(Total_Matrix_Localities_index==TRUE)
#Matrix_ages[Localities_above_zero]
#Matrix_subperiods[Localities_above_zero]

#Remove localities with less than n occurrences from occurrence matrices
Total_Matrix_Localities_index<-colSums(Total_Matrix_Localities)>Occurrence_Threshold #This number is n
Total_Matrix_Localities<-Total_Matrix_Localities[,colSums(Total_Matrix_Localities)>Occurrence_Threshold]#This number is n
Localities_above_threshold<-which(Total_Matrix_Localities_index==TRUE)
Ages_for_taxa<-Matrix_ages[Localities_above_threshold]
Subperiods_for_taxa<-Matrix_subperiods[Localities_above_threshold]

#Color localities based upon Age
Ages_for_color<-as.numeric(as.factor(Ages_for_taxa))

#Color localities based upon Subperiod
Subperiods_for_color<-as.numeric(as.factor(Subperiods_for_taxa))


#Quick and dirty mod for Age-based analyses. Need to figure easy way to remove uncertain ages 
Total_Matrix_Localities<-Total_Matrix_Localities[,c(1:5, 9:10)]
Subperiods_for_taxa<-Subperiods_for_taxa[c(1:5, 9:10)]
Subperiods_for_color<-as.numeric(as.factor(Subperiods_for_taxa))


#Convert Matrix to percentages
Total_Matrix_Localities_percentages<-Total_Matrix_Localities
for (k in 1:nrow(Total_Matrix_Localities)){
  for (p in 1:ncol(Total_Matrix_Localities)){
    Total_Matrix_Localities_percentages[k,p]<-(Total_Matrix_Localities[k,p]/sum(Total_Matrix_Localities[,p]))
    }
  }
Total_Matrix_Localities_percentages



#Run NMDS on Sites
library(vegan)
Total_Matrix_clustering<-t(Total_Matrix_Localities)
Total_Matrix_Localities_percentages_clustering<-t(Total_Matrix_Localities_percentages)
set.seed(2)
NMDS=metaMDS(Total_Matrix_clustering, k=2,trymax=1000)

#Plot stress to make sure is ok
stressplot(NMDS)

#Plot NMDS Results
colvec_for_stages <- c("red", "blue", "purple", "green")  
colvec_for_subperiods <- c("red", "blue", "green")


plot(NMDS,type="n")
ordihull(NMDS,Subperiods_for_color,draw="polygon",col=colvec_for_stages,label=F)
orditorp(NMDS,display="sites",cex=0.25, col = colvec_for_stages[Subperiods_for_color])

plot(NMDS,type="n")
points(NMDS,display="sites",cex=0.25, col = colvec_for_stages[Subperiods_for_color])
ordihull(NMDS,Subperiods_for_color,draw="polygon",col=colvec_for_stages,label=F)


#Run cluster analysis on Sites

simprof(Total_Matrix_Localities_percentages_clustering , num.expected=1000, num.simulated=999,
  method.cluster="average", method.distance="braycurtis", 
method.transform="identity", alpha=0.05,
sample.orientation="row", const=0,
silent=TRUE, increment=100)

distance_matrix <- vegdist(Total_Matrix_Localities_percentages_clustering, method = "bray") # distance matrix
clustering <- hclust(distance_matrix, method="average")
plot(clustering)

















#Calculates Total Number of Carbonate Occurrences in Each Stage. Last line spits out excel file of carbonate occurrences
Carbonate<-subset(New, New[,10] =="Carbonate")
Carbonate_Matrix<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
rownames(Carbonate_Matrix)<-Taxon_Names
colnames(Carbonate_Matrix)<-TimeBins
for (k in 1:length(Taxon_Names)){
  for (p in 1:length(TimeBins)){
    Carbonate_Taxon_Time<-subset(Carbonate, Carbonate[,4]==Taxon_Names[k] & Carbonate[,8]==TimeBins[p])
    Carbonate_Matrix[k,p]<-nrow(Carbonate_Taxon_Time)
  }
}
write.xlsx(Carbonate_Matrix, "/Users/jeffreythompson/Research/Substrate_Affinity/Matrices/Carbonate_Matrix_Stages.xlsx")

#Plot Stage Level Number of Occurrences in Carbonates for Each Family
Carbonate_Matrix_Plot<-Carbonate_Matrix
for (u in 1:nrow(Carbonate_Matrix_Plot)){
  for (x in 1:ncol(Carbonate_Matrix_Plot)){
if (Carbonate_Matrix_Plot[u,x]==0){
  Carbonate_Matrix_Plot[u,x]<-NA
    }
  }
}
geoscalePlot(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Carbonate_Matrix_Plot[1,], units=c( "Age", "Epoch"), type="b", lwd=2, data.lim=c(0, 80), label="Number of Occurrences", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(Carbonate_Matrix_Plot)){
  lines(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Carbonate_Matrix_Plot[x,], type="b", pch=x, col=x, lwd=2)
}
legend(320, 70, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)

#Calculates Total Number of Siliciclastic Occurrences in Each Stage. Last line spits out excel file of siliciclastic occurrences
Siliciclastic<-subset(New, New[,10] =="Siliciclastic")
Siliciclastic_Matrix<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
rownames(Siliciclastic_Matrix)<-Taxon_Names
colnames(Siliciclastic_Matrix)<-TimeBins
for (k in 1:length(Taxon_Names)){
  for (p in 1:length(TimeBins)){
    Siliciclastic_Taxon_Time<-subset(Siliciclastic, Siliciclastic[,4]==Taxon_Names[k] & Siliciclastic[,8]==TimeBins[p])
    Siliciclastic_Matrix[k,p]<-nrow(Siliciclastic_Taxon_Time)
  }
}
write.xlsx(Siliciclastic_Matrix, "/Users/jeffreythompson/Research/Substrate_Affinity/Matrices/Siliciclastic_Matrix_Stages.xlsx")

#Plot Stage Level Number of Occurrences in Siliciclastics for Each Family
Siliciclastic_Matrix_Plot<-Siliciclastic_Matrix
for (u in 1:nrow(Siliciclastic_Matrix_Plot)){
  for (x in 1:ncol(Siliciclastic_Matrix_Plot)){
    if (Siliciclastic_Matrix_Plot[u,x]==0){
      Siliciclastic_Matrix_Plot[u,x]<-NA
    }
  }
}
geoscalePlot(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Siliciclastic_Matrix_Plot[1,], units=c( "Age", "Epoch"), type="b", lwd=2, data.lim=c(0, 20), label="Number of Occurrences", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(Siliciclastic_Matrix_Plot)){
  lines(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Siliciclastic_Matrix_Plot[x,], type="b", pch=x, col=x, lwd=2)
}
legend(320, 15, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)

#Total Number of Fine Grained Occurrences in Each Stage. Last line cretes excel file of fine-grained occurrences in each stage
Fine_Grained<-subset(New, New[,11] =="Fine_Grained")
Fine_Matrix<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
rownames(Fine_Matrix)<-Taxon_Names
colnames(Fine_Matrix)<-TimeBins
for (k in 1:length(Taxon_Names)){
  for (p in 1:length(TimeBins)){
    Fine_Taxon_Time<-subset(Fine_Grained, Fine_Grained[,4]==Taxon_Names[k] & Fine_Grained[,8]==TimeBins[p])
    Fine_Matrix[k,p]<-nrow(Fine_Taxon_Time)
  }
}
write.xlsx(Fine_Matrix, "/Users/jeffreythompson/Research/Substrate_Affinity/Matrices/Fine_Grained_Matrix_Stages.xlsx")

#Plot Stage Level Number of Occurrences in Fine-Grained Sediments for Each Family
Fine_Matrix_Plot<-Fine_Matrix
for (u in 1:nrow(Fine_Matrix_Plot)){
  for (x in 1:ncol(Fine_Matrix_Plot)){
    if (Fine_Matrix_Plot[u,x]==0){
      Fine_Matrix_Plot[u,x]<-NA
    }
  }
}
geoscalePlot(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Fine_Matrix_Plot[1,], units=c( "Age", "Epoch"), type="b", lwd=2, data.lim=c(0, 70), label="Number of Occurrences", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(Fine_Matrix_Plot)){
  lines(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Fine_Matrix_Plot[x,], type="b", pch=x, col=x, lwd=2)
}
legend(320, 50, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)

#Total Number of Stage level collections used in analyses of Carb vs. Silic 
Tot_Stage_Carb_Silic<-Carbonate_Matrix+Siliciclastic_Matrix

#Total Number of Coarse Grained Occurrences in Each Stage. Last line cretes excel file of coarse-grained occurrences in each stage
Coarse_Grained<-subset(New, New[,11] =="Coarse_Grained")
Coarse_Matrix<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
rownames(Coarse_Matrix)<-Taxon_Names
colnames(Coarse_Matrix)<-TimeBins
for (k in 1:length(Taxon_Names)){
  for (p in 1:length(TimeBins)){
    Coarse_Taxon_Time<-subset(Coarse_Grained, Coarse_Grained[,4]==Taxon_Names[k] & Coarse_Grained[,8]==TimeBins[p])
    Coarse_Matrix[k,p]<-nrow(Coarse_Taxon_Time)
  }
}
write.xlsx(Coarse_Matrix, "/Users/jeffreythompson/Research/Substrate_Affinity/Matrices/Coarse_Grained_Matrix_Stages.xlsx")

#Plot Stage Level Number of Occurrences in Coarse-Grained Sediments for Each Family
Coarse_Matrix_Plot<-Coarse_Matrix
for (u in 1:nrow(Coarse_Matrix_Plot)){
  for (x in 1:ncol(Coarse_Matrix_Plot)){
    if (Coarse_Matrix_Plot[u,x]==0){
      Coarse_Matrix_Plot[u,x]<-NA
    }
  }
}
geoscalePlot(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Coarse_Matrix_Plot[1,], units=c( "Age", "Epoch"), type="b", lwd=2, data.lim=c(0, 20), label="Number of Occurrences", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(Coarse_Matrix_Plot)){
  lines(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Coarse_Matrix_Plot[x,], type="b", pch=x, col=x, lwd=2)
}
legend(320, 18, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)


#Total Number of Stage level collections used in analyses of coarse vs. fine
Tot_Stage_Fine_Coarse<-Fine_Matrix+Coarse_Matrix

#Subperiod Level Occurrences
#Total Number of Taxa in Each Subperiod
Total_Matrix_SubPeriod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
rownames(Total_Matrix_SubPeriod)<-Taxon_Names
colnames(Total_Matrix_SubPeriod)<-SubPeriods
for (k in 1:length(Taxon_Names)){
  for (p in 1:length(SubPeriods)){
    Total_Taxon_SubPeriods<-subset(New, New[,4]==Taxon_Names[k] & New[,19]==SubPeriods[p] )
    Total_Matrix_SubPeriod[k,p]<-nrow(Total_Taxon_SubPeriods)
  }
}

#Total Number of Carbonate Occurrences in Each Subperiod
Carbonate<-subset(New, New[,10] =="Carbonate")
Carbonate_Matrix_SubPeriod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
rownames(Carbonate_Matrix_SubPeriod)<-Taxon_Names
colnames(Carbonate_Matrix_SubPeriod)<-SubPeriods
for (k in 1:length(Taxon_Names)){
  for (p in 1:length(SubPeriods)){
    Carbonate_Taxon_Time<-subset(Carbonate, Carbonate[,4]==Taxon_Names[k] & Carbonate[,19]==SubPeriods[p])
    Carbonate_Matrix_SubPeriod[k,p]<-nrow(Carbonate_Taxon_Time)
  }
}
Carbonate_Matrix_SubPeriod_Plot<-Carbonate_Matrix_SubPeriod
for (u in 1:nrow(Carbonate_Matrix_SubPeriod_Plot)){
  for (x in 1:ncol(Carbonate_Matrix_SubPeriod_Plot)){
    if (Carbonate_Matrix_SubPeriod_Plot[u,x]==0){
      Carbonate_Matrix_SubPeriod_Plot[u,x]<-NA
    }
  }
}
geoscalePlot(c(341.05, 311.05), Carbonate_Matrix_SubPeriod_Plot[1,], units=c( "Epoch"), type="b", lwd=2, data.lim=c(0, 150), label="Number of Occurrences", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(Carbonate_Matrix_SubPeriod_Plot)){
  lines(c(341.05, 311.05), Carbonate_Matrix_SubPeriod_Plot[x,], type="b", pch=x, col=x, lwd=2)
}
legend(320, 150, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)


#Total Number of Siliciclastic Occurrences in Each Subperiod
Siliciclastic<-subset(New, New[,10] =="Siliciclastic")
Siliciclastic_Matrix_SubPeriod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
rownames(Siliciclastic_Matrix_SubPeriod)<-Taxon_Names
colnames(Siliciclastic_Matrix_SubPeriod)<-SubPeriods
for (k in 1:length(Taxon_Names)){
  for (p in 1:length(SubPeriods)){
    Siliciclastic_Taxon_Time<-subset(Siliciclastic, Siliciclastic[,4]==Taxon_Names[k] & Siliciclastic[,19]==SubPeriods[p])
    Siliciclastic_Matrix_SubPeriod[k,p]<-nrow(Siliciclastic_Taxon_Time)
  }
}
Siliciclastic_Matrix_SubPeriod_Plot<-Siliciclastic_Matrix_SubPeriod
for (u in 1:nrow(Siliciclastic_Matrix_SubPeriod_Plot)){
  for (x in 1:ncol(Siliciclastic_Matrix_SubPeriod_Plot)){
    if (Siliciclastic_Matrix_SubPeriod_Plot[u,x]==0){
      Siliciclastic_Matrix_SubPeriod_Plot[u,x]<-NA
    }
  }
}
geoscalePlot(c(341.05, 311.05), Siliciclastic_Matrix_SubPeriod_Plot[1,], units=c( "Epoch"), type="b", lwd=2, data.lim=c(0, 50), label="Number of Occurrences", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(Carbonate_Matrix_SubPeriod_Plot)){
  lines(c(341.05, 311.05), Siliciclastic_Matrix_SubPeriod_Plot[x,], type="b", pch=x, col=x, lwd=2)
}
legend(320, 50, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)

#Total number of occurrences used in carb vs. silic comparisons at the subperiod level
Tot_Subperiod_Carb_Silic<-Siliciclastic_Matrix_SubPeriod+Carbonate_Matrix_SubPeriod

#Total Number of Fine Grained Occurrences in Each Subperiod
Fine_Grained<-subset(New, New[,11] =="Fine_Grained")
Fine_Matrix_SubPeriod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
rownames(Fine_Matrix_SubPeriod)<-Taxon_Names
colnames(Fine_Matrix_SubPeriod)<-SubPeriods
for (k in 1:length(Taxon_Names)){
  for (p in 1:length(SubPeriods)){
    Fine_Taxon_Time<-subset(Fine_Grained, Fine_Grained[,4]==Taxon_Names[k] & Fine_Grained[,19]==SubPeriods[p])
    Fine_Matrix_SubPeriod[k,p]<-nrow(Fine_Taxon_Time)
  }
}
Fine_Matrix_SubPeriod_Plot<-Fine_Matrix_SubPeriod
for (u in 1:nrow(Fine_Matrix_SubPeriod_Plot)){
  for (x in 1:ncol(Fine_Matrix_SubPeriod_Plot)){
    if (Fine_Matrix_SubPeriod_Plot[u,x]==0){
      Fine_Matrix_SubPeriod_Plot[u,x]<-NA
    }
  }
}
geoscalePlot(c(341.05, 311.05), Fine_Matrix_SubPeriod_Plot[1,], units=c( "Epoch"), type="b", lwd=2, data.lim=c(0, 150), label="Number of Occurrences", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(Fine_Matrix_SubPeriod)){
  lines(c(341.05, 311.05), Fine_Matrix_SubPeriod_Plot[x,], type="b", pch=x, col=x, lwd=2)
}
legend(320, 150, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)

#Total Number of Coarse Grained Occurrences in Each Subperiod
Coarse_Grained<-subset(New, New[,11] =="Coarse_Grained")
Coarse_Matrix_SubPeriod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
rownames(Coarse_Matrix_SubPeriod)<-Taxon_Names
colnames(Coarse_Matrix_SubPeriod)<-SubPeriods
for (k in 1:length(Taxon_Names)){
  for (p in 1:length(SubPeriods)){
    Coarse_Taxon_Time<-subset(Coarse_Grained, Coarse_Grained[,4]==Taxon_Names[k] & Coarse_Grained[,19]==SubPeriods[p])
    Coarse_Matrix_SubPeriod[k,p]<-nrow(Coarse_Taxon_Time)
  }
}
Coarse_Matrix_SubPeriod_Plot<-Coarse_Matrix_SubPeriod
for (u in 1:nrow(Coarse_Matrix_SubPeriod_Plot)){
  for (x in 1:ncol(Coarse_Matrix_SubPeriod_Plot)){
    if (Coarse_Matrix_SubPeriod_Plot[u,x]==0){
      Coarse_Matrix_SubPeriod_Plot[u,x]<-NA
    }
  }
}
geoscalePlot(c(341.05, 311.05), Coarse_Matrix_SubPeriod_Plot[1,], units=c( "Epoch"), type="b", lwd=2, data.lim=c(0, 50), label="Number of Occurrences", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(Coarse_Matrix_SubPeriod_Plot)){
  lines(c(341.05, 311.05), Coarse_Matrix_SubPeriod_Plot[x,], type="b", pch=x, col=x, lwd=2)
}
legend(320, 50, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)

#Total number of occurrences used in fine vs. coarse comparisons at the subperiod level
Tot_Subperiod_Fine_Coarse<-Fine_Matrix_SubPeriod+Coarse_Matrix_SubPeriod

#Entire Carboniferous
#Populate matrices of all occurrences of each family.  
Taxon_Names<-c("Palaechinidae", "Archaeocidaridae", "Proterocidaridae", "Lepidesthidae", "Lepidocentridae")
for (p in 1:length(Taxon_Names)){
  Taxon_Matrix<-subset(New, New[,4]==Taxon_Names[p])
  assign(Taxon_Names[p], Taxon_Matrix)
}

#Creates Matrices for ocurrences of each family in Fine grained and coarse grained environments 
Grain_Size_Names<-c("Fine_Grained", "Coarse_Grained")
for (p in 1:length(Taxon_Names)){
  for (q in 1:length(Grain_Size_Names)){
    Taxon_Grain_Size<-subset(New, New[,4]==Taxon_Names[p] & New[,11]==Grain_Size_Names[q])
    assign(paste(Taxon_Names[p], Grain_Size_Names[q], sep="_"), Taxon_Grain_Size)
    Number<-nrow(Taxon_Grain_Size)
    assign(paste("N", Taxon_Names[p], Grain_Size_Names[q], sep="_"), Number)
  }
}

#Creates Matrices for ocurrences of each family in Carbonate and siliciclastic environments  
Carb_Silic_Names<-c("Carbonate", "Siliciclastic")
for (p in 1:length(Taxon_Names)){
  for (q in 1:length(Carb_Silic_Names)){
    Taxon_Carb_Silic<-subset(New, New[,4]==Taxon_Names[p] & New[,10]==Carb_Silic_Names[q])
    assign(paste(Taxon_Names[p], Carb_Silic_Names[q], sep="_"), Taxon_Carb_Silic)
    Number<-nrow(Taxon_Carb_Silic)
    assign(paste("N", Taxon_Names[p], Carb_Silic_Names[q], sep="_"), Number)
  }
}

#Creates vectors with number of occurrences of each clade in carbonate/siliciclastic and fine/coarse environments across entire Carboniferous. Used for total carboniferous calculations.
Numbers_Carb<-c(N_Palaechinidae_Carbonate, N_Archaeocidaridae_Carbonate, N_Proterocidaridae_Carbonate, N_Lepidesthidae_Carbonate, N_Lepidocentridae_Carbonate)
geoscalePlot(c(328.9), Numbers_Carb[1], units=c( "Period"), type="b", lwd=2, data.lim=c(0, 200), label="Number of Occurrences", erotate=0, direction="horizontal", abbrev=c("Subperiod"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:length(Numbers_Carb)){
  lines(c(328.9), Numbers_Carb[x], type="b", pch=x, col=x, lwd=2)
}
legend(320, 200, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)

Numbers_Silic<-c(N_Palaechinidae_Siliciclastic, N_Archaeocidaridae_Siliciclastic, N_Proterocidaridae_Siliciclastic, N_Lepidesthidae_Siliciclastic, N_Lepidocentridae_Siliciclastic)
geoscalePlot(c(328.9), Numbers_Silic[1], units=c( "Period"), type="b", lwd=2, data.lim=c(0, 50), label="Number of Occurrences", erotate=0, direction="horizontal", abbrev=c("Subperiod"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:length(Numbers_Silic)){
  lines(c(328.9), Numbers_Silic[x], type="b", pch=x, col=x, lwd=2)
}
legend(320, 50, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)

Numbers_Fine<-c(N_Palaechinidae_Fine_Grained, N_Archaeocidaridae_Fine_Grained, N_Proterocidaridae_Fine_Grained, N_Lepidesthidae_Fine_Grained, N_Lepidocentridae_Fine_Grained)
geoscalePlot(c(328.9), Numbers_Fine[1], units=c( "Period"), type="b", lwd=2, data.lim=c(0, 200), label="Number of Occurrences", erotate=0, direction="horizontal", abbrev=c("Subperiod"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:length(Numbers_Fine)){
  lines(c(328.9), Numbers_Fine[x], type="b", pch=x, col=x, lwd=2)
}
legend(320, 200, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)

Numbers_Coarse<-c(N_Palaechinidae_Coarse_Grained, N_Archaeocidaridae_Coarse_Grained, N_Proterocidaridae_Coarse_Grained, N_Lepidesthidae_Coarse_Grained, N_Lepidocentridae_Coarse_Grained)
geoscalePlot(c(328.9), Numbers_Coarse[1], units=c( "Period"), type="b", lwd=2, data.lim=c(0, 60), label="Number of Occurrences", erotate=0, direction="horizontal", abbrev=c("Subperiod"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:length(Numbers_Coarse)){
  lines(c(328.9), Numbers_Coarse[x], type="b", pch=x, col=x, lwd=2)
}
legend(320, 60, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)

Total_Fine<-sum(Numbers_Fine)
Total_Coarse<-sum(Numbers_Coarse)

Total_Carb<-sum(Numbers_Carb) #Remove?
Total_Silic<-sum(Numbers_Silic)#Remove?
Total_Proportion_Silic_Carb<-Total_Silic/(Total_Silic+Total_Carb)#Remove?
Total_Proportion_Carb_Silic<-Total_Carb/(Total_Silic+Total_Carb)#Remove?
Total_Proportion_Fine_Coarse<-Total_Fine/(Total_Fine+Total_Coarse)
Total_Proportion_Coarse_Fine<-Total_Coarse/(Total_Fine+Total_Coarse)

#Foote (2006, 2014) Calculations
#Stage Level Foote (2006, 2014) Affinities
Foote_Stage_Affinities_Carb_Silic<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
Foote_Stage_Affinities_Silic_Carb<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
rownames(Foote_Stage_Affinities_Carb_Silic)<-Taxon_Names
colnames(Foote_Stage_Affinities_Carb_Silic)<-TimeBins
rownames(Foote_Stage_Affinities_Silic_Carb)<-Taxon_Names
colnames(Foote_Stage_Affinities_Silic_Carb)<-TimeBins

for (u in 1:nrow(Carbonate_Matrix)){
  for (x in 1:ncol(Carbonate_Matrix)){
  Binom_Test_Carb_Silic<-binom.test(c(Carbonate_Matrix[u,x], Siliciclastic_Matrix[u,x]), p=Proportions_Stage_Carb_Silic[x], alternative=c("greater"))
  barplot(c(Carbonate_Matrix[u,x], Siliciclastic_Matrix[u,x]), main=c(Taxon_Names[u]), names=c("Carbonate", "Siliciclastic"), axes=TRUE, ylab="Number of Occurrences", xlab="Substrate")
  Foote_Stage_Affinities_Carb_Silic[u,x]<-Binom_Test_Carb_Silic$p.value
  if (Tot_Stage_Carb_Silic[u,x]<5){
    Foote_Stage_Affinities_Carb_Silic[u,x]<-"-"
  }
  write.xlsx(Foote_Stage_Affinities_Carb_Silic, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_Foote/Foote_Stage_Affinities_Carb_Silic.xlsx")
  
  Binom_Test_Silic_Carb<-binom.test(c(Siliciclastic_Matrix[u,x], Carbonate_Matrix[u,x]), p=Proportions_Stage_Silic_Carb[x], alternative=c("greater"))
  Foote_Stage_Affinities_Silic_Carb[u,x]<-Binom_Test_Silic_Carb$p.value
  if (Tot_Stage_Carb_Silic[u,x]<5){
    Foote_Stage_Affinities_Silic_Carb[u,x]<-"-"
  }
  write.xlsx(Foote_Stage_Affinities_Silic_Carb, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_Foote/Foote_Stage_Affinities_Silic_Carb.xlsx")
  }
}

#Subperiod Level Foote (2006, 2014) Affinities
Foote_Subperiod_Affinities_Carb_Silic<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
Foote_Subperiod_Affinities_Silic_Carb<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
rownames(Foote_Subperiod_Affinities_Carb_Silic)<-Taxon_Names
colnames(Foote_Subperiod_Affinities_Carb_Silic)<-SubPeriods
rownames(Foote_Subperiod_Affinities_Silic_Carb)<-Taxon_Names
colnames(Foote_Subperiod_Affinities_Silic_Carb)<-SubPeriods


for (u in 1:nrow(Carbonate_Matrix_SubPeriod)){
  for (x in 1:ncol(Carbonate_Matrix_SubPeriod)){
    Binom_Test_Carb_Silic<-binom.test(c(Carbonate_Matrix_SubPeriod[u,x], Siliciclastic_Matrix_SubPeriod[u,x]), p=Proportions_Subperiod_Carb_Silic[x], alternative=c("greater"))
    barplot(c(Carbonate_Matrix_SubPeriod[u,x], Siliciclastic_Matrix_SubPeriod[u,x]), main=c(Taxon_Names[u]), names=c("Carbonate", "Siliciclastic"), axes=TRUE, ylab="Number of Occurrences", xlab="Substrate")
    Foote_Subperiod_Affinities_Carb_Silic[u,x]<-Binom_Test_Carb_Silic$p.value
    if (Tot_Subperiod_Carb_Silic[u,x]<5){
      Foote_Subperiod_Affinities_Carb_Silic[u,x]<-"-"
    }
    write.xlsx(Foote_Subperiod_Affinities_Carb_Silic, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_Foote/Foote_Subperiod_Affinities_Carb_Silic.xlsx")
    
    Binom_Test_Silic_Carb<-binom.test(c(Siliciclastic_Matrix_SubPeriod[u,x], Carbonate_Matrix_SubPeriod[u,x]), p=Proportions_Subperiod_Silic_Carb[x], alternative=c("greater"))
    Foote_Subperiod_Affinities_Silic_Carb[u,x]<-Binom_Test_Silic_Carb$p.value
    if (Tot_Subperiod_Carb_Silic[u,x]<5){
      Foote_Subperiod_Affinities_Silic_Carb[u,x]<-"-"
    }
    write.xlsx(Foote_Subperiod_Affinities_Silic_Carb, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_Foote/Foote_Subperiod_Affinities_Silic_Carb.xlsx")
 }
}

#Entire Carboniferous Foote Affinities
Foote_Affinities_Carb_Silic_Carboniferous<-matrix(nrow=length(Taxon_Names), ncol=1, 0)
rownames(Foote_Affinities_Carb_Silic_Carboniferous)<-Taxon_Names
Foote_Affinities_Silic_Carb_Carboniferous<-matrix(nrow=length(Taxon_Names), ncol=1, 0)
rownames(Foote_Affinities_Silic_Carb_Carboniferous)<-Taxon_Names
Foote_Affinities_Coarse_Fine_Carboniferous<-matrix(nrow=length(Taxon_Names), ncol=1, 0)

for (x in 1:length(Numbers_Carb)){
  Binom_Test_Carb_Silic<-binom.test(c(Numbers_Carb[x], Numbers_Silic[x]), p=Proportion_Carb_Silic_Tot, alternative=c("greater"))
  barplot(c(Numbers_Carb[x], Numbers_Silic[x]), main=c(Taxon_Names[x]), names=c("Carbonate", "Siliciclastic"), axes=TRUE, ylab="Number of Occurrences", xlab="Substrate")
  Binom_Test_Silic_Carb<-binom.test(c(Numbers_Silic[x], Numbers_Carb[x]), p=Proportion_Silic_Carb_Tot, alternative=c("greater"))

  Foote_Affinities_Carb_Silic_Carboniferous[x]<-Binom_Test_Carb_Silic$p.value
  write.xlsx(Foote_Affinities_Carb_Silic_Carboniferous, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_Foote/Foote_Affinities_Carb_Silic_Carboniferous.xlsx")
  Foote_Affinities_Silic_Carb_Carboniferous[x]<-Binom_Test_Silic_Carb$p.value
  write.xlsx(Foote_Affinities_Silic_Carb_Carboniferous, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_Foote/Foote_Affinities_Silic_Carb_Carboniferous.xlsx")
  }

#Simpson and Harnik (2009) Bayesian affinity Calculations. Variable names are from Simpson and Harnik (2009)
#Stage Level Bayesian. Calculates PP's and plots PP of Carb at the stage level
Bayesian_Stage_Affinities_pH1E<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
Bayesian_Stage_Affinities_pH2E_second<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
rownames(Bayesian_Stage_Affinities_pH1E)<-Taxon_Names
colnames(Bayesian_Stage_Affinities_pH1E)<-TimeBins

for (u in 1:nrow(Carbonate_Matrix)){
  for (x in 1:ncol(Carbonate_Matrix)){
    n<-Carbonate_Matrix[u,x]+Siliciclastic_Matrix[u,x]
    pcarb<-Proportions_Stage_Carb_Silic[x]
    psilic<-Proportions_Stage_Silic_Carb[x]
    pH1=0.5
    pH2=0.5
    pEH1=choose(n,Carbonate_Matrix[u,x])*(pcarb^Carbonate_Matrix[u,x])*(1-pcarb)^(n-Carbonate_Matrix[u,x]) 
    pEH2=choose(n,Siliciclastic_Matrix[u,x])*(pcarb^Siliciclastic_Matrix[u,x])*(1-pcarb)^(n-Siliciclastic_Matrix[u,x])
    pH1E=(pEH1*pH1)/((pEH1*pH1)+(pEH2*pH2))
    Bayesian_Stage_Affinities_pH1E[u,x]<-pH1E
    if (Tot_Stage_Carb_Silic[u,x]<5){
     Bayesian_Stage_Affinities_pH1E[u,x]<-"-"
    }
    write.xlsx(Bayesian_Stage_Affinities_pH1E, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_Bayes/Bayesian_Affinities_Carb_Silic_Stages.xlsx")
    if (Proportions_Stage_Carb_Silic[x]<0.5){
      Bayesian_Stage_Affinities_pH1E[,x]<-"-"
    }
}
}

#Calculation of pH2E For Kasimovian
for (u in 1:nrow(Carbonate_Matrix)){
  for (x in 1:ncol(Carbonate_Matrix)){
    n<-Carbonate_Matrix[u,x]+Siliciclastic_Matrix[u,x]
    pcarb<-Proportions_Stage_Carb_Silic[x]
    psilic<-Proportions_Stage_Silic_Carb[x]
    pH1=0.5
    pH2=0.5
    pEH1=choose(n,Carbonate_Matrix[u,x])*(psilic^Carbonate_Matrix[u,x])*(1-psilic)^(n-Carbonate_Matrix[u,x]) 
    pEH2=choose(n,Siliciclastic_Matrix[u,x])*(psilic^Siliciclastic_Matrix[u,x])*(1-psilic)^(n-Siliciclastic_Matrix[u,x])
    pH2E=(pEH2*pH2)/((pEH1*pH1)+(pEH2*pH2))
  Bayesian_Stage_Affinities_pH2E_second[u,x]<-pH2E
    if (Tot_Stage_Carb_Silic[u,x]<5){
      Bayesian_Stage_Affinities_pH2E_second[u,x]<-"-"
    }
    else if (Proportions_Stage_Silic_Carb[x]<0.5){
      Bayesian_Stage_Affinities_pH2E_second[,x]<-"-"
    }
  }
}
Ones<-matrix(nrow=5, ncol=1, 1)
Bayesian_Stage_Affinities_pH1E[2:3,6]<-(Ones[2:3]-as.numeric(Bayesian_Stage_Affinities_pH2E_second[2:3,6]))
write.xlsx(Bayesian_Stage_Affinities_pH1E, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_Bayes/Bayesian_Affinities_Carb_Silic_Stages.xlsx")

geoscalePlot(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Bayesian_Stage_Affinities_pH1E[1,], units=c("Age", "Epoch"), type="b", lwd=2, data.lim=c(0, 1), label="Posterior Probability", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(Bayesian_Stage_Affinities_pH1E)){
  lines(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Bayesian_Stage_Affinities_pH1E[x,], type="b", pch=x, col=x, lwd=2)
}
legend(357, .35, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)
abline(.5,0, lty=2)

#Subperiod Level Bayesian Silic Carb. Calculates PP's and plots PP of Carb affinity at the subperiod level
Bayesian_Subperiod_Affinities_pH1E<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
Bayesian_Subperiod_Affinities_pH2E<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
rownames(Bayesian_Subperiod_Affinities_pH1E)<-Taxon_Names
colnames(Bayesian_Subperiod_Affinities_pH1E)<-TimeBins

for (u in 1:nrow(Carbonate_Matrix_SubPeriod)){
  for (x in 1:ncol(Carbonate_Matrix_SubPeriod)){
    n<-Carbonate_Matrix_SubPeriod[u,x]+Siliciclastic_Matrix_SubPeriod[u,x]
    pcarb<-Proportions_Subperiod_Carb_Silic[x]
    psilic<-Proportions_Subperiod_Silic_Carb[x]
    pH1=0.5
    pH2=0.5
    pEH1=choose(n,Carbonate_Matrix_SubPeriod[u,x])*(pcarb^Carbonate_Matrix_SubPeriod[u,x])*(1-pcarb)^(n-Carbonate_Matrix_SubPeriod[u,x]) 
    pEH2=choose(n,Siliciclastic_Matrix_SubPeriod[u,x])*(pcarb^Siliciclastic_Matrix_SubPeriod[u,x])*(1-pcarb)^(n-Siliciclastic_Matrix_SubPeriod[u,x])
    pH1E=(pEH1*pH1)/((pEH1*pH1)+(pEH2*pH2))
    Bayesian_Subperiod_Affinities_pH1E[u,x]<-pH1E
   if (Tot_Subperiod_Carb_Silic[u,x]<5){
     Bayesian_Subperiod_Affinities_pH1E[u,x]<-"-"
   }
    write.xlsx(Bayesian_Subperiod_Affinities_pH1E, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_Bayes/Bayesian_Affinities_Carb_Silic_Subperiod.xlsx")
   }
 }
geoscalePlot(c(341.05, 311.05), Bayesian_Subperiod_Affinities_pH1E[1,], units=c( "Epoch"), type="b", lwd=2, data.lim=c(0, 1), label="Posterior Probability", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(Bayesian_Subperiod_Affinities_pH1E)){
  lines(c(341.05, 311.05), Bayesian_Subperiod_Affinities_pH1E[x,], type="b", pch=x, col=x, lwd=2)
}
legend(357, .35, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)
abline(.5,0, lty=2)

#Period Level Bayesian Affinities
Bayesian_Affinities_Carb_Silic_Carboniferous<-matrix(nrow=length(Taxon_Names), ncol=1, 0)
for (x in 1:length(Numbers_Carb)){
n<-Numbers_Carb[x]+Numbers_Silic[x]
pcarb<-Proportion_Carb_Silic_Tot #Total Proportion of carbonate occurrences from all collections   
psilic<-Proportion_Silic_Carb_Tot #Total Proportion of siliciclastic occurrences from all collections 
pH1=0.5
pH2=0.5
pEH1=choose(n,Numbers_Carb[x])*(pcarb^Numbers_Carb[x])*(1-pcarb)^(n-Numbers_Carb[x]) 
pEH2=choose(n,Numbers_Silic[x])*(pcarb^Numbers_Silic[x])*(1-pcarb)^(n-Numbers_Silic[x])
pH1E=(pEH1*pH1)/((pEH1*pH1)+(pEH2*pH2))
Bayesian_Affinities_Carb_Silic_Carboniferous[x]<-pH1E
write.xlsx(Bayesian_Affinities_Carb_Silic_Carboniferous, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_Bayes/Bayesian_Affinities_Carb_Silic_Carboniferous.xlsx")
}
geoscalePlot(c(328.9), Bayesian_Affinities_Carb_Silic_Carboniferous[1], units=c( "Period"), type="b", lwd=2, data.lim=c(0, 1), label="Posterior Probability", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(Bayesian_Affinities_Carb_Silic_Carboniferous)){
  lines(c(328.9), Bayesian_Affinities_Carb_Silic_Carboniferous[x], type="b", pch=x, col=x, lwd=2)
}
legend(357, .95, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)
abline(.5,0, lty=2)

#Miller and Connolly (2001) Standardized Relative Affinity Calculations
# Miller and Connolly, 2001 Stage Level Carb Silic
Ac_Stage<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
As_Stage<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
rownames(Ac_Stage)<-Taxon_Names
colnames(Ac_Stage)<-TimeBins
rownames(As_Stage)<-Taxon_Names
colnames(As_Stage)<-TimeBins

for (u in 1:nrow(Carbonate_Matrix)){
  for (x in 1:ncol(Carbonate_Matrix)){
    Ac_Stage[u,x]<-Carbonate_Matrix[u,x]/(Carbonate_Matrix[u,x]+Siliciclastic_Matrix[u,x])
    As_Stage[u,x]<-Siliciclastic_Matrix[u,x]/(Carbonate_Matrix[u,x]+Siliciclastic_Matrix[u,x])
  }
}

Ac_Stage_Plot<-Ac_Stage
for (u in 1:nrow(Ac_Stage_Plot)){
  for (x in 1:ncol(Ac_Stage_Plot)){
    if (Ac_Stage_Plot[u,x]=="0"){
      Ac_Stage_Plot[u,x]<-NA
    }
  }
}

geoscalePlot(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Ac_Stage_Plot[1,], units=c( "Age", "Epoch"), type="b", lwd=2, data.lim=c(0, 1), label=expression(paste("A"["c"])), erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(Ac_Stage_Plot)){
  lines(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Ac_Stage_Plot[x,], type="b", pch=x, col=x, lwd=2)
}
legend(320, 70, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)



#Random Affinities

RNDc_Final_Stage<-array(0, dim=c(length(Taxon_Names), ncol=length(TimeBins), 1000))
dimnames(RNDc_Final_Stage)[[1]]<-Taxon_Names
dimnames(RNDc_Final_Stage)[[2]]<-TimeBins

RNDs_Final_Stage<-array(0, dim=c(length(Taxon_Names), ncol=length(TimeBins), 1000))
dimnames(RNDs_Final_Stage)[[1]]<-Taxon_Names
dimnames(RNDs_Final_Stage)[[2]]<-TimeBins

RNDc_Final_SD_Stage<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
RNDs_Final_SD_Stage<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
rownames(RNDc_Final_SD_Stage)<-Taxon_Names
colnames(RNDs_Final_SD_Stage)<-TimeBins
rownames(RNDc_Final_SD_Stage)<-Taxon_Names
colnames(RNDs_Final_SD_Stage)<-TimeBins

Atotc<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
Atots<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
rownames(Atotc)<-Taxon_Names
colnames(Atots)<-TimeBins
rownames(Atotc)<-Taxon_Names
colnames(Atots)<-TimeBins

SRA_Carb_Stage<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
SRA_Silic_Stage<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
rownames(SRA_Carb_Stage)<-Taxon_Names
colnames(SRA_Carb_Stage)<-TimeBins
rownames(SRA_Silic_Stage)<-Taxon_Names
colnames(SRA_Silic_Stage)<-TimeBins

    for (u in 1:nrow(Carbonate_Matrix)){
      for (x in 1:ncol(Carbonate_Matrix)){
        for (p in 1:1000){
        
      if (Carbonate_Matrix[u,x]==0 & Siliciclastic_Matrix[u,x]==0) next
         CarbSilic<-subset(New, New[,8]==TimeBins[x])
         CarbSilic<-subset(CarbSilic, CarbSilic[,10]=="Carbonate"| CarbSilic[,10]=="Siliciclastic")
         Rand_Sample_Carb<-CarbSilic[sample(nrow(CarbSilic), size=(Carbonate_Matrix[u,x]+Siliciclastic_Matrix[u,x]), replace=TRUE),]
      
          Rand_Sample_Carb<-as.matrix(Rand_Sample_Carb)
          if (ncol(Rand_Sample_Carb)==1) Rand_Sample_Carb<-t(Rand_Sample_Carb)
      
          RNDc<-matrix(0, nrow = (Carbonate_Matrix[u,x]+Siliciclastic_Matrix[u,x]), ncol = ncol(New))
          RNDs<-matrix(0, nrow = (Carbonate_Matrix[u,x]+Siliciclastic_Matrix[u,x]), ncol = ncol(New))     

          for (n in 1:(Carbonate_Matrix[u,x]+Siliciclastic_Matrix[u,x])){
            if (Rand_Sample_Carb[n,10]=="Carbonate"){
              RNDc[n,]<-as.matrix(Rand_Sample_Carb[n,])
            }
          }
          
          for (n in 1:(Carbonate_Matrix[u,x]+Siliciclastic_Matrix[u,x])){
            if (Rand_Sample_Carb[n,10]=="Siliciclastic"){
              RNDs[n,]<-as.matrix(Rand_Sample_Carb[n,])
            }
          }
       
RNDc<-RNDc[apply(RNDc!=0,1, any),,drop=FALSE] #This and the below line of code just remove 0 rows from the RND affinities matrices
RNDs<-RNDs[apply(RNDs!=0,1, any),,drop=FALSE]

Ac_Rand<-nrow(RNDc)/(nrow(RNDc)+nrow(RNDs))
As_Rand<-nrow(RNDs)/(nrow(RNDc)+nrow(RNDs))

RNDc_Final_Stage[u,x,p]<-Ac_Rand
RNDs_Final_Stage[u,x,p]<-As_Rand
}

RNDc_Final_SD_Stage[u,x]<-sd(RNDc_Final_Stage[u,x,])
RNDs_Final_SD_Stage[u,x]<-sd(RNDs_Final_Stage[u,x,])

Atotc[u,x]<-mean(RNDc_Final_Stage[u,x,])
Atots[u,x]<-mean(RNDs_Final_Stage[u,x,])

SRA_Carb_Stage[u,x]<-(Ac_Stage[u,x]-Atotc[u,x])/RNDc_Final_SD_Stage[u,x]
if (Tot_Stage_Carb_Silic[u,x]<5){
  SRA_Carb_Stage[u,x]<-"-"
}
write.xlsx(SRA_Carb_Stage, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_SRA/SRA_Carb_Stage.xlsx")
SRA_Silic_Stage[u,x]<-(As_Stage[u,x]-Atots[u,x])/RNDs_Final_SD_Stage[u,x]
if (Tot_Stage_Carb_Silic[u,x]<5){
  SRA_Silic_Stage[u,x]<-"-"
}
write.xlsx(SRA_Silic_Stage, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_SRA/SRA_Silic_Stage.xlsx")
}
}
geoscalePlot(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), SRA_Carb_Stage[1,], units=c("Age", "Epoch"), type="b", lwd=2, data.lim=c(-3, 2), label="Standardized Relative Affinity", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,305.75), cex.ts=.9)
for (x in 2:nrow(SRA_Carb_Stage)){
  lines(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), SRA_Carb_Stage[x,], type="b", pch=x, col=x, lwd=2)
}
legend(357, -1.4, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)
abline(0,0, lty=2)
      
geoscalePlot(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), SRA_Silic_Stage[1,], units=c("Age", "Epoch"), type="b", lwd=2, data.lim=c(-2, 3), label="Standardized Relative Affinity", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(SRA_Silic_Stage)){
  lines(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), SRA_Silic_Stage[x,], type="b", pch=x, col=x, lwd=2)
}
legend(357, 3, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)
abline(0,0, lty=2)

# Miller and Connolly, 2001 Carb Silic Subperiod level
Ac_Subperiod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
As_Subperiod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
rownames(Ac_Subperiod)<-Taxon_Names
colnames(Ac_Subperiod)<-SubPeriods
rownames(As_Subperiod)<-Taxon_Names
colnames(As_Subperiod)<-SubPeriods

for (u in 1:nrow(Carbonate_Matrix_SubPeriod)){
  for (x in 1:ncol(Carbonate_Matrix_SubPeriod)){
    Ac_Subperiod[u,x]<-Carbonate_Matrix_SubPeriod[u,x]/(Carbonate_Matrix_SubPeriod[u,x]+Siliciclastic_Matrix_SubPeriod[u,x])
    As_Subperiod[u,x]<-Siliciclastic_Matrix_SubPeriod[u,x]/(Carbonate_Matrix_SubPeriod[u,x]+Siliciclastic_Matrix_SubPeriod[u,x])
  }
}


Ac_Subperiod_Plot<-Ac_Subperiod
for (u in 1:nrow(Ac_Subperiod_Plot)){
  for (x in 1:ncol(Ac_Subperiod_Plot)){
    if (Ac_Subperiod_Plot[u,x]=="0"){
      Ac_Subperiod_Plot[u,x]<-NA
    }
  }
}

geoscalePlot(c(341.05, 311.05), Ac_Subperiod_Plot[1,], units=c( "Epoch"), type="b", lwd=2, data.lim=c(0, 1), label=expression(paste("A"["c"])), erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(Ac_Subperiod_Plot)){
  lines(c(341.05, 311.05), Ac_Subperiod_Plot[x,], type="b", pch=x, col=x, lwd=2)
}
legend(320, 50, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)



#Random Affinities Carb Silic Subperiod Level
RNDc_Final_Subperiod<-array(0, dim=c(length(Taxon_Names), ncol=length(SubPeriods), 1000))
dimnames(RNDc_Final_Subperiod)[[1]]<-Taxon_Names
dimnames(RNDc_Final_Subperiod)[[2]]<-SubPeriods

RNDs_Final_Subperiod<-array(0, dim=c(length(Taxon_Names), ncol=length(SubPeriods), 1000))
dimnames(RNDs_Final_Subperiod)[[1]]<-Taxon_Names
dimnames(RNDs_Final_Subperiod)[[2]]<-SubPeriods

RNDc_Final_SD_Subperiod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
RNDs_Final_SD_Subperiod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
rownames(RNDc_Final_SD_Subperiod)<-Taxon_Names
colnames(RNDs_Final_SD_Subperiod)<-SubPeriods
rownames(RNDc_Final_SD_Subperiod)<-Taxon_Names
colnames(RNDs_Final_SD_Subperiod)<-SubPeriods

Atotc_Subperiod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
Atots_Subperiod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
rownames(Atotc_Subperiod)<-Taxon_Names
colnames(Atotc_Subperiod)<-SubPeriods
rownames(Atots_Subperiod)<-Taxon_Names
colnames(Atots_Subperiod)<-SubPeriods

SRA_Carb_Subperiod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
SRA_Silic_Subperiod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
rownames(SRA_Carb_Subperiod)<-Taxon_Names
colnames(SRA_Carb_Subperiod)<-SubPeriods
rownames(SRA_Silic_Subperiod)<-Taxon_Names
colnames(SRA_Silic_Subperiod)<-SubPeriods

for (u in 1:nrow(Carbonate_Matrix_SubPeriod)){
  for (x in 1:ncol(Carbonate_Matrix_SubPeriod)){
    for (p in 1:1000){
      
      if (Carbonate_Matrix_SubPeriod[u,x]==0 & Siliciclastic_Matrix_SubPeriod[u,x]==0) next
      CarbSilic<-subset(New, New[,19]==SubPeriods[x])
      CarbSilic<-subset(CarbSilic, CarbSilic[,10]=="Carbonate"| CarbSilic[,10]=="Siliciclastic")
      Rand_Sample_Carb<-CarbSilic[sample(nrow(CarbSilic), size=(Carbonate_Matrix_SubPeriod[u,x]+Siliciclastic_Matrix_SubPeriod[u,x]), replace=TRUE),]
      
      
      Rand_Sample_Carb<-as.matrix(Rand_Sample_Carb)
      if (ncol(Rand_Sample_Carb)==1) Rand_Sample_Carb<-t(Rand_Sample_Carb)
      
      RNDc<-matrix(0, nrow = (Carbonate_Matrix_SubPeriod[u,x]+Siliciclastic_Matrix_SubPeriod[u,x]), ncol = ncol(New))
      RNDs<-matrix(0, nrow = (Carbonate_Matrix_SubPeriod[u,x]+Siliciclastic_Matrix_SubPeriod[u,x]), ncol = ncol(New))     
      
      for (n in 1:(Carbonate_Matrix_SubPeriod[u,x]+Siliciclastic_Matrix_SubPeriod[u,x])){
        if (Rand_Sample_Carb[n,10]=="Carbonate"){
          RNDc[n,]<-as.matrix(Rand_Sample_Carb[n,])
        }
      }
      
      for (n in 1:(Carbonate_Matrix_SubPeriod[u,x]+Siliciclastic_Matrix_SubPeriod[u,x])){
        if (Rand_Sample_Carb[n,10]=="Siliciclastic"){ #Altered this so both pull from 1 sample
          RNDs[n,]<-as.matrix(Rand_Sample_Carb[n,])
        }
      }
      
      
      RNDc<-RNDc[apply(RNDc!=0,1, any),,drop=FALSE]
      RNDs<-RNDs[apply(RNDs!=0,1, any),,drop=FALSE]
      
      Ac_Rand<-nrow(RNDc)/(nrow(RNDc)+nrow(RNDs))
      As_Rand<-nrow(RNDs)/(nrow(RNDc)+nrow(RNDs))
      
      RNDc_Final_Subperiod[u,x,p]<-Ac_Rand
      RNDs_Final_Subperiod[u,x,p]<-As_Rand
      
      
    }
    RNDc_Final_SD_Subperiod[u,x]<-sd(RNDc_Final_Subperiod[u,x,])
    RNDs_Final_SD_Subperiod[u,x]<-sd(RNDs_Final_Subperiod[u,x,])
    
    Atotc_Subperiod[u,x]<-mean(RNDc_Final_Subperiod[u,x,])
    Atots_Subperiod[u,x]<-mean(RNDs_Final_Subperiod[u,x,])
    
    SRA_Carb_Subperiod[u,x]<-(Ac_Subperiod[u,x]-Atotc_Subperiod[u,x])/RNDc_Final_SD_Subperiod[u,x]
    if (Tot_Subperiod_Carb_Silic[u,x]<5){
      SRA_Carb_Subperiod[u,x]<-"-"
    }
    write.xlsx(SRA_Carb_Subperiod, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_SRA/SRA_Carb_Subperiod.xlsx")
    SRA_Silic_Subperiod[u,x]<-(As_Subperiod[u,x]-Atots_Subperiod[u,x])/RNDs_Final_SD_Subperiod[u,x]
    if (Tot_Subperiod_Carb_Silic[u,x]<5){
      SRA_Silic_Subperiod[u,x]<-"-"
    }
    write.xlsx(SRA_Silic_Subperiod, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_SRA/SRA_Silic_Subperiod.xlsx")
  }
}

geoscalePlot(c(341.05, 311.05), SRA_Carb_Subperiod[1,], units=c("Epoch"), type="b", lwd=2, data.lim=c(-3, 2), label="Standardized Relative Affinity", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(SRA_Carb_Subperiod)){
  lines(c(341.05, 311.05), SRA_Carb_Subperiod[x,], type="b", pch=x, col=x, lwd=2)
}
legend(320, -1.5, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)
abline(0,0, lty=2)

geoscalePlot(c(341.05, 311.05), SRA_Silic_Subperiod[1,], units=c( "Epoch"), type="b", lwd=2, data.lim=c(-2, 3), label="Standardized Relative Affinity", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(SRA_Silic_Subperiod)){
  lines(c(341.05, 311.05), SRA_Silic_Subperiod[x,], type="b", pch=x, col=x, lwd=2)
}
legend(320, 3, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)
abline(0,0, lty=2)

# Miller and Connolly, 2001 Stage Level Coarse Fine
Afine_Stage<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
Acoarse_Stage<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
rownames(Afine_Stage)<-Taxon_Names
colnames(Afine_Stage)<-TimeBins
rownames(Acoarse_Stage)<-Taxon_Names
colnames(Acoarse_Stage)<-TimeBins

for (u in 1:nrow(Fine_Matrix)){
  for (x in 1:ncol(Fine_Matrix)){
    Afine_Stage[u,x]<-Fine_Matrix[u,x]/(Fine_Matrix[u,x]+Coarse_Matrix[u,x])
    Acoarse_Stage[u,x]<-Coarse_Matrix[u,x]/(Fine_Matrix[u,x]+Coarse_Matrix[u,x])
  }
}


Afine_Stage_Plot<-Afine_Stage
for (u in 1:nrow(Afine_Stage_Plot)){
  for (x in 1:ncol(Afine_Stage_Plot)){
    if (Afine_Stage_Plot[u,x]=="0"){
      Afine_Stage_Plot[u,x]<-NA
    }
  }
}

geoscalePlot(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Afine_Stage_Plot[1,], units=c( "Age", "Epoch"), type="b", lwd=2, data.lim=c(0, 1), label=expression(paste("A"["fine"])), erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(Afine_Stage_Plot)){
  lines(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), Afine_Stage_Plot[x,], type="b", pch=x, col=x, lwd=2)
}
legend(320, 70, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)





#Random Affinities Coarse Fine Stage Level
RNDfine_Final_Stage<-array(0, dim=c(length(Taxon_Names), ncol=length(TimeBins), 1000))
dimnames(RNDfine_Final_Stage)[[1]]<-Taxon_Names
dimnames(RNDfine_Final_Stage)[[2]]<-TimeBins

RNDcoarse_Final_Stage<-array(0, dim=c(length(Taxon_Names), ncol=length(TimeBins), 1000))
dimnames(RNDcoarse_Final_Stage)[[1]]<-Taxon_Names
dimnames(RNDcoarse_Final_Stage)[[2]]<-TimeBins

RNDfine_Final_SD_Stage<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
RNDcoarse_Final_SD_Stage<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
rownames(RNDfine_Final_SD_Stage)<-Taxon_Names
colnames(RNDcoarse_Final_SD_Stage)<-TimeBins
rownames(RNDfine_Final_SD_Stage)<-Taxon_Names
colnames(RNDfine_Final_SD_Stage)<-TimeBins

Atotfine_Stage<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
Atotcoarse_Stage<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
rownames(Atotfine_Stage)<-Taxon_Names
colnames(Atotfine_Stage)<-TimeBins
rownames(Atotcoarse_Stage)<-Taxon_Names
colnames(Atotcoarse_Stage)<-TimeBins

SRA_Fine_Stage<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
SRA_Coarse_Stage<-matrix(nrow=length(Taxon_Names), ncol=length(TimeBins), 0)
rownames(SRA_Fine_Stage)<-Taxon_Names
colnames(SRA_Fine_Stage)<-TimeBins
rownames(SRA_Coarse_Stage)<-Taxon_Names
colnames(SRA_Coarse_Stage)<-TimeBins

for (u in 1:nrow(Fine_Matrix)){
  for (x in 1:ncol(Fine_Matrix)){
    for (p in 1:1000){
      
      if (Fine_Matrix[u,x]==0 & Coarse_Matrix[u,x]==0) next
     FineCoarse<-subset(New, New[,8]==TimeBins[x])
     FineCoarse<-subset(FineCoarse, FineCoarse[,11]=="Fine_Grained"| FineCoarse[,11]=="Coarse_Grained")
     Rand_Sample_Fine<-FineCoarse[sample(nrow(FineCoarse), size=(Fine_Matrix[u,x]+Coarse_Matrix[u,x]), replace=TRUE),]
     
      Rand_Sample_Fine<-as.matrix(Rand_Sample_Fine)
      if (ncol(Rand_Sample_Fine)==1) Rand_Sample_Fine<-t(Rand_Sample_Fine)
      
      RNDfine<-matrix(0, nrow = (Fine_Matrix[u,x]+Coarse_Matrix[u,x]), ncol = ncol(New))
      RNDcoarse<-matrix(0, nrow = (Fine_Matrix[u,x]+Coarse_Matrix[u,x]), ncol = ncol(New))     
      
      for (n in 1:(Fine_Matrix[u,x]+Coarse_Matrix[u,x])){
        if (Rand_Sample_Fine[n,11]=="Fine_Grained"){
          RNDfine[n,]<-as.matrix(Rand_Sample_Fine[n,])
        }
      }
      
      for (n in 1:(Fine_Matrix[u,x]+Coarse_Matrix[u,x])){
        if (Rand_Sample_Fine[n,11]=="Coarse_Grained"){ 
          RNDcoarse[n,]<-as.matrix(Rand_Sample_Fine[n,])
        }
      }
      
      RNDfine<-RNDfine[apply(RNDfine!=0,1, any),,drop=FALSE]
      RNDcoarse<-RNDcoarse[apply(RNDcoarse!=0,1, any),,drop=FALSE]
      
      Afine_Rand<-nrow(RNDfine)/(nrow(RNDfine)+nrow(RNDcoarse))
      Acoarse_Rand<-nrow(RNDcoarse)/(nrow(RNDfine)+nrow(RNDcoarse))
      
      RNDfine_Final_Stage[u,x,p]<-Afine_Rand
      RNDcoarse_Final_Stage[u,x,p]<-Acoarse_Rand
      
    }
    RNDfine_Final_SD_Stage[u,x]<-sd(RNDfine_Final_Stage[u,x,])
    RNDcoarse_Final_SD_Stage[u,x]<-sd(RNDcoarse_Final_Stage[u,x,])
    
    Atotfine_Stage[u,x]<-mean(RNDfine_Final_Stage[u,x,])
    Atotcoarse_Stage[u,x]<-mean(RNDcoarse_Final_Stage[u,x,])
    
    SRA_Fine_Stage[u,x]<-(Afine_Stage[u,x]-Atotfine_Stage[u,x])/RNDfine_Final_SD_Stage[u,x]
    if (Tot_Stage_Fine_Coarse[u,x]<5){
      SRA_Fine_Stage[u,x]<-"-"
    }
    write.xlsx(SRA_Fine_Stage, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_SRA/SRA_Fine_Stage.xlsx")
    SRA_Coarse_Stage[u,x]<-(Acoarse_Stage[u,x]-Atotcoarse_Stage[u,x])/RNDcoarse_Final_SD_Stage[u,x]
    if (Tot_Stage_Fine_Coarse[u,x]<5){
      SRA_Coarse_Stage[u,x]<-"-"
    }
    write.xlsx(SRA_Coarse_Stage, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_SRA/SRA_Coarse_Stage.xlsx")
  }
}
geoscalePlot(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), SRA_Fine_Stage[1,], units=c("Age", "Epoch"), type="b", lwd=2, data.lim=c(-2, 2), label="Standardized Relative Affinity", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(SRA_Fine_Stage)){
  lines(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), SRA_Fine_Stage[x,], type="b", pch=x, col=x, lwd=2)
}
legend(316.75, 2, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=14)
abline(0,0, lty=2)

geoscalePlot(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), SRA_Coarse_Stage[1,], units=c("Age", "Epoch"), type="b", lwd=2, data.lim=c(-2, 2), label="Standardized Relative Affinity", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(SRA_Coarse_Stage)){
  lines(c(352.8, 338.8, 327.05, 319.2, 311.1, 305.35, 301.3), SRA_Coarse_Stage[x,], type="b", pch=x, col=x, lwd=2)
}
legend(316.75, -.5, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=14)
abline(0,0, lty=2)

# Miller and Connolly, 2001 Subperiod Level Coarse Fine
Afine_Subperiod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
Acoarse_Subperiod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
rownames(Afine_Subperiod)<-Taxon_Names
colnames(Afine_Subperiod)<-SubPeriods
rownames(Acoarse_Subperiod)<-Taxon_Names
colnames(Acoarse_Subperiod)<-SubPeriods

for (u in 1:nrow(Fine_Matrix_SubPeriod)){
  for (x in 1:ncol(Fine_Matrix_SubPeriod)){
    Afine_Subperiod[u,x]<-Fine_Matrix_SubPeriod[u,x]/(Fine_Matrix_SubPeriod[u,x]+Coarse_Matrix_SubPeriod[u,x])
    Acoarse_Subperiod[u,x]<-Coarse_Matrix_SubPeriod[u,x]/(Fine_Matrix_SubPeriod[u,x]+Coarse_Matrix_SubPeriod[u,x])
  }
}

Afine_Subperiod_Plot<-Afine_Subperiod
for (u in 1:nrow(Afine_Subperiod_Plot)){
  for (x in 1:ncol(Afine_Subperiod_Plot)){
    if (Afine_Subperiod_Plot[u,x]=="0"){
      Afine_Subperiod_Plot[u,x]<-NA
    }
  }
}

geoscalePlot(c(341.05, 311.05), Afine_Subperiod_Plot[1,], units=c( "Epoch"), type="b", lwd=2, data.lim=c(0, 1), label=expression(paste("A"["fine"])), erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(Afine_Subperiod_Plot)){
  lines(c(341.05, 311.05), Afine_Subperiod_Plot[x,], type="b", pch=x, col=x, lwd=2)
}
legend(320, 50, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)





#Random Affinities Coarse Fine Subperiod Level
RNDfine_Final_Subperiod<-array(0, dim=c(length(Taxon_Names), ncol=length(SubPeriods), 1000))
dimnames(RNDfine_Final_Subperiod)[[1]]<-Taxon_Names
dimnames(RNDfine_Final_Subperiod)[[2]]<-SubPeriods

RNDcoarse_Final_Subperiod<-array(0, dim=c(length(Taxon_Names), ncol=length(SubPeriods), 1000))
dimnames(RNDcoarse_Final_Subperiod)[[1]]<-Taxon_Names
dimnames(RNDcoarse_Final_Subperiod)[[2]]<-SubPeriods

RNDfine_Final_SD_Subperiod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
RNDcoarse_Final_SD_Subperiod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
rownames(RNDfine_Final_SD_Subperiod)<-Taxon_Names
colnames(RNDfine_Final_SD_Subperiod)<-SubPeriods
rownames(RNDcoarse_Final_SD_Subperiod)<-Taxon_Names
colnames(RNDcoarse_Final_SD_Subperiod)<-SubPeriods

Atotfine_Subperiod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
Atotcoarse_Subperiod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
rownames(Atotfine_Subperiod)<-Taxon_Names
colnames(Atotfine_Subperiod)<-SubPeriods
rownames(Atotcoarse_Subperiod)<-Taxon_Names
colnames(Atotcoarse_Subperiod)<-SubPeriods

SRA_Fine_Subperiod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
SRA_Coarse_Subperiod<-matrix(nrow=length(Taxon_Names), ncol=length(SubPeriods), 0)
rownames(SRA_Fine_Subperiod)<-Taxon_Names
colnames(SRA_Fine_Subperiod)<-SubPeriods
rownames(SRA_Coarse_Subperiod)<-Taxon_Names
colnames(SRA_Coarse_Subperiod)<-SubPeriods

for (u in 1:nrow(Fine_Matrix_SubPeriod)){
  for (x in 1:ncol(Fine_Matrix_SubPeriod)){
    for (p in 1:1000){
      
      if (Fine_Matrix_SubPeriod[u,x]==0 & Coarse_Matrix_SubPeriod[u,x]==0) next
     FineCoarse<-subset(New, New[,19]==SubPeriods[x])
     FineCoarse<-subset(FineCoarse, FineCoarse[,11]=="Fine_Grained"| FineCoarse[,11]=="Coarse_Grained")
     Rand_Sample_Fine<-FineCoarse[sample(nrow(FineCoarse), size=(Fine_Matrix_SubPeriod[u,x]+Coarse_Matrix_SubPeriod[u,x]), replace=TRUE),]
      
      Rand_Sample_Fine<-as.matrix(Rand_Sample_Fine)
      if (ncol(Rand_Sample_Fine)==1) Rand_Sample_Fine<-t(Rand_Sample_Fine)
      
      RNDfine<-matrix(0, nrow = (Fine_Matrix_SubPeriod[u,x]+Coarse_Matrix_SubPeriod[u,x]), ncol = ncol(New))
      RNDcoarse<-matrix(0, nrow = (Fine_Matrix_SubPeriod[u,x]+Coarse_Matrix_SubPeriod[u,x]), ncol = ncol(New))     
      
      for (n in 1:(Fine_Matrix_SubPeriod[u,x]+Coarse_Matrix_SubPeriod[u,x])){
        if (Rand_Sample_Fine[n,11]=="Fine_Grained"){
          RNDfine[n,]<-as.matrix(Rand_Sample_Fine[n,])
        }
      }
      
      for (n in 1:(Fine_Matrix_SubPeriod[u,x]+Coarse_Matrix_SubPeriod[u,x])){
        if (Rand_Sample_Fine[n,11]=="Coarse_Grained"){
          RNDcoarse[n,]<-as.matrix(Rand_Sample_Fine[n,])
        }
      }
      
      
      RNDfine<-RNDfine[apply(RNDfine!=0,1, any),,drop=FALSE]
      RNDcoarse<-RNDcoarse[apply(RNDcoarse!=0,1, any),,drop=FALSE]
      
      Afine_Rand<-nrow(RNDfine)/(nrow(RNDfine)+nrow(RNDcoarse))
      Acoarse_Rand<-nrow(RNDcoarse)/(nrow(RNDfine)+nrow(RNDcoarse))
      
      RNDfine_Final_Subperiod[u,x,p]<-Afine_Rand
      RNDcoarse_Final_Subperiod[u,x,p]<-Acoarse_Rand
      
      
    }
    RNDfine_Final_SD_Subperiod[u,x]<-sd(RNDfine_Final_Subperiod[u,x,])
    RNDcoarse_Final_SD_Subperiod[u,x]<-sd(RNDcoarse_Final_Subperiod[u,x,])
    
    Atotfine_Subperiod[u,x]<-mean(RNDfine_Final_Subperiod[u,x,])
    Atotcoarse_Subperiod[u,x]<-mean(RNDcoarse_Final_Subperiod[u,x,])
    
    SRA_Fine_Subperiod[u,x]<-(Afine_Subperiod[u,x]-Atotfine_Subperiod[u,x])/RNDfine_Final_SD_Subperiod[u,x]
    if (Tot_Subperiod_Fine_Coarse[u,x]<5){
      SRA_Fine_Subperiod[u,x]<-"-"
    }
    write.xlsx(SRA_Fine_Subperiod, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_SRA/SRA_Fine_Subperiod.xlsx")
    SRA_Coarse_Subperiod[u,x]<-(Acoarse_Subperiod[u,x]-Atotcoarse_Subperiod[u,x])/RNDcoarse_Final_SD_Subperiod[u,x]
    if (Tot_Subperiod_Fine_Coarse[u,x]<5){
      SRA_Coarse_Subperiod[u,x]<-"-"
    }
    write.xlsx(SRA_Coarse_Subperiod, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_SRA/SRA_Coarse_Subperiod.xlsx")
    
  }
}
geoscalePlot(c(341.05, 311.05), SRA_Fine_Subperiod[1,], units=c( "Epoch"), type="b", lwd=2, data.lim=c(-2, 2), label="Standardized Relative Affinity", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(SRA_Fine_Subperiod)){
  lines(c(341.05, 311.05), SRA_Fine_Subperiod[x,], type="b", pch=x, col=x, lwd=2)
}
legend(320, -.75, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)
abline(0,0, lty=2)

geoscalePlot(c(341.05, 311.05), SRA_Coarse_Subperiod[1,], units=c( "Epoch"), type="b", lwd=2, data.lim=c(-2, 2), label="Standardized Relative Affinity", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:nrow(SRA_Coarse_Subperiod)){
  lines(c(341.05, 311.05), SRA_Coarse_Subperiod[x,], type="b", pch=x, col=x, lwd=2)
}
legend(320, 2, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)
abline(0,0, lty=2)


#Calculates carb siliciclastic affinities across the entire Carbonferous for all clades following methods of Miller and Connolly, 2001
RNDc_Final<-matrix(0, nrow=1000, ncol=1)
RNDs_Final<-matrix(0, nrow=1000, ncol=1)

Ac_Carboniferous<-matrix(0, nrow=1, ncol=5)
As_Carboniferous<-matrix(0, nrow=1, ncol=5)
colnames(Ac_Carboniferous)<-Taxon_Names
colnames(As_Carboniferous)<-Taxon_Names

for (x in 1:length(Numbers_Carb)){
  Ac_Carboniferous[x]<-Numbers_Carb[x]/(Numbers_Carb[x]+Numbers_Silic[x])
  As_Carboniferous[x]<-Numbers_Silic[x]/(Numbers_Carb[x]+Numbers_Silic[x])
}

geoscalePlot(c(328.9), Ac_Carboniferous[1], units=c( "Period"), type="b", lwd=2, data.lim=c(0, 1), label=expression(paste("A"["c"])), erotate=0, direction="horizontal", abbrev=c("Subperiod"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:length(Ac_Carboniferous)){
  lines(c(328.9), Ac_Carboniferous[x], type="b", pch=x, col=x, lwd=2)
}
legend(320, 200, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)


RNDc_Carboniferous_Final_SD<-matrix(0, ncol=length(Numbers_Carb))
RNDs_Carboniferous_Final_SD<-matrix(0, ncol=length(Numbers_Carb))
colnames(RNDc_Carboniferous_Final_SD)<-Taxon_Names
colnames(RNDs_Carboniferous_Final_SD)<-Taxon_Names


Atotc_Carboniferous<-matrix(0, ncol=length(Numbers_Carb))
Atots_Carboniferous<-matrix(0, ncol=length(Numbers_Carb))
colnames(Atotc_Carboniferous)<-Taxon_Names
colnames(Atots_Carboniferous)<-Taxon_Names

SRA_Carb_Carboniferous<-matrix(0, ncol=length(Numbers_Carb))
SRA_Silic_Carboniferous<-matrix(0, ncol=length(Numbers_Carb))
colnames(SRA_Carb_Carboniferous)<-Taxon_Names
colnames(SRA_Silic_Carboniferous)<-Taxon_Names

#Caucluate Random Affinities at the Period Level
for (x in 1:length(Numbers_Carb)){
for (p in 1:1000){
  
  CarbSilic<-subset(New, New[,10]=="Carbonate"| New[,10]=="Siliciclastic")
  Rand_Sample_Carb<-CarbSilic[sample(nrow(CarbSilic), size=(Numbers_Carb[x]+Numbers_Silic[x]), replace=TRUE),]

RNDc<-matrix(0, nrow =(Numbers_Carb[x]+Numbers_Silic[x]), ncol = ncol(New))
RNDs<-matrix(0, nrow =(Numbers_Carb[x]+Numbers_Silic[x]), ncol = ncol(New))

  for (n in 1:(Numbers_Carb[x]+Numbers_Silic[x])){
    if (Rand_Sample_Carb[n,10]=="Carbonate"){
      RNDc[n,]<-as.matrix(Rand_Sample_Carb[n,])
    }
  }

for (n in 1:(Numbers_Carb[x]+Numbers_Silic[x])){
  if (Rand_Sample_Carb[n,10]=="Siliciclastic"){
    RNDs[n,]<-as.matrix(Rand_Sample_Carb[n,])
  }
}


RNDc<-RNDc[apply(RNDc!=0,1, any),,drop=FALSE]
RNDs<-RNDs[apply(RNDs!=0,1, any),,drop=FALSE]

Ac_Rand<-nrow(RNDc)/(nrow(RNDc)+nrow(RNDs))
As_Rand<-nrow(RNDs)/(nrow(RNDc)+nrow(RNDs))

RNDc_Final[p]<-Ac_Rand
RNDs_Final[p]<-As_Rand
}

RNDc_Carboniferous_Final_SD[x]<-sd(RNDc_Final)
RNDs_Carboniferous_Final_SD[x]<-sd(RNDs_Final)

Atotc_Carboniferous[x]<-mean(RNDc_Final)
Atots_Carboniferous[x]<-mean(RNDs_Final)

SRA_Carb_Carboniferous[x]<-(Ac_Carboniferous[x]-Atotc_Carboniferous[x])/RNDc_Carboniferous_Final_SD[x]
write.xlsx(SRA_Carb_Carboniferous, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_SRA/SRA_Carb_Carboniferous.xlsx")
SRA_Silic_Carboniferous[x]<-(As_Carboniferous[x]-Atots_Carboniferous[x])/RNDs_Carboniferous_Final_SD[x]
write.xlsx(SRA_Silic_Carboniferous, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_SRA/SRA_Silic_Carboniferous.xlsx")
}

geoscalePlot(c(328.9), SRA_Carb_Carboniferous[1], units=c( "Period"), type="b", lwd=2, data.lim=c(-3, 3), label="Standardized Relative Affinity", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:ncol(SRA_Carb_Carboniferous)){
  lines(c(328.9), SRA_Carb_Carboniferous[x], type="b", pch=x, col=x, lwd=2)
}
legend(357, -.75, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)
abline(0,0, lty=2)

geoscalePlot(c(328.9), SRA_Silic_Carboniferous[1], units=c( "Period"), type="b", lwd=2, data.lim=c(-3, 3), label="Standardized Relative Affinity", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:ncol(SRA_Silic_Carboniferous)){
  lines(c(328.9), SRA_Silic_Carboniferous[x], type="b", pch=x, col=x, lwd=2)
}
legend(357, -.75, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)
abline(0,0, lty=2)

#Calculates environmental affinities Coarse Fine across the entire Carbonferous for all clades following methods of Miller and Connolly, 2001

RNDfine_Final<-matrix(0, nrow=1000, ncol=1)
RNDcoarse_Final<-matrix(0, nrow=1000, ncol=1)

Afine_Carboniferous<-matrix(0, nrow=1, ncol=5)
Acoarse_Carboniferous<-matrix(0, nrow=1, ncol=5)
colnames(Afine_Carboniferous)<-Taxon_Names
colnames(Acoarse_Carboniferous)<-Taxon_Names

for (x in 1:length(Numbers_Carb)){
  Afine_Carboniferous[x]<-Numbers_Fine[x]/(Numbers_Fine[x]+Numbers_Coarse[x])
  Acoarse_Carboniferous[x]<-Numbers_Coarse[x]/(Numbers_Fine[x]+Numbers_Coarse[x])
}

geoscalePlot(c(328.9), Afine_Carboniferous[1], units=c( "Period"), type="b", lwd=2, data.lim=c(0, 1), label=expression(paste("A"["fine"])), erotate=0, direction="horizontal", abbrev=c("Subperiod"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:length(Afine_Carboniferous)){
  lines(c(328.9), Afine_Carboniferous[x], type="b", pch=x, col=x, lwd=2)
}
legend(320, 200, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)




RNDfine_Carboniferous_Final_SD<-matrix(0, ncol=length(Numbers_Fine))
RNDcoarse_Carboniferous_Final_SD<-matrix(0, ncol=length(Numbers_Fine))
colnames(RNDfine_Carboniferous_Final_SD)<-Taxon_Names
colnames(RNDcoarse_Carboniferous_Final_SD)<-Taxon_Names

Atotfine_Carboniferous<-matrix(0, ncol=length(Numbers_Fine))
Atotcoarse_Carboniferous<-matrix(0, ncol=length(Numbers_Fine))
colnames(Atotfine_Carboniferous)<-Taxon_Names
colnames(Atotcoarse_Carboniferous)<-Taxon_Names

SRA_Fine_Carboniferous<-matrix(0, ncol=length(Numbers_Fine))
SRA_Coarse_Carboniferous<-matrix(0, ncol=length(Numbers_Fine))
colnames(SRA_Fine_Carboniferous)<-Taxon_Names
colnames(SRA_Coarse_Carboniferous)<-Taxon_Names

#Caucluate Random Affinities For Whole Carboniferous for fine coarse

for (x in 1:length(Numbers_Fine)){
  for (p in 1:1000){
    
    FineCoarse<-subset(New, New[,11]=="Fine_Grained"| New[,11]=="Coarse_Grained")
    Rand_Sample_Fine<-FineCoarse[sample(nrow(FineCoarse), size=(Numbers_Fine[x]+Numbers_Coarse[x]), replace=TRUE),]
  
    
    RNDfine<-matrix(0, nrow =(Numbers_Fine[x]+Numbers_Coarse[x]), ncol = ncol(New))
    RNDcoarse<-matrix(0, nrow =(Numbers_Fine[x]+Numbers_Coarse[x]), ncol = ncol(New))
    
    for (n in 1:(Numbers_Fine[x]+Numbers_Coarse[x])){
      if (Rand_Sample_Fine[n,11]=="Fine_Grained"){
        RNDfine[n,]<-as.matrix(Rand_Sample_Fine[n,])
      }
    }
    
    for (n in 1:(Numbers_Fine[x]+Numbers_Coarse[x])){
      if (Rand_Sample_Fine[n,11]=="Coarse_Grained"){
        RNDcoarse[n,]<-as.matrix(Rand_Sample_Fine[n,])
      }
    }
    
    
    RNDfine<-RNDfine[apply(RNDfine!=0,1, any),,drop=FALSE]
    RNDcoarse<-RNDcoarse[apply(RNDcoarse!=0,1, any),,drop=FALSE]
    
    Afine_Rand<-nrow(RNDfine)/(nrow(RNDfine)+nrow(RNDcoarse))
    Acoarse_Rand<-nrow(RNDcoarse)/(nrow(RNDfine)+nrow(RNDcoarse))
    
    RNDfine_Final[p]<-Afine_Rand
    RNDcoarse_Final[p]<-Acoarse_Rand
  }
  
  RNDfine_Carboniferous_Final_SD[x]<-sd(RNDfine_Final)
  RNDcoarse_Carboniferous_Final_SD[x]<-sd(RNDcoarse_Final)
  
  Atotfine_Carboniferous[x]<-mean(RNDfine_Final)
  Atotcoarse_Carboniferous[x]<-mean(RNDcoarse_Final)
  
  SRA_Fine_Carboniferous[x]<-(Afine_Carboniferous[x]-Atotfine_Carboniferous[x])/RNDfine_Carboniferous_Final_SD[x]
  write.xlsx(SRA_Fine_Carboniferous, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_SRA/SRA_Fine_Carboniferous.xlsx")
  SRA_Coarse_Carboniferous[x]<-(Acoarse_Carboniferous[x]-Atotcoarse_Carboniferous[x])/RNDcoarse_Carboniferous_Final_SD[x]
  write.xlsx(SRA_Coarse_Carboniferous, "/Users/jeffreythompson/Research/Substrate_Affinity/Results_SRA/SRA_Coarse_Carboniferous.xlsx")
}

geoscalePlot(c(328.9), SRA_Fine_Carboniferous[1], units=c( "Period"), type="b", lwd=2, data.lim=c(-3, 3), label="Standardized Relative Affinity", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:ncol(SRA_Fine_Carboniferous)){
  lines(c(328.9), SRA_Fine_Carboniferous[x], type="b", pch=x, col=x, lwd=2)
}
legend(357, -.75, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)
abline(0,0, lty=2)

geoscalePlot(c(328.9), SRA_Coarse_Carboniferous[1], units=c( "Period"), type="b", lwd=2, data.lim=c(-3, 3), label="Standardized Relative Affinity", erotate=0, direction="horizontal", abbrev=c("Stage"), age.lim=c(357.1,301.1), cex.ts=.9)
for (x in 2:ncol(SRA_Coarse_Carboniferous)){
  lines(c(328.9), SRA_Coarse_Carboniferous[x], type="b", pch=x, col=x, lwd=2)
}
legend(357, -.75, Taxon_Names, cex=.6, pch=c(1:5), col=c(1:5), text.width=15)
abline(0,0, lty=2)
#End of Code