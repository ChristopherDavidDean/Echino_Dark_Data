
#The following three lines install packages necessary for data manipulation and exporting data
install.packages("dplyr")
install.packages("geoscale")
install.packages("xlsx")
install.packages("vegan")
Collections<-read.csv("/Users/jeffreythompson/Research/Substrate_Affinity/Submission_1/Supplemental_Table_2.csv") #Imports Carboniferous PaleobioDB collections

#This imports the raw database of echinoid occurrences from the Supplemental_Table_1 file
setwd("/Users/jeffreythompson/Documents/Research/Palaeozoic_Echinoid_Diversity")

Total<-read.delim("Carboniferous_Occurrences.txt")


# New is the primary working database for all analyses performed in the manuscript.
# The next few lines of code remove duplicate occurrences such that only unique occurrences remain for analyases
New<-distinct(Total, Species, Genus, Family, Locality, Formation, Lithology, Grain.Size, .keep_all=TRUE)
Try<-New[!duplicated(New[c(3, 6, 7, 4, 11, 10)]),]
Replicates_1<-Try[Try[,3]=="?" & duplicated(Try[,c(6, 7, 4, 10)]),]
Try<-anti_join(Try, Replicates_1, by=c("Genus", "Locality", "Formation", "Family", "Lithology", "Grain.Size"))
Replicates_2<-Try[Try[,11]=="?" & duplicated(Try[,c(3, 6, 7, 4, 10)]),] 
New<-anti_join(Try, Replicates_2, by=c("Genus", "Locality", "Formation", "Family", "Lithology", "Grain.Size")) 
New<-as.matrix(New)


#Creates vectors for names of stages and Subperiods
TimeBins<-c("Tournaisian", "Visean", "Serpukhovian", "Bashkirian", "Moscovian", "Kasimovian", "Gzhelian")
SubPeriods<-c("Mississippian", "Pennsylvanian")
Taxon_Names<-c("Palaechinidae", "Archaeocidaridae", "Proterocidaridae", "Lepidesthidae", "Lepidocentridae")






#This imports the raw database of echinoid occurrences from the Supplemental_Table_1 file
setwd("/Users/jeffreythompson/Documents/Research/Palaeozoic_Echinoid_Diversity")

Total<-read.delim("Carboniferous_Occurrences.txt")




New<-data.frame(New)




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




