
#Install packages "dplyr", "geoscale", "xlsx", and "vegan" necessary for analysis
#also need to set working directory
setwd("/Users/jeffreythompson/Documents/Research/Palaeozoic_Echinoid_Diversity")
install.packages("dplyr")
library(dplyr)
install.packages("ggplot2")
library(ggplot2)
install.packages("geoscale")
library(geoscale)
install.packages("rJava")
library(rJava)
install.packages("xlsx")
library(xlsx)
install.packages("vegan")
library(vegan)
library(readxl)
install.packages("ggplot2")
library(ggplot2)
install.packages("RColorBrewer")
library(RColorBrewer)



Dataset<- read.xlsx("Palaeozoic_Occurrences_JRT_Updating_New.xlsx", 1)




#Want to create a matrix for preservation score alongside genus, family, and TimeBins
Genus_Names<-unique(Dataset$Genus)
Time_Bins<-unique(Dataset$Age)
Preservation_Scores<-unique(Dataset$Preservation.Score)
Family_Names<-unique(Dataset$Family)


#Create matrices for abundances of each preservational type per genus
genus_pres_score_matrix<-matrix(nrow = length(Genus_Names), ncol = length(preservation_scores),0)
rownames(genus_pres_score_matrix)<-Genus_Names
colnames(genus_pres_score_matrix)<-Preservation_Scores

genus_pres_score_matrix_PA<-matrix(nrow = length(Genus_Names), ncol = length(preservation_scores),0)
rownames(genus_pres_score_matrix_PA)<-Genus_Names
colnames(genus_pres_score_matrix_PA)<-Preservation_Scores


for (b in 1:length(Genus_Names)){
  for (c in 1:length(Preservation_Scores)){
    Total_Taxon_Stage<-subset(Dataset, Dataset$Genus==Genus_Names[b] & Dataset$Preservation.Score==Preservation_Scores[c] )
    genus_pres_score_matrix[b,c]<-nrow(Total_Taxon_Stage)
    if (genus_pres_score_matrix[b,c]>=1){
      genus_pres_score_matrix_PA[b,c]<-1
    }
  }
}


#Plot showing abundance of preservation score in each genus (preservation score on y & abundance of score on x, with different genus stacked)

yep_G_A<-colorRampPalette(brewer.pal(12, "Paired")) (16)
sum(genus_pres_score_matrix[,5])

Table<-table(Dataset$Genus, Dataset$Preservation.Score)

barplot_abundance_pres_score_genus_A<-barplot(Table, col= yep_G_A, 
                                              las=2, cex.axis= 0.75, cex.names= 0.5, horiz=T, xlim = c(0,950),
                                              main = "Abundance of preservation scores per genus", 
                                              xlab = "Total frequency of preservation score", ylab = "Preservation Score")
#text(y= barplot_abundance_pres_score_genus_A, x=sums_G_A, labels = sums_G_A, pos=4, cex=0.7, col = "red")
legend(x=800, y= 9.8, cex = 0.2, legend = rownames(Table), pch = 16, col = yep_G_A)


#Plot showing abundance of preservation score in each genus (genus on x & abundance of preservation score on y, with the different preservation scores stacked)
yep_G_B<-colorRampPalette(brewer.pal(12, "Paired")) (8)
transformed_g_ps_matrix<-t(genus_pres_score_matrix)

Table<-table(Dataset$Preservation.Score, Dataset$Genus)

barplot_abundance_pres_score_genus_B<-barplot(as.matrix(transformed_g_ps_matrix[,1:ncol(transformed_g_ps_matrix)],), col= yep_G_B, las=2, cex.axis= 0.75, cex.names= 0.52, horiz=T, xlim = c(0,800),
                                              main = "Abundance of preservation scores per genus", xlab = "Abundance of preservation Score")
legend(x=200, y= 45,cex = 0.5, legend = rownames(Table), pch = 16, col = yep_G_B)
title(ylab = "Genus", line = 0, cex.lab=1)





#Create matrices for abundances of each preservational type per family
family_pres_score_matrix<-matrix(nrow = length(Family_Names), ncol = length(Preservation_Scores),0)
rownames(family_pres_score_matrix)<-Family_Names
colnames(family_pres_score_matrix)<-Preservation_Scores

family_pres_score_matrix_PA<-matrix(nrow = length(Family_Names), ncol = length(Preservation_Scores),0)
rownames(family_pres_score_matrix_PA)<-Family_Names
colnames(family_pres_score_matrix_PA)<-Preservation_Scores
family_pres_score_matrix_PA

for (e in 1:length(Family_Names)){
  for (c in 1:length(Preservation_Scores)){
    Total_Taxon_Stage<-subset(Dataset[], Dataset$Family==Family_Names[e] & Dataset$Preservation.Score==Preservation_Scores[c] )
    family_pres_score_matrix[e,c]<-nrow(Total_Taxon_Stage)
    if (family_pres_score_matrix[e,c]>=1){
      family_pres_score_matrix_PA[e,c]<-1
    }
  }
}



#Want to plot the family_pres_score_matrix relationship
#Plot showing abundance of preservation score in each family (preservation score on y & abundance of score on x, with different family stacked)
yep_F<-colorRampPalette(brewer.pal(12, "Paired")) (15)
Table<-table(Dataset$Family, Dataset$Preservation.Score)

barplot_abundance_pres_score_family<-barplot(Table, col= yep_F, 
                                             las=2, cex.axis= 0.75, cex.names= 0.5, horiz=T, xlim = c(0,950),
                                             main = "Abundance of preservation scores per family", 
                                             xlab = "Total frequency of preservational score", 
                                             ylab = "Preservational State")
text(y= barplot_abundance_pres_score_family, x=sums_F, labels = sums_F, pos=4, cex=0.8, col = "red")
legend(x=780, y= 6, cex = 0.5, legend = rownames(Table), pch = 16, col = yep_F)

#Plot showing abundance of preservation score in each family (family on x & abundance of preservation score on y, with the different preservation scores stacked)
yep_F_B<-colorRampPalette(brewer.pal(12, "Paired")) (7)
Table<-table(Dataset$Preservation.Score, Dataset$Family)
barplot_abundance_pres_score_family_B<-barplot(Table, col= yep_F_B, 
                                               las=2, cex.axis= 0.75, cex.names= 0.46, horiz=T, xlim = c(0,1200),
                                               main = "Abundance of preservational types per family", 
                                               xlab = "Number of specimens per preservational type")
legend(x=800, y= 10,cex = 0.5, legend = rownames(Table), pch = 16, col = yep_F_B)
title(ylab = "Family", line = 0, cex.lab=1)






