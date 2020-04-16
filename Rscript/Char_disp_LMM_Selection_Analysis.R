#######################################################
### R Code for publication in prep
### Title :  Evidence for character convergence, but not displacement in root traits

# Sara Colom and Regina Baucom #


#                   Main procedures:
#         2018 competition field experiment 
#   Linear mixed models for fixed and random effects
#           Principal Componant Analysis 
#           Calculation of phenotypic distance
#   Testing Main Prediction of Character displacement
#           Selection analysis on PC's        


#######################################################
### *Note that this code does not show all preliminary#
### analysis                                        ###
#######################################################

#######################################################
### Getting started								
###################



# Read in main libraries and save aesthetics for ggplot2

library(sjstats)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(readxl)
library(emmeans)
library(factoextra)
library(FactoMineR)
library(ggpubr)
library(kableExtra)
library(DT)
library(plotly)
library(corrplot)
library(tm)
library(tidyverse)
library(httr)
library(janitor)
library(missMDA)
library(ggpubr)
library(randomForest)
library(rdist)

#Save theme text setting
Tx<-theme(axis.text.x = element_text(face="bold",  
                                     size=8),
          axis.text.y = element_text(face="bold", 
                                     size=8))+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

Tx2<-theme(axis.text.y = element_text(face="bold", 
                                      size=20),
           axis.title.y = element_text(face="bold", 
                                       size=20))+
  theme(axis.text.x = element_text(vjust = 1, hjust=1,face="bold",angle=0,size=20),
        axis.title.x = element_text(angle=0,size=20,face="bold"),
        plot.title=element_text(size=25,hjust=0))


TxWhite<-theme(axis.text.x = element_text(face="bold",size=20,color="white"),
               axis.text.y = element_text(face="bold", 
                                          size=20,color="white"))+
  theme(axis.text.x = element_text(vjust = 1, hjust=1))

# Black themed background
BlackTheme<-theme(
  # get rid of panel grids
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Change plot and panel background
  plot.background=element_rect(fill = "black"),
  panel.background = element_rect(fill = 'black'))

FocusTraits=c("trans.SKL_NODES", "trans.DIA_STM", "trans.DIA_STM_SIMPLE", 
              "trans.AREA", "trans.AVG_DENSITY", "trans.TD_AVG", "trans.WIDTH_MAX", 
              "trans.D10", "trans.D20", "trans.D30", "trans.D40", "trans.D50", 
              "trans.D60", "trans.D70", "trans.D80", "trans.D90", "trans.RDISTR_X", 
              "trans.RDISTR_Y", "trans.SKL_DEPTH", "trans.SKL_WIDTH", "trans.RTP_COUNT", 
              "trans.STA_RANGE", "trans.STA_MAX", "trans.RTA_RANGE", "trans.RTA_MAX", 
              "trans.ADVT_COUNT", "trans.BASAL_COUNT", "trans.ADVT_ANG", "trans.BASAL_ANG", 
              "trans.HYP_DIA", "trans.TAP_DIA", "trans.MAX_DIA_90", "trans.DROP_50"
)


# Read in data

  # Fitness data
  Fit<-read.csv("~/Desktop/CharacterDisplacementRootTraits/CleanData/FitPA4.csv")
  
  # Root trait data
  tot2<-read.csv("~/Desktop/CharacterDisplacementRootTraits/CleanData/totPA4.csv")
  
  # Size
  size<-read.csv("~/Desktop/CharacterDisplacementRootTraits/CleanData/SizeData.csv")
  
  
# Sample sizes
  table(Fit$Trt,Fit$Species)
  table(tot2$Trt,tot2$Species)

  
# Create vector with list of variables to exclude from performing normalizing transformation
  EXCLUDE<-c("Position", "Species","EXP_CODE", "Image_Name","IMAGE_ID", "CIR_RATIO", "X_PIXEL", "Full.Path","Name","Type","Date","URL","Last.Updated","Folder_Name","UniqId.x","File_copy","UniqId.y",
             "Y_PIXEL", "X_SCALE", "Y_SCALE", "COMP_TIME", "Comment", "Block", 
             "Order", "ML", "GerminationDate", "Comment1", "Dead_plant", "DeathCause", 
             "RootsHarvested", "SeedsCounted", "Trt", "Comp", "Combos", "Population"
  )
  
# Load libraries used to perform transformation on root traits
  
  library(rcompanion)
  library(MASS)
  
# Transform all of the root traits of interest--i.e., variables excluded from list above  
  predictors<-tot2[-which(names(tot2)%in%EXCLUDE)]
  library(caret)
  trans = preProcess(predictors, 
                     c("BoxCox", "center", "scale")) # Box-Cox transformation, zero mean centered and 1 SD scaled
  predictorsTrans = data.frame(
    trans = predict(trans, predictors))

  
  
  ##################################################  
# Principal Componant Analysis on traits of interest
  ##################################################   
  
# Important to have non-NA values for 'alone' treatment 'Combos' variable in 
#  competition so it is not removed when 'na.omit' is applied

# Raw data
tot2$Combos<-as.character(tot2$Combos) # Make character
tot2[which(tot2$Trt=="Alone"),"Combos"]<-"none" # Replace NA with 'none' in 'Combos' column
tot2$Combos<-as.factor(tot2$Combos) # Make factor

# Transformed data
totTrans$Combos<-as.character(totTrans$Combos) # Make character
totTrans[which(totTrans$Trt=="Alone"),"Combos"]<-"none" # Replace NA with 'none' in 'Combos' column
totTrans$Combos<-as.factor(totTrans$Combos) # Make factor

# Create a data set of FULL observations only for PCA
Dirt<-(na.omit(totTrans[c('ML','Species','Trt','Block','Combos','Position',FocusTraits)])) # Omit rows with missing information

Qually<-which(names(Dirt)%in%c("Species","Trt","Block","ML","Position","Combos")) # Save qualitative variables as vector

ObsRoot.pca<-PCA(Dirt,quali.sup=c(Qually),scale.unit=F,graph=T) # Run PCA

# Incorporate PC values as traits/columns 
Dirt$PC1<-as.vector(ObsRoot.pca$ind$coord[,1])
Dirt$PC2<-as.vector(ObsRoot.pca$ind$coord[,2])
Dirt$PC3<-as.vector(ObsRoot.pca$ind$coord[,3])
Dirt$PC4<-as.vector(ObsRoot.pca$ind$coord[,4])

#  Subset for traits of interest from transformed root data 
FocusTraits<-names(totTrans[1:33])
DirtSub<-totTrans[which(names(totTrans)%in%c(FocusTraits,"Species","Trt","Block","ML","Position","Combos"))]

DirtSub$Combos<-as.character(DirtSub$Combos)
DirtSub[which(DirtSub$Trt=="Alone"),]$Combos<-"None"
DirtSub$Combos<-as.factor(DirtSub$Combos)



##################################################  
#         Linear Mixed Model ANOVAS
##################################################   

# Final model--reported in Table 2.
# Evaluating across species

# PC 1
PC1model<-lmer(PC1~Species+Block+Trt+(1|ML),Dirt)

anova(PC1model) 
ranova(PC1model)

# PC 2
PC2model<-lmer(PC2~Species+Block+Trt+(1|ML),Dirt)
anova(PC2model)
ranova(PC2model) 

#PC 3
PC3model<-lmer(PC3~Species+Block+Trt+(1|ML),Dirt)

anova(PC3model) 
ranova(PC3model)

#PC 4
PC4model<-lmer(PC4~Species+Block+Trt+(1|ML),Dirt)

anova(PC4model)
ranova(PC4model) 

# Final model--referenced in main text
# Evaluating within species

DirtP<-droplevels(Dirt%>%filter(Species=="Ip")) # Subset for I. purpurea
DirtH<-droplevels(Dirt%>%filter(Species!="Ip")) # Subset for I. purpurea

# PC 1
PC1model<-lmer(PC1~Block+Trt+(1|ML),DirtP)

anova(PC1model) # Treatment effect and NO block effect
ranova(PC1model) # No maternal line effect

# PC 2
PC2model<-lmer(PC2~Block+Trt+(1|ML),DirtP)

anova(PC2model)
ranova(PC2model)


#PC 3
PC3model<-lmer(PC3~Block+Trt+(1|ML),DirtP)
# Comment: singular fit
anova(PC3model) # Species and Block effects
ranova(PC3model)

#PC 4
PC4model<-lmer(PC4~Block+Trt+(1|ML),DirtP)

anova(PC4model)
ranova(PC4model) # Significant maternal line effect for PC 4



##################################################  
#       Remove Block Effects from PC's
##################################################   

# Perform linear regression of each trait and seed number onto block and work with residuals

Data1<-totTrans
Data1[1:33]<-NA 

for(i in 1:33) {
  Residual<-lm(totTrans[,i]~Block,totTrans,na.action="na.exclude")$residuals
  Data1[names(Residual),i]<-Residual
}

BlkRmv<-Data1



##################################################  
#      PCA of Dirt Focus Traits w/t block*
##################################################   


## Use data set w Block effects removed

BlkRmv$Combos<-as.character(BlkRmv$Combos)
BlkRmv[which(BlkRmv$Trt=="Alone"),"Combos"]<-"none"
BlkRmv$Combos<-as.factor(BlkRmv$Combos)

BlkRmvFull<-na.omit(BlkRmv) # Subset full data set for PCA

Qually<-which(names(BlkRmvFull)%in%c("Species","Population","Trt","Block","ML","Position","Combos")) # Qualitative variables

res.pca<-PCA(BlkRmvFull,quali.sup=c(Qually),scale.unit=T,graph=F)

BlkRmvFull$Species<-as.factor(BlkRmvFull$Species)

# PC 1 v PC 2 plot

p<-fviz_pca_ind(res.pca, label="none", habillage=BlkRmvFull$Species,
                addEllipses=TRUE, ellipse.level=0.95)+
  theme_bw()

pca<-p + scale_color_manual(values=c('#999999','#E69F00'))+
  theme_classic()+Tx+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  ggtitle("")+
  xlab("PCA 1 (22.5%)")+
  ylab("PCA 2 (20.0%)")+
  Tx2
pca



## PCA 3 v 2 plot

q<-fviz_pca_ind(res.pca, label="none", habillage=BlkRmvFull$Species,
                addEllipses=TRUE, ellipse.level=0.95,axes = c(2,3))+
  theme_bw()
pcaQ<-q + scale_color_manual(values=c('#999999','#E69F00'))+
  theme_classic()+Tx+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  ggtitle("")+
  xlab("PCA 3 (14.1 %)")+
  ylab("PCA 4 (10.5 %)")

pcaQ

# Scree Polot --- Supplimentary Fig. 1

eigenvalues <- res.pca$eig
barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        main = "",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eigenvalues), eigenvalues[, 2], 
      type="b", pch=19, col = "red")+
  
  fviz_screeplot(res.pca, ncp=10,barcolor="gold",barfill="gold",line="white")

## Plot the % Contribution of each trait on each PC, Fig. 3

Contrib<-data.frame(res.pca$var$contrib)
RowNames<-row.names(Contrib)
Contrib$Trait<-gsub("trans.","",RowNames)

# Use position=position_dodge()
q<-ggplot(data=Contrib, aes(reorder(Trait, -Dim.1),Dim.1)) +
  geom_bar(stat="identity", position=position_dodge(),fill="darkgrey")+
  theme_classic()+
  xlab("")+
  ylab("")+
  ggtitle("")

q<-q +Tx+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size=12),axis.text.y = element_text(size=15),plot.title = element_text(hjust = 0.5,size=20))+
  ggtitle("PC1 (Root topology)")

r<-ggplot(data=Contrib, aes(reorder(Trait, -Dim.2), y=Dim.2)) +
  geom_bar(stat="identity", position=position_dodge(),fill="darkgrey")+
  theme_classic()+
  xlab("")+
  ylab("")+
  ggtitle("")

r<-r +Tx+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size=12),axis.text.y = element_text(size=15),plot.title = element_text(hjust = 0.5,size=20))+
  ggtitle("PC2 (Root architecture)")


s<-ggplot(data=Contrib, aes(reorder(Trait, -Dim.3), y=Dim.3)) +
  geom_bar(stat="identity", position=position_dodge(),fill="darkgrey")+
  theme_classic()+
  xlab("")+
  ylab("")+
  ggtitle("")

s<-s +Tx+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size=12),axis.text.y = element_text(size=15),plot.title = element_text(hjust = 0.5,size=20))+
  ggtitle("PC3 (Root size)")

t<-ggplot(data=Contrib, aes(reorder(Trait, -Dim.4), y=Dim.4)) +
  geom_bar(stat="identity", position=position_dodge(),fill="darkgrey")+
  theme_classic()+
  xlab("")+
  ylab("")+
  ggtitle("")

t<-t +Tx+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size=12),axis.text.y = element_text(size=15),plot.title = element_text(hjust = 0.5,size=20))+
  ggtitle("PC4 (Root morphology)")

# Libraries used to make it a pannel figure
library(grid)
library(cowplot)
library(gridExtra)

# Create elements for the figure pannel.

# Common y title
y.grob <- textGrob("Loading Score", 
                   gp=gpar(col="black", fontsize=25), rot=90)
# Common x title
x.grob <- textGrob("Traits", 
                   gp=gpar(col="black", fontsize=25), rot=0)

plot<-plot_grid(q,r,s,t, align='vh', vjust=1, scale = 1,labels = "AUTO",label_size = 18,label_x=c(0.15),label_y=c(0.95))

grid.arrange(arrangeGrob(plot, left = y.grob,bottom=x.grob))



##################################################  
#      Estimate Phenotypic Distance btw PC's
##################################################   


# Calculate disimilarity matrix
library(rdist)
dim(res.pca$ind$coord)

# Separate by species
Index<-which(BlkRmvFull$Species%in%"Ihed") # Index the Ihed rows

Ihed.PCA1_2<-res.pca$ind$coord[Index,1:2]
dim(Ihed.PCA1_2)

Ip.PCA1_2<-res.pca$ind$coord[-Index,1:2]
dim(Ip.PCA1_2)


## Calculate euclidean distance for each combination pairing:
# 1) Subset between Ipurp and Ihed
# 2) Merge by Position
# 3) Calculate euclidean distance between first two dimensions, first 3 and first four

TraitsAll<-BlkRmvFull

# 1) Subset between Ipurp and Ihed

TraitsAll$PCA1<-res.pca$ind$coord[,1]
TraitsAll$PCA2<-res.pca$ind$coord[,2]
TraitsAll$PCA3<-res.pca$ind$coord[,3]
TraitsAll$PCA4<-res.pca$ind$coord[,4]

Ipurp<-droplevels(subset(TraitsAll,Species=="Ip"))
Ipurp<-droplevels(subset(Ipurp,Trt=="Inter"))[c("Species","Trt","Block","ML","Position","Combos","PCA1","PCA2","PCA3","PCA4")]
Ihed<-droplevels(subset(TraitsAll,Species!="Ip"))[c("Species","Trt","Block","ML","Position","Combos","PCA1","PCA2","PCA3","PCA4")]

# Merge by Position
colnames(Ihed)[c(1,4,7,8,9,10)]<-paste(colnames(Ihed)[c(1,4,7,8,9,10)],"_Competitor",sep="") # Modify competitor's column names
Both<-merge(Ipurp,Ihed)


# Loop to calculate euclidean distance of PCA1 only
Both$PhenDist_PCA1<-0
for(i in 1:nrow(Both)){
  Focal<- Both[i,7]
  Comp<- Both[i,13]
  Both[i,]$PhenDist_PCA1<-cdist(Focal,Comp)
}

# Loop to calculate euclidean distance of PCA2 only
Both$PhenDist_PCA2<-0
for(i in 1:nrow(Both)){
  Focal<- Both[i,8]
  Comp<- Both[i,14]
  Both[i,]$PhenDist_PCA2<-cdist(Focal,Comp)
}

# Loop to calculate euclidean distance of PCA3 only
Both$PhenDist_PCA3<-0
for(i in 1:nrow(Both)){
  Focal<- Both[i,9]
  Comp<- Both[i,15]
  Both[i,]$PhenDist_PCA3<-cdist(Focal,Comp)
}

# Loop to calculate euclidean distance of PCA4 only
Both$PhenDist_PCA4<-0
for(i in 1:nrow(Both)){
  Focal<- Both[i,10]
  Comp<- Both[i,16]
  Both[i,]$PhenDist_PCA4<-cdist(Focal,Comp)
}


#######
#     Creating componant of Figure panel--Fig. 4
######

PC2<-ggplot(TraitsAll)+
  geom_histogram(aes(PCA2,fill=Species),alpha=.5)+
  #geom_density(aes(PCA2,fill=Species),alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")+
  theme_classic()+
  ylab("")+
  xlab("PCA 2")+
  Tx2

PC1<-ggplot(TraitsAll)+
  geom_histogram(aes(PCA1,fill=Species))+
  #geom_density(aes(PCA2,fill=Species),alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")+
  theme_classic()+
  Tx2+
  xlab("PCA 1")+
  ylab("Count")


PC3<-ggplot(TraitsAll)+
  geom_histogram(aes(PCA3,fill=Species),alpha=.5)+
  #geom_density(aes(PCA2,fill=Species),alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")+
  theme_classic()+
  ylab("")+
  xlab("PCA 3")+
  Tx2

PC4<-ggplot(TraitsAll)+
  geom_histogram(aes(PCA4,fill=Species),alpha=.5)+
  #geom_density(aes(PCA2,fill=Species),alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")+
  theme_classic()+
  ylab("")+
  xlab("PCA 4")+
  Tx2


blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )


ggarrange(PC2, blankPlot, pca, PC1, 
          ncol=2, nrow=2, widths=c(2, 2), heights=c(2, 3),legend="none")

ggarrange(PC4, blankPlot, pcaQ, PC3, 
          ncol=2, nrow=2, widths=c(2, 2), heights=c(2, 3),legend="none")

ggarrange(pca,pcaQ, pcal,PC1,PC2,PC4, 
          ncol=3, nrow=2,legend="none",labels = "AUTO",label.x = 0.3,label.y = 1)



##################################################  
#      Treatment effect on Fitness 
##################################################   

# (1) Open leaf number (plant size data) combined w Fitness data
# (2) Perform linear mixed model ANOVA reported in main text

Size=read.csv("~/Desktop/CharacterDisplacementRootTraits/CleanData/SizeData.csv")
Fit2=merge(Size,Fit)
Fit2$Block=as.factor(Fit2$Block) # fix str

# Full model (Reported)
SN1<-lmer(SeedNumber~Trt+Block+Leaf.Number+Block:Trt+(1|ML),Fit2)
anova(SN1)

emmeans(SN1, ~Trt) # Averaged over treatment 

##################################################  
#      Estimate fitness, rel. fitness
##################################################   


SdMnSpeciesTrt<- aggregate(SeedNumber ~ Species+Trt, Fit2, mean) # Mean residual. NOTE I am using the version with block effects removed, only
names(SdMnSpeciesTrt)[3]<- "MeanSdNm" #Rename coloumn for mean seed number
# Merge average fitness of species and treatment 
Fitmean<-merge(Fit2,SdMnSpeciesTrt,by=c("Species","Trt"))

# Calculate relative fitness as the observed seed number by the total mean seed number of that species and treatment.
Fitmean$Rel_Fit<-Fitmean$SeedNumber/Fitmean$MeanSdNm

Fitmean$Combos<-as.character(Fitmean$Combos)
Fitmean[which(Fitmean$Trt=="Alone"),"Combos"]<-"none"
Fitmean$Combos<-as.factor(Fitmean$Combos)

# Extract residuals of Block and size
Fitmean$SeedNumberResid<-NA
SeedResiduals<-(lm(Rel_Fit~Block+Leaf.Number,Fitmean))$residuals
Fitmean[names(SeedResiduals),"SeedNumberResid"]<-SeedResiduals

# Selection on phenotypic distance

# Average phenotypic distance and relative fitness by treatment, maternal line and combination type

PhenoDistPCA1.Average<-aggregate(PhenDist_PCA1~ML+Combos,Both,FUN=mean)
PhenoDistPCA2.Average<-aggregate(PhenDist_PCA2~ML+Combos,Both,FUN=mean)
PhenoDistPCA3.Average<-aggregate(PhenDist_PCA3~ML+Combos,Both,FUN=mean)
PhenoDistPCA4.Average<-aggregate(PhenDist_PCA4~ML+Combos,Both,FUN=mean)

# Combine and save the phenotypic distances
AveragePhenoDist<-(unique(cbind(PhenoDistPCA1.Average,PhenoDistPCA2.Average,PhenoDistPCA4.Average)))

AveragePhenoDist<-AveragePhenoDist[, !duplicated(colnames(AveragePhenoDist))] 

# Remove NA's in the combos for alone
Fitmean$Combos<-as.character(Fitmean$Combos)
Fitmean[which(Fitmean$Trt=="Alone"),]$Combos<-"none"
Fitmean$Combos<-as.factor(Fitmean$Combos)

RelFitMean<-aggregate(SeedResiduals~ML+Combos+Trt,Fitmean,mean)

AveragePhenoDist<-merge(RelFitMean,AveragePhenoDist)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#               Examine evidence for character displacement 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# PCA 1
summary(lm(SeedResiduals~PhenDist_PCA1,AveragePhenoDist))

# PCA 2
summary(lm(SeedResiduals~PhenDist_PCA2,AveragePhenoDist)) # Evidence of a *Negative* linear effect for PC2 (Root architecture)
summary(lm(SeedResiduals~log(PhenDist_PCA2),AveragePhenoDist)) # Evidence of a *Negative* linear effect for PC2 (Root architecture)

# PCA 4
summary(lm(SeedResiduals~PhenDist_PCA4,AveragePhenoDist))

# Plot phenotypic distance regression


# Create a generic plotting function
PlotRegressionPD<-function(Data,mapping,Method=lm,Se=F,xTitle,colorpoint="red",colorline="gold",...){
  ggplot()+
    geom_point(data=Data,mapping,color=colorpoint,...)+
    theme_classic()+
    geom_smooth(data=Data,mapping,method=Method,se=Se,color= colorline,linetype="dashed")+
    xlab(xTitle)+
    ylab("Relative fitness")
}



PlotRegressionPD(Data=AveragePhenoDist,alpha=0.5,mapping=aes(PhenDist_PCA2,SeedResiduals),xTitle="",colorpoint = "black",colorline="black",size=5)+
  ylab("Relative fitness")+
  theme(axis.text.x = element_text(face="bold",  
                                   size=25),
        axis.text.y = element_text(face="bold", 
                                   size=25))+
  theme(axis.text.x = element_text( vjust = 1, hjust=1))+
  theme(axis.title = element_text( size=30))+
  xlab( expression(paste("Root architecture PC2", " (",sqrt(PC[paste("I.purpurea")]^{2} - PC[paste("I.hederacea")]^{2}),")")))



########################################################
# Perform linear regression of PCs onto relative fitness
########################################################

# Subset for I purpurea species
IpPA4<-droplevels(TraitsAll%>%filter(Species=="Ip"))

# Calculate mean of traits by treatment and maternal line

pcFamilyMeans<-aggregate(list(IpPA4[c("PCA1","PCA2","PCA3","PCA4")]),by=list(IpPA4$Trt,IpPA4$ML),FUN=mean) #

colnames(pcFamilyMeans)[1:2]<-c("Trt","ML")

for(i in c("PCA1","PCA2","PCA3","PCA4")){
  Quadratic=pcFamilyMeans[i]*pcFamilyMeans[i]
  pcFamilyMeans=cbind(pcFamilyMeans,Quadratic)
}

dim(pcFamilyMeans)

colnames(pcFamilyMeans)[7:10]<-paste(c("PCA1","PCA2","PCA3","PCA4"),"2",sep="_")

RelFitMean<-aggregate(SeedNumberResid~ML+Trt,Fitmean,mean)

pcFamilyMeans<-merge(pcFamilyMeans,RelFitMean)
PCAall=pcFamilyMeans

# Write combined data to clean data folder
# write.csv(pcFamilyMeans,"pcFamilyMeans.csv",row.names = F)

# Subset by treatment

PCAalone<-droplevels(pcFamilyMeans%>%filter(Trt=="Alone"))
PCAcomp<-droplevels(pcFamilyMeans%>%filter(Trt!="Alone"))


PC1.res.alone<-summary(lm(SeedNumberResid~PCA1,PCAalone))
PC2.res.alone<-summary(lm(SeedNumberResid~PCA2,PCAalone))
PC3.res.alone<-summary(lm(SeedNumberResid~PCA3,PCAalone))
PC4.res.alone<-summary(lm(SeedNumberResid~PCA4,PCAalone)) # Evidence for positive selection on PC4 (Marginally Significant)


PC1.res.comp<-summary(lm(SeedNumberResid~PCA1,PCAcomp)) 
PC2.res.comp<-summary(lm(SeedNumberResid~PCA2,PCAcomp))
PC3.res.comp<-summary(lm(SeedNumberResid~PCA3,PCAcomp))
PC4.res.comp<-summary(lm(SeedNumberResid~PCA4,PCAcomp)) # Evidence for negative selection on PC4


## ANCOVA for PCA4 (Root morphology)

anova(lm(SeedNumberResid~PCA4,pcFamilyMeans))

# Output

#          Df Sum Sq  Mean Sq F value Pr(>F)
#PCA4       1 0.0069 0.006904  0.0699 0.7923

## Quadtratic seleciton on PCs between combinations

PCAcomp$PCA3_2<-PCAcomp$PCA3*PCAcomp$PCA3
PCAcomp$PCA1_2<-PCAcomp$PCA1*PCAcomp$PCA1
PCAcomp$PCA4_2<-PCAcomp$PCA4*PCAcomp$PCA4

PCAcomp$PCA2_2<-PCAcomp$PCA2*PCAcomp$PCA2


summary(lm(SeedNumberResid~PCA1_2+PCA1,PCAcomp)) 
summary(lm(SeedNumberResid~PCA2_2+PCA2,PCAcomp))
summary(lm(SeedNumberResid~PCA3_2+PCA3,PCAcomp))
summary(lm(SeedNumberResid~PCA4_2+PCA4,PCAcomp))


summary(lm(SeedNumberResid~PCA1_2+PCA1,PCAalone))
summary(lm(SeedNumberResid~PCA2_2+PCA2,PCAalone))
summary(lm(SeedNumberResid~PCA3_2+PCA3,PCAalone))
summary(lm(SeedNumberResid~PCA4_2+PCA4,PCAalone))

 # No evidence of quadratic selection

# Plot selection on PC4 (root morphology)

TraitsAllFit<-merge(TraitsAll,RelFitMean)

TraitsAllFit<-merge(TraitsAll,RelFitMean2)

Regress<-ggplot(pcFamilyMeans)+
  geom_point(aes(PCA4, SeedNumberResid,color=Trt),size=5,alpha=0.5)+
  scale_color_manual(values=c("red","black"))+
  stat_smooth(aes(PCA4, SeedNumberResid,color=Trt),alpha=0.5,method="lm", formula=y~x,se=F,fullrange = T)+
  #scale_linetype_manual(values=c("twodash", "solid"))+
  theme_classic()+
  ylab("Relative fitness")+
  xlab("Root morphology (PC4)")+
  theme(axis.text.x = element_text(  
    size=20),
    axis.text.y = element_text( 
      size=20),axis.title.y = element_text(size=20),
    axis.title.x = element_text( 
      size=20))+
  theme(axis.text.x = element_text(vjust = 1, hjust=1))+
  theme(axis.text.x = element_text(vjust = 1, hjust=1))+
  labs(color = "Treatment")+
  theme(legend.text=element_text(size=12),legend.title=element_text(size=12),legend.position = 'top',
        legend.direction = "horizontal")

# Extra


## Selection analysis on a few main traits

#(1) Subset for traits of interest
#(2) Obtain family mean of traits (average ea trait by ml, treatment and species)
#(3) Average trait by class (e.g. average D% traits, architecture traits etc)

#### Average root traits by maternal line, treatment and combination ####

# The following for loop calculates the means of root traits by ML x Trt x Species for subset of traits

D=grep("trans.D[0-9]{2}",names(BlkRmvFull))
architecture=grep("STA_MAX|RDISTR_X",names(BlkRmvFull))
morphology=grep("TD_AVG|HYP_DIA|TAP_DIA|MAX_DIA",names(BlkRmvFull))


trait<-c(D,architecture,morphology)
TraitNames=names(BlkRmvFull[trait])

# Aggregate 
for(i in trait) {
  if(i==trait[1]){
    df5<-aggregate(BlkRmv[,i] ~ML+Trt, BlkRmv, mean)
    colnames(df5)[3]<-paste(colnames(BlkRmv)[i])
  }
  else{
    df10<-aggregate(BlkRmv[,i] ~ML+Trt, BlkRmv, mean) 
    colnames(df10)[3]<-paste(colnames(BlkRmv)[i])
    df5<-merge(df10,df5,by=c("ML","Trt"))
  }
}


head(df5)



#  Combine maternal line means and fitness means
RelFitMean2<-aggregate(SeedNumberResid~ML+Trt,RelFitMean,mean)
df5Fitness<-merge(df5,RelFitMean2)


Index<-which(names(df5Fitness)%in% c('ML','Trt','Combos'))
Alone<-df5Fitness%>%filter(Trt=="Alone")
Comp<-df5Fitness%>%filter(Trt!="Alone")

# Selection gradient
summary(lm(SeedNumberResid ~trans.STA_MAX+trans.RDISTR_X+trans.D90+trans.TD_AVG+trans.MAX_DIA_90, Alone))
summary(lm(SeedNumberResid ~trans.STA_MAX+trans.RDISTR_X+trans.D90+trans.TD_AVG+trans.MAX_DIA_90, Comp))

anova(lm(SeedNumberResid ~trans.STA_MAX+trans.RDISTR_X+trans.D90+trans.TD_AVG,df5Fitness))


#####################################
# Post-Hoc Linear mixed model on single traits
#####################################

# Soil root tissue angle
STAmodel=lmer(trans.STA_MAX~Block+Trt+(1|ML),DirtP)

anova(STAmodel)
ranova(STAmodel)

# Horizontal root distribution
Xmodel=lmer(trans.RDISTR_X~Block+Trt+(1|ML),DirtP)

anova(Xmodel)
ranova(Xmodel)

# Tip diameter
TDmodel=lmer(trans.TD_AVG~Block+Trt+(1|ML),DirtP)
anova(TDmodel)
ranova(TDmodel)

# Hyp diameter
HYPmodel=lmer(trans.HYP_DIA~Block+Trt+(1|ML),DirtP)
anova(HYPmodel)
ranova(HYPmodel)

# Maximum diameter
Maxmodel=lmer(trans.MAX_DIA_90~Block+Trt+(1|ML),DirtP)
anova(Maxmodel)
ranova(Maxmodel)

# Maximum diameter
Maxmodel=lmer(trans.MAX_DIA_90~Block+Trt+(1|ML),DirtP)
anova(Maxmodel)
ranova(Maxmodel)

# Run for loop to save results of LMM on the traits that were significant from PC back regression
Vars<-c("trans.RDISTR_X","trans.MAX_DIA_90","trans.STA_MAX","trans.HYP_DIA","trans.TD_AVG","Trt","Block","ML")
Else=which(names(DirtP)%in%Vars)
D=grep("trans.D[0-9]{2}",names(DirtP))

SubsetDirt=DirtP[c(Else,D)]
dim(SubsetDirt)

Response=as.list(names(SubsetDirt[4:17]))

for(i in 1:14){
  df=SubsetDirt
  
  if(i==1){
    Trait=paste(Response[i])
    ModelFixed=data.frame(anova(lmer(paste(Response[i],"~ Block+Trt+(1|ML)"), data=df)))
    ModelRand=data.frame(ranova(lmer(paste(Response[i],"~ Block+Trt+(1|ML)"), data=df)))
    
    
    ModelRand$Effect=row.names(ModelRand)
    ModelFixed$Effect=row.names(ModelFixed)
    ModelRand$Trait=paste(Trait)
    ModelFixed$Trait=paste(Trait)
    print("Fixed effects estimated")
  }
  else{
    Trait=paste(Response[i])
    ModelFixed1=data.frame(anova(lmer(paste(Response[i],"~ Block+Trt+(1|ML)"), data=df)))
    ModelRand1=data.frame(ranova(lmer(paste(Response[i],"~ Block+Trt+(1|ML)"), data=df)))
    
    ModelRand1$Effect=row.names(ModelRand1)
    ModelFixed1$Effect=row.names(ModelFixed1)
    ModelRand1$Trait=paste(Trait)
    ModelFixed1$Trait=paste(Trait)
    ModelRand=rbind(ModelRand,ModelRand1)
    ModelFixed=rbind(ModelFixed,ModelFixed1)
    print("Random effects estimated")
  }
}

# Remove the estimates for intercept in the random effects results
ModelRand<-ModelRand[grep("ML",ModelRand$Effect),]
ModelRand$Effect<-"Mat. Line"

# Examine results
DT::datatable(ModelFixed)
DT::datatable(ModelRand)

# Round to the first three digits
#ModelFixed[1:6]=round(ModelFixed[1:6],3)
#ModelRand[1:6]=round(ModelRand[1:6],3)
