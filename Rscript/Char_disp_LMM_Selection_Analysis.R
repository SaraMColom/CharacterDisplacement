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
library(ggpubr)
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
setwd("C:/Users/Sara Colom/Desktop/character-displacement/CharacterDisplacement")

  # Fitness data
  Fit<-read.csv("CleanData/FitPA4.csv")
  
  # Root trait data
  tot2<-read.csv("CleanData/totPA4.csv")
  
  # Size
  size<-read.csv("CleanData/SizeData.csv")
  
  
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
totTrans<-predictorsTrans
totTrans$Species<-as.factor(tot2$Species)
totTrans$Population<-as.factor(tot2$Population)
totTrans$Block<-as.factor(tot2$Block) 
totTrans$ML<-as.factor(tot2$ML) 
totTrans$Trt<-as.factor(tot2$Trt)
totTrans$Position<-as.factor(tot2$Position)
totTrans$Combos<-as.factor(tot2$Combos)
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

str(Dirt)
Dirt$Block<-as.factor(Dirt$Block)

##################################################  
#         Linear Mixed Model ANOVAS
##################################################   



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

anova(PC2model) # Marginal treatment effect
ranova(PC2model)

#PC 3
PC3model<-lmer(PC3~Block+Trt+(1|ML),DirtP)
# Comment: singular fit
anova(PC3model) 
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

########################################################
#### PCA ####
########################################################



# FOCAL TRAITS
FocusTraits=c("trans.SKL_NODES", "trans.DIA_STM", "trans.DIA_STM_SIMPLE",
              "trans.AREA", "trans.AVG_DENSITY", "trans.TD_AVG", "trans.WIDTH_MAX",
              "trans.D10", "trans.D20", "trans.D30", "trans.D40", "trans.D50",
              "trans.D60", "trans.D70", "trans.D80", "trans.D90", "trans.RDISTR_X",
              "trans.RDISTR_Y", "trans.SKL_DEPTH", "trans.SKL_WIDTH", "trans.RTP_COUNT",
              "trans.STA_RANGE", "trans.STA_MAX", "trans.RTA_RANGE", "trans.RTA_MAX",
              "trans.ADVT_COUNT", "trans.BASAL_COUNT", "trans.ADVT_ANG", "trans.BASAL_ANG",
              "trans.HYP_DIA", "trans.TAP_DIA", "trans.MAX_DIA_90", "trans.DROP_50"
)


# Make sure that 'Alone' treatment is not excluded from PCA analysis
BlkRmv$Combos<-as.character(BlkRmv$Combos)
BlkRmv[which(BlkRmv$Trt=="Alone"),"Combos"]<-"none"
BlkRmv$Combos<-as.factor(BlkRmv$Combos)

BlkRmvFull<-na.omit(BlkRmv)
#BlkRmvFull <-BlkRmvFull[c(FocusTraits[-which(FocusTraits%in%c("trans.D10","trans.D20","trans.DIA_STM_SIMPLE"))],"Species","Population","Trt","Block","ML","Position","Combos")]
Qually<-which(names(BlkRmvFull)%in%c("Species","Population","Trt","Block","ML","Position","Combos"))
res.pca<-PCA(BlkRmvFull,quali.sup=c(Qually),scale.unit=T,graph=F)
BlkRmvFull$Species<-as.factor(BlkRmvFull$Species)

# PC 3 v PC 2 plot

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




## PCA 3 v 4 plot

q<-fviz_pca_ind(res.pca, label="none", habillage=BlkRmvFull$Species,
                addEllipses=TRUE, ellipse.level=0.95,axes = c(3,4))+
  theme_bw()
pcaQ<-q + scale_color_manual(values=c('#999999','#E69F00'))+
  theme_classic()+Tx+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  ggtitle("")+
  xlab("PCA 3 (14.1 %)")+
  ylab("PCA 4 (10.5 %)")

pcaQ



# Scree Polot --- Supplimentary Fig. 2

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

## Plot the % Contribution of each trait on each PC (Extra)

RootCodes <- read.csv("C:/Users/Sara Colom/Desktop/RootCodes.csv", encoding="UTF-8")
Contrib<-data.frame(res.pca$var$contrib)
RowNames<-row.names(Contrib)
Contrib$Trait<-gsub("trans.","",RowNames)

# Order both by trait names so they are correctly replaced
Contrib<-Contrib[order(Contrib$Trait),]
RootCodes<-RootCodes[order(RootCodes$TraitCode.),]

Contrib[which(Contrib$Trait%in%RootCodes$TraitCode.),"Trait"] <- paste(RootCodes[which(RootCodes$TraitCode.%in%Contrib$Trait),1])


#######     Fig. 3 plotting % contribution per dimension

# Use position=position_dodge()
q<-ggplot(data=Contrib, aes(reorder(Trait, -Dim.1),Dim.1)) +
  geom_bar(stat="identity", position=position_dodge(),fill="darkgrey")+
  theme_classic()+
  xlab("")+
  ylab("")+
  ggtitle("")

q<-q +Tx+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=11,color="black"),axis.text.y = element_text(size=15,color="black"),plot.title = element_text(hjust = 0.5,size=20))+
  ggtitle("PC1 (Root topology)")

r<-ggplot(data=Contrib, aes(reorder(Trait, -Dim.2), y=Dim.2)) +
  geom_bar(stat="identity", position=position_dodge(),fill="darkgrey")+
  theme_classic()+
  xlab("")+
  ylab("")+
  ggtitle("")

r<-r +Tx+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=11,color="black"),axis.text.y = element_text(size=15,color="black"),plot.title = element_text(hjust = 0.5,size=20))+
  ggtitle("PC2 (Root architecture)")


s<-ggplot(data=Contrib, aes(reorder(Trait, -Dim.3), y=Dim.3)) +
  geom_bar(stat="identity", position=position_dodge(),fill="darkgrey")+
  theme_classic()+
  xlab("")+
  ylab("")+
  ggtitle("")

s<-s +Tx+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=11,color="black"),axis.text.y = element_text(size=15,color="black"),plot.title = element_text(hjust = 0.5,size=20))+
  ggtitle("PC3 (Root size)")

t<-ggplot(data=Contrib, aes(reorder(Trait, -Dim.4), y=Dim.4)) +
  geom_bar(stat="identity", position=position_dodge(),fill="darkgrey")+
  theme_classic()+
  xlab("")+
  ylab("")+
  ggtitle("")

t<-t +Tx+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=11,color="black"),axis.text.y = element_text(size=15,color="black"),plot.title = element_text(hjust = 0.5,size=20))+
  ggtitle("PC4 (Root morphology)")

# Explorepercent contribution of ea trait to ea dimension --> report the % for the axis it contributes the most
DT::datatable(Contrib)

# Plot the correlation between traits and PC
corrplot(res.pca$var$cos2[,1:4], is.corr=FALSE,tl.col="black")

# Libraries used to make it a pannel figure
library(grid)
library(cowplot)
library(gridExtra)

# Create elements for the figure pannel.

# Common y title
y.grob <- textGrob("Contribution %", 
                   gp=gpar(col="black", fontsize=25), rot=90)
# Common x title
x.grob <- textGrob("Traits", 
                   gp=gpar(col="black", fontsize=25), rot=0)

plot<-plot_grid(q,r,s,t, align='vh', vjust=1, scale = 1,labels = "AUTO",label_size = 18,label_x=c(0.09),label_y=c(0.95))

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



##################################################  
#      Treatment effect on Fitness 
##################################################   

# (1) Open leaf number (plant size data) combined w Fitness data
# (2) Perform linear mixed model ANOVA reported in main text

Fit2=merge(size,Fit)
Fit2$Block=as.factor(Fit2$Block) # fix str

# Full model (Reported)
SN1<-lmer(SeedNumber~Trt+Block+Leaf.Number+Block:Trt+(1|ML),Fit2)
anova(SN1)

summary(SN1)

emmeans(SN1, ~Trt) # Averaged over treatment 

FitMeansTrt<-data.frame(emmeans(SN1, ~Trt))


##################################################  
#      Estimate fitness, rel. fitness
##################################################   


SdMnSpeciesTrt<- aggregate(SeedNumber ~ Species + Trt, Fit2, mean) # Mean residual. NOTE I am using the version with block effects removed, only
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
SeedResiduals<-(lm(Rel_Fit~Block+Leaf.Number,Fitmean))$residuals # Removed Size and Blaock here
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

# Add standard error of each phenotypic distance measured
library(plotrix)

PhenoDistStandErr<-aggregate(list(Both[grep("PhenDist",names(Both))]),by=list(Both$Combos,Both$ML),FUN=std.error) #

colnames(PhenoDistStandErr)[1:2]=c("Combos","ML")
colnames(PhenoDistStandErr)[grep("Phen",names(PhenoDistStandErr))]=paste(colnames(PhenoDistStandErr)[grep("Phen",names(PhenoDistStandErr))],"SE",sep="_")

# Combine phenotypic distance standard error and the mean values of phenotypic distances
CombPhenStdErrMean<-merge(PhenoDistStandErr,AveragePhenoDist,by=c("ML","Combos"))




# Remove NA's in the combos for alone
Fitmean$Combos<-as.character(Fitmean$Combos)
Fitmean[which(Fitmean$Trt=="Alone"),]$Combos<-"none"
Fitmean$Combos<-as.factor(Fitmean$Combos)

RelFitMean<-aggregate(SeedResiduals~ML+Combos+Trt,Fitmean,mean)

AveragePhenoDist<-merge(RelFitMean,AveragePhenoDist)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#               Examine evidence for character displacement 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Save an object listing the combinations with more than one bio rep 
Keep<-row.names(table(Ipurp$Combos)[(table(Ipurp$Combos)>1)])

# PCA 1
summary(lm(SeedResiduals~PhenDist_PCA1,AveragePhenoDist%>%filter(Combos%in%Keep)))

# PCA 2
summary(lm(SeedResiduals~PhenDist_PCA2,AveragePhenoDist%>%filter(Combos%in%Keep))) # Evidence of a *Negative* linear effect for PC2 (Root architecture)

  # Non linear examination w logistic curve
summary(lm(SeedResiduals~log(PhenDist_PCA2),AveragePhenoDist%>%filter(Combos%in%Keep))) # Evidence of a *Negative* linear effect for PC2 (Root architecture)

# PCA 4
summary(lm(SeedResiduals~PhenDist_PCA4,AveragePhenoDist%>%filter(Combos%in%Keep)))

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

PlotRegressionPD(Data=AveragePhenoDist,alpha=0.5,mapping=aes(log(PhenDist_PCA2),SeedResiduals),xTitle="",colorpoint = "black",colorline="black",size=5)+
  ylab("Relative fitness")+
  theme(axis.text.x = element_text(face="bold",  
                                   size=25),
        axis.text.y = element_text(face="bold", 
                                   size=25))+
  theme(axis.text.x = element_text( vjust = 1, hjust=1))+
  theme(axis.title = element_text( size=30))+
  xlab( expression(paste("Root architecture PC2", " (",sqrt(PC[paste("I.purpurea")]^{2} - PC[paste("I.hederacea")]^{2}),")")))


PlotRegressionPD(Data=AveragePhenoDist,alpha=0.5,mapping=aes(PhenDist_PCA4,SeedResiduals),xTitle="",colorpoint = "black",colorline="black",size=5)+
  ylab("Relative fitness")+
  theme(axis.text.x = element_text(face="bold",  
                                   size=25),
        axis.text.y = element_text(face="bold", 
                                   size=25))+
  theme(axis.text.x = element_text( vjust = 1, hjust=1))+
  theme(axis.title = element_text( size=30))+
  xlab( expression(paste("Root morphology PC4", " (",sqrt(PC[paste("I.purpurea")]^{2} - PC[paste("I.hederacea")]^{2}),")")))


# Plot the logistic model w root architecture

ggplot() +
  geom_point() +
  geom_smooth(data=AveragePhenoDist, aes(PhenDist_PCA2,SeedResiduals),method = lm, formula = y ~ log(x) , se =TRUE,fullrange=TRUE)



########################################################
#           Estimate Family means of traits
########################################################

# Subset for I purpurea species
IpPA4<-droplevels(TraitsAll%>%filter(Species=="Ip"))

# Calculate mean of traits by treatment and maternal line

pcFamilyMeans<-aggregate(list(IpPA4[c("PCA1","PCA2","PCA3","PCA4")]),by=list(IpPA4$Trt,IpPA4$ML),FUN=mean) #
pcFamilySE<-aggregate(list(IpPA4[c("PCA1","PCA2","PCA3","PCA4")]),by=list(IpPA4$Trt,IpPA4$ML),FUN=std.error) #

colnames(pcFamilyMeans)[1:2]<-c("Trt","ML")
colnames(pcFamilySE)[1:2]<-c("Trt","ML")
colnames(pcFamilySE)[grep("PC",names(pcFamilySE))]<-paste(colnames(pcFamilySE)[grep("PC",names(pcFamilySE))],"SE",sep="_")


for(i in c("PCA1","PCA2","PCA3","PCA4")){
  Quadratic=pcFamilyMeans[i]*pcFamilyMeans[i]
  pcFamilyMeans=cbind(pcFamilyMeans,Quadratic)
}

dim(pcFamilyMeans)

colnames(pcFamilyMeans)[7:10]<-paste(c("PCA1","PCA2","PCA3","PCA4"),"2",sep="_")

RelFitMean2<-aggregate(SeedNumberResid~ML+Trt,Fitmean,mean) # Average fitness by maternal line and treatment only

pcFamilyMeans<-merge(pcFamilyMeans,RelFitMean2)
PCAall=pcFamilyMeans



###################### ###################### ###################### ###################### ############
###############       Plotting phenotypic distance with variation/Standard Error            ############
###################### ###################### ###################### ###################### ############

CombPhenStdErrMean<-merge(RelFitMean,CombPhenStdErrMean)

# Combine phenotypic distance and fitness and plot w standard error


 # Lets color in the outliers
CombPhenStdErrMean$Outlier<-ifelse(CombPhenStdErrMean$PhenDist_PCA2>8,"Yes1","No")
CombPhenStdErrMean[CombPhenStdErrMean$PhenDist_PCA2>5&CombPhenStdErrMean$PhenDist_PCA2<8,"Outlier"]<-"Yes2"
Out=CombPhenStdErrMean[which(CombPhenStdErrMean$Outlier!="No"),"Combos"]

Fitmean$Outlier<-ifelse(Fitmean$Combos%in%Out[1],"Yes1","No")
Fitmean[which(Fitmean$Combos%in%Out[2]),"Outlier"]<-"Yes2"


ggplot(data=CombPhenStdErrMean,aes(Combos,PhenDist_PCA2,fill=Outlier))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=PhenDist_PCA2-PhenDist_PCA2_SE,ymax=PhenDist_PCA2+PhenDist_PCA2_SE),width=.4,                    # Width of the error bars
                position=position_dodge(.9))+
  theme_classic()+
  scale_fill_manual(values = c("grey","blue","gold"))+
  theme(axis.text.x=element_text(angle=90,size=12,vjust=0.5,color="black"),axis.text.y=element_text(size=12,color="black"),legend.position = "none",axis.title = element_text(size=20))+
  ylab("Phenotypic Distance PC2\n")+
  xlab("\nCombination pairing")


ggplot(data=Fitmean,aes(Combos,SeedNumberResid,fill=Outlier))+
  geom_boxplot()+
   theme_classic()+
  scale_fill_manual(values = c("grey","blue","gold"))+
  theme(axis.text.x=element_text(angle=90,size=12,vjust=0.5,color="black"),axis.text.y=element_text(size=12,color="black"),legend.position = "none",axis.title = element_text(size=20))+
  ylab("Relative Fitness\n")+
  xlab("\nCombination pairing")


# Count the number of Bio rep per maternal line by maternaline species combinations (i.e. Combos)

SampleNumber<-
  Both%>%
 count(Combos)

SampleNumberCombos<-data.frame(SampleNumber)


SampleNumber<-
 Fitmean%>%
  count(Combos)

SampleNumberCombosFit<-data.frame(SampleNumber)



################    Plotting confidence intervals around regression

# Combine raw phenotypic distance and average relative fitness to visualize variation 
BothFit=merge(RelFitMean,Both)
AveragePhenoDist<-merge(AveragePhenoDist,CombPhenStdErrMean[c("Combos","Outlier")])

BothFit<-merge(BothFit,AveragePhenoDist[c("Combos","Outlier")])


# Error around relative fitness
FitStandError<-aggregate(list(Fitmean[c("SeedNumberResid")]),by=list(Fitmean$Combos,Fitmean$ML),FUN=std.error)
colnames(FitStandError)[1:3]<-c("Combos","ML","FitSE")


AveragePhenoDist2<-merge(AveragePhenoDist,CombPhenStdErrMean[c("Combos","Outlier","PhenDist_PCA2_SE")])
AveragePhenoDist2<-merge(AveragePhenoDist2,FitStandError)

# Plot CI based on raw points
RegressError1<-ggplot()+
  #geom_point(data=BothFit,aes(PhenDist_PCA2,SeedResiduals,color=Outlier),size=5,alpha=0.2)+
  geom_point(data=AveragePhenoDist2,aes(PhenDist_PCA2,SeedResiduals,color=Outlier),size=5,alpha=0.5)+
  geom_smooth(data=AveragePhenoDist2,aes(PhenDist_PCA2,SeedResiduals),color="red",method="lm",size=1.2,se=F,linetype="dashed",fullrange=TRUE)+
  geom_errorbar(data=AveragePhenoDist2,aes(x=PhenDist_PCA2,ymin=SeedResiduals-FitSE,ymax=SeedResiduals+FitSE))+
  geom_errorbarh(data=AveragePhenoDist2,aes(y=SeedResiduals,xmin=PhenDist_PCA2-PhenDist_PCA2_SE,xmax=PhenDist_PCA2+PhenDist_PCA2_SE))+
  #geom_abline(slope=-0.0444,intercept = -0.0311,linetype="dashed",color="orange",size=1.5)+
  #geom_abline(slope=0.0019,intercept = 0.1071,linetype="dashed",color="orange",size=1.5)+
  #geom_smooth(data=AveragePhenoDist, aes(PhenDist_PCA2,SeedResiduals),method = lm, formula = y ~ log(x) , se =TRUE,fullrange=TRUE,color="blue")+
  theme_classic()+
  scale_color_manual(values = c("black","blue","gold"))+
  theme_classic()+
  ylab("Relative Fitness")+
  xlab( expression(paste("Root architecture PC2", " (",sqrt(PC[paste("I.purpurea")]^{2} - PC[paste("I.hederacea")]^{2}),")")))+
ggtitle(" ")+
  theme(axis.text.x = element_text(  
    size=20),
    axis.text.y = element_text(
      size=20),axis.title.y = element_text(size=20),
    axis.title.x = element_text(
      size=20))+
  theme(axis.text.x = element_text(vjust = 1, hjust=1,color="black"))+
  theme(axis.text.y = element_text(color="black"))+ 
  theme(plot.title=element_text(size=25,hjust=0.5))+
  theme(legend.position = "none")


RegressError1



########################################################
# Perform linear regression of PCs onto relative fitness
########################################################

# Selection Analysis
# Subset by treatment

PCAalone<-droplevels(pcFamilyMeans%>%filter(Trt=="Alone"))
PCAcomp<-droplevels(pcFamilyMeans%>%filter(Trt!="Alone"))


PC1.res.alone<-summary(lm(SeedNumberResid~PCA1,PCAalone))
PC1.res.alone
PC2.res.alone<-summary(lm(SeedNumberResid~PCA2,PCAalone))
PC2.res.alone
PC3.res.alone<-summary(lm(SeedNumberResid~PCA3,PCAalone))
PC3.res.alone # Root size--No evidence of selection. Not in the table.
PC4.res.alone<-summary(lm(SeedNumberResid~PCA4,PCAalone))
PC4.res.alone

PC1.res.comp<-summary(lm(SeedNumberResid~PCA1,PCAcomp)) 
PC1.res.comp
PC2.res.comp<-summary(lm(SeedNumberResid~PCA2,PCAcomp))
PC2.res.comp
PC3.res.comp<-summary(lm(SeedNumberResid~PCA3,PCAcomp))
PC3.res.comp # Root size--No evidence of selection. Not in the table.
PC4.res.comp<-summary(lm(SeedNumberResid~PCA4,PCAcomp%>%filter(PCA4>-1))) # Remove outlier w
PC4.res.comp # NOT reported
PC4.res.comp<-summary(lm(SeedNumberResid~PCA4,PCAcomp)) # w/o Remove outlier 
PC4.res.comp # Reported

## ANCOVA for PCA4 (Root morphology)

anova(lm(SeedNumberResid~PCA4*Trt,pcFamilyMeans))

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


########################################################
      # Plot selection on PC4 (root morphology)
########################################################


TraitsAllFit<-merge(TraitsAll,RelFitMean)
TraitsAllFit<-merge(TraitsAll,RelFitMean2)


IpPA4b<-merge(IpPA4,RelFitMean2)

pcFamilyMeans<-merge(pcFamilyMeans,pcFamilySE)
FitStandError2<-aggregate(list(Fitmean[c( "SeedNumberResid")]),by=list(Fitmean$Trt,Fitmean$ML),FUN=std.error) #Estimate Standard Error

colnames(FitStandError2)<-c("Trt","ML","Fit_SE")
pcFamilyMeans<-merge(FitStandError2,pcFamilyMeans)


#### Fig. 4 A and B

RegressAlone<-ggplot(pcFamilyMeans%>%filter(Trt=="Alone"))+
  geom_point(aes(PCA4, SeedNumberResid),size=5,alpha=0.5)+
  #geom_errorbar(aes(x=PCA4,ymin=SeedNumberResid-Fit_SE,ymax=SeedNumberResid+Fit_SE))+ # Uncomment to add error bars
  #geom_errorbarh(aes(y=SeedNumberResid,xmin=PCA4-PCA4_SE,xmax=PCA4+PCA4_SE))+ # Uncomment to add error bars
  scale_color_manual(values=c("red","black"))+
  stat_smooth(aes(PCA4, SeedNumberResid),alpha=0.5,method="lm", formula=y~x,se=F,fullrange = T,color="black")+
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

RegressAlone


RegressComp<-ggplot(data=pcFamilyMeans%>%filter(Trt!="Alone"))+
  geom_point(aes(PCA4, SeedNumberResid),size=5,alpha=0.5)+
  #geom_errorbar(aes(x=PCA4,ymin=SeedNumberResid-Fit_SE,ymax=SeedNumberResid+Fit_SE))+ # Uncomment to add error bars
  #geom_errorbarh(aes(y=SeedNumberResid,xmin=PCA4-PCA4_SE,xmax=PCA4+PCA4_SE))+ # Uncomment to add error bars
  stat_smooth(data=pcFamilyMeans%>%filter(Trt!="Alone"),aes(PCA4, SeedNumberResid),alpha=0.5,method="lm", formula=y~x,se=F,fullrange = T,color="black")+
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


RegressComp




#####################################################
  # Post-Hoc Linear mixed model on single traits
#####################################################


# Run for loop to save results of LMM on the traits that were significant from PC back regression
Vars<-c("trans.SKL_NODES", "trans.AVG_DENSITY", "trans.TD_AVG", "trans.RTP_COUNT", 
        "trans.STA_RANGE", "trans.STA_MAX", "trans.RTA_RANGE", "trans.ADVT_COUNT", 
        "trans.HYP_DIA", "trans.TAP_DIA", "trans.MAX_DIA_90","Species","ML","Trt","Block","Combos")


SubsetDirt=DirtP[c(Vars)]
dim(SubsetDirt)


df=SubsetDirt
df$Block<-as.factor(df$Block)

Response=as.list(names(SubsetDirt[1:11]))

for(i in 1:11){
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

# Does competitor impact basal root number?
df$Comp=sub(".*\\-","",df$Combos)


#####################################################
    # Post Hoc character displacement test
#####################################################

# Index traits of interest


Ipurp<-droplevels(subset(TraitsAll,Species=="Ip"))
Ipurp<-droplevels(subset(Ipurp,Trt=="Inter"))[c(Vars,"Combos","Position")]
Ihed<-droplevels(subset(TraitsAll,Species!="Ip"))[c(Vars,"Combos","Position" )]



# Merge by Position
colnames(Ihed)[1:15]<-paste(colnames(Ihed)[1:15],"_Competitor",sep="") # Modify competitor's column names
Both2<-merge(Ipurp,Ihed)

# Calculate euclidean distances

# Define a distance function


f<-function(x)
{
  return(mean(dist(x)))
}

# Estimates ALL pairwise* Distances 

Both2$PhenDist_SKL=apply(Both2[,grep("SKL_NODES",names(Both2))],1,f) # SKL NODES

Both2$PhenDist_Dstm=apply(Both2[,grep("DIA_STM",names(Both2))],1,f) # DIA STM

Both2$PhenDist_AvDen=apply(Both2[,grep("AVG_DENSITY",names(Both2))],1,f)

Both2$PhenDist_TDav=apply(Both2[,grep("TD_AVG",names(Both2))],1,f)

Both2$PhenDist_RTPcount=apply(Both2[,grep("RTP_COUNT",names(Both2))],1,f)

Both2$PhenDist_STA=apply(Both2[,grep("STA_RANGE",names(Both2))],1,f)

Both2$PhenDist_RTA=apply(Both2[,grep("RTA_RANGE",names(Both2))],1,f)

Both2$PhenDist_ADVTcount=apply(Both2[,grep("ADVT_COUNT",names(Both2))],1,f)

Both2$PhenDist_STAmax=apply(Both2[,grep("STA_MAX",names(Both2))],1,f)

Both2$PhenDist_TAPdia=apply(Both2[,grep("TAP_DIA",names(Both2))],1,f)

Both2$PhenDist_Hypdia=apply(Both2[,grep("HYP_DIA",names(Both2))],1,f)

Both2$PhenDist_90dia=apply(Both2[,grep("MAX_DIA_90",names(Both2))],1,f)



#Both2SE <-aggregate(list(Both2[c( "PhenDist_BASALcount")]),by=list(Both2$Combos),FUN=std.error) #Estimate Standard Error
#colnames(Both2SE)[1]<-c("Combos")
#colnames(Both2SE)[grep("Phen",names(Both2SE))]<-paste(colnames(Both2SE)[grep("Phen",names(Both2SE))],"SE",sep="_")



DistanceMeans<-aggregate(list(Both2[grep("PhenDist",names(Both2))]),by=list(Both2$Trt,Both2$ML,Both2$Combos),FUN=mean) # Mean phenotypic distance

colnames(DistanceMeans)[1:3]=c("Trt","ML","Combos") # Rename coloumns



# Merge to fitness averaged at level of combo

PhenDistFit=merge(RelFitMean,DistanceMeans)
PhenDistFit<-droplevels(PhenDistFit%>%filter(Trt=="Inter"))


#Both2a<-merge(PhenDistFit,Both2SE,by="Combos")
#Both2a<-unique(droplevels(Both2a))


# Lets color in the outliers
#Both2a$Outlier<-ifelse(Both2a$PhenDist_BASALcount>2.4,"Yes","No")
#Out2=Both2a[which(Both2a$Outlier=="Yes"),"Combos"]

#Fitmean$Outlier2<-ifelse(Fitmean$Combos%in%Out2[1],"Yes","No")




# Perform linear regression
summary(lm(SeedResiduals~PhenDist_SKL,PhenDistFit))
#summary(lm(SeedResiduals~PhenDist_Dstm,PhenDistFit))
summary(lm(SeedResiduals~PhenDist_AvDen,PhenDistFit))
summary(lm(SeedResiduals~PhenDist_TDav,PhenDistFit))
summary(lm(SeedResiduals~PhenDist_RTPcount,PhenDistFit))
summary(lm(SeedResiduals~PhenDist_STA,PhenDistFit)) 
summary(lm(SeedResiduals~PhenDist_RTA,PhenDistFit))
summary(lm(SeedResiduals~PhenDist_ADVTcount,PhenDistFit))
summary(lm(SeedResiduals~PhenDist_STAmax,PhenDistFit))# 
summary(lm(SeedResiduals~PhenDist_Hypdia,PhenDistFit))# 
summary(lm(SeedResiduals~PhenDist_90dia,PhenDistFit))# 


summary(lm(SeedResiduals~PhenDist_TAPdia,PhenDistFit))



##################################################################################
      # Post=Hoc Analysis for Combination effect on root architecture
##################################################################################

# PC 2--at the level of a combination BLOCK removed


TraitsAll$Comp<-sub(".*\\-","",TraitsAll$Combos)
PC2model2<-lm(PCA2~Trt+Comp+Combos,TraitsAll) # Marginal competitor effect, significant Ml by Ml effect
anova(PC2model2)

PC1model2<-lm(PCA1~Trt+Comp+Combos,TraitsAll) # Marginal competitor effect, significant Ml by Ml effect
anova(PC1model2) # Marginal combination effect

PC3model2<-lm(PCA3~Trt+Comp+Combos,TraitsAll) # Marginal competitor effect, significant Ml by Ml effect
anova(PC3model2) # Marginal combination effect


PC4model2<-lm(PCA4~Trt+Comp+Combos,TraitsAll) # Marginal competitor effect, significant Ml by Ml effect
anova(PC4model2) # Significant treatment effect and combination effect


#########################################
    # Fig. S5 Additional (I.hed plot)
#########################################



sizeIhed<-size%>%
  filter(Species=="Ihed")%>%
  group_by(ML,Combos)%>%
  summarise("MeanSize"=mean(Leaf.Number),"SE_Size"=std.error(Leaf.Number))

PhenDistIhed<-Both%>%
  #filter(Species=="Ihed")%>%
  group_by(Combos)%>%
  summarise("MeanPhenDistPC2"=mean(PhenDist_PCA2),"SE_PhenDist"=std.error(PhenDist_PCA2))

HedData<-merge(sizeIhed,PhenDistIhed)


summary(lm(MeanSize~log(MeanPhenDistPC2),HedData%>%filter(MeanSize<20))) # marginally significant when we filter our the 
# size greater than 20; significant LOG relationship

hedModel<-lm(PCA2~Combos,droplevels(TraitsAll%>%filter(Species=="Ihed")))
anova(hedModel) # Not significant

ggplot(data=HedData%>%filter(MeanSize<20),aes(MeanPhenDistPC2,MeanSize))+
  geom_point()+
  geom_errorbar(aes(ymin=MeanSize-SE_Size,ymax=MeanSize+SE_Size))+
  geom_errorbarh(aes(xmin=MeanPhenDistPC2-SE_PhenDist,xmax=MeanPhenDistPC2+SE_PhenDist))+
  theme_classic()+
  scale_fill_manual(values = c("grey","blue","gold"))+
  theme(axis.text.x=element_text(angle=90,size=12,vjust=0.5,color="black"),axis.text.y=element_text(size=12,color="black"),legend.position = "none",axis.title = element_text(size=20))+
  ylab("Plant Size\n")+
  xlab("\nPhenotypic Distance Root Architecture (PC2)")



#########################################
# Plot phenotypic distance measure of root architecture vs raw value of trait
#########################################

plot(Both$PCA2,Both$PhenDist_PCA2)
abline(Both$PCA2,Both$PhenDist_PCA2)

summary(lm(PhenDist_PCA2 ~ PCA2, Both))







