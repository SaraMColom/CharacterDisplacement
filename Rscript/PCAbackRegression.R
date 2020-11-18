# Load libraries

library(factoextra)
library(FactoMineR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

########################################################
### Set up
########################################################
setwd("../CleanData/")


# Read in data
Fit<-read.csv("FitPA4.csv")
# Root trait data
tot2<-read.csv("totPA4.csv")
# Size
size<-read.csv("SizeData.csv",na.strings = c("."," "))


# Correct data structure
size$Block=as.factor(size$Block)
Fit$Block=as.factor(Fit$Block)
tot2$Block=as.factor(tot2$Block)

# Subset for I. purpurea
Purp<-droplevels(subset(size,size$Species=="Ip"&size$Population=="PA4"))


### TRANSFORM Root trait data, box-cox, center, scale

# Exlude these columns from transformation

EXCLUDE<-c("Position", "Species","EXP_CODE", "Image_Name","IMAGE_ID", "CIR_RATIO", "X_PIXEL", "Full.Path","Name","Type","Date","URL","Last.Updated","Folder_Name","UniqId.x","File_copy","UniqId.y",
           "Y_PIXEL", "X_SCALE", "Y_SCALE", "COMP_TIME", "Comment", "Block",
           "Order", "ML", "GerminationDate", "Comment1", "Dead_plant", "DeathCause",
           "RootsHarvested", "SeedsCounted", "Trt", "Comp", "Combos", "Population"
) 

########################################################
# Transform data
########################################################

predictors<-tot2[-which(names(tot2)%in%EXCLUDE)]
library(caret)
trans = preProcess(predictors,
                   c("BoxCox", "center", "scale"))
predictorsTrans = data.frame(
  trans = predict(trans, predictors))
#Incorporate factor variables to the transformed data set
totTrans<-predictorsTrans
totTrans$Species<-tot2$Species
totTrans$Population<-tot2$Population
totTrans$Block<-tot2$Block
totTrans$ML<-tot2$ML
totTrans$Trt<-tot2$Trt
totTrans$Position<-tot2$Position
totTrans$Combos<-tot2$Combos
Focus<-paste("trans.",names(tot2[,c(1:33)]),sep="")
# Make factor variables factors
totTrans[c("Population","Species","ML","Combos","Block","Trt")]<-lapply(totTrans[c("Population","Species","ML","Combos","Block","Trt")],as.factor)


# Remove Block effects
Data1<-totTrans
Data1[1:33]<-NA
for(i in 1:33) {
  Residual<-lm(totTrans[,i]~Block,totTrans,na.action="na.exclude")$residuals
  Data1[names(Residual),i]<-Residual
}
BlkRmv<-Data1 # Focal traits with block effect removed.


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

p<-fviz_pca_ind(res.pca, label="none", habillage=BlkRmvFull$Species,
                addEllipses=TRUE, ellipse.level=0.95)+
  theme_bw()

pca<-p + scale_color_manual(values=c('#999999','#E69F00'))+
  theme_classic()+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  ggtitle("")

pca

# Plot the top 30 contributing individuals
fviz_pca_biplot(res.pca, label="var",col.var="contrib",repel=T,
                select.var = list(contrib = 10))+
  scale_color_gradient2(low="green", mid="red",
                        high="red", midpoint=96) +
  theme_minimal()

# Plot the ten top traits contributing to ea PC for the top 4 PC's
# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10) # PC1
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10) # PC2
fviz_contrib(res.pca, choice = "var", axes = 3, top = 10) # PC3
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 4, top = 10) # PC4


# Estimate the EigenVectors of each trait
library(Rcmdr)

# Extract the loading score of each trait


Loads<-sweep(res.pca$var$coord,2,sqrt(res.pca$eig[1:ncol(res.pca$var$coord),1]),FUN="/") 

LoadsDf<-data.frame(Loads)
LoadsDf$Traits<-gsub("trans.","",x=row.names(Loads))

PC1Loads<-ggplot(LoadsDf)+
  geom_bar(aes(Traits,Dim.1),stat="identity")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90,color="black",size=15),plot.title = element_text(hjust=0.5,size=20))+
  ylab(" ")+
  xlab(" ")+
  ggtitle("PC1")

PC2Loads<-ggplot(LoadsDf)+
  geom_bar(aes(Traits,Dim.2),stat="identity")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90,color="black",size=15),plot.title = element_text(hjust=0.5,size=20))+
  ylab(" ")+
  xlab(" ")+
  ggtitle("PC2")


PC3Loads<-ggplot(LoadsDf)+
  geom_bar(aes(Traits,Dim.3),stat="identity")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90,color="black",size=15),plot.title = element_text(hjust=0.5,size=20))+
  ylab(" ")+
  xlab(" ")+
  ggtitle("PC3")

PC4Loads<-ggplot(LoadsDf)+
  geom_bar(aes(Traits,Dim.4),stat="identity")+
    theme_classic()+
  theme(axis.text.x=element_text(angle=90,color="black",size=15),plot.title = element_text(hjust=0.5,size=20))+
  ylab(" ")+
  xlab(" ")+
  ggtitle("PC4")

library(grid)


Grid=cowplot::plot_grid(PC1Loads,PC2Loads,PC3Loads,PC4Loads,nrow=2,ncol=2,align="hv")

x.grob<-textGrob("Traits",gp=gpar(col=
                                  "black",fontsize=25),rot=0)

y.grob<-textGrob("Loading Score",gp=gpar(col="black",fontsize=25),rot=90)


gridExtra::grid.arrange(gridExtra::arrangeGrob(Grid,bottom=x.grob,left=y.grob,padding=unit(1, units='in'),nrow=1))

# Making a copy of the traits data with block effects removed

TraitsAll<-BlkRmvFull

# 1) Subset between Ipurp and Ihed

TraitsAll$PCA1<-res.pca$ind$coord[,1]
TraitsAll$PCA2<-res.pca$ind$coord[,2]
TraitsAll$PCA3<-res.pca$ind$coord[,3]
TraitsAll$PCA4<-res.pca$ind$coord[,4]

########################################################
# Calculate relative fitness
########################################################

# Reduce Purp to combine to leaf fitness data
Purp.Red<-Purp[c("Position","Leaf.Number","Species")]
colnames(Purp.Red)[3]<-"ID"

Fit2<-merge(Purp.Red,Fit)

SdMnSpeciesTrt<- aggregate(SeedNumber ~ Species+Trt, Fit2, mean) # Mean residual
names(SdMnSpeciesTrt)[3]<- "MeanSdNm" #Rename coloumn for mean seed number

# Merge average fitness by species and treatment 
Fitmean<-merge(Fit2,SdMnSpeciesTrt,by=c("Species","Trt"))

# Calculate relative fitness as the observed seed number by the total mean seed number of that species and treatment.
Fitmean$Rel_Fit<-Fitmean$SeedNumber/Fitmean$MeanSdNm

Fitmean$Combos<-as.character(Fitmean$Combos)
Fitmean[which(Fitmean$Trt=="Alone"),"Combos"]<-"none"
Fitmean$Combos<-as.factor(Fitmean$Combos)

# Extract residuals of block and size
Fitmean$SeedNumberResid<-NA
SeedResiduals<-(lm(Rel_Fit~Block+Leaf.Number,Fitmean))$residuals
Fitmean[names(SeedResiduals),"SeedNumberResid"]<-SeedResiduals


# All looks good. Use SeedNumberResid
plot(Fitmean$SeedNumber,Fitmean$Rel_Fit)
plot(Fitmean$SeedNumberResid,Fitmean$Rel_Fit)

RelFitMean<-aggregate(SeedNumberResid~ML+Trt,Fitmean,mean) # Average relative fitness by treatment and maternal line


########################################################
# Perform linear regression of PCs onto relative fitness
########################################################

# Obtain the family and treatment means of individual PC scores

# Subset for I purpurea species
IpPA4<-droplevels(TraitsAll%>%filter(Species=="Ip"))

# Calculate mean of traits by treatment and maternal line

pcFamilyMeans<-aggregate(list(IpPA4[c("PCA1","PCA2","PCA3","PCA4")]),by=list(IpPA4$Trt,IpPA4$ML),FUN=mean) #

colnames(pcFamilyMeans)[1:2]<-c("Trt","ML") # Rename columns


# Combine relative fitness and family means of the PC's

pcFamilyMeans<-merge(pcFamilyMeans,RelFitMean)
PCAall=pcFamilyMeans

pcFamilyMeans=pcFamilyMeans


PCAalone<-droplevels(pcFamilyMeans%>%filter(Trt=="Alone"))
PCAcomp<-droplevels(pcFamilyMeans%>%filter(Trt!="Alone"))

# Estimate the selection gradient for each treatment and PC trait
PC1.res.alone<-summary(lm(SeedNumberResid~PCA1,PCAalone))
PC2.res.alone<-summary(lm(SeedNumberResid~PCA2,PCAalone))
PC3.res.alone<-summary(lm(SeedNumberResid~PCA3,PCAalone)) 
PC4.res.alone<-summary(lm(SeedNumberResid~PCA4,PCAalone))


PC1.res.comp<-summary(lm(SeedNumberResid~PCA1,PCAcomp)) 
PC2.res.comp<-summary(lm(SeedNumberResid~PCA2,PCAcomp))
PC3.res.comp<-summary(lm(SeedNumberResid~PCA3,PCAcomp))  
#PC4.res.comp<-summary(lm(SeedNumberResid~PCA4,PCAcomp%>%filter(PCA4>-1)))
PC4.res.comp<-summary(lm(SeedNumberResid~PCA4,PCAcomp))


# Subset the linear coefficients and standar errors from selection analysis

AloneResults<-list(PC1.res.alone,PC2.res.alone,PC3.res.alone,PC4.res.alone)
CompResults<-list(PC1.res.comp,PC2.res.comp,PC3.res.comp,PC4.res.comp)



Empty=list()
for (i in AloneResults){
  for (j in c(1,2,4)){
coefficient<-(AloneResults[[j]])$coefficients[2,1]
Empty[j]<-coefficient
  }
}

SelGradAlone<-do.call('rbind',Empty)

Empty=list()
for (i in CompResults){
  for (j in c(1,2,4)){
    coefficient<-(CompResults[[j]])$coefficients[2,1]
    Empty[j]<-coefficient
  }
}

SelGradCompetition<-do.call('rbind',Empty)


# Repeat steps for standard error

Empty=list()
for (i in AloneResults){
  for (j in c(1,2,4)){
    coefficient<-(AloneResults[[j]])$coefficients[2,2]
    Empty[j]<-coefficient
  }
}

StErAlone<-do.call('rbind',Empty)

Empty=list()
for (i in CompResults){
  for (j in c(1,2,4)){
    coefficient<-(CompResults[[j]])$coefficients[2,2]
    Empty[j]<-coefficient
  }
}

StErCompetition<-do.call('rbind',Empty)



########################################################
# Do back regression
########################################################

# Eigenvectors multiplied by the vector of selection gradient
EigenVectors<-Loads[,c(1,2,4)]

# For alone treatment
SelectionCoef.Alone.BR<-EigenVectors %*% SelGradAlone

# For alone treatment
SelectionCoef.Comp.BR<-EigenVectors %*% SelGradCompetition

# Confidence Intervals

# 1. get the square value of eigenvectors
EigenVecSq<-EigenVectors*EigenVectors

# 2. get the square value of standard errors
StErrSqAlone<-StErAlone*StErAlone
StErrSqComp<-StErCompetition*StErCompetition

# Matrix multiplication between standard error squared and eigenvectors squared, and then square root the output
AloneSE.BR<-sqrt((EigenVecSq%*%StErrSqAlone))
CompSE.BR<-sqrt((EigenVecSq%*%StErrSqComp))


# Alone
AloneBR.Res<-data.frame(cbind(SelectionCoef.Alone.BR,AloneSE.BR))
colnames(AloneBR.Res)<-c("Beta","SE")
AloneBR.Res$tValue=AloneBR.Res$Beta/AloneBR.Res$SE

#Competition
CompBR.Res<-data.frame(cbind(SelectionCoef.Comp.BR,CompSE.BR))
colnames(CompBR.Res)<-c("Beta","SE")
CompBR.Res$tValue=CompBR.Res$Beta/CompBR.Res$SE


AloneBR.Res$Trt="Alone"
CompBR.Res$Trt="Comp"

AloneBR.Res$Traits<-row.names(AloneBR.Res)
CompBR.Res$Traits<-row.names(CompBR.Res)

Total.BR.Res<-rbind(AloneBR.Res,CompBR.Res)

#write.csv(Total.BR.Res,"ResultsBackRegression.csv",row.names = F)

### Estimate significance with confidence intervals of 2 SE

## 2 SE
Total.BR.Res$SE2=Total.BR.Res$SE*1.96
Total.BR.Res$Upp=Total.BR.Res$Beta+Total.BR.Res$SE2
Total.BR.Res$Low=Total.BR.Res$Beta-Total.BR.Res$SE2

DT::datatable(Total.BR.Res)
DT::datatable(Total.BR.Res[-which(Total.BR.Res$Upp>0 & Total.BR.Res$Low<0),])
Significant<-Total.BR.Res[-which(Total.BR.Res$Upp>0 & Total.BR.Res$Low<0),]

Save<-Total.BR.Res[which(Total.BR.Res$Traits%in%Significant$Traits),]
Save$Traits=row.names(Save)
Save$Traits<-gsub("1","",Save$Traits)
DT::datatable(Save[c(1,2,4,7,8)])


