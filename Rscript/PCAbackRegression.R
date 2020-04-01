library(factoextra)
library(FactoMineR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

########################################################
### Set up
########################################################



# Read in data
Fit<-read.csv("~/Desktop/CharacterDisplacementRootTraits/CleanData/FitPA4.csv")# Root trait data
tot2<-read.csv("~/Desktop/CharacterDisplacementRootTraits/CleanData/totPA4.csv")
# Size/Fitness data
size<-read.csv("~/Desktop/CharacterDisplacementRootTraits/CleanData/SizeFitData.csv",na.strings=c("."," "))


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

# Select the top 30 contributing individuals
fviz_pca_biplot(res.pca, label="var",col.var="contrib",repel=T,
                select.var = list(contrib = 10))+
  scale_color_gradient2(low="green", mid="red",
                        high="red", midpoint=96) +
  theme_minimal()


# Alternative eigen value and vector
Res.cor<-cor(BlkRmvFull[-Qually])
traitsNames<-row.names(Res.cor)
eigen<-eigen(Res.cor)
eigen$values


# EigenVectors
res.pca1<-prcomp(BlkRmvFull[-Qually],center=F) # the prcomp funciton provideds eigenvectors as output

Loads<-data.frame(res.pca1$rotation)
# write.csv(Loads,"LoadingScoresPCA.csv",row.names = F)


# Obtain the family and treatment means of individual PC scores

TraitsAll<-BlkRmvFull

# 1) Subset between Ipurp and Ihed

TraitsAll$PCA1<-res.pca$ind$coord[,1]
TraitsAll$PCA2<-res.pca$ind$coord[,2]
TraitsAll$PCA3<-res.pca$ind$coord[,3]
TraitsAll$PCA4<-res.pca$ind$coord[,4]

aggregate(PCA1~Species,TraitsAll,mean)
aggregate(PCA2~Species,TraitsAll,mean)
aggregate(PCA3~Species,TraitsAll,mean)
aggregate(PCA4~Species,TraitsAll,mean)

aggregate(PCA1~Species,TraitsAll,sd)
aggregate(PCA2~Species,TraitsAll,sd)
aggregate(PCA3~Species,TraitsAll,sd)
aggregate(PCA4~Species,TraitsAll,sd)


########################################################
# Perform linear regression of PCs onto relative fitness
########################################################


pcFamilyMeans<-read.csv("~/Desktop/CharacterDisplacementRootTraits/CleanData/pcFamilyMeans.csv")

PCAalone<-droplevels(pcFamilyMeans%>%filter(Trt=="Alone"))
PCAcomp<-droplevels(pcFamilyMeans%>%filter(Trt!="Alone"))


PC1.res.alone<-summary(lm(SeedNumberResid~PCA1,PCAalone))
PC2.res.alone<-summary(lm(SeedNumberResid~PCA2,PCAalone))
PC3.res.alone<-summary(lm(SeedNumberResid~PCA3,PCAalone))
PC4.res.alone<-summary(lm(SeedNumberResid~PCA4,PCAalone))


PC1.res.comp<-summary(lm(SeedNumberResid~PCA1,PCAcomp))
PC2.res.comp<-summary(lm(SeedNumberResid~PCA2,PCAcomp))
PC3.res.comp<-summary(lm(SeedNumberResid~PCA3,PCAcomp))
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
EigenVectors<-as.matrix(Loads[c(1,2,4)])

# For alone treatment
SelectionCoef.Alone.BR<-EigenVectors %*% SelGradAlone

# For alone treatment
SelectionCoef.Comp.BR<-EigenVectors %*% SelGradCompetition

#   Estimate beta and Standard Error

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

#Competition
CompBR.Res<-data.frame(cbind(SelectionCoef.Comp.BR,CompSE.BR))
colnames(CompBR.Res)<-c("Beta","SE")



AloneBR.Res$Trt="Alone"
CompBR.Res$Trt="Comp"

AloneBR.Res$Traits<-row.names(AloneBR.Res)
CompBR.Res$Traits<-row.names(CompBR.Res)

Total.BR.Res<-rbind(AloneBR.Res,CompBR.Res)


### Estimate significance with confidence intervals of 2 SE

## 2 SE
Total.BR.Res$SE2=Total.BR.Res$SE*1.96
Total.BR.Res$Upp=Total.BR.Res$Beta+Total.BR.Res$SE2
Total.BR.Res$Low=Total.BR.Res$Beta-Total.BR.Res$SE2

DT::datatable(Total.BR.Res)
DT::datatable(Total.BR.Res[-which(Total.BR.Res$Upp>0 & Total.BR.Res$Low<0),])


