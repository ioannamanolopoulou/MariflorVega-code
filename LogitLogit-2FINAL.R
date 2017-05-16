#Things to do 
#Model Selection - Stepwise using AIC
library("pROC")


PresBack = function(AllData) #P/BG logistic regression
{
  ## Define Ntilde which has 0/1 as N, rather than -1/1
  AllData$Ntilde = AllData$N 
  AllData$Ntilde[AllData$Ntilde==-1] = 0
  ##
  
  #take a random sample of size N rather than all background points. 
  AllData[AllData$Ntilde==0,] = AllData[sample(which(AllData$Ntilde==0),sum(AllData$Ntilde==0),replace=TRUE),]
  #fit logistic regression
  lm1 = glm(Ntilde ~ elev + grad + flora, family=binomial(link="logit"),data=AllData) 
  #  lm1 = glm(Ntilde ~ elev + grad + flora + x + I(x^2) + y + I(y^2),family=binomial,data=AllData) 
  #show coefficients
  cat(lm1$coefficients)
  #print(summary(lm1))
  return(lm1)
}


PresAbs = function(AllData) #bernoulli GLM with complementary log log link - same as IPP model for PA data
{
  lm1 = glm(N ~ elev + grad,family=binomial(link="cloglog"),data=AllData) #IPP model loglog link
  #  lm1 = glm(N ~ elev + grad + x + I(x^2) + y + I(y^2),family=binomial(link="cloglog"),data=AllData)
  cat(lm1$coefficients)
  #  print(summary(lm1))
  return(lm1)
}

#par(mfrow=c(3,2))
beta.initial = c(0,0,0)
library(plotly)
library(spatstat)
load("AllData.RData")

nrow(PresenceAbsence)
set.seed(23948)
test.rows<-sample(1:nrow(PresenceAbsence),2131)
PresAbsTest<-PresenceAbsence[test.rows,]
PresAbsTrain<-PresenceAbsence[-test.rows,]

#plot the points in the test area
AllData=PresAbsTest
plot(AllData$x[AllData$N==1],AllData$y[AllData$N==1],col='red',cex=0.5,pch=20,xlab="x",ylab="y",main="red:presence,\n blue:absence")
points(AllData$x[AllData$N==0],AllData$y[AllData$N==0],col='blue',cex=0.5,pch=20)


#SAVE BG PIXELS FOR PLOTTING
cat("RUNNING UNBIASED (COMPLETE) DATA")
# 3a
cat("\nRunning Pres/BG, unbiased\n")
AllData = PresenceBackground

#### save pixel background image for plotting
totalrows = length(unique(AllData$y))
totalcols = length(unique(AllData$x))
map.predict = AllData[AllData$N==-1,]

#totalrows = length(unique(PresAbsTest$y))
#totalcols = length(unique(PresAbsTest$x))
#map.predict = PresAbsTest

#####


#PRESENCE-ABSENCE MODEL
#UNBIASED- GOLD STANDARD
cat("\nRunning Pres/Abs, unbiased\n")
#IPP MODEL (cloglog bernoulli glm)
AllData = PresAbsTrain
lm1 = PresAbs(AllData)
lm1$coefficients
pred<-predict(lm1,PresAbsTest, type="response")
range(pred)
auc(PresAbsTest$N,pred)
#0.6179
####

#### plot heatmap
PredictSpatial = matrix(predict(lm1,map.predict, type="response"),nrow = totalrows, ncol = totalcols)
print(plot_ly(z = PredictSpatial, type = "heatmap"))
#"true heatmap"
######

#### plot points
plot(AllData$x[AllData$N==0],AllData$y[AllData$N==0],col='red',cex=0.5,pch=20,xlab="x",ylab="y",main="red:absence,\n blue:presence")
points(AllData$x[AllData$N==1],AllData$y[AllData$N==1],col='blue',cex=0.5,pch=20)
#"true occurances/ absences - true dist"
######

#PRESENCE_BG IPP MODELS
#### fit point process directly to points versus background
cat("\nPoint process\n")

#ALL BG POINTS
fitPP<-ppm(bei, ~ elev + grad + flora,data=bei.extra,Poisson()) 
cat(fitPP$coef)
plot(fitPP)
auc(fitPP)

#obs = 0.6138, theo = 0.5953


fitPPuc<-ppm(bei,~elev+grad, data=bei.extra, Poisson())
fitPPuc
cat(fitPPuc$coef)
auc(fitPPuc)
#obs= 0.6138, theo = 0.5952

#######




#PRESENCE-BG LOGISTIC REGRESSION
cat("\n\nRUNNING BIASED DATA, BUT BIAS INDEPENDENT OF PROCESS")

#INDEPT BIAS
cat("\nRunning Pres/BG, biased\n")
AllData = PresenceBackground.biased1
lm1 = PresBack(AllData)

pred<-predict(lm1, PresAbsTest, type="response")
?predict.glm
range(pred) #0.0311 - #0.3717
auc(PresAbsTest$N,pred)
#0.5767

PredictSpatial = matrix(predict(lm1, map.predict, type="response"),nrow = totalrows, ncol = totalcols)
print(plot_ly(z = PredictSpatial, type = "heatmap"))

plot(AllData$x[AllData$N==-1],AllData$y[AllData$N==-1],col='red',cex=0.5,pch=20,xlab="x",ylab="y",main="red:background,\n blue:presence")
points(AllData$x[AllData$N==1],AllData$y[AllData$N==1],col='blue',cex=0.5,pch=20)


#COMPLEX BIAS
cat("\nRunning pres/BG, very biased\n")
AllData = PresenceBackground.biased2
lm1 = PresBack(AllData)

pred<-predict(lm1,PresAbsTest, type="response")
range(pred)
roc(PresAbsTest$N,pred)
#0.6214

PredictSpatial = matrix(predict(lm1, map.predict, type="response"),nrow = totalrows, ncol = totalcols)
print(plot_ly(z = PredictSpatial, type = "heatmap"))

plot(AllData$x[AllData$N==-1],AllData$y[AllData$N==-1],col='red',cex=0.5,pch=20,xlab="x",ylab="y",main="red:background,\n blue:presence")
points(AllData$x[AllData$N==1],AllData$y[AllData$N==1],col='blue',cex=0.5,pch=20)


















##########################

# 4b
cat("\nRunning Pres/Abs, biased\n")
AllData = PresenceAbsence.biased1

lm1 = PresAbs(AllData)
set.seed(21029)
for (i in 1:k){
  train<-AllData[folds$subsets[folds$which !=i],]
  test<-AllData[folds$subsets[folds$which ==i],]
  
  lm1<-PresAbs(train)
  
  pred<-predict(lm1,test, type="response")
  
}

roc(test$N,pred)





###
PredictSpatial = matrix(predict(lm1,map.predict),nrow = totalrows, ncol = totalcols)
print(plot_ly(z = PredictSpatial, type = "heatmap"))
###

###
plot(AllData$x[AllData$N==0],AllData$y[AllData$N==0],col='red',cex=0.5,pch=20,xlab="x",ylab="y",main="red:absence,\n blue:presence")
points(AllData$x[AllData$N==1],AllData$y[AllData$N==1],col='blue',cex=0.5,pch=20)
###

#cat("\nPoint process\n")
#fitPP<-ppm(bei, ~ elev + grad + flora,data=bei.extra,Poisson()) 
#cat(fitPP$coef)

cat("\n\nRUNNING BADLY BIASED DATA")
# 3b
load("AllData.RData")


##
PredictSpatial = matrix(predict(lm1,map.predict),nrow = totalrows, ncol = totalcols)
print(plot_ly(z = PredictSpatial, type = "heatmap"))
##

###
plot(AllData$x[AllData$N==-1],AllData$y[AllData$N==-1],col='red',cex=0.5,pch=20,xlab="x",ylab="y",main="red:background,\n blue:presence")
points(AllData$x[AllData$N==1],AllData$y[AllData$N==1],col='blue',cex=0.5,pch=20)

# 4b
cat("\nRunning pres/Abs, very biased\n")
AllData = PresenceAbsence.biased2
lm1 = PresAbs(AllData)

##
PredictSpatial = matrix(predict(lm1,map.predict),nrow = totalrows, ncol = totalcols)
print(plot_ly(z = PredictSpatial, type = "heatmap"))
##

##
plot(AllData$x[AllData$N==0],AllData$y[AllData$N==0],col='red',cex=0.5,pch=20,xlab="x",ylab="y",main="red:absence,\n blue:presence")
points(AllData$x[AllData$N==1],AllData$y[AllData$N==1],col='blue',cex=0.5,pch=20)
##

#cat("\nPoint process\n")
#fitPP<-ppm(bei, ~ elev + grad + flora,data=bei.extra,Poisson()) 
#cat(fitPP$coef)


#install.packages("devtools")
#library(devtools)
#install_github("wfithian/multispeciesPP")
library(multispeciesPP)
mod<-multispeciesPP(~x+y+elev+flora,~grad+flora, PA=PresenceAbsence, BG=PresenceBackground,PO=PresenceOnly)
PresenceOnly<-subset(PresenceBackground, PresenceBackground$N==1)
?multispeciesPP

#########
#ten fold cross validation
#install.packages("cvTools")
#library(cvTools)
#k<-10
#folds<-cvFolds(NROW(AllData),K=k)
#library(pROC)
#set.seed(20399)
#for (i in 1:k){
#train<-AllData[folds$subsets[folds$which !=i],]
#test<-AllData[folds$subsets[folds$which ==i],]

#lm1<-PresAbs(train)

#pred<-predict(lm1,newdata=test, type="response")

#}

#roc(test$N,pred)
#test
#PresAbsTest<-merge(test,PresenceAbsence, by.x=c("x","y","elev","grad","flora"), by.y=c("x","y","elev","grad","flora" ),sort=F,all.x=T)
#PresAbsTest<-nrow(unique(PresAbsTest$N.x))
#auc(pred,PresAbsTest$N.y)
#0.6144
#why so low
##########
