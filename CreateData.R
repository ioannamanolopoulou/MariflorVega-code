library(data.table)
library(spatstat)

setwd('~/Dropbox/Teaching/Student Projects/UG projects/Nattanit Srisamrual/')


##### IGNORE FROM HERE #####
gridwidth = 5

bei.extra$flora = bei.extra$elev
bei.extra$flora$v = matrix(runif(prod(dim(bei.extra$flora$v))),nrow=nrow(bei.extra$flora$v))

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

dt<- data.table(x=bei$x, y=bei$y)
dt[, xc:= as.numeric(as.character(cut(x, breaks=seq(-gridwidth/2,1000 + gridwidth/2,by=gridwidth), labels=seq(0,1000,by=gridwidth))))]
dt[, yc:= as.numeric(as.character(cut(y, breaks=seq(-gridwidth/2,500 + gridwidth/2,by=gridwidth), labels=seq(0,500,by=gridwidth))))]
dt<- dt[, list(N=length(x)),by=c('xc','yc')]

dg<- as.data.table(expand.grid(xc=seq(0,1000,by=gridwidth) , yc=seq(0,500,by=gridwidth)))
dg[, xi:= xc/gridwidth+1]
dg[, yi:= yc/gridwidth+1]
tmp<- dg[, list(grad= bei.extra$grad$v[yi, xi]), by=c('xc','yc')]
dg<- merge(dg, tmp, by=c('xc','yc'))
tmp<- dg[, list(elev= bei.extra$elev$v[yi, xi]), by=c('xc','yc')]
dg<- merge(dg, tmp, by=c('xc','yc'))
tmp<- dg[, list(flora= bei.extra$flora$v[yi, xi]), by=c('xc','yc')]
dg<- merge(dg, tmp, by=c('xc','yc'))

dt<- merge(dg, dt, by=c('xc','yc'), all.x=1)
set(dt, dt[, which(is.na(N))], 'N', 0)

dt$xi = NULL
dt$yi = NULL
##### TO HERE #####

##dt now contains xc,yc,grad,elev,flora,N, where xc and yc are now the mid-points of the
# pixel where each observation belongs. 

####################################################################################
####################################################################################
####################################################################################
##### COMPLETE OBSERVATION #######
####################################################################################
####################################################################################
####################################################################################

####PRESENCE-ABSENCE

nobs = sum(dt$N) + sum(dt$N==0)
#nobs = sum(dt$N) + sum(dt$N==0) + length(dt$N)
PresenceAbsence = array(0,dim=c(nobs,6))
counter = 1
for(i in 1:length(dt$N))
{
 # PresenceAbsence[counter,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],-1)
#  counter = counter + 1
  
  if(dt$N[i] == 0)
  {
    PresenceAbsence[counter,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],dt$flora[i],0)
    counter = counter + 1
  }
  if(dt$N[i] > 0)
  {
    for(j in counter:(counter + dt$N[i] -1))
    {
      PresenceAbsence[j,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],dt$flora[i],1)
    }
    counter = counter + dt$N[i]
  }
}

PresenceAbsence = as.data.frame(PresenceAbsence)
names(PresenceAbsence) = c("x","y","grad","elev","flora","N")

############
####PRESENCE-BACKGROUND
############

nobs = sum(dt$N) + length(dt$N)
PresenceBackground = array(0,dim=c(nobs,6))
counter = 1
for(i in 1:length(dt$N))
{
  
  PresenceBackground[counter,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],dt$flora[i],-1)
  counter = counter + 1
  
  if(dt$N[i] > 0)
  {
    for(j in counter:(counter + dt$N[i] -1))
    {
      PresenceBackground[j,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],dt$flora[i],1)
    }
    
    counter = counter + dt$N[i] 
  }
}

PresenceBackground = as.data.frame(PresenceBackground)
names(PresenceBackground) = c("x","y","grad","elev","flora","N")

#nbackground = sum(AllData$N)
#emptycells = which(AllData$N==0)
#deleteempties = sample(emptycells,length(emptycells)-nbackground)
#AllData = AllData[-deleteempties,] 
#AllData = as.data.frame(AllData)

####################################################################################
####################################################################################
####################################################################################
##### BIAS FUNCTION INDEPENDENT OF GRAD, ELEV #######
####################################################################################
####################################################################################
####################################################################################

# Now define observation probability to subsample only some of the gridpoints
subsample.p = exp(-1+dt$flora)
which.observe = rbinom(length(subsample.p),1,subsample.p)

############
####PRESENCE-ABSENCE
############

obs = sum(dt$N[which.observe==1]) + sum(dt$N[which.observe==1]==0) 
#nobs = sum(dt$N[which.observe==1]) + sum(dt$N[which.observe==1]==0) + length(dt$N)
PresenceAbsence.biased1= array(0,dim=c(nobs,6))
counter = 1
for(i in 1:length(dt$N))
{
 # PresenceAbsence.biased1[counter,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],-1)
#  counter = counter + 1
  
  if(dt$N[i] == 0 && which.observe[i]==1)
  {
    PresenceAbsence.biased1[counter,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],dt$flora[i],0)
    counter = counter + 1
  }
  if(dt$N[i] > 0  && which.observe[i]==1)
  {
    for(j in counter:(counter + dt$N[i] -1))
    {
      PresenceAbsence.biased1[j,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],dt$flora[i],1)
    }
    counter = counter + dt$N[i]
  }
}

PresenceAbsence.biased1 = as.data.frame(PresenceAbsence.biased1)
names(PresenceAbsence.biased1) = c("x","y","grad","elev","flora","N")

############
############
####PRESENCE-BACKGROUND
############

nobs = sum(dt$N[which.observe==1]) + length(dt$N)
PresenceBackground.biased1 = array(0,dim=c(nobs,6))
counter = 1
for(i in 1:length(dt$N))
{
  PresenceBackground.biased1[counter,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],dt$flora[i],-1)
  counter = counter + 1
  
  if(dt$N[i] > 0 && which.observe[i]==1)
  {
    for(j in counter:(counter + dt$N[i] -1))
    {
      PresenceBackground.biased1[j,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],dt$flora[i],1)
    }  
    counter = counter + dt$N[i] 
  }
}

PresenceBackground.biased1 = as.data.frame(PresenceBackground.biased1)
names(PresenceBackground.biased1) = c("x","y","grad","elev","flora","N")

############
############
####PRESENCE-BACKGROUND BOTH SAME BIAS
############

nobs = sum(dt$N[which.observe==1]) + sum(dt$N[which.observe==1]==0) 
#nobs = sum(dt$N[which.observe==1]) + sum(dt$N[which.observe==1]==0) + length(dt$N)
PresenceBackground.bothbiased1= array(0,dim=c(nobs,6))
counter = 1
for(i in 1:length(dt$N))
{
  if(which.observe[i]==1)
  {
    PresenceBackground.bothbiased1[counter,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],dt$flora[i],-1)
    counter = counter + 1
  }
  if(dt$N[i] > 0  && which.observe[i]==1)
  {
    for(j in counter:(counter + dt$N[i] -1))
    {
      PresenceBackground.bothbiased1[j,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],dt$flora[i],1)
    }
    counter = counter + dt$N[i]
  }
}

PresenceBackground.bothbiased1 = as.data.frame(PresenceBackground.bothbiased1)
names(PresenceBackground.bothbiased1) = c("x","y","grad","elev","flora","N")

####################################################################################
####################################################################################
####################################################################################
##### ANOTHER BIAS FUNCTION #######
####################################################################################
####################################################################################
####################################################################################


# Now define logistic observation function to subsample only some of the gridpoints
subsample.p = exp(-1+3*dt$grad)

which.observe = rbinom(length(subsample.p),1,subsample.p)

#########
####PRESENCE-ABSENCE
########

nobs = sum(dt$N[which.observe==1]) + sum(dt$N[which.observe==1]==0) 
#nobs = sum(dt$N[which.observe==1]) + sum(dt$N[which.observe==1]==0) + length(dt$N)
PresenceAbsence.biased2= array(0,dim=c(nobs,6))
counter = 1
for(i in 1:length(dt$N))
{
 # PresenceAbsence.biased2[counter,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],-1)
#  counter = counter + 1
  
  if(dt$N[i] == 0 && which.observe[i]==1)
  {
    PresenceAbsence.biased2[counter,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],dt$flora[i],0)
    counter = counter + 1
  }
  if(dt$N[i] > 0  && which.observe[i]==1)
  {
    for(j in counter:(counter + dt$N[i] -1))
    {
      PresenceAbsence.biased2[j,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],dt$flora[i],1)
    }
    counter = counter + dt$N[i]
  }
}

PresenceAbsence.biased2 = as.data.frame(PresenceAbsence.biased2)
names(PresenceAbsence.biased2) = c("x","y","grad","elev","flora","N")

############
####PRESENCE-BACKGROUND
###########

nobs = sum(dt$N[which.observe==1]) + length(dt$N)
PresenceBackground.biased2 = array(0,dim=c(nobs,6))
counter = 1
for(i in 1:length(dt$N))
{
  PresenceBackground.biased2[counter,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],dt$flora[i],-1)
  counter = counter + 1
  
  if(dt$N[i] > 0 && which.observe[i]==1)
  {
    for(j in counter:(counter + dt$N[i] -1))
    {
      PresenceBackground.biased2[j,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],dt$flora[i],1)
    }  
    counter = counter + dt$N[i] 
  }
}

PresenceBackground.biased2 = as.data.frame(PresenceBackground.biased2)
names(PresenceBackground.biased2) = c("x","y","grad","elev","flora","N")


##############
###PRESENCE-BACKGROUND WITH SAME BIAS
##############

obs = sum(dt$N[which.observe==1]) + sum(dt$N[which.observe==1]==0) 
#nobs = sum(dt$N[which.observe==1]) + sum(dt$N[which.observe==1]==0) + length(dt$N)
PresenceBackground.bothbiased2= array(0,dim=c(nobs,6))
counter = 1
for(i in 1:length(dt$N))
{
  if(which.observe[i]==1)
  {
    PresenceBackground.bothbiased2[counter,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],dt$flora[i],-1)
    counter = counter + 1
  }
  if(dt$N[i] > 0  && which.observe[i]==1)
  {
    for(j in counter:(counter + dt$N[i] -1))
    {
      PresenceBackground.bothbiased2[j,] = c(dt$xc[i],dt$yc[i],dt$grad[i],dt$elev[i],dt$flora[i],1)
    }
    counter = counter + dt$N[i]
  }
}

PresenceBackground.bothbiased2 = as.data.frame(PresenceBackground.bothbiased2)
names(PresenceBackground.bothbiased2) = c("x","y","grad","elev","flora","N")

####################################################################################
####################################################################################
####################################################################################


save(PresenceAbsence,PresenceBackground,PresenceAbsence.biased1,PresenceBackground.biased1,PresenceBackground.bothbiased1,PresenceAbsence.biased2,PresenceBackground.biased2,PresenceBackground.bothbiased2,file="AllData.RData")
