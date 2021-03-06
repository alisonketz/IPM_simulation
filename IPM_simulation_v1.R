
###
### Dam Walsh AAH IPM simulation code
###

setwd("~/Documents/IPM/IPM_simulation")

### preliminaries

library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(ggplot2)
library(coda)
library(knitr)
library(RcppGSL)
library(lubridate)

# setwd("/media/dwalsh/ADATA UFD/regularization_simulations")
# setwd("D:/regularization_simulations")

cores<-detectCores()

##############################################################################################################
###
## hazgen
#  Function to generate hunting (Oct-Dec) and non-hunting (Jan-Sept)
#  #hazards
#  Assume total hazard =hunting + nonhunting during Oct-Dec
#  Hunting hazard=ratioh2nh*nonhunting hazard
#  Arguments prob=annual survival prob, 
#  ratioh2nh = ratio of hunting/nonhunting hazards
###
##############################################################################################################


#hunting season begins 2017-09-16,2018-09-15
yday("2017-09-16")
yday("2017-09-29")#Date of hunting season begins

#
#

haz.gen<-function(prob,ratioh2nh){
  overall.haz <- -log(prob)
  haz.nonhunt <- ((1/(272+(365-272)*ratioh2nh))*overall.haz) #nonhunting hazard - daily
  haz.hunt <- haz.nonhunt*ratioh2nh #hunting hazard - daily
  return(list(hunt = haz.hunt, nonhunt = haz.nonhunt))
}

Sage <- c(0.4,0.5,0.55,0.65,0.7) #annual survival by age class(0,1,..,3,4+)

haz.young<-haz.gen(Sage[1],15) #annual survival young-of-the year

haz.adult1<-haz.gen(Sage[2],20) #annual survival probability -
## adult 1 yrs
haz.adult2<-haz.gen(Sage[3],15)

haz.adult3<-haz.gen(Sage[4],10)

haz.old <-haz.gen(Sage[5],2)  #annual survival probability
## adult >4yrs

hazmale <- rbind(c(haz.young$nonhunt,haz.adult1$nonhunt,haz.adult2$nonhunt,
                   haz.adult3$nonhunt, haz.old$nonhunt),
                 c(haz.young$hunt,haz.adult1$hunt,haz.adult2$hunt,
                   haz.adult3$hunt, haz.old$hunt))
hazfemale<-hazmale
Shazf<-hazfemale  #Survival for hunting and non-hunting periods

for(i in 1:ncol(hazfemale)){
  Shazf[1,i]<-exp(-hazfemale[1,i]*272)
  Shazf[2,i]<-exp(-hazfemale[2,i]*(365-272))
}

fecundity<-c(0,0,0.95,0.95,0.95) #birth prob for yr(0,1,2-3,>4)
nfawnprob<-c(0.05,0.5,0.4,0.05) #number fawns probabilities

Ninit.male<-c(20380,5661,2441,1538,868) #initial pop. sizes

Ninit.female<-c(18511,12714,8601,9723,5626) #initial pop. sizes
nyears<-5

#True parameter values
truew <-sum(mean(fecundity[3:length(fecundity)])*nfawnprob*0:3)

propMatrix<-diag(Sage)
propMatrix<-rbind(Sage*(fecundity)*truew*0.5,propMatrix)

projyr<-10 #how many years to project
trueN<-matrix(0, projyr, ncol(propMatrix))
trueN[1,]<-Ninit.female


for(i in 2:projyr){
  tempN<-propMatrix%*%trueN[i-1,]
  trueN[i,]<-tempN[-length(tempN)]
  trueN[i,ncol(trueN)]<-trueN[i,ncol(trueN)]+tempN[length(tempN)]
  
}
truetotalpop<-rowSums(trueN)

truer<-0.85

## datagen
sourceCpp("popsim.cpp")
datainput<-nxtyr(Ninit.male,Ninit.female,hazmale,hazfemale,fecundity,nfawnprob,0.5,truer,0.1,3,10,cores)

# nxtyr(Nmale, Nfemale, hazmale,hazfemale,rprob, fawnprob, sexpr, report,saprob,maxfawn, years, cores){
#   //Function to generate age-specific population and harvest sizes
#   //Arguments:
#   //Nmale - vector of current number of male individuals in each class of interest
#   //Nfemale - vector of current number of female individuals in each class of interest
#   //haz - matrix of hunting and non-hunting hazards for each class for each sex
#   //rprob - age-specific birth probability
#   //fawnprob - probability of the number of fawns for each doe [0:maxfawn]
#   //sexpr - prob fawn is male
#   //report - prob of reporting harvest
#   //saprob - prob harvested animal sexed and aged
#   //maxfawn - maximum number of fawns per doe
#   //years - number of years to run process
#   //cores - number of threads for parallel processing via openMP





## simsetup
#Function for generating beta distribution from quantiles
quantile2beta<-function(mu,c25,c975){
  alpha.temp<-seq(1,500,by=0.01)   #possible alpha values
  beta.temp<-(alpha.temp*(1-mu))/mu #possible beta values
  objective<-rep(0,length(alpha.temp)) #container for objective function values
  
  for(i in 1:length(alpha.temp)){
    qtemp<-qbeta(c(0.025,0.975),alpha.temp[i],beta.temp[i])
    objective[i]<-sum((qtemp[2]-c975)^2,(qtemp[1]-c25)^2)  #sums of squares
  }
  index<-which(objective==min(objective))  #minimize
  return(c(alpha.temp[index],beta.temp[index]))
}
