
### IPM simulation based on
### Dam Walsh's AAH IPM simulation code
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

cores=detectCores()

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
yday("2017-09-16") #259
yday("2018-09-15") #258
yday("2017-09-29")#Date of hunting season begins, #272

hunt.start=yday("2018-09-15")
#
#

haz.gen<-function(prob,ratioh2nh,hunt.start=yday("2018-09-15")){
  overall.haz = -log(prob)
  haz.nonhunt = ((1/(hunt.start+(365-hunt.start+7)*ratioh2nh))*overall.haz) #nonhunting hazard - daily
  haz.hunt = haz.nonhunt*ratioh2nh #hunting hazard - daily
  return(list(hunt = haz.hunt, nonhunt = haz.nonhunt))
}


Sage = c(0.5,seq(.7,.81,by=.02)) #annual survival by age class(0,1,2,...,5,6+)
num_ageclass=length(Sage) #must equal num_ageclass
# num_ageclass = 7

hazmale=matrix(NA,2,num_ageclass)#first row is nonhunt, secondrow is hunt hazard

for(i in 1:(num_ageclass-1)){
  hazmale[,i]=unlist(haz.gen(Sage[i],sample(10:20,1,replace=T)))[2:1]
}
hazmale[,num_ageclass]=unlist(haz.gen(Sage[i],2))[2:1]
hazfemale=hazmale

Shazf=hazfemale  #Survival for hunting and non-hunting periods

for(i in 1:num_ageclass){
  Shazf[1,i]=exp(-hazfemale[1,i]*hunt.start)
  Shazf[2,i]=exp(-hazfemale[2,i]*(365-hunt.start+7))
}

fecundity=c(0,0,rep(.95,num_ageclass-2)) #birth prob for yr(0,1,2-5,6+)

#number of fawns a doe can have is {0,1,2,3}, the probability of which of these outcomes is
nfawnprob=c(0.05,0.5,0.4,0.05) #number fawns probabilities


#adding 2 more age classes to initial population sizes...
# Ninit.male=c(20380,5661,2441,1538,868) #initial pop. sizes
# Ninit.female=c(18511,12714,8601,9723,5626) #initial pop. sizes
# for(i in 2:length(Ninit.male)){
#   cat((Ninit.male[i-1] - Ninit.male[i])/Ninit.male[i-1],"\n")
#   #cat(exp(-2*(Ninit.male[i-1] - Ninit.male[i])/Ninit.male[i-1]),"\n")
# }
# # plot(exp(-5*seq(0,1,.01)))
# 20280-20280*.7222
# 868-868*.4 # = 521
# 530-530*.45 # = 292
# 
# for(i in 2:length(Ninit.female)){
#   cat((Ninit.female[i-1] - Ninit.female[i])/Ninit.female[i-1],"\n")
#   #cat(exp(-2*(Ninit.male[i-1] - Ninit.male[i])/Ninit.male[i-1]),"\n")
# }
# 5626-5626*.35#[1] 3656.9
# 3660-3660*.38#[1] 2269.2


Ninit.male=c(20380,5661,2441,1538,868,529,282) #initial pop. sizes

Ninit.female=c(18511,12714,8601,9723,5626,3660,2269) #initial pop. sizes

nyears=2018-2002+1
nyears

plot(Ninit.male)
plot(Ninit.female)
#True parameter values
truew =sum(mean(fecundity[3:length(fecundity)])*nfawnprob*0:3)

propMatrix=diag(Sage)
propMatrix=rbind(Sage*(fecundity)*truew*0.5,propMatrix)

projyr=nyears #how many years to project
trueN=matrix(0, projyr, ncol(propMatrix))
trueN[1,]=Ninit.female


for(i in 2:projyr){
  tempN=propMatrix%*%trueN[i-1,]
  trueN[i,]=tempN[-length(tempN)]
  trueN[i,ncol(trueN)]=trueN[i,ncol(trueN)]+tempN[length(tempN)]
  
}
truetotalpop=rowSums(trueN)

truer=0.95

## datagen
sourceCpp("popsim.cpp")

sex.prob = .5
sex.aged.prob=.1
max.fawn=3

data.sim = nxtyr(Ninit.male,Ninit.female,hazmale,hazfemale,fecundity,nfawnprob,sex.prob,truer,sex.aged.prob,max.fawn,nyears,cores)

data.sim

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
quantile2beta=function(mu,c25,c975){
  alpha.temp=seq(1,500,by=0.01)   #possible alpha values
  beta.temp=(alpha.temp*(1-mu))/mu #possible beta values
  objective=rep(0,length(alpha.temp)) #container for objective function values
  
  for(i in 1:length(alpha.temp)){
    qtemp=qbeta(c(0.025,0.975),alpha.temp[i],beta.temp[i])
    objective[i]=sum((qtemp[2]-c975)^2,(qtemp[1]-c25)^2)  #sums of squares
  }
  index=which(objective==min(objective))  #minimize
  return(c(alpha.temp[index],beta.temp[index]))
}


# Npostmale = //post-hunt male population at end of biological year (5/14)
# Npostfemale= //post-hunt female population at end of biological year (5/14)
# harvestm=//harvested male population
# harvestf=//harvested female population
# oharvest= //reported harvest - first column=antlered/ second column = antlerless
# saharvm=//reported and sex/age -males
# saharvf=//reported and sex/age -females

data.sim$Npostmale
data.sim$Npostmale
data.sim$harvestm
data.sim$harvestf


###
### Fitting the AAH IPM model 
###



modelcode = nimbleCode({
  
  ###
  ### Prior models
  ###
  
  #Recruitment
  #unisex fawns per female
  rho~dbeta(10,10)
  
  #Reporting Rates
  #different for antlered and antlerless groups
  for (j in 1:2) {
    cll.report[j]~dnorm(p.cll.report[j],tau.cll.report)
  }
  
  #change in variable to probability scale
  for (j in 1:2) {
    for (a in 1:n.ageclass) {
      report[j,a]<-exp(-exp(cll.report[ant.idx[j,a]]))
    }
  }
  
  #Hunting Season Survival
  for (k in 1:2) {
    
    #shrinkage parameter for cloglog HS
    tau.cll.s.hunt[k]~dgamma(hs.alpha[k],hs.beta[k])
    
    #long-term mean cloglog HS
    mu.cll.s.hunt[k]~dnorm(p.cll.s.hunt[k], 2)
    
    for (t in 1:(T+1)) {
      cll.s.hunt[t,k]~dnorm(mu.cll.s.hunt[k],tau.cll.s.hunt[k])
    }
  }
  
  # for(i in 1:3){
  #   LHR[i]~dnorm(0,0.00001)
  # } 					#offset for differences among age classes
  # LHR[4]=0
  
  for(t in 1:(T+1)){
    for(j in 1:2) {
      for(a in 1:A){
        s.hunt[y,j,a]<-exp(-exp(cll.s.hunt[y,ant.idx[j,a]]))#+LHR[hr.idx[j,a]]
      }
    }
  } 	#change in variable to probability scale
  
  #Survival Outside the Hunting Season
  #shrinkage parameter for cloglog NS
  tau.cll.s.nat~dgamma(s.nat.alpha,s.nat.beta)
  
  #long-term mean cloglog NS
  #change in variable to probability scale
  mu.cll.s.nat~dnorm(pr.cll.s.nat, 2)
  for (t in 1:T) {
    cll.s.nat[t]~dnorm(mu.cll.s.nat,tau.cll.s.nat)
    s.nat[t]<-exp(-exp(cll.s.nat[t]))
  }
  
  ###
  ### Population Process Model  - Population Projection
  ###
  
  #N.p temporarily holds the projected age class
  for(t in 2:(T+1)){
    for(j in 1:2){
      for(a in 1:n.ageclass){
        N.project[t,j,a]<-N[t-1,j,a]*s.hunt[t-1,j,a]*s.nat[t-1]
      }
    }
    
    for(j in 1:2) {
      #max age class
      N.project.maxclass[t,j]<-N.project[t,j,(n.ageclass-1)]+N.project[t,j,n.ageclass]
      N[t,j,n.ageclass]<-N.p.Aclass[t,j]
      for(a in 2:(n.ageclass-1)){
        N[t,j,a]<-N.project[t,j,(a-1)]
      }
      N.p.1[t,j]<-(sum(N.p[t,1,1:n.ageclass]))*fec  						#females * unisex fawns per female
      N[t,j,1]<-N.p.1[t,j]
    }
  }
  
  
  ###
  ### Likelihood for the data
  ###
  
  for (t in 1:T){
    for (j in 1:2) {
      for (a in 1:n.ageclass.obs){
        hunt[t,j,a]<-N[t,j,a]*(1-s.hunt[t,j,a])*report[j,a]
      }
    }
    mu.O[t,1]<-sum(hunt[t,1,1:A])+hunt[t,2,1]
    mu.O[t,2]<-sum(hunt[t,2,2:A])
    
    #harvest data by antlered or antlerless group
    for (k in 1:2) {
      tau.O[t,k]<-1/mu.O[t,k]
      O[t,k]~dnorm(mu.O[t,k], tau.O[t,k])
    }
    
    p.noant[t,1]<-hunt[t,1,1]/mu.O[t,1]
    p.noant[t,2]<-hunt[t,1,2]/mu.O[t,1]
    p.noant[t,3]<-hunt[t,1,3]/mu.O[t,1]
    p.noant[t,4]<-hunt[t,1,4]/mu.O[t,1]
    p.noant[t,5]<-sum(hunt[t,1,5:6])/mu.O[t,1]
    p.noant[t,6]<-sum(hunt[t,1,7:9])/mu.O[t,1]
    p.noant[t,7]<-sum(hunt[t,1,10:12])/mu.O[t,1]
    p.noant[t,8]<-hunt[t,1,13]/mu.O[t,1]
    p.noant[t,9]<-1-sum(P1[t,1:8])
    
    p.ant[t,1]<-hunt[t,2,2]/mu.O[t,2]
    p.ant[t,2]<-hunt[t,2,3]/mu.O[t,2]
    p.ant[t,3]<-hunt[t,2,4]/mu.O[t,2]
    p.ant[t,4]<-sum(hunt[t,2,5:6])/mu.O[t,2]
    p.ant[t,5]<-sum(hunt[t,2,7:9])/mu.O[t,2]
    p.ant[t,6]<-sum(hunt[t,2,10:12])/mu.O[t,2]
    p.ant[t,7]<-hunt[t,2,13]/mu.O[t,2]
    
    
    #aged antlerless harvest data
    for(i in 1:8){
      aged.noant[t,i]<-C[t,1,i]
    }
    aged.noant[t,9]<-C[t,2,1]
    aged.noant[t,1:9]~dmulti(p.noant[t,1:9],C.total[t,1])
    
    #aged antlered harvest data
    for(i in 2:8){A2[t,(i-1)]<-C[t,2,i]}
    aged.ant[t,1:7]~dmulti(p.ant[t,1:7],C.total[t,2])
  }
  
  #Auxiliary Likelihood Construction
  for(t in 1:T){
    cll.a.s.nat[t]~dnorm(cll.s.nat[t],tau.cll.ant.s.nat[t])
  }
  
})



nimData <- list(

  )

nimConsts <-list(
  n.ageclass = num_ageclass
  )


#specify initial values
nimInits <- list(
  rho=beta(1,1,1),
  
)

modelout<-nimbleModel(code= modelcode,
                      name="aah_ipm_sim",
                      constants = nimConsts,
                      data = nimData,
                      inits = nimInits)
