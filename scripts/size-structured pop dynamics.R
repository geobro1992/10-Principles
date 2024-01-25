##########################
# PLETHODON SALAMANDER PROJECTIONS

rm(list=ls(all=TRUE))
library(reshape)
#library(ggplot2)
library(plyr)
library(fields)

#Model parameters:
Npop=10  #initial pop size (not important)
Tmax=500 #max years
Tfishing=250 #when we start fishing
Trecovery=500 #stop fishing 

fishatmat=T # this switches which type of fishing (and graphing) we choose. If true, it graphs all life history graphs and age-selective fishing

#Life history parameters      	
recruitsize=1
Amax=25 #max age
amat=6 #age at maturation - this varies depending on the matration function chosen

Linf=5.5 #max size
k=.3   # Growth function param
fishingsize = .75*Linf #size-selective fishery threshold size 


M=0.18   #natural mortality coef  

c=7.49   #egg slope 
d=-23.6   #egg intercept

N0 = 100*Linf  # intitial number of eggs     
alpha = 3         # BH recruitment

q=1 #maturation ogive steepness (1 = very steep)

mu_f=0.2 #fishing mortality
harvest=c(0, mu_f)  

fishB=rep(0, length(harvest))
fishY=rep(0, length(harvest))
fishP=rep(0, length(harvest))
Recruits=rep(0, length(harvest))


#set up matrices
N=matrix(nrow=Amax, ncol=Tmax, data=0)

L=matrix(nrow=Amax, ncol=1, data=1)
mu=matrix(nrow=Amax, ncol=1, data=0)
E=matrix(nrow=Amax, ncol=Tmax, data=0)
S=matrix(nrow=1, ncol=Tmax, data=0)
P=matrix(nrow=1, ncol=Tmax)
B=matrix(nrow=1, ncol=Tmax)
Y=matrix(nrow=1, ncol=Tmax) 

#set up vectors to hold age-specific variables 
Pmale=rep(0, Amax) #probability of sex change
fert=rep(0, Tmax) #fert probability
pmat=rep(0, Amax) #maturation probability
eggs=rep(0, Amax) #fecundity
mu=rep(0, Amax) #natural mortality
select=rep(0, Amax) #fishery selectivity

N[1:10, 1]=Npop #initial population size
L[1]=recruitsize #larval size at recruitment



#First: Define the age-specific growth, 

for (a in 1:(Amax-1)) {
  
  #VBGE
  L[a+1]= Linf*(1-exp(-k)) + (L[a])*exp(-k) #length at age
  
}#end 1st age loop


# Second: define age-specific mortality and maturation, and fecundity

for (a in 1:(Amax)) {
  
  #natural mortality
  mu[a]= 1-(1/(1+Linf/L[a]*(exp(k)-1)))^(M/k)  #natural mortality function of Length
  
  #maturation ogive
  pmat[a]=ifelse( a < amat, 0, 1)     #knife edge maturation
  
  #FECUNDITY
  eggs[a]= c*L[a] + d #fecundity is a function of individual mass-at-age 
  
} #end 2nd age loop


#Third: define fishing selectivity    
for (a in 1:(Amax)) {
  #JUST AFTER AGE AT MATURATION
  if(fishatmat==1) {
    select[a] = ifelse( a < amat, 0, 1)  #knife edge gear selectivity right after maturation  
  } else {
    #JUST AFTER SIZE AT MATURATION
    #select[a]= ifelse( L[a] < (Lmat), 0, 1)  #gear selectivity same for all fish > maturation
    ##SIZE THRESHOLD
    select[a] = ifelse(L[a]< fishingsize, 0, 1)     
  }
} #end 3rd age loop


#Now for these age specific rates, simulate pop dynamics through time       
for(t in 1:(Tmax-1)) {
  
  
  E[,t]=N[, t]*pmat*eggs  #assuming spawning occurs between 1 t and the next
  
  P[t]= sum(E[,t])  #assuming fertilization is 100%
  
  beta= (alpha*P[t]-1)/(N0*P[t])  
  
  N[1,t+1]= alpha*P[t]/(1+beta*P[t]) #this is the N_0 class that is born and enters in the next time step.. 
  
  for (age in 1:(Amax-1)) {
    
    Fishing= if (Tfishing < t & Trecovery > t) {
      Fishing = select[age]*mu_f } else {
        Fishing = 0
      }
    
    N[age+1,t+1] = N[age, t]*exp(-mu[age]-Fishing)
    B[t]=sum(N[,t])
    Y[t]=sum((N[,t]) * (1-exp(-mu[age]-Fishing)) *(Fishing/(mu[age]+Fishing)))
    
    
    #assuming mortality happens after spawning
    
  } #end second age loop
  
} #end t loop

sums=colSums(N[amat:Amax,]) #calculate spawning popsize at each time


LEP_unfished=sum(E[,Tfishing-2])
LEP_fished=sum(E[,Tfishing+10])
SPR.sal=LEP_fished/LEP_unfished

pd.sal = round(100*(1 - (sums[Tfishing+100] / sums[Tfishing-2])))
sums.sal = sums


par(bg=NA)
plot(100:Tmax, sums[100:length(sums)], type="l", lwd = 10, col = "#d95f02", axes = F, xlab="", ylab="")

dev.copy(png,'salamander.png')
dev.off()


plot(1:Tmax, sums[1:length(sums)], type="l", lwd = 10, col = "#d95f02", xlab="Time", ylab="Abundance", ylim = c(0, 1000))

dev.copy(png,'salamander_supp.png')
dev.off()


########################
## juvenile mortality
#####################

##########################
# PLETHODON SALAMANDER PROJECTIONS

rm(list=ls(all=TRUE))
library(reshape)
#library(ggplot2)
library(plyr)
library(fields)

#Model parameters:
Npop=10  #initial pop size (not important)
Tmax=500 #max years
Tfishing=250 #when we start fishing
Trecovery=500 #stop fishing 

fishatmat=T # this switches which type of fishing (and graphing) we choose. If true, it graphs all life history graphs and age-selective fishing

#Life history parameters      	
recruitsize=1
Amax=25 #max age
amat=6 #age at maturation - this varies depending on the matration function chosen

Linf=5.5 #max size
k=.3   # Growth function param
fishingsize = .75*Linf #size-selective fishery threshold size 


M=0.18   #natural mortality coef  

c=7.49   #egg slope 
d=-23.6   #egg intercept

N0 = 100*Linf  # intitial number of eggs     
alpha = 3         # BH recruitment

q=1 #maturation ogive steepness (1 = very steep)

nest.success=0.2 #fishing mortality
harvest=c(0, mu_f)  

fishB=rep(0, length(harvest))
fishY=rep(0, length(harvest))
fishP=rep(0, length(harvest))
Recruits=rep(0, length(harvest))


#set up matrices
N=matrix(nrow=Amax, ncol=Tmax, data=0)

L=matrix(nrow=Amax, ncol=1, data=1)
mu=matrix(nrow=Amax, ncol=1, data=0)
E=matrix(nrow=Amax, ncol=Tmax, data=0)
S=matrix(nrow=1, ncol=Tmax, data=0)
P=matrix(nrow=1, ncol=Tmax)
B=matrix(nrow=1, ncol=Tmax)
Y=matrix(nrow=1, ncol=Tmax) 

#set up vectors to hold age-specific variables 
Pmale=rep(0, Amax) #probability of sex change
fert=rep(0, Tmax) #fert probability
pmat=rep(0, Amax) #maturation probability
eggs=rep(0, Amax) #fecundity
mu=rep(0, Amax) #natural mortality
select=rep(0, Amax) #fishery selectivity

N[1:10, 1]=Npop #initial population size
L[1]=recruitsize #larval size at recruitment



#First: Define the age-specific growth, 

for (a in 1:(Amax-1)) {
  
  #VBGE
  L[a+1]= Linf*(1-exp(-k)) + (L[a])*exp(-k) #length at age
  
}#end 1st age loop


# Second: define age-specific mortality and maturation, and fecundity

for (a in 1:(Amax)) {
  
  #natural mortality
  mu[a]= 1-(1/(1+Linf/L[a]*(exp(k)-1)))^(M/k)  #natural mortality function of Length
  
  #maturation ogive
  pmat[a]=ifelse( a < amat, 0, 1)     #knife edge maturation
  
  #FECUNDITY
  eggs[a]= c*L[a] + d #fecundity is a function of individual mass-at-age 
  
} #end 2nd age loop


#Third: define fishing selectivity    
for (a in 1:(Amax)) {
  #JUST AFTER AGE AT MATURATION
  if(fishatmat==1) {
    select[a] = ifelse( a >= amat, 0, 1)  #knife edge gear selectivity right after maturation  
  } else {
    #JUST AFTER SIZE AT MATURATION
    #select[a]= ifelse( L[a] < (Lmat), 0, 1)  #gear selectivity same for all fish > maturation
    ##SIZE THRESHOLD
    select[a] = ifelse(L[a] > fishingsize, 0, 1)     
  }
} #end 3rd age loop


#Now for these age specific rates, simulate pop dynamics through time       
for(t in 1:(Tmax-1)) {
  
  
  E[,t]=N[, t]*pmat*eggs  #assuming spawning occurs between 1 t and the next
  
  P[t] = if (Tfishing < t & Trecovery > t) { 
    P[t] = sum(E[,t])*nest.success} else {
      P[t] = sum(E[,t])
    }
  
  beta= (alpha*P[t]-1)/(N0*P[t])  
  
  N[1,t+1]= alpha*P[t]/(1+beta*P[t]) #this is the N_0 class that is born and enters in the next time step.. 
  
  for (age in 1:(Amax-1)) {
    
    N[age+1,t+1] = N[age, t]*exp(-mu[age])
    B[t]=sum(N[,t])
    
    
    #assuming mortality happens after spawning
    
  } #end second age loop
  
} #end t loop

sums=colSums(N[amat:Amax,]) #calculate spawning popsize at each time


LEP_unfished=sum(E[,Tfishing-2])
LEP_fished=sum(E[,Tfishing+10])
SPR.sal=LEP_fished/LEP_unfished

pd.sal = round(100*(1 - (sums[Tfishing+100] / sums[Tfishing-2])))
sums.sal = sums

par(bg=NA)
plot(100:Tmax, sums[100:length(sums)], type="l", lwd = 10, col = "#d95f02", axes = F, xlab="", ylab="")

dev.copy(png,'salamander2.png')
dev.off()