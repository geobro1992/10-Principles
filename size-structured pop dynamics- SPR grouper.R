rm(list=ls(all=TRUE))
library(reshape)
 #library(ggplot2)
 library(plyr)
library(fields)

#Model parameters:
Npop=100  #initial pop size (not important)
Tmax=500 #max years
Tfishing=200 #when we start fishing
Trecovery=500 #stop fishing 
  
  fishatmat=F # this switches which type of fishing (and graphing) we choose. If true, it graphs all life history graphs and age-selective fishing
  #Life history parameters      	
   recruitsize=4
    Amax=40 #max age

 	amat=5.5 #age at maturation - this varies depending on the matration function chosen
	Linf=99 #max size
	 fishingsize = .75*Linf #size-selective fishery threshold size 
        
       k=.2   # Growth function param
	  kappa=0.18   #natural mortality coef  
	 
	 eps=0.3 #natural mortality exp
	  
	  c= 2   #egg prod coef 
	  d=3 #egg &sperm prod exp
	  g = 7000000#sperm prod coef
	  v = 1 #mass-at-size coef
	  wexp =  3  #mass-at-size exp
	  h=1 #sperm fertilization effectiveness
 	    N0 = 100000*Linf       
  	    alpha = .8
  	
  	    q=.8 #maturation ogive steepness (1 = very steep)
  
   mu_f=0.2 #fishing mortality
  #harvest = seq(0, 0.8, by = 0.1)
  harvest=c(0, mu_f)  
   
  fishB=rep(0, length(harvest))
  fishY=rep(0, length(harvest))
  fishP=rep(0, length(harvest))
 Recruits=rep(0, length(harvest))
    
    
    #set up matrices
N=matrix(nrow=Amax, ncol=Tmax, data=0)
 
W=matrix(nrow=Amax, ncol=1, data=1)
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

N[1, 1]=Npop #initial population size
L[1]=recruitsize #larval size at recruitment
 
 
  
#First: Define the age-specific growth, mortality, and maturation function
  
  for (a in 1:(Amax-1)) {
  	
  	
#GROWTH 
		 
	L[a+1]= Linf*(1-exp(-k)) + (L[a])*exp(-k) #length at age
     W[a]= v*L[a]^wexp #weight at age
     
#natural mortality
   mu[a]= kappa + eps/W[a]  #natural mortality function of Length (Lorenzen curve)
     
 
 #SIZE-DEP MATURATION
    #  pmat[a]= 1/(1+exp(-q*(L[a]-Lmat)))   #probability an indvidual of age a is mature - logistic function of length  


###  AGE-DEP MATURATION:
       pmat[a]= 1/(1+exp(-q*(a-amat))) 
     #pmat[a]=ifelse( a < amat, 0, 1)     #knife edge maturation
     
     
 #FISHING SELECTIVITY    
     #JUST AFTER AGE AT MATURATION
     if(fishatmat==1) {
   select[a] = ifelse( a < amat+1, 0, 1)  #knife edge gear selectivity right after maturation  
        } else {
      #JUST AFTER SIZE AT MATURATION
       #select[a]= ifelse( L[a] < (Lmat), 0, 1)  #gear selectivity same for all fish > maturation
      ##SIZE THRESHOLD
      	select[a] = ifelse(L[a]< fishingsize, 0, 1)     
          }
 #FECUNDITY
       eggs[a]= c*L[a]^wexp #fecundity is a function of individual mass-at-age 
 
   } #end 1st age loop
# plot(mu~W, type="l")


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
 	B[t]=sum(N[,t]*W)
 	Y[t]=sum((N[,t]*W) * (1-exp(-mu[age]-Fishing)) *(Fishing/(mu[age]+Fishing)))
 	
 	
 	 #assuming mortality happens after spawning
		
    } #end second age loop
    
   
       } #end t loop
     sums=colSums(N[amat:Amax,]) #calculate spawning popsize at each time
 
    
  LEP_unfished=sum(E[,Tfishing-2])
  LEP_fished=sum(E[,Tfishing+10])
  SPR=LEP_fished/LEP_unfished
 
 if(fishatmat==1) {
 
 quartz()
 par(mfrow=c(3,2))

   plot(1:Amax, L, type="l", main="von Bert growth", ylab="Size", xlab="Age")
    legend("bottomright", legend=c( paste("Linf = ", Linf), paste("von Bert k = ", k)))

   plot(1:(Amax-1), eggs[-Amax]*pmat[-Amax], type="l", main="Age-specific Fecundity", ylab="Number of female progeny", xlab= "Age")
    legend("bottomright", legend=c(paste("50% mature at age", amat)))

   
  P_= 1:N0
  R=alpha*P_/(1+beta*P_)
   plot(P_, R, type="l", main="Beverton Holt Recruitment", xlab="Eggs(t)", ylab="Recruits(t+1)")
    legend("bottomright", legend=c(paste("alpha = ", round(alpha,6)), paste("beta = ", round(beta, 6))))

   
  
   plot(sums[1:(Tmax-1)], N[amat,2:Tmax], main="Stock Recruitment Curve", xlab="Spawners", ylab="1st time spawners", type="l", lwd=1.5)
  #abline(coef=c(0,1),  lwd=.5)
   
    legend("bottomright", legend=c(paste("M = ", kappa)))
  
   
   plot(1:Tmax, sums, type="l", main = "Fishing (F = 0.2) one year after maturity", xlab="Time", ylab="Abundance")
   legend("bottomright", legend=c( paste("Fishing mature fish only: SPR = ", round(SPR, 2))))
  } else {
  
  plot(1:Tmax, sums, type="l", main = "Fishing (F = 0.2) limit at 3/4 Linf ", xlab="Time", ylab="Abundance")
  #legend("bottomright", legend=c( paste("Fishing based on size: SPR = ", round(SPR, 2))))


   }
  