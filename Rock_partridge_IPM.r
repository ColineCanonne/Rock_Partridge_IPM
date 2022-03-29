######## ----------   ----------   ----------   JAGS code  ----------   ----------   ----------   ########
model{
##########################
# IPM built with 4 blocks: SURVIVAL, FECUNDITY, NETWORK TREND and  POPULATION MODEL 
# Y years: from 2000 to 2020
# R Spring counting sites  
# R2 Summer counting sites (using pointing dogs)
##########################

# **** SURVIVAL  ****

### Priors & constraints
ISJ ~ dunif(0,1) 				#initial state
capture_l ~ dnorm(0,0.01) 			#slope of capture effect 
recap ~ dunif(0,1) 				#recapture probability (to change emitters)
ms_l ~ dnorm(0,0.01)              		#mean monthly survival on logit scale
age_l ~ dnorm(0,0.01)           			#age effect on logit scale
sigma_s ~ dunif(0.1,3) 			#sd of random year effect 
tau_s <- pow(sigma_s, -2) 
# Monthly survival 
for (t in 1:Y){								
  eps.s[t] ~ dnorm(0,tau_s)                             		#random year effet
  logit(msad[t])  <- ms_l  + eps.s[t]                    		#adult monthly survival
  logit(msadC[t])  <- ms_l  + capture_l + eps.s[t]    		#first month after capture
  logit(msj[t]) <- ms_l  + age_l + eps.s[t]             		#juvenile monthly survival
  logit(msjC[t]) <-  ms_l + age_l + capture_l + eps.s[t] #first month after capture monthly survival
}#t
# Annual survival from May to May 
for (t in 1:Y){
  ysa[t] <- pow(msad[t],12)	
  ysj[t] <- pow(msad[t],5)*pow(msj[t],7)
}#t

### CMR Matrix (detection p=1)
# Initial states
I.STATE[1,1] <- 0
I.STATE[1,2] <- 0
I.STATE[1,3] <- ISJ
I.STATE[1,4] <- 1-ISJ
I.STATE[1,5] <- 0
I.STATE[1,6] <- 0 
# Transitions 
for (o in 1:(n.occasions-1)){
  PHI[1,o,1] <- msj[year[o]]*(1-age_tr[o])*(1-recap) # age_tr=1 in March, otherwise 0 : change age class
  PHI[1,o,2] <- msj[year[o]]*(age_tr[o])
  PHI[1,o,3] <- msj[year[o]]*recap
  PHI[1,o,4] <- msj[year[o]]*recap
  PHI[1,o,5] <- 1-msj[year[o]]
  PHI[1,o,6] <- 0
  PHI[2,o,1] <- 0
  PHI[2,o,2] <- msad[year[o]]*(1-recap)
  PHI[2,o,3] <- 0
  PHI[2,o,4] <- msad[year[o]]*recap
  PHI[2,o,5] <- 1-msad[year[o]]
  PHI[2,o,6] <- 0
  PHI[3,o,1] <- msjC[year[o]]
  PHI[3,o,2] <- 0
  PHI[3,o,3] <- 0
  PHI[3,o,4] <- 0
  PHI[3,o,5] <- 1-msjC[year[o]]
  PHI[3,o,6] <- 0
  PHI[4,o,1] <- 0
  PHI[4,o,2] <- msadC[year[o]]
  PHI[4,o,3] <- 0
  PHI[4,o,4] <- 0
  PHI[4,o,5] <- 1-msadC[year[o]]
  PHI[4,o,6] <- 0
  PHI[5,o,1] <- 0
  PHI[5,o,2] <- 0 
  PHI[5,o,3] <- 0
  PHI[5,o,4] <- 0
  PHI[5,o,5] <- 0
  PHI[5,o,6] <- 1
    PHI[6,o,1] <- 0
  PHI[6,o,2] <- 0
  PHI[6,o,3] <- 0
  PHI[6,o,4] <- 0
  PHI[6,o,5] <- 0
  PHI[6,o,6] <- 1
}#occ

### CMR Likelihood 
for (h in 1:nch){ 
 	z[h,f[h]] ~ dcat(I.STATE[1,])
}#h individuals
for (h in 1:nch){  
	for (o in (f[h]:(l[h]-1))){
		z[h,o+1] ~ dcat(PHI[z[h,o],o,])
	}#occ
}#h
	
# **** FECUNDITY  ****

### Priors & constraints
bp1 ~ dunif(0,1)     				#probability of initiating at least 1 clutch
bp2 ~ dunif(0,1)     				#probability of 1 vs 2 clutches | at least one clutch
ns ~ dunif(0,1)		          		#rearing success until August
bs_mean ~ dunif(2.5,10)	     		#mean clutch size
sigma_bs_year ~ dunif(0.1,3) 		#sd of clutch size random year effect
tau_bs_year <- pow(sigma_bs_year,-2) 
sigma_bs_site ~ dunif(0.1,3) 		 	#sd of clutch size random site effect
tau_bs_site <- pow(sigma_bs_site,-2) 
sigma_bs_res ~ dunif(0.1,3) 			#sd of clutch size residual variability
tau_bs_res <- pow(sigma_bs_res,-2) 

### Likelihood
for (i in 1:nb_bp1){
  clutch[i] ~ dbin(bp1,1)                  
}#i
for (i in 1:nb_bp2){
  double_clutch[i] ~ dbin(bp2,1)           
}#i
for (i in 1:nb_ns){
	success[i] ~ dbin(ns,1)                
}#i
for (g in 1:R2){
  site[g] ~ dnorm(0,tau_bs_site)
   for (t in fc_dog[g]:lc_dog[g]){
    mu[g,t] <- bs_mean  + site[g]  + eps.f[t]
	SR[g,t] ~ dnorm(mu[g,t],tau_bs_ts)            
   }#t
}#g
for (t in 1:Y){
    eps.f[t] ~ dnorm(0,tau_bs_year)                  
	bs[t] <- bs_mean + eps.f[t]                     
	N_c[t] <- bp1*(bp2*2*ns*bs[t] +(1-bp2)*1*ns*bs[t] )   # Nb chicks/female at the end of August
}#t

# **** NETWORK TREND ****
### Priors & constraints
for (i in 1:R){ 
	Nm[i,fc[i]] ~ dunif(0, 100)		#initial population size on spring counting sites
	  for (t in 1:(fc[i]-1)){
		log(Nm[i,t]) <- 0       		#set up NA not in likelihood 
    }#t
}#g
sigma_proc ~ dunif(0.1,5) 			#sd random site effect around network trend 
tau_proc <- pow(sigma_proc, -2) 

### Likelihood	
for (i in 1:R){
  for (t in fc[i]:(Y-1)) {
    log(Nm[i, t + 1]) <- log(Nm[i, t]) + r[t]
	}#t
  for (t in fc[i]:Y) {
    logy[i, t] ~ dnorm(log(Nm[i, t] + 1), tau_proc)   
  }#t
}#i

# **** POPULATION MODEL (from May to May)  ****
#Age class initial population sizes
Nf[1,1] <- 430
Nf[2,1] <- 570

for (t in 1:(Y-1)){
  Nf_aug[t] ~ dbin(pow(msad[t],3),Nf[1,t] + Nf[2,t]) # Nb of females end of August year t
  N_c_aug[t] ~ dpois(Nf_aug[t]*N_c[t])                    # Nb of chicks end of August year t
  Nf_c_aug[t] ~ dbin(0.5, N_c_aug[t] )                     # Nb of female chicks end of August year t
  Nf[1,t+1] ~ dbin(pow(msad[t],2)*pow(msj[t],7),Nf_c_aug[t]) # female chicks surviving until May t+1
  Nf[2,t+1] ~ dbin(ysad[t], Nf[1,t] + Nf[2,t])     # sub-adults and adults surviving until May t+1
  r[t] <- log(Nf_tot[t+1]/(Nf_tot[t]+1))
}#t

for (t in 1:Y){ 
    Nf_tot[t] <- Nf[1,t] + Nf[2,t]
}#t

}#model

inits <- function(){  # Inits from last model
  list(bp1 = 0.88
       ,bp2 = 0.55
       ,ns = 0.5
       ,ms_l = 2.62
       ,bs_mean = 4.5
       ,sigma_s = 0.25
       ,sigma_bs_site = 1.52
       ,sigma_bs_ts = 2.93
       ,sigma_proc = 0.32
  )} #*
