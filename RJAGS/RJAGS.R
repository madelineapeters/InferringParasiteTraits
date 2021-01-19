##########################################
## Packages
##########################################
library(rjags)
library(dplyr)
library(tidyr)
library(dclone)
library(ggplot2)

##########################################
## Data
##########################################
# Read in original data
# Create summary dataframe with parasitaemia timeseries for acutely infected mice, separated by stage
# Create ring-stage dataframe, which has log-transformed ring-stage parasitaemia
# Create late-stage dataframe, which has log-transformed late-stage parasitaemia
# Store log-transformed parasitaemias into a data list
# Also specify step-size h, number of steps (d) and number of time points (e) in data list
# Add hyperparameters for priors

RFD=read.csv(paste(getwd(),"Real_Fucking_Data.csv",sep="/"))
names(RFD)[5:7] = c("ddr","dds","ddt") #donor cell-donor parasites rings and late stages
names(RFD)[9:11] = c("drr","drs","drt") #donor cell-recipient parasites rings and late stages
names(RFD)[13:15] = c("rdr","rds","rdt") #recipient cell-donor parasites rings and late stages
names(RFD)[17:19] = c("rrr","rrs","rrt") #recipient cell-recipient parasites rings and late stages
RFD = RFD[,2:21]

RFD=RFD %>% 
  mutate(., totD = DcDp+DcRp+DcNp) %>% 
  mutate(., totR = RcDp+RcRp+RcNp) %>% 
  mutate(.,Dpara = 100*DcDp/totD) %>% 
  mutate(.,Rpara = 100*RcDp/totR) %>% 
  mutate(., DparaR = 100*DcRp/totD) %>% 
  mutate(., RparaR = 100*RcRp/totR)

jags.dataA = RFD[51:100,] %>% 
  filter(.,Time<=24) %>% 
  #mutate(., Ring = (Dpara*ddr/(ddr+dds+ddt))) %>% 
  #mutate(., Late = (Dpara*(dds+ddt)/(ddr+dds+ddt))) %>% 
  mutate(., lnR = (Dpara*ddr/(ddr+dds+ddt))) %>% 
  mutate(., lnL = (Dpara*(dds+ddt)/(ddr+dds+ddt))) %>% 
  select(.,Mouse,Time,lnR,lnL)

jags.dataN = RFD[101:150,] %>% 
  filter(.,Time<=24) %>% 
  #mutate(., Ring = (Dpara*ddr/(ddr+dds+ddt))) %>% 
  #mutate(., Late = (Dpara*(dds+ddt)/(ddr+dds+ddt))) %>% 
  mutate(., lnR = (Dpara*ddr/(ddr+dds+ddt))) %>% 
  mutate(., lnL = (Dpara*(dds+ddt)/(ddr+dds+ddt))) %>% 
  select(.,Mouse,Time,lnR,lnL)

jags.spread.RA = jags.dataA %>% 
  select(.,Time,Mouse,lnR) %>% 
  spread(.,Mouse,lnR) %>%
  select(.,-Time)

jags.spread.LA = jags.dataA %>% 
  select(.,Time,Mouse,lnL) %>% 
  spread(.,Mouse,lnL) %>%
  select(.,-Time)

jags.spread.RN = jags.dataN %>% 
  select(.,Time,Mouse,lnR) %>% 
  spread(.,Mouse,lnR) %>%
  select(.,-Time)

jags.spread.LN = jags.dataN %>% 
  select(.,Time,Mouse,lnL) %>% 
  spread(.,Mouse,lnL) %>%
  select(.,-Time)

#For the hierarchical parasite maturation model
mixed.data.list = list(
  Ring = bind_cols(jags.spread.RA,jags.spread.RN),
  Late = bind_cols(jags.spread.LA,jags.spread.LN),
  time.int = c(2, 5, 8, 11, 15, 19, 22, 25),
  h = 1/24,
  d = 25,
  e = 8,
  nR = 5,
  nL = 10)
mixed.initials.norm = list(
  lambdaA = 0.65,
  kA = 0.3,
  kN = 0.3,
  c.R = 0.89,
  tau.c = 0.1,
  P0 = 5,
  mu = 0.32,
  sigma = 0.26,
  sd = 0.1)
red.mixed.initials.norm = list(
  kA = 0.3,
  kN = 0.3,
  c.mu = 0.89,
  tau.c = 0.1,
  P0 = 5,
  mu = 0.32,
  sigma = 0.26,
  sd = 0.1)

#For the pooled parasite maturation model
pooled.data.list = list(
  RingA = jags.spread.RA,
  LateA = jags.spread.LA,
  RingN = jags.spread.RN,
  LateN = jags.spread.LN,
  time.int = c(2, 5, 8, 11, 15, 19, 22, 25),
  h = 1/24,
  d = 25,
  e = 8,
  np = 5)
pooled.initials.norm = list(
  lambdaA = 0.65,
  kA = 0.44,
  kN = 0.3,
  c = 0.89,
  P0 = 5,
  mu = 0.32,
  sigma = 0.26,
  sd = 0.1)
red.pooled.initials.norm = list(
  kA = 0.44,
  kN = 0.3,
  c = 0.89,
  P0 = 5,
  mu = 0.32,
  sigma = 0.26,
  sd = 0.1)

##########################################
## Posteriors
##########################################
####Full model, hierarchical####
full.hier.jags = jags.model('Full.Mixed.StandClear.txt',
                            data = mixed.data.list,
                            inits = mixed.initials.norm,
                            n.chains = 3,
                            n.adapt = 20000)

update(full.hier.jags,n.iter=20000)

full.hier.jags.out = coda.samples(model=full.hier.jags,
                                  variable.names=c('lambdaA','kA','kN','c.mu','tau.c','P0','mu','sigma','sd'),
                                  n.iter=1000,thin=1)

plot(full.hier.jags.out)


full.hier.dev = dic.samples(full.hier.jags, 
                            variable.names=c("sd","lambdaA","kA","kN","mu","c.mu","tau.c","sigma","P0"), 
                            n.iter=50000)

####Full model, pooled####
full.pooled.jags = jags.model('Full.Pooled.StandClear.txt',
                              data = pooled.data.list,
                              inits = pooled.initials.norm,
                              n.chains = 3,
                              n.adapt = 20000)    

update(full.pooled.jags,20000)

full.pooled.jags.out = coda.samples(full.pooled.jags,
                                    c('sd','kA','kN','lambdaA','c','mu','sigma','P0'),
                                    1000,thin=1)  

plot(full.pooled.jags.out)

full.pooled.dev = dic.samples(full.pooled.jags,variable.names=c('sd','lambdaA','kA','kN','c','mu','sigma','P0'),n.iter=50000)

####Reduced model, hierarchical####
red.hier.jags = jags.model('Red.Mixed.StandClear.txt',
                           data = mixed.data.list,
                           inits = red.mixed.initials.norm,
                           n.chains = 3,
                           n.adapt = 20000)

update(red.hier.jags,n.iter=20000)

red.hier.jags.out = coda.samples(model=red.hier.jags,
                                 variable.names=c('lambdaA','kA','kN','c.mu','tau.c','P0','mu','sigma','sd'),
                                 n.iter=1000,thin=1)

plot(red.hier.jags.out)


red.hier.dev = dic.samples(red.hier.jags, 
                           variable.names=c("sd","kA","kN","mu","c.mu","tau.c","sigma","P0"), 
                           n.iter=50000)

####Reduced model, pooled####
red.pooled.jags = jags.model('Red.Pooled.StandClear.txt',
                             data = pooled.data.list,
                             inits = red.pooled.initials.norm,
                             n.chains = 3,
                             n.adapt = 20000)    

update(red.pooled.jags,20000)

red.pooled.jags.out = coda.samples(red.pooled.jags,
                                   c('sd','kA','kN','c','mu','sigma','P0'),
                                   1000)  

plot(red.pooled.jags.out)

red.pooled.dev = dic.samples(red.pooled.jags,variable.names=c('sd','kA','kN','c','mu','sigma','P0'),n.iter=50000)

##########################################
## Data cloning
##########################################
######Data cloning for pooled full model#########
full.pooled.model = custommodel("model {
                                 
                                 #Fixed effects
                                 kN ~ dunif(0,0.35)
                                 kA ~ dnorm(0.32,1/pow(0.1,2))
                                 lambdaN = 1
                                 lambdaA ~ dnorm(0.65,1/pow(0.07,2))
                                 c ~ dunif(0,1)
                                 P0 ~ dnorm(5.25,1/pow(0.2,2))
                                 mu ~ dnorm(0.32, 1/pow(0.05,2))
                                 sigma ~ dnorm(0.26,1/pow(0.1,2))
                                 sd ~ dunif(0.01,0.3)
                                 P0.adj = pnorm(1,mu,1/pow(sigma,2)) - pnorm(0,mu,1/pow(sigma,2))
                                 
                                 ################################################
                                 ## ODEs (Euler method)
                                 ################################################
                                 
                                 for (k in 1:K){
                                 
                                 #Initial conditions
                                 RA[1,k] <- (P0/P0.adj)*(pnorm(0.5,mu,1/pow(sigma,2))-pnorm(0,mu,1/pow(sigma,2))) 	# Ring-stage
                                 LA[1,k] <- (P0/P0.adj)*(pnorm(1,mu,1/pow(sigma,2))-pnorm(0.5,mu,1/pow(sigma,2)))	# Late-stage
                                 RN[1,k] <- (P0/P0.adj)*(pnorm(0.5,mu,1/pow(sigma,2))-pnorm(0,mu,1/pow(sigma,2))) 	# Ring-stage
                                 LN[1,k] <- (P0/P0.adj)*(pnorm(1,mu,1/pow(sigma,2))-pnorm(0.5,mu,1/pow(sigma,2)))	# Late-stage
                                 t[1,k] <- 0
                                 switch[1,k] <- 1
                                 
                                 for (g in 2:d){	# For each timestep; T = total number of time steps (25, 1+24)
                                 
                                 t[g,k] <- g/d
                                 switchA[g,k] <- ifelse(t[g,k]*lambdaA<0.5, 1, 0)
                                 switchN[g,k] <- ifelse(t[g,k]*lambdaN<0.5, 1, 0)
                                 
                                 # Changes from time t-1 to time t (timestep h = 1 hour)
                                 dRA[g-1,k] <-  kA*P0*dnorm(1-t[g,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c*t[g,k]) - switchA[g,k]*P0*dnorm(0.5-t[g,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c*t[g,k]) - c*RA[g-1,k] - (1-switchA[g,k])*kA*P0*dnorm(1.5-t[g,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c*t[g,k])
                                 dLA[g-1,k] <- switchA[g,k]*P0*dnorm(0.5-t[g,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c*t[g,k]) - P0*dnorm(1-t[g,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c*t[g,k]) - c*LA[g-1,k] + (1-switchA[g,k])*kA*P0*dnorm(1.5-t[g,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c*t[g,k])
                                 
                                 dRN[g-1,k] <-  kN*P0*dnorm(1-t[g,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c*t[g,k]) - switchN[g,k]*P0*dnorm(0.5-t[g,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c*t[g,k]) - c*RN[g-1,k] - (1-switchN[g,k])*kN*P0*dnorm(1.5-t[g,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c*t[g,k])
                                 dLN[g-1,k] <- switchN[g,k]*P0*dnorm(0.5-t[g,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c*t[g,k]) - P0*dnorm(1-t[g,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c*t[g,k]) - c*LN[g-1,k] + (1-switchN[g,k])*kN*P0*dnorm(1.5-t[g,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c*t[g,k])
                                 
                                 # Update variables (Euler, timestep h = 1)
                                 RA[g,k] <- max(0.001, RA[g-1,k] + dRA[g-1,k]*h)
                                 LA[g,k] <- max(0.001, LA[g-1,k] + dLA[g-1,k]*h)
                                 
                                 RN[g,k] <- max(0.001, RN[g-1,k] + dRN[g-1,k]*h)
                                 LN[g,k] <- max(0.001, LN[g-1,k] + dLN[g-1,k]*h)
                                 
                                 } #end time t
                                 
                                 
                                 ################################################
                                 ## Likelihood
                                 ################################################
                                 
                                 for(t in 1:e){ #t indexes the time intervals of interest
                                 for (z in 1:(np)){ # a indexes for individual mice
                                 
                                 LateA[t,z,k] ~ dnorm(((LA[time.int[t],k])),1/pow(sd,2))
                                 RingA[t,z,k] ~ dnorm(((RA[time.int[t],k])),1/pow(sd,2))
                                 
                                 LateN[t,z,k] ~ dnorm(((LN[time.int[t],k])),1/pow(sd,2))
                                 RingN[t,z,k] ~ dnorm(((RN[time.int[t],k])),1/pow(sd,2))
                                 
                                 }
                                 }
                                 
                                 }
                                 }")

set.seed(300)

e = 8
np = 5

pooled.data.list = list(
  RingA = dcdim(array(as.matrix(jags.spread.RA),c(e,np,1))),
  LateA = dcdim(array(as.matrix(jags.spread.LA),c(e,np,1))),
  RingN = dcdim(array(as.matrix(jags.spread.RN),c(e,np,1))),
  LateN = dcdim(array(as.matrix(jags.spread.LN),c(e,np,1))),
  time.int = c(2, 5, 8, 11, 15, 19, 22, 25),
  h = 1/24,
  d = 25,
  e = 8,
  np = 5,
  K = 1)
full.pooled.inits = list(
  lambdaA = 0.65,
  kA = 0.3,
  kN = 0.3,
  c = 0.89,
  P0 = 5,
  mu = 0.32,
  sigma = 0.26,
  sd = 0.1)
ifun = function(model, n.clones) {
  dclone(list(Y = dcdim(full.pooled.inits)),n.clones)
}
full.pooled.dcfit = dc.fit(data = pooled.data.list, params = c('lambdaA','kA','kN','c','mu','P0','sigma','sd'), 
                            model = full.pooled.model, inits = full.pooled.inits,
                            n.clones = c(1,2,4,6,8), unchanged = c("time.int","h","d","e","np"), multiply = "K",
                            n.iter = 10000)
summary(full.pooled.dcfit)
plot(full.pooled.dcfit)
dctable(full.pooled.dcfit)
plot(dctable(full.pooled.dcfit))
dcdiag(full.pooled.dcfit)
plot(dcdiag(full.pooled.dcfit))
pairs(full.pooled.dcfit)

DC1.pooled<-dctable(full.pooled.dcfit)
DC2.pooled<-dcdiag(full.pooled.dcfit)

plot(DC1.pooled,which=1:length(DC1.pooled),type=c('var'))
plot(DC1.pooled,which=1:length(DC1.pooled),type=c('all'))

######Data cloning for hierarchical full model#########
full.hier.model <- custommodel("model {
                               
                               for (k in 1:K){
                               for (z in 1:nL){
                               c[z,k] ~ dnorm(c.mu,1/pow(tau.c,2))
                               }
                               }
                               
                               
                               #Parameters for random effect
                               c.mu ~ dunif(0,2)
                               tau.c ~ dunif(0.01,0.2)
                               
                               #Fixed effects
                               kN ~ dunif(0,0.3)
                               kA ~ dunif(0.25,0.65)
                               lambdaN = 1
                               lambdaA ~ dnorm(0.65,1/pow(0.07,2))
                               mu ~ dunif(0,1)
                               P0 ~ dnorm(5.25,1/pow(0.2,2))
                               sigma ~ dnorm(0.26,1/pow(0.1,2))
                               sd ~ dunif(0.01,0.6)
                               
                               ################################################
                               ## ODEs (Euler method)
                               ################################################
                               
                               for (k in 1:K){
                               ########Acute mice
                               for (z in 1:nR){
                               P0.adj[z,k] = pnorm(1,mu,1/pow(sigma,2)) - pnorm(0,mu,1/pow(sigma,2))
                               #Initial conditions
                               R[1,z,k] <- (P0/P0.adj[z,k])*(pnorm(0.5,mu,1/pow(sigma,2))-pnorm(0,mu,1/pow(sigma,2))) 	# Ring-stage
                               L[1,z,k] <- (P0/P0.adj[z,k])*(pnorm(1,mu,1/pow(sigma,2))-pnorm(0.5,mu,1/pow(sigma,2)))	# Late-stage
                               t[1,z,k] <- 0
                               switchA[1,z,k] <- 1
                               
                               for (g in 2:d){	# For each timestep; T = total number of time steps (25, 1+24)
                               
                               t[g,z,k] <- g/d
                               switchA[g,z,k] <- ifelse(t[g,z,k]*lambdaA<0.5, 1, 0)
                               
                               # changes from time t-1 to time t (timestep h = 1 hour)
                               dR[g-1,z,k] <-  kA*P0*dnorm(1-t[g,z,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c[z,k]*t[g,z,k]) - switchA[g,z,k]*P0*dnorm(0.5-t[g,z,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c[z,k]*t[g,z,k]) - c[z,k]*R[g-1,z,k] - (1-switchA[g,z,k])*kA*P0*dnorm(1.5-t[g,z,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c[z,k]*t[g,z,k])
                               dL[g-1,z,k] <- switchA[g,z,k]*P0*dnorm(0.5-t[g,z,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c[z,k]*t[g,z,k]) - P0*dnorm(1-t[g,z,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c[z,k]*t[g,z,k]) - c[z,k]*L[g-1,z,k] + (1-switchA[g,z,k])*kA*P0*dnorm(1.5-t[g,z,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c[z,k]*t[g,z,k])
                               
                               # Update variables (Euler, timestep h = 1)
                               R[g,z,k] <- max(0, R[g-1,z,k] + dR[g-1,z,k]*h)
                               L[g,z,k] <- max(0, L[g-1,z,k] + dL[g-1,z,k]*h)
                               
                               } #end time t
                               }
                               
                               for (z in (nR+1):nL){
                               P0.adj[z,k] = pnorm(1,mu,1/pow(sigma,2)) - pnorm(0,mu,1/pow(sigma,2))
                               #Initial conditions
                               R[1,z,k] <- (P0/P0.adj[z,k])*(pnorm(0.5,mu,1/pow(sigma,2))-pnorm(0,mu,1/pow(sigma,2))) 	# Ring-stage
                               L[1,z,k] <- (P0/P0.adj[z,k])*(pnorm(1,mu,1/pow(sigma,2))-pnorm(0.5,mu,1/pow(sigma,2)))	# Late-stage
                               t[1,z,k] <- 0
                               switchN[1,z,k] <- 1
                               
                               for (g in 2:d){	# For eac[z,k]h timestep; T = total number of time steps (25, 1+24)
                               
                               t[g,z,k] <- g/d
                               switchN[g,z,k] <- ifelse(t[g,z,k]*lambdaN<0.5, 1, 0)
                               
                               # c[z,k]hanges from time t-1 to time t (timestep h = 1 hour)
                               dR[g-1,z,k] <-  kN*P0*dnorm(1-t[g,z,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c[z,k]*t[g,z,k]) - switchN[g,z,k]*P0*dnorm(0.5-t[g,z,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c[z,k]*t[g,z,k]) - c[z,k]*R[g-1,z,k] - (1-switchN[g,z,k])*kN*P0*dnorm(1.5-t[g,z,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c[z,k]*t[g,z,k])
                               dL[g-1,z,k] <- switchN[g,z,k]*P0*dnorm(0.5-t[g,z,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c[z,k]*t[g,z,k]) - P0*dnorm(1-t[g,z,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c[z,k]*t[g,z,k]) - c[z,k]*L[g-1,z,k] + (1-switchN[g,z,k])*kN*P0*dnorm(1.5-t[g,z,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c[z,k]*t[g,z,k])
                               
                               # Update variables (Euler, timestep h = 1)
                               R[g,z,k] <- max(0.001, R[g-1,z,k] + dR[g-1,z,k]*h)
                               L[g,z,k] <- max(0.001, L[g-1,z,k] + dL[g-1,z,k]*h)
                               
                               } #end time t
                               }
                               ################################################
                               ## Likelihood
                               ################################################
                               
                               for(t in 1:e){ #t indexes the time intervalAs of interest
                               for (z in 1:nR){ #z indexes for individual mice
                               
                               Late[t,z,k] ~ dnorm(((L[time.int[t],z,k])),1/pow(sd,2))
                               Ring[t,z,k] ~ dnorm(((R[time.int[t],z,k])),1/pow(sd,2))
                               
                               }
                               for (z in (nR+1):nL){ #z indexes for individual mice
                               
                               Late[t,z,k] ~ dnorm(((L[time.int[t],z,k])),1/pow(sd,2))
                               Ring[t,z,k] ~ dnorm(((R[time.int[t],z,k])),1/pow(sd,2))
                               
                               }
                               }
                               
                               } #end K loop
                               }")

set.seed(300)

e = 8
nL = 10

full.data.list <- list(
  Ring = dcdim(array(as.matrix(bind_cols(jags.spread.RA,jags.spread.RN)),c(e,nL,1))),
  Late = dcdim(array(as.matrix(bind_cols(jags.spread.LA,jags.spread.LN)),c(e,nL,1))),
  time.int = c(2, 5, 8, 11, 15, 19, 22, 25),
  h = 1/24,
  d = 25,
  e = 8,
  nR = 5,
  nL = 10,
  K = 1)
full.hier.inits <- list(
  lambdaA = 0.65,
  kA = 0.3,
  kN = 0.3,
  c.mu = 0.89,
  tau.c = 0.1,
  P0 = 5,
  mu = 0.32,
  sigma = 0.26,
  sd = 0.1)
ifun <- function(model, n.clones) {
  dclone(list(Y = dcdim(full.hier.inits)),n.clones)
}
full.hier.dcfit <- dc.fit(data = full.data.list, params = c('lambdaA','kA','kN','c.mu','tau.c','mu','P0','sigma','sd'), 
                          model = full.hier.model, inits = full.hier.inits,
                          n.clones = c(1,2,4,6,8), unchanged = c("time.int","h","d","e","nR","nL"), multiply = "K",
                          n.iter = 10000)
summary(full.hier.dcfit)
plot(full.hier.dcfit)
dctable(full.hier.dcfit)
plot(dctable(full.hier.dcfit))
dcdiag(full.hier.dcfit)
plot(dcdiag(full.hier.dcfit))
pairs(full.hier.dcfit)

DC1.full<-dctable(full.hier.dcfit)
DC2.full<-dcdiag(full.hier.dcfit)

plot(DC1.full,which=1:length(DC1.full),type=c('var'))
plot(DC1.full,which=1:length(DC1.full),type=c('all'))

######Data cloning for pooled reduced model#########
red.pooled.model <- custommodel("model {
                                
                                #Parameters for random effect
                                
                                #Fixed effects
                                kN ~ dunif(0,0.35)
                                kA ~ dnorm(0.32,1/pow(0.1,2))
                                lambdaN = 1
                                lambdaA = 1
                                c ~ dunif(0.5,1)
                                P0 ~ dnorm(5.25,1/pow(0.2,2))
                                mu ~ dnorm(0.32, 1/pow(0.05,2))
                                sigma ~ dnorm(0.26,1/pow(0.1,2))
                                sd ~ dunif(0.01,0.3)
                                P0.adj = pnorm(1,mu,1/pow(sigma,2)) - pnorm(0,mu,1/pow(sigma,2))
                                
                                ################################################
                                ## ODEs (Euler method)
                                ################################################
                                
                                for (k in 1:K){
                                
                                #Initial conditions
                                RA[1,k] <- (P0/P0.adj)*(pnorm(0.5,mu,1/pow(sigma,2))-pnorm(0,mu,1/pow(sigma,2))) 	# Ring-stage
                                LA[1,k] <- (P0/P0.adj)*(pnorm(1,mu,1/pow(sigma,2))-pnorm(0.5,mu,1/pow(sigma,2)))	# Late-stage
                                RN[1,k] <- (P0/P0.adj)*(pnorm(0.5,mu,1/pow(sigma,2))-pnorm(0,mu,1/pow(sigma,2))) 	# Ring-stage
                                LN[1,k] <- (P0/P0.adj)*(pnorm(1,mu,1/pow(sigma,2))-pnorm(0.5,mu,1/pow(sigma,2)))	# Late-stage
                                t[1,k] <- 0
                                switch[1,k] <- 1
                                
                                for (g in 2:d){	# For each timestep; T = total number of time steps (25, 1+24)
                                
                                t[g,k] <- g/d
                                switchA[g,k] <- ifelse(t[g,k]*lambdaA<0.5, 1, 0)
                                switchN[g,k] <- ifelse(t[g,k]*lambdaN<0.5, 1, 0)
                                
                                # Changes from time t-1 to time t (timestep h = 1 hour)
                                dRA[g-1,k] <-  kA*P0*dnorm(1-t[g,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c*t[g,k]) - switchA[g,k]*P0*dnorm(0.5-t[g,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c*t[g,k]) - c*RA[g-1,k] - (1-switchA[g,k])*kA*P0*dnorm(1.5-t[g,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c*t[g,k])
                                dLA[g-1,k] <- switchA[g,k]*P0*dnorm(0.5-t[g,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c*t[g,k]) - P0*dnorm(1-t[g,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c*t[g,k]) - c*LA[g-1,k] + (1-switchA[g,k])*kA*P0*dnorm(1.5-t[g,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c*t[g,k])
                                
                                dRN[g-1,k] <-  kN*P0*dnorm(1-t[g,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c*t[g,k]) - switchN[g,k]*P0*dnorm(0.5-t[g,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c*t[g,k]) - c*RN[g-1,k] - (1-switchN[g,k])*kN*P0*dnorm(1.5-t[g,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c*t[g,k])
                                dLN[g-1,k] <- switchN[g,k]*P0*dnorm(0.5-t[g,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c*t[g,k]) - P0*dnorm(1-t[g,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c*t[g,k]) - c*LN[g-1,k] + (1-switchN[g,k])*kN*P0*dnorm(1.5-t[g,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c*t[g,k])
                                
                                # Update variables (Euler, timestep h = 1)
                                RA[g,k] <- max(0.001, RA[g-1,k] + dRA[g-1,k]*h)
                                LA[g,k] <- max(0.001, LA[g-1,k] + dLA[g-1,k]*h)
                                
                                RN[g,k] <- max(0.001, RN[g-1,k] + dRN[g-1,k]*h)
                                LN[g,k] <- max(0.001, LN[g-1,k] + dLN[g-1,k]*h)
                                
                                } #end time t
                                
                                
                                ################################################
                                ## Likelihood
                                ################################################
                                
                                for(t in 1:e){ #t indexes the time intervals of interest
                                for (z in 1:(np)){ # a indexes for individual mice
                                
                                LateA[t,z,k] ~ dnorm(((LA[time.int[t],k])),1/pow(sd,2))
                                RingA[t,z,k] ~ dnorm(((RA[time.int[t],k])),1/pow(sd,2))
                                
                                LateN[t,z,k] ~ dnorm(((LN[time.int[t],k])),1/pow(sd,2))
                                RingN[t,z,k] ~ dnorm(((RN[time.int[t],k])),1/pow(sd,2))
                                
                                }
                                }
                                
                                }
                                }")

set.seed(300)

e = 8
np = 5

pooled.data.list <- list(
  RingA = dcdim(array(as.matrix(jags.spread.RA),c(e,np,1))),
  LateA = dcdim(array(as.matrix(jags.spread.LA),c(e,np,1))),
  RingN = dcdim(array(as.matrix(jags.spread.RN),c(e,np,1))),
  LateN = dcdim(array(as.matrix(jags.spread.LN),c(e,np,1))),
  time.int = c(2, 5, 8, 11, 15, 19, 22, 25),
  h = 1/24,
  d = 25,
  e = 8,
  np = 5,
  K = 1)
red.pooled.inits <- list(
  kA = 0.3,
  kN = 0.3,
  c = 0.89,
  P0 = 5,
  mu = 0.32,
  sigma = 0.26,
  sd = 0.1)
ifun <- function(model, n.clones) {
  dclone(list(Y = dcdim(red.pooled.inits)),n.clones)
}
red.pooled.dcfit <- dc.fit(data = pooled.data.list, params = c('kA','kN','c','mu','P0','sigma','sd'), 
                           model = red.pooled.model, inits = red.pooled.inits,
                           n.clones = c(1,2,4,6,8), unchanged = c("time.int","h","d","e","np"), multiply = "K",
                           n.iter = 10000)
summary(red.pooled.dcfit)
plot(red.pooled.dcfit)
dctable(red.pooled.dcfit)
plot(dctable(red.pooled.dcfit))
dcdiag(red.pooled.dcfit)
plot(dcdiag(red.pooled.dcfit))
pairs(red.pooled.dcfit)

DC1.pooled<-dctable(red.pooled.dcfit)
DC2.pooled<-dcdiag(red.pooled.dcfit)

plot(DC1.pooled,which=1:length(DC1.pooled),type=c('var'))
plot(DC1.pooled,which=1:length(DC1.pooled),type=c('all'))

######Data cloning for hierarchical reduced model#########

red.hier.model <- custommodel("model {
                              
                              for (k in 1:K){
                              for (z in 1:nL){
                              c[z,k] ~ dnorm(c.mu,1/pow(tau.c,2))
                              }
                              }
                              
                              
                              #Parameters for random effect
                              c.mu ~ dunif(0,1)
                              tau.c ~ dunif(0.01,0.2)
                              
                              #Fixed effects
                              kN ~ dunif(0,0.3)
                              kA ~ dunif(0.25,0.65)
                              lambdaN = 1
                              lambdaA = 1
                              mu ~ dunif(0,1)
                              P0 ~ dnorm(5.25,1/pow(0.2,2))
                              sigma ~ dnorm(0.26,1/pow(0.1,2))
                              sd ~ dunif(0.01,0.6)
                              
                              ################################################
                              ## ODEs (Euler method)
                              ################################################
                              
                              for (k in 1:K){
                              ########Acute mice
                              for (z in 1:nR){
                              P0.adj[z,k] = pnorm(1,mu,1/pow(sigma,2)) - pnorm(0,mu,1/pow(sigma,2))
                              #Initial conditions
                              R[1,z,k] <- (P0/P0.adj[z,k])*(pnorm(0.5,mu,1/pow(sigma,2))-pnorm(0,mu,1/pow(sigma,2))) 	# Ring-stage
                              L[1,z,k] <- (P0/P0.adj[z,k])*(pnorm(1,mu,1/pow(sigma,2))-pnorm(0.5,mu,1/pow(sigma,2)))	# Late-stage
                              t[1,z,k] <- 0
                              switchA[1,z,k] <- 1
                              
                              for (g in 2:d){	# For each timestep; T = total number of time steps (25, 1+24)
                              
                              t[g,z,k] <- g/d
                              switchA[g,z,k] <- ifelse(t[g,z,k]*lambdaA<0.5, 1, 0)
                              
                              # Changes from time t-1 to time t (timestep h = 1 hour)
                              dR[g-1,z,k] <-  kA*P0*dnorm(1-t[g,z,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c[z,k]*t[g,z,k]) - switchA[g,z,k]*P0*dnorm(0.5-t[g,z,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c[z,k]*t[g,z,k]) - c[z,k]*R[g-1,z,k] - (1-switchA[g,z,k])*kA*P0*dnorm(1.5-t[g,z,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c[z,k]*t[g,z,k])
                              dL[g-1,z,k] <- switchA[g,z,k]*P0*dnorm(0.5-t[g,z,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c[z,k]*t[g,z,k]) - P0*dnorm(1-t[g,z,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c[z,k]*t[g,z,k]) - c[z,k]*L[g-1,z,k] + (1-switchA[g,z,k])*kA*P0*dnorm(1.5-t[g,z,k]*lambdaA, mu, 1/pow(sigma,2))*lambdaA*exp(-c[z,k]*t[g,z,k])
                              
                              # Update variables (Euler, timestep h = 1)
                              R[g,z,k] <- max(0, R[g-1,z,k] + dR[g-1,z,k]*h)
                              L[g,z,k] <- max(0, L[g-1,z,k] + dL[g-1,z,k]*h)
                              
                              } #end time t
                              }
                              
                              for (z in (nR+1):nL){
                              P0.adj[z,k] = pnorm(1,mu,1/pow(sigma,2)) - pnorm(0,mu,1/pow(sigma,2))
                              #Initial conditions
                              R[1,z,k] <- (P0/P0.adj[z,k])*(pnorm(0.5,mu,1/pow(sigma,2))-pnorm(0,mu,1/pow(sigma,2))) 	# Ring-stage
                              L[1,z,k] <- (P0/P0.adj[z,k])*(pnorm(1,mu,1/pow(sigma,2))-pnorm(0.5,mu,1/pow(sigma,2)))	# Late-stage
                              t[1,z,k] <- 0
                              switchN[1,z,k] <- 1
                              
                              for (g in 2:d){	# For each timestep; T = total number of time steps (25, 1+24)
                              
                              t[g,z,k] <- g/d
                              switchN[g,z,k] <- ifelse(t[g,z,k]*lambdaN<0.5, 1, 0)
                              
                              # Changes from time t-1 to time t (timestep h = 1 hour)
                              dR[g-1,z,k] <-  kN*P0*dnorm(1-t[g,z,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c[z,k]*t[g,z,k]) - switchN[g,z,k]*P0*dnorm(0.5-t[g,z,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c[z,k]*t[g,z,k]) - c[z,k]*R[g-1,z,k] - (1-switchN[g,z,k])*kN*P0*dnorm(1.5-t[g,z,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c[z,k]*t[g,z,k])
                              dL[g-1,z,k] <- switchN[g,z,k]*P0*dnorm(0.5-t[g,z,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c[z,k]*t[g,z,k]) - P0*dnorm(1-t[g,z,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c[z,k]*t[g,z,k]) - c[z,k]*L[g-1,z,k] + (1-switchN[g,z,k])*kN*P0*dnorm(1.5-t[g,z,k]*lambdaN, mu, 1/pow(sigma,2))*lambdaN*exp(-c[z,k]*t[g,z,k])
                              
                              # Update variables (Euler, timestep h = 1)
                              R[g,z,k] <- max(0.001, R[g-1,z,k] + dR[g-1,z,k]*h)
                              L[g,z,k] <- max(0.001, L[g-1,z,k] + dL[g-1,z,k]*h)
                              
                              } #end time t
                              }
                              ################################################
                              ## Likelihood
                              ################################################
                              
                              for(t in 1:e){ #t indexes the time intervalAs of interest
                              for (z in 1:nR){ #z indexes for individual mice
                              
                              Late[t,z,k] ~ dnorm(((L[time.int[t],z,k])),1/pow(sd,2))
                              Ring[t,z,k] ~ dnorm(((R[time.int[t],z,k])),1/pow(sd,2))
                              
                              }
                              for (z in (nR+1):nL){ #z indexes for individual mice
                              
                              Late[t,z,k] ~ dnorm(((L[time.int[t],z,k])),1/pow(sd,2))
                              Ring[t,z,k] ~ dnorm(((R[time.int[t],z,k])),1/pow(sd,2))
                              
                              }
                              }
                              
                              } #end K loop
                              }")

set.seed(150)

e = 8
nL = 10

hier.data.list <- list(
  Ring = dcdim(array(as.matrix(bind_cols(jags.spread.RA,jags.spread.RN)),c(e,nL,1))),
  Late = dcdim(array(as.matrix(bind_cols(jags.spread.LA,jags.spread.LN)),c(e,nL,1))),
  time.int = c(2, 5, 8, 11, 15, 19, 22, 25),
  h = 1/24,
  d = 25,
  e = 8,
  nR = 5,
  nL = 10,
  K = 1)
red.hier.inits <- list(
  kA = 0.3,
  kN = 0.3,
  c.mu = 0.89,
  tau.c = 0.1,
  P0 = 5,
  mu = 0.32,
  sigma = 0.26,
  sd = 0.1)
ifun <- function(model, n.clones) {
  dclone(list(Y = dcdim(red.hier.inits)),n.clones)
}
red.hier.dcfit <- dc.fit(data = hier.data.list, params = c('kA','kN','c.mu','tau.c','mu','P0','sigma','sd'), 
                         model = red.hier.model, inits = red.hier.inits,
                         n.clones = c(1,2,4,6,8), unchanged = c("time.int","h","d","e","nR","nL"), multiply = "K",
                         n.iter = 10000)
summary(red.hier.dcfit)
plot(red.hier.dcfit)
dctable(red.hier.dcfit)
plot(dctable(red.hier.dcfit))
dcdiag(red.hier.dcfit)
plot(dcdiag(red.hier.dcfit))
pairs(red.hier.dcfit)

DC1.full<-dctable(red.hier.dcfit)
DC2.full<-dcdiag(red.hier.dcfit)

plot(DC1.full,which=1:length(DC1.full),type=c('var'))
plot(DC1.full,which=1:length(DC1.full),type=c('all'))


##########################################
## Plotting
##########################################
#Basic dataframe of best fit for full model
para.full = summary(full.pooled.dcfit)$statistics[,1]
#para.full = summary(dcfit.full)$statistics[,1]
lambdaN = 1
lambdaA = para.full[5]
kA = para.full[1]
kN = para.full[2]
c = para.full[4]
mu = para.full[6]
sigma = para.full[8]
P0 = para.full[3]
full.sd =  para.full[7]

full.para.est = c(lambdaA,lambdaN,kA,kN,c,mu,sigma,P0)
full.paraA = full.para.est[c(1,3,5:8)]
full.paraN = full.para.est[c(2,4:8)]

full.plot.A = PM.to.DDE.full(full.paraA,timedata,ADcDp,'plot',0)
full.plot.N = PM.to.DDE.full(full.paraN,timedata,NDcDp,'plot',0)

#Basic dataframe of best fit for reduced model
para.red = summary(red.pooled.dcfit)$statistics[,1]
#para.red = summary(dcfit.pooled)$statistics[,1]
lambdaN = 1
lambdaA = 1
kA = para.red[1]
kN = para.red[2]
c = para.red[4]
mu = para.red[5]
sigma = para.red[7]
P0 = para.red[3]
red.sd =  para.red[6]

red.para.est = c(lambdaA,lambdaN,kA,kN,c,mu,sigma,P0)
red.paraA = red.para.est[c(1,3,5:8)]
red.paraN = red.para.est[c(2,4:8)]

red.plot.A = PM.to.DDE.full(red.paraA,timedata,ADcDp,'plot',0)
red.plot.N = PM.to.DDE.full(red.paraN,timedata,NDcDp,'plot',0)
##################
ADcDp.adj = ADcDp %>% 
  filter(.,Time<=24) %>% 
  mutate(., Late = (Spec.Para*(Troph+Schiz)/(Ring+Troph+Schiz))) %>% 
  mutate(., Ring = (Spec.Para*Ring/(Ring+Troph+Schiz))) %>%
  select(.,Mouse,Time,Ring,Late)
NDcDp.adj = NDcDp %>% 
  filter(.,Time<=24) %>% 
  mutate(., Late = (Spec.Para*(Troph+Schiz)/(Ring+Troph+Schiz))) %>% 
  mutate(., Ring = (Spec.Para*Ring/(Ring+Troph+Schiz))) %>% 
  select(.,Mouse,Time,Ring,Late)
#Make dataframe for plotting empirical data
ADcDp.adj$Host = "Acutely infected"; NDcDp.adj$Host = "Naive"
names(ADcDp.adj)[3:4] = c('Ring-stage','Late-stage')
names(NDcDp.adj)[3:4] = c('Ring-stage','Late-stage')
DcDp.adj = bind_rows(ADcDp.adj,NDcDp.adj) %>% gather(.,"Stage","Parasitaemia",3:4)

#Make dataframe for plotting best fits
red.plot.A$Host = "Acutely infected"; red.plot.N$Host = "Naive"
full.plot.A$Host = "Acutely infected"; full.plot.N$Host = "Naive"
names(red.plot.A)[2:3]  = c('Ring-stage','Late-stage')
names(red.plot.N)[2:3]  = c('Ring-stage','Late-stage')
names(full.plot.A)[2:3]  = c('Ring-stage','Late-stage')
names(full.plot.N)[2:3]  = c('Ring-stage','Late-stage')

red.model.fits = bind_rows(red.plot.A,red.plot.N) %>% gather(.,"Stage","Parasitaemia",2:3)
full.model.fits = bind_rows(full.plot.A,full.plot.N) %>% gather(.,"Stage","Parasitaemia",2:3)

dat_text1 <- data.frame(
  label = c(paste(expression(lambda[A])), "", "", ""),
  Host   = c("Acutely infected","Naive","Acutely infected","Naive"),
  Stage = c("Ring-stage","Ring-stage","Late-stage","Late-stage"),
  x     = rep(17,4),
  y     = rep(3.3,4)
)
dat_text2 <- data.frame(
  label = c(paste("=",round(para.full[5],2),sep=" "), "", "", ""),
  Host   = c("Acutely infected","Naive","Acutely infected","Naive"),
  Stage = c("Ring-stage","Ring-stage","Late-stage","Late-stage"),
  x     = rep(21,4),
  y     = rep(3.3,4)
)
dat_text3 <- data.frame(
  label = c(paste(expression(lambda[A])), "", "", ""),
  Host   = c("Acutely infected","Naive","Acutely infected","Naive"),
  Stage = c("Ring-stage","Ring-stage","Late-stage","Late-stage"),
  x     = rep(17,4),
  y     = rep(2.9,4)
)
dat_text4 <- data.frame(
  label = c(paste("=","1.00",sep=" "), "", "", ""),
  Host   = c("Acutely infected","Naive","Acutely infected","Naive"),
  Stage = c("Ring-stage","Ring-stage","Late-stage","Late-stage"),
  x     = rep(21,4),
  y     = rep(2.9,4)
)

dat_text5 <- data.frame(
  label = c("A", "c)", "b)","d)"),
  Host   = c("Acutely infected","Naive"),
  Stage = c("Ring-stage","Ring-stage","Late-stage","Late-stage"),
  x     = rep(1,4),
  y     = rep(3.8,4)
)

#Order facetting labels
DcDp.adj$Stage_f = factor(DcDp.adj$Stage, levels=c('Ring-stage','Late-stage'))
red.model.fits$Stage_f = factor(red.model.fits$Stage, levels=c('Ring-stage','Late-stage'))
full.model.fits$Stage_f = factor(full.model.fits$Stage, levels=c('Ring-stage','Late-stage'))

dat_text1$Stage_f = factor(dat_text1$Stage, levels=c('Ring-stage','Late-stage'))
dat_text2$Stage_f = factor(dat_text2$Stage, levels=c('Ring-stage','Late-stage'))
dat_text3$Stage_f = factor(dat_text3$Stage, levels=c('Ring-stage','Late-stage'))
dat_text4$Stage_f = factor(dat_text4$Stage, levels=c('Ring-stage','Late-stage'))
dat_text5$Stage_f = factor(dat_text5$Stage, levels=c('Ring-stage','Late-stage'))

#Make plot without text annotations
h = ggplot()+
  geom_ribbon(data=full.model.fits,aes(x=Time,ymax=Parasitaemia+full.sd,ymin=pmax(0,Parasitaemia-full.sd)),fill="seagreen",alpha=0.25)+
  geom_ribbon(data=red.model.fits,aes(x=Time,ymax=Parasitaemia+pooled.sd,ymin=pmax(0,Parasitaemia-red.sd)),fill="steelblue1",alpha=0.25)+
  geom_line(data=DcDp.adj,aes(x=Time,y=Parasitaemia,group=Mouse))+
  geom_point(data=DcDp.adj,aes(x=Time,y=Parasitaemia,group=Mouse),shape=21,size=4)+
  geom_line(data=full.model.fits,aes(x=Time,y=Parasitaemia,col='seagreen'),size=1)+
  geom_line(data=red.model.fits,aes(x=Time,y=Parasitaemia,col='steelblue1'),size=1)+
  theme_bw()+
  xlab('Hour post-transfusion')+ylab('Stage-structured parasitaemia (% RBC)')+
  facet_grid(Host~Stage_f)+
  theme(text = element_text(size=12),legend.position="none",strip.text = element_text(size=12))+ylim(0,4)+
  scale_colour_manual(values=c('seagreen','steelblue1'))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank())+
  scale_x_continuous(breaks=c(1,seq(4,24,by=4)))

#Add text annotations
i = h + geom_text(data=dat_text3,aes(x = x, y = y, label = label),parse=T,col='steelblue1',size=6)+
  geom_text(data=dat_text4,aes(x = x, y = y, label = label),col='steelblue1',size=6)+
  geom_text(data=dat_text1,aes(x = x, y = y, label = label),parse=T,col='seagreen',size=6)+
  geom_text(data=dat_text2,aes(x = x, y = y, label = label),col='seagreen',size=6)+
  geom_text(data=dat_text5,aes(x = x, y = y, label = label),size=6)
