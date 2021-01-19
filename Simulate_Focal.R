model.select = "full" #"full" for fitting lambdaA, "reduced" for setting lambdaA to 1

##############################
# Define function for cluster
##############################
optim.fit.fun = function(i,model,data.A,data.N,xc.select){
  
  .libPaths("~/Rpkgs")
  library(dplyr)
  library(R.utils)
  
  #Define time series of interest
  timedata = c(1,4,6,10,14,18,21,24)
  
  #Define parasite maturation piecewise function
  mat.piecewise = function(x,t,lambdaA,Beta,c,x.c,mu,sigma,P0) { 
    P<-if (((x>=t*lambdaA) & (x<x.c))) {P0*dnorm(x-t*lambdaA,mean=mu,sd=sigma)} else if (((x>=t*lambdaA) & (x>=x.c) & ((x-t*lambdaA)<=x.c))) {
      P0*dnorm(x-t*lambdaA,mean=mu,sd=sigma)*exp(-(c/lambdaA)*(x-x.c))} else if (((x>=t*lambdaA) & (x>=x.c) & ((x-t*lambdaA)>x.c))) {
        P0*dnorm(x-t*lambdaA,mean=mu,sd=sigma)*exp(-c*t)} else if (((x<t*lambdaA) & (x<x.c) & ((t*lambdaA-x)>=(1-x.c)))) {
          Beta*P0*dnorm(1+x-t*lambdaA,mean=mu,sd=sigma)*exp(-(1*c/lambdaA)*(1-x.c))} else if (((x<t*lambdaA) & (x<x.c) & ((t*lambdaA-x)<(1-x.c)))) {
            Beta*P0*dnorm(1+x-t*lambdaA,mean=mu,sd=sigma)*exp(-(c/lambdaA)*(t*lambdaA-x-(1-1)*x.c))} else if (((x<t*lambdaA) & (x>=x.c) & ((t*lambdaA-x)>=(1-x.c)))) {
              Beta*P0*dnorm(1+x-t*lambdaA,mean=mu,sd=sigma)*exp(-(1*c/lambdaA)*(1-x.c))*exp(-(c/lambdaA)*(x-x.c))} else if (((x<t*lambdaA) & (x>=x.c) & ((t*lambdaA-x)<(1-x.c)))) {
                Beta*P0*dnorm(1+x-t*lambdaA,mean=mu,sd=sigma)*exp(-(c/lambdaA)*(t*lambdaA-1*x.c))}
    
    return(P)
  } #end function specifying maturation function
  
  #Define function that calculates fit of parasite maturation model to simulated data (paradf2); returns fit for one host type
  PM.to.DDE.full = function(paravec,timedata,paradf2,xc){
    
    out=1000
    
    tmp = try({
      
      Ringfun=function(t){
        integrand=Vectorize(function(x){
          mat.piecewise(x,t=t,lambdaA=paravec[1],Beta=paravec[2],c=paravec[3],x.c=xc,mu=paravec[4],sigma=paravec[5],P0=paravec[6])
        })
        integrated=integrate(integrand,lower=0,upper=0.5)$value
        return(integrated)
      }
      
      Latefun=function(t){
        integrand=Vectorize(function(x){
          mat.piecewise(x,t=t,lambdaA=paravec[1],Beta=paravec[2],c=paravec[3],x.c=xc,mu=paravec[4],sigma=paravec[5],P0=paravec[6])
        })
        integrated=integrate(integrand,lower=0.5,upper=1)$value
        return(integrated)
      }
      
      ###############
      Rings.I=sapply(timedata/24,Ringfun)
      Late.I=sapply(timedata/24,Latefun)
      
      I.ring.fit = sqrt((paradf2$Ring-Rings.I)^2)
      I.late.fit = sqrt((paradf2$Late-Late.I)^2)
      
      out = sum(I.ring.fit,I.late.fit)
      
    })
    if(inherits(try, "try-error")) {
      recoverFromFailure()
    }
    
    return(out)
    
  } #end function
  
  if (model == "reduced"){
    
    #Define joint fitting functions for use with optim
    joint.PM.fit.red = function(para,data.A,data.N,xc,bias){
      
      para = c(sapply(para[c(1:3)],FUN=exp),para[4],sapply(para[c(5:6)],FUN=exp))
      paraA = c(1,unlist(para[c(2:6)]))
      paraN = c(1,unlist(para[c(1,3:6)]))
      
      Acute.fit = PM.to.DDE.full(paraA,timedata,data.A,xc)
      Naive.fit = PM.to.DDE.full(paraN,timedata,data.N,xc)
      
      fit = Acute.fit + Naive.fit
      return(fit)
      
    }
    
    #Set upper and lower parameter bounds for fitting the reduced parasite maturation model
    theta_min = c(kN=0,kA=0,c=0,mu=-2,sigma=0.05,P0=4)
    theta_max = c(kN=1,kA=1,c=3,mu=2,sigma=3,P0=20)
    
    optim.con = NA
    while (!(optim.con %in% c(0,1))){
      
      #Generate initial guess of optim
      optim.guess = sapply(1:length(theta_min),FUN=function(x){runif(1,theta_min[x],theta_max[x])})
      optim.guess[c(1:3,5:6)] = sapply(c(1:3,5:6),FUN=function(x){log(optim.guess[x])})
      
      #Run optim for reduced model
      optim.ans = NULL
      tryCatch({
        optim.ans = withTimeout({
          optim(optim.guess,joint.PM.fit.red,data.A=data.A,data.N=data.N,xc=xc.select,bias=bias.select,method=c("Nelder-Mead"),control=c(maxit=1000))
        }, timeout = 300, onTimeout = "silent")
      }, TimeoutException = function(ex) {
        message("Timeout. Skipping.")
      })
      
      optim.con = optim.ans$convergence
      
    }
    
    while (optim.con == 1){
      
      optim.guess = optim.ans$par
      optim.ans = optim(optim.guess,joint.PM.fit.red,data.A=data.A,data.N=data.N,xc=xc.select,bias=bias.select,method=c("Nelder-Mead"),control=c(maxit=1000))
      optim.con = optim.ans$convergence
      
    }
    
    optim.guess = c(log(1),optim.ans$par)
    optim.val = optim.ans$value
    optim.con = optim.ans$convergence
    optim.list = c(optim.guess,optim.val,optim.con,model) 
    
  } else {
    
    #Define joint fitting functions for use with optim
    joint.PM.fit.full = function(para,data.A,data.N,xc,bias){
      
      if(exp(para[1])>1){fit = 1000} else {
        para = c(sapply(para[c(1:4)],FUN=exp),para[5],sapply(para[c(6:7)],FUN=exp))
        paraA = unlist(para[c(1,3:7)])
        paraN = c(1,unlist(para[c(2,4:7)]))
        
        Acute.fit = PM.to.DDE.full(paraA,timedata,data.A,xc)
        Naive.fit = PM.to.DDE.full(paraN,timedata,data.N,xc)
        fit = Acute.fit + Naive.fit
      }

      return(fit)
      
    }
    
    #Set upper and lower parameter bounds for fitting the full parasite maturation model
    theta_min = c(lambdaA=0.5,kN=0,kA=0,c=0,mu=-2,sigma=0.05,P0=4)
    theta_max = c(lambdaA=1,kN=1,kA=1,c=3,mu=2,sigma=3,P0=20)
    
    optim.con = NA
    while (!(optim.con %in% c(0,1))){
      
      #Generate initial guess of optim
      optim.guess = sapply(1:length(theta_min),FUN=function(x){runif(1,theta_min[x],theta_max[x])})
      optim.guess[c(1:4,6:7)] = sapply(c(1:4,6:7),FUN=function(x){log(optim.guess[x])})
      
      #Run optim for reduced model
      optim.ans = NULL
      tryCatch({
        optim.ans = withTimeout({
          optim(optim.guess,joint.PM.fit.full,data.A=data.A,data.N=data.N,xc=xc.select,bias=bias.select,method=c("Nelder-Mead"),control=c(maxit=1000))
        }, timeout = 300, onTimeout = "silent")
      }, TimeoutException = function(ex) {
        message("Timeout. Skipping.")
      })
      
      optim.con = optim.ans$convergence
      
    }
    
    while (optim.con == 1){
      
      optim.guess = optim.ans$par
      optim.ans = optim(optim.guess,joint.PM.fit.full,data.A=data.A,data.N=data.N,xc=xc.select,bias=bias.select,method=c("Nelder-Mead"),control=c(maxit=1000))
      optim.con = optim.ans$convergence
      
    }
    
    optim.guess = c(optim.ans$par)
    optim.val = optim.ans$value
    optim.con = optim.ans$convergence
    optim.list = c(optim.guess,optim.val,optim.con,model) 
  }
  
  return(optim.list)
  
}

############################
# Parameters
############################
#Specify parameters used to make data frame of parameter combinations for focal data sets
omega = 9 
fracA = 0.036 #Khoury et al. 2017
fracN = 0.031 #Khoury et al. 2017
rho = 150 #Cromer et al. 2006
pRA = 0.0045 #Khoury et al. 2014 SI
pRN = 0.031
pRD = 0.01
kappaA = omega*(fracA*(rho*pRD+(1-pRD)))/((fracA*(rho*pRD+(1-pRD)))+(1-fracA)*(rho*pRA+(1-pRA)))
kappaN = omega*(fracN*(rho*pRD+(1-pRD)))/((fracN*(rho*pRD+(1-pRD)))+(1-fracN)*(rho*pRN+(1-pRN)))
gammaA = 0.9 #Khoury et al. 2015
gammaN = 0.45 #Khoury et al. 2015
sigma1 = 1.1
sigma2 = 1.6
epsilon.list = seq(0,1,0.25)
phi.list = c(0,1,2,3)

#Create data frames for 104 unique parameter combinations
para.list = as.data.frame(matrix(nrow=104,ncol=10))
names(para.list) = c("kappaN","kappaA","gammaN","gammaA","sigma1","sigma2","pR","rho","epsilon","phi")

para.list$kappaN = kappaN
para.list$kappaA = kappaA
para.list$gammaN = gammaN
para.list$gammaA = gammaA
para.list$sigma1 = sigma1
para.list$sigma2 = sigma2
para.list$pR = pRD
para.list$rho = rho
para.list$epsilon = c(sapply(1:length(epsilon.list),FUN=function(i){
  rep(epsilon.list[i],4)
}),rep(NA,24))
para.list$phi = phi.list

para.list$bias = c(rep("full",20),rep("stages",20),rep("full",20),rep("stages",20),rep("sequestration",8),rep("parasitemia",8),rep("none",8))
para.list$xc = c(rep(0,40),rep(0.5,40),rep(0,4),rep(0.5,4),rep(0,4),rep(0.5,4),rep(0,4),rep(0.5,4))

############################
# Run cluster
############################
.libPaths("~/Rpkgs")
library(dplyr)
library(PBSddesolve)
library(parallel)

optim.joint = data.frame()
starting.i = nrow(optim.joint)+1

for (i in starting.i:104){
  
  ############################
  # Simulated data set-up
  ############################
  
  #Choose parameters from data frame specifying parameter combinations
  paraset = para.list[i,1:10]
  paraA = paraset[c(2,4:10)] #subset parameters for acutely infected host
  paraN = paraset[c(1,3,5:10)] #subset parameters for naive host
  #Specify bias to be applied to simulated data set
  bias.select = para.list$bias[i] 
  #Specify xc
  xc.select = as.double(para.list$xc[i])
  
  #Function to create simulated data set with specified measurement bias applied; returns "Time", "Ring" and "Late" columns.
  DDE.opt = function(parameters,bias){
    
    out = 1000
    
    parameters = unlist(parameters)
    
    ###DDE
    DDE.fun = function (t, y, parms){
      
      #Set parameters
      mu = parms[1]
      P0 = parms[2] #initial number of donor pRBCs
      
      #Fit parameters
      gamma = parms[3] #clearance rate
      kappa = parms[4] #invasion ratio (# new donor cells infected per burst infected donor cell)
      phi = parms[5] #sequestration rate
      sigma1 = parms[6] #shape 1 term for Beta distribution of parasite age
      sigma2 = parms[7] #shape 2 term
      rho = parms[8]
      
      P0H = parms[9]
      
      #Variables to track (this is where R sets these equal to the initial values or updates them)
      #Uninfected cells
      UN = y[1] #normocytes
      UR = y[2] #reticulocytes
      #Ring-stage cells
      RN = y[3]
      RR = y[4]
      #Troph-stage cells
      LTN = y[5]
      LTR = y[6]
      #Schizont-stage cells
      LSN = y[7]
      LSR = y[8]
      #Sequestered cells
      SQN = y[9]
      SQR = y[10]
      #Original infected cells
      PR = y[11]
      PT = y[12]
      PS = y[13]
      PSQ = y[14]
      #Host-parasite infected cells
      HR = y[15]
      HLT = y[16]
      HLS = y[17]
      HSQ = y[18]
      
      #Past time values
      if (t > 0.5){lag0.5 = pastvalue(t-0.5)}
      if (t > 0.75){lag0.75 = pastvalue(t-0.75)}
      
      #Transfused population survival functions
      
      if (t <= 0.25){
        ST = exp(-gamma*t)
        SM = exp(-(gamma+phi)*t)
        SQ = (phi/(phi+gamma))*(1-SM)
      }
      if ((t > 0.25)&(t <= 0.5)){
        ST = exp(-gamma*0.25)
        SM = exp(-(phi*0.25+gamma*t))
        SQ = (phi/(phi+gamma))*(1-exp(-(phi*0.25+gamma*0.25)))*exp(-gamma*(t-0.25))
      }
      if ((t > 0.5)){
        ST = exp(-gamma*0.25)
        SM = exp(-(phi*0.25+gamma*0.5))
        SQ = (phi/(phi+gamma))*(1-exp(-(phi*0.25+gamma*0.25)))*exp(-gamma*0.25)
      }
      
      #Transfused population differential equations
      dPRdt = - pulsebeta(P0,sigma1,sigma2,0.5-t)
      dPTdt = pulsebeta(P0,sigma1,sigma2,0.5-t) - gamma*PT - pulsebeta(P0,sigma1,sigma2,0.75-t)*ST
      dPSdt = pulsebeta(P0,sigma1,sigma2,0.75-t)*ST - gamma*PS - phi*PS - pulsebeta(P0,sigma1,sigma2,1-t)*SM
      dPSQdt = phi*PS - pulsebeta(P0,sigma1,sigma2,1-t)*SQ
      
      #Uninfected population differential equations
      dUNdt = - kappa*UN/(UN+rho*UR)*pulsebeta(P0,sigma1,sigma2,1-t)*(SQ+SM) - kappa*UN/(UN+rho*UR)*pulsebeta(P0H,1,1,1-t)
      dURdt = - kappa*rho*UR/(UN+rho*UR)*pulsebeta(P0,sigma1,sigma2,1-t)*(SQ+SM) - kappa*rho*UR/(UN+rho*UR)*pulsebeta(P0H,1,1,1-t)
      
      #Newly infected population differential equations
      dRNdt = kappa*UN/(UN+rho*UR)*pulsebeta(P0,sigma1,sigma2,1-t)*(SQ+SM)
      dRRdt = kappa*rho*UR/(UN+rho*UR)*pulsebeta(P0,sigma1,sigma2,1-t)*(SQ+SM) 
      dLTNdt = 0
      dLTRdt = 0
      dLSNdt = 0
      dLSRdt = 0
      dSQNdt = 0
      dSQRdt = 0
      
      #Donor cells infected with host parasites
      dHRdt = kappa*UN/(UN+rho*UR)*pulsebeta(P0H,1,1,1-t) + kappa*rho*UR/(UN+rho*UR)*pulsebeta(P0H,1,1,1-t)
      dHLTdt = 0
      dHLSdt = 0
      dHSQdt = 0
      
      if ((t>0.5)&(t<=0.75)){
        
        dRNdt = dRNdt - kappa*lag0.5[1]/(lag0.5[1]+rho*lag0.5[2])*pulsebeta(P0,sigma1,sigma2,1.5-t)*((phi/(phi+gamma))*(1-exp(-phi*(t-0.5)))+exp(-phi*(t-0.5))*exp(-gamma*(t-0.5))) 
        
        dRRdt = dRRdt - kappa*rho*lag0.5[2]/(lag0.5[1]+rho*lag0.5[2])*pulsebeta(P0,sigma1,sigma2,1.5-t)*((phi/(phi+gamma))*(1-exp(-phi*(t-0.5)))+exp(-phi*(t-0.5))*exp(-gamma*(t-0.5)))
        
        dLTNdt = kappa*lag0.5[1]/(lag0.5[1]+rho*lag0.5[2])*pulsebeta(P0,sigma1,sigma2,1.5-t)*((phi/(phi+gamma))*(1-exp(-phi*(t-0.5)))+exp(-phi*(t-0.5))*exp(-gamma*(t-0.5))) - gamma*LTN
        
        dLTRdt = kappa*rho*lag0.5[2]/(lag0.5[1]+rho*lag0.5[2])*pulsebeta(P0,sigma1,sigma2,1.5-t)*((phi/(phi+gamma))*(1-exp(-phi*(t-0.5)))+exp(-phi*(t-0.5))*exp(-gamma*(t-0.5))) - gamma*LTR
        
        dHRdt = dHRdt - kappa*lag0.5[1]/(lag0.5[1]+rho*lag0.5[2])*pulsebeta(P0H,1,1,1.5-t) - kappa*rho*lag0.5[2]/(lag0.5[1]+rho*lag0.5[2])*pulsebeta(P0H,1,1,1.5-t)
        
        dHLTdt = kappa*lag0.5[1]/(lag0.5[1]+rho*lag0.5[2])*pulsebeta(P0H,1,1,1.5-t) + kappa*rho*lag0.5[2]/(lag0.5[1]+rho*lag0.5[2])*pulsebeta(P0H,1,1,1.5-t) - gamma*HLT
        
      }
      if (t>0.75){
        
        dRNdt = dRNdt - kappa*lag0.5[1]/(lag0.5[1]+rho*lag0.5[2])*pulsebeta(P0,sigma1,sigma2,1.5-t)*((phi/(phi+gamma))*(1-exp(-phi*0.25))+exp(-phi*0.25)*exp(-gamma*0.25))*exp(-gamma*(t-0.75))
        
        dRRdt = dRRdt - kappa*rho*lag0.5[2]/(lag0.5[1]+rho*lag0.5[2])*pulsebeta(P0,sigma1,sigma2,1.5-t)*((phi/(phi+gamma))*(1-exp(-phi*0.25))+exp(-phi*0.25)*exp(-gamma*0.25))*exp(-gamma*(t-0.75))
        
        dLTNdt = kappa*lag0.5[1]/(lag0.5[1]+rho*lag0.5[2])*pulsebeta(P0,sigma1,sigma2,1.5-t)*((phi/(phi+gamma))*(1-exp(-phi*0.25))+exp(-phi*0.25)*exp(-gamma*0.25))*exp(-gamma*(t-0.75)) - 
          gamma*LTN - kappa*lag0.75[1]/(lag0.75[1]+rho*lag0.75[2])*pulsebeta(P0,sigma1,sigma2,1.75-t)*((phi/(phi+gamma))*(1-exp(-phi*(t-0.75)))+exp(-phi*(t-0.75))*exp(-gamma*(t-0.75)))*exp(-0.25*gamma)
        
        dLTRdt = kappa*rho*lag0.5[2]/(lag0.5[1]+rho*lag0.5[2])*pulsebeta(P0,sigma1,sigma2,1.5-t)*((phi/(phi+gamma))*(1-exp(-phi*0.25))+exp(-phi*0.25)*exp(-gamma*0.25))*exp(-gamma*(t-0.75)) - 
          gamma*LTR - kappa*rho*lag0.75[2]/(lag0.75[1]+rho*lag0.75[2])*pulsebeta(P0,sigma1,sigma2,1.75-t)*((phi/(phi+gamma))*(1-exp(-phi*(t-0.75)))+exp(-phi*(t-0.75))*exp(-gamma*(t-0.75)))*exp(-0.25*gamma)
        
        dLSNdt = kappa*lag0.75[1]/(lag0.75[1]+rho*lag0.75[2])*pulsebeta(P0,sigma1,sigma2,1.75-t)*((phi/(phi+gamma))*(1-exp(-phi*(t-0.75)))+exp(-phi*(t-0.75))*exp(-gamma*(t-0.75)))*exp(-0.25*gamma) - phi*LSN - gamma*LSN
        
        dLSRdt = kappa*rho*lag0.75[2]/(lag0.75[1]+rho*lag0.75[2])*pulsebeta(P0,sigma1,sigma2,1.75-t)*((phi/(phi+gamma))*(1-exp(-phi*(t-0.75)))+exp(-phi*(t-0.75))*exp(-gamma*(t-0.75)))*exp(-0.25*gamma) - phi*LSR - gamma*LSR
        
        dSQNdt = phi*LSN
        
        dSQRdt = phi*LSR
        
        dHRdt = dHRdt - kappa*lag0.5[1]/(lag0.5[1]+rho*lag0.5[2])*pulsebeta(P0H,1,1,1.5-t) - kappa*rho*lag0.5[2]/(lag0.5[1]+rho*lag0.5[2])*pulsebeta(P0H,1,1,1.5-t)
        
        dHLTdt = kappa*lag0.5[1]/(lag0.5[1]+rho*lag0.5[2])*pulsebeta(P0H,1,1,1.5-t) + kappa*rho*lag0.5[2]/(lag0.5[1]+rho*lag0.5[2])*pulsebeta(P0H,1,1,1.5-t) - gamma*HLT -
          kappa*lag0.75[1]/(lag0.75[1]+rho*lag0.75[2])*pulsebeta(P0H,1,1,1.75-t)*exp(-gamma*0.25) - kappa*rho*lag0.75[2]/(lag0.75[1]+rho*lag0.75[2])*pulsebeta(P0H,1,1,1.75-t)*exp(-gamma*0.25)
        
        dHLSdt = kappa*lag0.75[1]/(lag0.75[1]+rho*lag0.75[2])*pulsebeta(P0H,1,1,1.75-t)*exp(-gamma*0.25) + kappa*rho*lag0.75[2]/(lag0.75[1]+rho*lag0.75[2])*pulsebeta(P0H,1,1,1.75-t)*exp(-gamma*0.25) - gamma*HLS - phi*HLS
        
        dHSQdt = phi*HLS
        
      }
      
      return(c(dUNdt, dURdt,dRNdt,dRRdt,dLTNdt,dLTRdt,dLSNdt,dLSRdt,dSQNdt,dSQRdt,dPRdt,dPTdt,dPSdt,dPSQdt,dHRdt,dHLTdt,dHLSdt,dHSQdt))
      
    } # end of age-structured initial infection function # end of age-structured initial infection function
    
    timedata = c(1 , 4 , 7 ,10 ,14 ,18 ,21 ,24)
    
    tmp = try({
      
      ###Parameters
      #Proportion reticulocytes
      pRD = parameters[5]
      #Reticulocyte preference
      rho = parameters[6]
      
      #Acute donor invasion rate
      kappa.eff = parameters[1]
      
      #Clearance rate
      gamma = parameters[2] 
      
      #Sequestration rate
      phi = parameters[8]
      
      #Beta distribution parameters (initial parasite age distribution)
      sigma1 = parameters[3]
      sigma2 = parameters[4]
      
      #Natural death
      mu = 0.025
      
      #Original parasitaemia
      P0 = 50
      P0H = 0#parameters[9]
      
      #Rate describing probability of correct reticulocyte parasite age classification according to exponetial distribution
      epsilon = parameters[7]
      
      parametersA = unname(unlist((c(mu,P0,gamma,kappa.eff,phi,sigma1,sigma2,rho,P0H))))
      
      ###Functions
      pulsebeta = function(initialI, sigma1, sigma2, time){
        res = rep(NA, length(time))
        for (num in 1:length(time)){
          res[num] = initialI*(dbeta(time[num], sigma1, sigma2))
        }
        return(res)
      }
      
      ###DDE set-up
      N0 = (1000-P0)*(1-pRD)
      R0 = (1000-P0)*pRD
      
      cbeta = function(t){
        dbeta(t,sigma1,sigma2)
      }
      PR0 = integrate(cbeta,lower=0,upper=0.5)[[1]]
      PT0 = integrate(cbeta,lower=0.5,upper=0.75)[[1]]
      PS0 = integrate(cbeta,lower=0.75,upper=1)[[1]]
      yinit = c(y1 = N0, #initital state variable values
                y2 = R0,
                y3 = 0,
                y4 = 0,
                y5 = 0,
                y6 = 0,
                y7 = 0,
                y8 = 0,
                y9 = 0,
                y10 = 0,
                y11 = P0*PR0,
                y12 = P0*PT0,
                y13 = P0*PS0,
                y14 = 0,
                y15 = 0,
                y16 = 0,
                y17 = 0,
                y18 = 0)
      times = seq(0,1,0.01)
      
      ###Run DDE
      ddeA = dde(y = yinit, times = times, func = DDE.fun, parms = parametersA, hbsize = 50000)
      
      ###Clean dde output
      out.A = as.data.frame(ddeA)
      names(out.A) = c('Time','UN','UR','RN',"RR","LTN","LTR","LSN","LSR","SQN","SQR","PR","PT","PS","PQ","HR","HLT","HLS","HSQ")
      
      if (bias == "full"){
        
        out.A = out.A %>% 
          mutate(.,Total = UN+UR+RN+RR+LTN+LTR+LSN+LSR+PR+PT+PS+HR+HLT+HLS) %>% 
          mutate(.,Ring = RN+RR+PR) %>% 
          mutate(.,Late = LTN+LTR+LSN+LSR+PT+PS)
        
        out.A = out.A %>% mutate(.,Ring.adj = Ring - epsilon*(RR)) %>% mutate(.,Late.adj = Late + epsilon*(RR))
        
        dataA = out.A %>% 
          select(.,Time,Ring.adj,Late.adj,Total) %>%
          mutate(.,Ring = 100*Ring.adj/Total) %>%
          mutate(.,Late = 100*Late.adj/Total) %>%
          select(.,Time,Ring,Late)
        
      } else if (bias == "sequestration"){ #only sequestration bias
        
        out.A = out.A %>% 
          mutate(.,Total = UN+UR+RN+RR+LTN+LTR+LSN+LSR+PR+PT+PS+HR+HLT+HLS) %>% 
          mutate(.,Ring = RN+RR+PR) %>% 
          mutate(.,Late = LTN+LTR+LSN+LSR+PT+PS)
        
        out.A = out.A %>% mutate(.,Ring.adj = Ring) %>% mutate(.,Late.adj = Late)
        
        dataA = out.A %>% 
          select(.,Time,Ring.adj,Late.adj,Total) %>%
          mutate(.,Ring = 100*Ring.adj/1000) %>%
          mutate(.,Late = 100*Late.adj/1000) %>%
          select(.,Time,Ring,Late)
        
      } else if (bias == "stages"){ #only stage-structure bias
        
        out.A = out.A %>% 
          mutate(.,Total = UN+UR+RN+RR+LTN+LTR+LSN+LSR+PR+PT+PS+HR+HLT+HLS+SQN+SQR+PQ) %>% 
          mutate(.,Ring = RN+RR+PR) %>% 
          mutate(.,Late = LTN+LTR+LSN+LSR+PT+PS+PQ+SQN+SQR)
        
        out.A = out.A %>% mutate(.,Ring.adj = Ring - epsilon*(RR)) %>% mutate(.,Late.adj = Late + epsilon*(RR))
        
        dataA = out.A %>% 
          select(.,Time,Ring.adj,Late.adj,Total) %>%
          mutate(.,Ring = 100*Ring.adj/1000) %>%
          mutate(.,Late = 100*Late.adj/1000) %>%
          select(.,Time,Ring,Late)
        
      } else if (bias == "parasitaemia") {
        
        out.A = out.A %>% 
          mutate(.,Total = UN+UR+RN+RR+LTN+LTR+LSN+LSR+PR+PT+PS+HR+HLT+HLS+SQN+SQR+PQ) %>% 
          mutate(.,Ring = RN+RR+PR) %>% 
          mutate(.,Late = LTN+LTR+LSN+LSR+PT+PS+PQ+SQN+SQR)
        
        out.A = out.A %>% mutate(.,Ring.adj = Ring) %>% mutate(.,Late.adj = Late)
        
        dataA = out.A %>% 
          select(.,Time,Ring.adj,Late.adj,Total) %>%
          mutate(.,Ring = 100*Ring.adj/Total) %>%
          mutate(.,Late = 100*Late.adj/Total) %>%
          select(.,Time,Ring,Late)
        
      } else if (bias == "none") {
        
        out.A = out.A %>% 
          mutate(.,Total = UN+UR+RN+RR+LTN+LTR+LSN+LSR+PR+PT+PS+HR+HLT+HLS+SQN+SQR+PQ) %>% 
          mutate(.,Ring = RN+RR+PR) %>% 
          mutate(.,Late = LTN+LTR+LSN+LSR+PT+PS+PQ+SQN+SQR)
        
        out.A = out.A %>% mutate(.,Ring.adj = Ring) %>% mutate(.,Late.adj = Late)
        
        dataA = out.A %>% 
          select(.,Time,Ring.adj,Late.adj,Total) %>%
          mutate(.,Ring = 100*Ring.adj/1000) %>%
          mutate(.,Late = 100*Late.adj/1000) %>%
          select(.,Time,Ring,Late)
        
      } else {
        print(paste("Incorrect bias type. Options include 'full', 'sequestration', 'stages', parasitaemia', 'none'.",sep=""))
      }
      
      out = dataA[c(4 , 16 , 29  ,42 , 58 , 75 , 88 ,100),]
      out$Time = timedata
      
    })
    if(inherits(try, "try-error")) {
      recoverFromFailure()
    }
    return(out)
}

  #Create simulated data sets using parameter subsets for each host type and specified bias
  data.A = DDE.opt(paraA,bias.select)
  data.N = DDE.opt(paraN,bias.select)

  ############################
  # Run cluster f*nction
  ############################
  cl = makeCluster(40)
  ID.list = 1:1000
  
  detest = clusterApplyLB(cl, ID.list, optim.fit.fun, model=model.select, data.A, data.N, xc.select)
  stopCluster(cl)
  optim.df = as.data.frame(matrix(unlist(detest), nrow = length(ID.list), ncol = round(length(unlist(detest)) / length(ID.list)), byrow = T))
  
  optim.df = as.data.frame(sapply(1:9,FUN=function(x){
    as.numeric(as.vector(optim.df[,x]))
  }))
  
  optim.best = optim.df %>% top_n(-1,V8) %>% mutate(.,model=model.select)
  names(optim.best) = c("lambdaA","kN","kA","c","mu","sigma","P0","fit","con","model")
  optim.joint = bind_rows(optim.joint,optim.best)
  
  write.csv(optim.joint,paste("Focal_",model.select,"csv",sep="."),row.names=FALSE)
}

