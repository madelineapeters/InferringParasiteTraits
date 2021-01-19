#Load required packages
suppressMessages(library(tidyverse))
suppressMessages(library(reghelper))
suppressMessages(library(ggpattern))

#### Data wrangling ####
#Read in fits
inoc.full = read.csv("~/Desktop/Mideo.lab2/Parasite.maturation/Dryad/Inoculum/Inoculum_full.csv",stringsAsFactors=FALSE) %>% select(.,1:9)
names(inoc.full) = paste(names(inoc.full),"full",sep=".")
inoc.red = read.csv("~/Desktop/Mideo.lab2/Parasite.maturation/Dryad/Inoculum/Inoculum_reduced.csv",stringsAsFactors=FALSE) %>% select(.,1:9)
names(inoc.red) = paste(names(inoc.red),"red",sep=".")

#Transform variables from log values and bind columns.
for (r in 1:nrow(inoc.red)){
  
  row.vals = inoc.red[r,]
  
  new.vals = c(sapply(row.vals[1:4],exp),
               row.vals[5],
               sapply(row.vals[6:7],exp),
               row.vals[8:9]) 
  
  inoc.red[r,] = new.vals
}
for (r in 1:nrow(inoc.full)){
  
  row.vals = inoc.full[r,]
  
  new.vals = c(sapply(row.vals[1:4],exp),
               row.vals[5],
               sapply(row.vals[6:7],exp),
               row.vals[8:9]) 
  
  inoc.full[r,] = new.vals
}

inoc.joint = bind_cols(inoc.full,inoc.red)

#Calcuate F-test for each row and store resulting p-value.
## Define function for F-test (from Greischar et al. 2016 supplementary code)
getF = function(n, nparms1,sse1,nparms2,sse2){
  # function to return p-value associated with F statistic
  # which follows an F distribution with df1 = p-q, df2 = n-p
  # first determine which model has more parameters/
  # smaller sum squared error:
  if(nparms1>nparms2 & sse2>=sse1){
    sseOmega = sse1
    p = nparms1
    sseomega = sse2
    q = nparms2
    Fstatistic = ((sseomega-sseOmega)/(p-q))/(sseOmega/(n-p))
    pVal = pf(q=Fstatistic, df1=(p-q), df2=(n-p), lower.tail=FALSE)
  } # end test for model1 is more complex
  if(nparms2>nparms1 & sse1>=sse2){
    sseOmega = sse2
    p = nparms2
    sseomega = sse1
    q = nparms1
    Fstatistic = ((sseomega-sseOmega)/(p-q))/(sseOmega/(n-p))
    pVal = pf(q=Fstatistic, df1=(p-q), df2=(n-p), lower.tail=FALSE)
  } # end test for model2 is more complex
  if(nparms2>nparms1 & sse2 > sse1){
    print("Full model not sufficiently better than reduced model.")
    pVal = 1
  } # if model2 has more parameters but doesn't fit better, that's an error
  if(nparms1>nparms2 & sse1 > sse2){
    print("Full model not sufficiently better than reduced model.")
    pVal = 1
  } # if model2 has more parameters but doesn't fit better, that's an error
  return(pVal)
}

#Calculate p-value for F-test
inoc.joint$p = sapply(1:nrow(inoc.joint),FUN=function(r){
  ss1 = inoc.joint$fit.red[r]
  ss2 = inoc.joint$fit.full[r]
  p = getF(32, 6,ss1,7,ss2)
  return(p)
})

#Specify parameters used to make data frame of parameter combinations for focal data sets
omega = 9
fracA = 0.036
fracN = 0.031
rho = 150
pRA = 0.0045
pRN = 0.031
pRD = 0.01
kappaA = omega*(fracA*(rho*pRD+(1-pRD)))/((fracA*(rho*pRD+(1-pRD)))+(1-fracA)*(rho*pRA+(1-pRA)))
kappaN = omega*(fracN*(rho*pRD+(1-pRD)))/((fracN*(rho*pRD+(1-pRD)))+(1-fracN)*(rho*pRN+(1-pRN)))
gammaA = 0.9
gammaN = 0.45
sigma1 = 1
sigma2 = 1.1

#Add true parameter values
para.list = as.data.frame(matrix(nrow=55,ncol=10))
names(para.list) = c("kappaN","kappaA","gammaN","gammaA","sigma1","sigma2","pR","rho","epsilon","phi")

para.list$kappaN = kappaN
para.list$kappaA = kappaA
para.list$gammaN = gammaN
para.list$gammaA = gammaA
para.list$sigma1 = sigma1
para.list$sigma2 = sigma2
para.list$pR = pRD
para.list$rho = rho
para.list$epsilon = seq(0,1, by=0.25)
para.list$phi = 0
para.list$frac = c(rep(0,5),rep(0.031,5),rep(0.1,5),rep(0.2,5),rep(0.3,5),rep(0.4,5),rep(0.5,5),rep(0.6,5),rep(0.7,5),rep(0.8,5),rep(0.9,5)) 
para.list$bias = "stages"
para.list$xc = 0.5

inoc.joint = bind_cols(inoc.joint,para.list)

#Create data frame that contains best model for each row (i.e., full OR reduced parasite maturation model) based on p-value and correct bias labels. Create "Score" column for logistic regression later on.
inoc.optim.full =  filter(inoc.joint,p<0.05) %>% 
  mutate(.,MSS = fit.red - fit.full) %>% 
  select(.,-contains("red"))
inoc.optim.red =  filter(inoc.joint,p>=0.05) %>% 
  mutate(.,MSS = fit.red - fit.full) %>% 
  select(.,-contains("full"))
names.list = c(paste(c("lambdaA","kN","kA","c","mu","sigma","P0"),"fit",sep="."),"Value","Convergence","p",names(para.list)[1:8],"epsilon","phi","frac","bias","xc","MSS")
names(inoc.optim.full) = names.list
names(inoc.optim.red) = names.list

inoc.best = bind_rows(inoc.optim.full,inoc.optim.red)
inoc.best$index = factor(inoc.best$frac,levels=c(0,0.031,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),labels=c(0,1,2,3,4,5,6,7,8,9,10))
inoc.best$index = as.integer(inoc.best$index)
inoc.best = filter(inoc.best,frac<0.9)

#### Plotting ####

#Create heatmap figure showing estimated developmental rate as a function of stage misclassification and the fraction of reticulocytes in the transfused parasite cohort.
lambdaA.scale_fill = scale_fill_gradient2(low="#F0E442",mid="#009E73",midpoint=31,high="#CC79A7",name="Cycle \nlength \n(hours)",limits=c(24,38),breaks=c(24,32,38),labels=c("24","32","38"))

epsilon.frac.lambdaA = ggplot()+
  geom_tile(data=inoc.best,aes(x=epsilon,y=index,fill=24/lambdaA.fit))+
  geom_text(data=filter(inoc.best,lambdaA.fit!=1),aes(x=epsilon,y=index,label=round(24/lambdaA.fit,0)),size=5)+
  theme_classic()+
  ylab("% reticulocytes in \ntransfused parasite cohort")+
  xlab(expression("Stage misclassification"~epsilon~"(%)"))+
  scale_size_continuous(guide=FALSE)+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11),labels=c("0","3","10","20","30","40","50","60","70","80","90"))+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+lambdaA.scale_fill+
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size=12),
    legend.text = element_text(size=12),
    axis.title = element_text(size=16),
    legend.title = element_text(size=15)
  )

#Create heatmap figure showing estimated clearance and replication rates as a function of stage misclassification and the fraction of reticulocytes in the transfused parasite cohort.
c.scale_fill = scale_fill_gradient2(low="#F0E442",mid="#009E73",midpoint=0.45,high="#CC79A7",name=expression(italic(c)),limits=c(0,0.9),breaks=c(0,0.45,0.9),labels=c("0","0.45","0.9"))

epsilon.frac.c = ggplot()+
  geom_tile(data=filter(inoc.best,index<11),aes(x=epsilon,y=index,fill=c.fit))+
  geom_text(data=filter(inoc.best,index<11),aes(x=epsilon,y=index,label=round(c.fit,2)),size=5)+
  theme_classic()+
  ylab("% reticulocytes in \ntransfused parasite cohort")+
  xlab(expression("Stage misclassification"~epsilon~"(%)"))+
  scale_size_continuous(guide=FALSE)+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11),labels=c("0","3","10","20","30","40","50","60","70","80","90"))+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+c.scale_fill+
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size=12),
    legend.text = element_text(size=12),
    axis.title = element_text(size=16),
    legend.title = element_text(size=15)
  )

kA.scale_fill = scale_fill_gradient2(low="#F0E442",mid="#009E73",midpoint=0.25,high="#CC79A7",name=expression(k[A]),limits=c(0,0.5),breaks=c(0,0.25,0.5),labels=c("0","0.25","0.5"))

epsilon.frac.kA = ggplot()+
  geom_tile(data=filter(inoc.best,index<11),aes(x=epsilon,y=index,fill=kA.fit))+
  geom_text(data=filter(inoc.best,index<11),aes(x=epsilon,y=index,label=round(kA.fit,2)),size=5)+
  theme_classic()+
  ylab("% reticulocytes in \ntransfused parasite cohort")+
  xlab(expression("Stage misclassification"~epsilon~"(%)"))+
  scale_size_continuous(guide=FALSE)+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10),labels=c("0","3","10","20","30","40","50","60","70","80"))+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+kA.scale_fill+
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size=12),
    legend.text = element_text(size=12),
    axis.title = element_text(size=16),
    legend.title = element_text(size=15)
  )