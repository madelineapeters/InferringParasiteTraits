#Load required packages
suppressMessages(library(tidyverse))
suppressMessages(library(reghelper))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))

#### Data wrangling ####
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

#Create data frames for 20 unique parameter combinations
variable.list = c("kappaN","kappaA","gammaN","gammaA","alpha","beta","pR","rho","epsilon","phi","frac","bias","xc")
variable = variable.list[c]

para.list = as.data.frame(matrix(nrow=50,ncol=13))
names(para.list) = c("kappaN","kappaA","gammaN","gammaA","alpha","beta","pR","rho","epsilon","phi","frac","bias","xc")

para.list$kappaN = kappaN
para.list$kappaA = kappaA
para.list$gammaN = gammaN
para.list$gammaA = gammaA
para.list$alpha = sigma1
para.list$beta = sigma2
para.list$pR = pRD
para.list$rho = rho
para.list$epsilon = NA
para.list$phi = 0
para.list$frac = NA
para.list$bias = "stages"
para.list$xc = 0.5

#Specify varied parameter (col 1 = epsilon, col 2 = frac retics)
epsilon.frac.opt = matrix(c(1,0.5,0.5,0.5,0.25,0.8),nrow=3,byrow=T)

#Create empty data frames for full and reduced model results
full.df = data.frame()
red.df = data.frame()

#Iterate over parameter combinations
for (set in 1:nrow(epsilon.frac.opt)){
  
  #Read in data frames for full model
  fit.list = list.files("~/Sensitivity",pattern=paste("Sensitivity_epsilon_",epsilon.frac.opt[set,1],"_frac_",epsilon.frac.opt[set,2],".*.full.csv",sep=""))
  
  sens.df = data.frame()
  
  for (i in 1:length(fit.list)){
    
    fit.name = fit.list[i]
    #get parameter
    para = strsplit(fit.name,paste("Sensitivity_epsilon_",epsilon.frac.opt[set,1],"_frac_",epsilon.frac.opt[set,2],".",sep=""))[[1]][2] %>% strsplit(.,".full.csv") %>% unlist
    
    #read in fits
    fit.df = read.csv(paste("~/Sensitivity",fit.name,sep="/"),stringsAsFactors=FALSE)
    
    if (set == 1){
      c = which(para == c("kappaN","kappaA","gammaN","gammaA","sigma1","sigma2","pR","rho","epsilon","phi","frac","bias","xc"))
      
      para.set = read.csv(paste("~/Sensitivity/para_epsilon_",epsilon.frac.opt[set,1],"_frac_",epsilon.frac.opt[set,2],".",para,".csv",sep=""))
      lower = c(kappaN,kappaA,gammaN,0,sigma1,sigma2,pR,rho,epsilon,0)[c]
      upper = c(1,1,3,0.9,2,sigma2,0.03,300,epsilon,3)[c]
      fit.df = fit.df[para.set[,c]>=lower&para.set[,c]<=upper,]
      para.set = para.set[para.set[,c]>=lower&para.set[,c]<=upper,]
    } else {
    
    #get parameter index
    c = which(para == c("kappaN","kappaA","gammaN","gammaA","sigma1","sigma2","pR","rho","epsilon","phi","frac","bias","xc"))
    
    para.set = para.list
    para.set$epsilon = epsilon.frac.opt[set,1]
    para.set$frac = epsilon.frac.opt[set,2]
    
    lower = c(kappaN,kappaA,gammaN,0,sigma1,sigma2,pR,rho,epsilon,0)[c]
    upper = c(1,1,3,0.9,2,sigma2,0.03,300,epsilon,3)[c]
    
    para.set[,c] = seq(lower,upper,length.out=50)[1:nrow(fit.df)]
    
    }
    
    para.df = select(para.set,c[1])
    names(para.df) = "value"
    para.df$parameter = para
    
    fit.df = bind_cols(para.df,fit.df)
    if (set == 1){fit.df = sample_n(fit.df,50,replace=FALSE)}
    
    names(fit.df)[-c(1,2,15)] = paste(names(fit.df)[-c(1,2,15)],"full",sep=".")
    sens.df = bind_rows(sens.df,fit.df)
    sens.df$Set = set
  }
  full.df = bind_rows(full.df,sens.df)
  
  #Read in data frames for reduced model
  fit.list = list.files("~/Sensitivity",pattern=paste("Sensitivity_epsilon_",epsilon.frac.opt[set,1],"_frac_",epsilon.frac.opt[set,2],".*.reduced.csv",sep=""))
  
  sens.df = data.frame()
  
  for (i in 1:length(fit.list)){
    
    fit.name = fit.list[i]
    #get parameter
    para = strsplit(fit.name,paste("Sensitivity_epsilon_",epsilon.frac.opt[set,1],"_frac_",epsilon.frac.opt[set,2],".",sep=""))[[1]][2] %>% strsplit(.,".reduced.csv") %>% unlist
    
    #read in fits
    fit.df = read.csv(paste("~/Sensitivity",fit.name,sep="/"),stringsAsFactors=FALSE)
    names(fit.df)[-13] = paste(names(fit.df)[-13],"reduced",sep=".")
    fit.df$parameter = para
    fit.df$Set = set
    sens.df = bind_rows(sens.df,fit.df)
    
  }
  red.df = bind_rows(red.df,sens.df) %>% distinct()
}

#Transform variables from log values and bind columns.
for (r in 1:nrow(full.df)){
  
  row.vals = full.df[r,]
  
  new.vals = c(row.vals[1:2],sapply(row.vals[3:6],exp),
               row.vals[7],
               sapply(row.vals[8:9],exp),
               row.vals[10:16]) 
  
  full.df[r,] = new.vals
}

for (r in 1:nrow(red.df)){
  
  row.vals = red.df[r,]
  
  new.vals = c(sapply(row.vals[1:4],exp),
               row.vals[5],
               sapply(row.vals[6:7],exp),
               row.vals[8:15]) 
  
  red.df[r,] = new.vals
}

full.df$parameter = factor(full.df$parameter,levels=c("kappaA","kappaN","gammaA","gammaN","sigma1","sigma2","phi","rho","pR"),labels=c("kappa[A]","kappa[N]","gamma[A]","gamma[N]","sigma[1]","sigma[2]","phi","rho","d[R]"))
red.df$parameter = factor(red.df$parameter,levels=c("kappaA","kappaN","gammaA","gammaN","sigma1","sigma2","phi","rho","pR"),labels=c("kappa[A]","kappa[N]","gamma[A]","gamma[N]","sigma[1]","sigma[2]","phi","rho","d[R]"))

sens.joint = left_join(full.df,red.df,by=c("parameter","index","Set"))

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
    print("Error: full model not sufficiently better than reduced model.")
    pVal = 1
  } # if model2 has more parameters but doesn't fit better, that's an error
  if(nparms1>nparms2 & sse1 > sse2){
    print("Error: full model not sufficiently better than reduced model.")
    pVal = 1
  } # if model2 has more parameters but doesn't fit better, that's an error
  return(pVal)
}
sens.joint$p = sapply(1:nrow(sens.joint),FUN=function(r){
  ss1 = sens.joint$fit.reduced[r]
  ss2 = sens.joint$fit.full[r]
  p = getF(32, 6,ss1,7,ss2)
  return(p)
})
sens.joint$p[is.na(sens.joint$p)]=1

#Create data frame that contains best model for each row (i.e., full OR reduced parasite maturation model) based on p-value and correct bias labels.
sens.optim.full =  filter(sens.joint,p<0.05) %>% 
  mutate(.,Index = index) %>% mutate(.,BestModel = "full") %>% mutate(.,set=Set) %>% select(.,-contains("red"),-index,-Set)
sens.optim.red =  filter(sens.joint,p>=0.05) %>% 
  mutate(.,Index = index) %>% mutate(.,BestModel = "reduced") %>% mutate(.,set=Set) %>% select(.,-contains("full"),-index,-Set)
names.list = c("Value","Parameter",paste(c("lambdaA","kN","kA","c","mu","sigma","P0"),"fit",sep="."),"Fit","Convergence","Model","Epsilon","Frac","p","Index","BestModel","Set")
names(sens.optim.full) = names.list
names(sens.optim.red) = names.list

sens.best = bind_rows(sens.optim.full,sens.optim.red)
sens.best$BestModel = factor(sens.best$BestModel,levels=c("full","reduced"),labels=c("Full model \n(delayed maturation)","Reduced \nmodel"))

#### Plotting ####
#Plot developmental rate as a function of first shape parameter and parasite stage structure at various shape parameter values (Figure S4)
alpha.df = sens.best %>% filter(.,Parameter=="sigma[1]")
alpha.df$FracRing = sapply(1:nrow(alpha.df),FUN=function(r){
  sigma1 = alpha.df$Value[r]
  sigma2 = 1.1
  pbeta(0.5,sigma1,sigma2)
})

alpha.cycle = ggplot()+geom_point(data=filter(alpha.df,Value>0.5,Set==1),aes(x=Value,y=24/lambdaA.fit,col=Model),size=2)+
  ylab("Estimated cycle length (hours)")+
  scale_x_continuous(name=expression(sigma[1]),sec.axis=sec_axis(~ pbeta(0.5,.,1.1),breaks=c(pbeta(0.5,0.5,1.1),pbeta(0.5,1,1.1),pbeta(0.5,1.5,1.1),pbeta(0.5,2,1.1),pbeta(0.5,2.5,1.1),pbeta(0.5,3,1.1),pbeta(0.5,3.5,1.1),pbeta(0.5,3.78,1)),name="Percent ring-stages",labels=c("73%","53%","39%","28%","20%","14%","10%","7%")))+
  scale_y_continuous(breaks=c(24,28,32,36,40))+
  scale_color_manual(values=c("black","darkgrey"))+
  theme_bw()+
  theme(
    strip.background=element_blank(),
    strip.text = element_text(size=14),
    axis.title = element_text(size=15),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = c(0.75,0.9),
    legend.title = element_blank(),
    legend.background = element_rect(colour="black"),
    legend.text = element_text(size=12),
    axis.text= element_text(size=13)
  )

structure.df = as.data.frame(matrix(c(
  dbeta(seq(0,1,0.01),0.5,1.1),
  dbeta(seq(0,1,0.01),1,1.1),
  dbeta(seq(0,1,0.01),1.5,1.1),
  dbeta(seq(0,1,0.01),2,1.1),
  dbeta(seq(0,1,0.01),3,1.1)
)),ncol=1)

structure.df$alpha = c(rep(0.5,101),
                       rep(1,101),
                       rep(1.5,101),
                       rep(2,101),
                       rep(3,101))
structure.df$age = seq(0,1,0.01)

label.df = structure.df %>% slice(.,c(101,202,303,404,505,606))
label.df$label = factor(label.df$alpha,levels=c(0.5,1,1.5,2,3),labels=c(expression(sigma[1] == 0.5),expression(sigma[1] == 1),expression(sigma[1] == 1.5),expression(sigma[1] == 2),expression(sigma[1] == 3)))

structure.facet = ggplot()+
  geom_line(data=structure.df,aes(x=age,y=10*V1,col=factor(alpha)),size=2)+
  ylab("")+
  geom_text(data=label.df,aes(x=0.65,y=40,label=label),size=5,parse=TRUE)+
  facet_grid(~alpha)+
  theme_bw()+
  theme(
    axis.title.y = element_text(size=15),
    axis.text.x=element_blank(),
    axis.text.y = element_text(colour="white"),
    axis.ticks.y = element_line(colour="white"),
    axis.title.x=element_blank(),
    panel.grid=element_blank(),
    axis.ticks.x=element_blank(),
    legend.position="none",
    strip.text=element_blank(),
    strip.background=element_blank(),
    plot.title=element_text(hjust=0.5,size=15)
  )+scale_color_viridis_d()+ggtitle("Parasite stage structure")

grid.arrange(structure.facet,alpha.cycle,nrow=2,heights=c(0.25,1))

#Plot developmental rate as a function of data-generating model parameters (Figure 5)
sens.best = bind_rows(sens.optim.full,sens.optim.red)
sens.best$BestModel = factor(sens.best$BestModel,levels=c("full","reduced"),labels=c("Full model \n(delayed maturation)","Reduced \nmodel"))

sens.best$Set = factor(sens.best$Set,levels=c(1,2,3),labels=c("1.0 / 50*'%'","0.5 / 50*'%'","0.25 / 80*'%'"))

p = ggplot()+
  geom_point(data=filter(sens.best,Parameter=="sigma[1]",Value>0.5,BestModel=="Full model \n(delayed maturation)"),aes(x=Value,y=24/lambdaA.fit),size=1.5)+
  geom_point(data=filter(sens.best,Parameter=="sigma[1]",Value>0.5,BestModel!="Full model \n(delayed maturation)"),aes(x=Value,y=24/lambdaA.fit),col="grey",size=1.5)+
  geom_point(data=filter(sens.best,Parameter!="sigma[1]",BestModel=="Full model \n(delayed maturation)"),aes(x=Value,y=24/lambdaA.fit),size=1.5)+
  geom_point(data=filter(sens.best,Parameter!="sigma[1]",BestModel!="Full model \n(delayed maturation)"),aes(x=Value,y=24/lambdaA.fit),col="grey",size=1.5)+
  facet_grid(Set~Parameter,labeller=label_parsed,scales="free_x")+
  theme_bw()+
  ylab("Estimated cycle length (hours)")+
  xlab("Parameter value")+
  scale_y_continuous(breaks=c(24,36,48),limits=c(22,50))+
  scale_alpha_manual(values=c(1,1))+
  scale_x_continuous(breaks=equal_breaks(n=2,s=0.05),expand=c(0.05,0))+
  theme(
    strip.background=element_blank(),
    strip.text = element_text(size=14),
    strip.text.x = element_text(color="#009E73"),
    strip.text.y = element_text(color="#56B4E9"),
    axis.title = element_blank(),
    axis.text = element_text(size=12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.title = element_text(size=13,hjust=0.5),
    legend.position = "bottom",
    legend.background = element_rect(colour="black"),
    #legend.title = element_text(size=13),
    legend.text = element_text(size=12,hjust=0),
    panel.spacing.x = unit(1.5, "lines"),
    panel.spacing.y = unit(1.5, "lines"),
    plot.title = element_text(size=13,hjust=0),
    legend.direction = "horizontal"
  )

title = ggdraw() + 
  draw_label(
    "Adjusted DDE model parameter",
    #fontface = 'bold',
    x = 0.5,
    hjust = 0.5,
    size=11
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 0),
    plot.background = element_rect(fill="#009E73",color=NA)
  )
title2 = ggdraw() + 
  draw_label(
    parse(text="Stage~misclassification~epsilon~'/'~'%'~Reticulocytes"),
    #fontface = 'bold',
    x = 0.5,
    hjust = 0.5,angle=-90,size=11
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(7, 0, 0, 0),
    plot.background = element_rect(fill="#56B4E9",color=NA)
  )

x.text = ggdraw() + 
  draw_label(
    "Parameter value",
    #fontface = 'bold',
    x = 0.5,
    hjust = 0.5,
    size=15
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 0),
    #plot.background = element_rect(fill="#009E73",color=NA)
  )
y.text = ggdraw() + 
  draw_label(
    parse(text="Estimated~cycle~length~(hours)"),
    #fontface = 'bold',
    x = 0.5,
    hjust = 0.5,angle=90,size=15
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(7, 0, 0, 0),
    #plot.background = element_rect(fill="#56B4E9",color=NA)
  )

plot_top = plot_grid(NULL,title,NULL,rel_widths=c(0.03,1,0),nrow=1)

plot_right = plot_grid(NULL,title2,NULL,rel_heights=c(0.08,1,0.08),nrow=3)

plot_stack = plot_grid(plot_top,p,x.text,rel_heights=c(0.08,1,0.04),rel_widths=c(1),nrow=3)

plot_row = plot_grid(y.text,plot_stack,plot_right,rel_widths=c(0.04,1,0.045),nrow=1)

plot_row

