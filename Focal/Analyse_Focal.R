#Load required packages
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(RColorBrewer))
suppressMessages(library(reghelper))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))

##### Data wrangling #####

#Read in fits
focal.full = read.csv("~/Focal/Focal_full.csv",stringsAsFactors=FALSE)
names(focal.full) = paste(names(focal.full),"full",sep=".")
focal.full$index = 1:nrow(focal.full)

focal.red = read.csv("~/Focal/Focal_reduced.csv",stringsAsFactors=FALSE) %>% select(.,-bias,-xc,-phi.1,-f)
names(focal.red) = paste(names(focal.red),"red",sep=".")
focal.red$index = 1:nrow(focal.red)

#Transform variables from log values and bind columns.
for (r in 1:nrow(focal.red)){
  
  row.vals = focal.red[r,]
  
  new.vals = c(sapply(row.vals[1:4],exp),
               row.vals[5],
               sapply(row.vals[6:7],exp),
               row.vals[8:11]) 
  
  focal.red[r,] = new.vals
}
for (r in 1:nrow(focal.full)){
  
  row.vals = focal.full[r,]
  
  new.vals = c(sapply(row.vals[1:4],exp),
               row.vals[5],
               sapply(row.vals[6:7],exp),
               row.vals[8:11]) 
  
  focal.full[r,] = new.vals
}
focal.joint = full_join(focal.full,focal.red,by=c("index"))

#Calcuate F-test for each row and store resulting p-value.
#Define function for F-test (from Greischar et al. 2016 supplementary code)
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
focal.joint$p = sapply(1:nrow(focal.joint),FUN=function(r){
  ss1 = focal.joint$fit.red[r]
  ss2 = focal.joint$fit.full[r]
  p = getF(32, 6,ss1,7,ss2)
  return(p)
})

#Add true parameter values
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
epsilon.list = seq(0,1,0.25)
phi.list = c(0,1,2,3)

#Create data frames for unique parameter combinations
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
para.list$epsilon = c(rep(c(sapply(1:length(epsilon.list),FUN=function(i){
  rep(epsilon.list[i],4)
})),4),rep(NA,24))
para.list$phi = phi.list

para.list$bias = c(rep("full",20),rep("stages",20),rep("full",20),rep("stages",20),rep("sequestration",8),rep("parasitemia",8),rep("none",8))
para.list$xc = c(rep(0,40),rep(0.5,40),rep(0,4),rep(0.5,4),rep(0,4),rep(0.5,4),rep(0,4),rep(0.5,4))
para.list$index = 1:nrow(para.list)

focal.joint =left_join(focal.joint,para.list,by=c("index"))

#Create data frame that contains best model for each row (i.e., full OR reduced parasite maturation model) based on p-value and correct bias labels. Split overall data frame into two: one for each clearance structure.
focal.optim.full =  filter(focal.joint,p<0.05) %>% 
  select(.,-contains("red"),-index)
focal.optim.red =  filter(focal.joint,p>=0.05) %>% 
  select(.,-contains("full"),-index)
names.list = c(paste(c("lambdaA","kN","kA","c","mu","sigma","P0"),"fit",sep="."),"Value","Convergence","Model","p",names(para.list)[1:8],"epsilon","phi","Bias","xc")
names(focal.optim.full) = names.list
names(focal.optim.red) = names.list

focal.best = bind_rows(focal.optim.full,focal.optim.red)

levels.list = c("full","sequestration","stages","parasitemia","none")
labels.list = c("All biases","Absolute abundance of circulating parasites",
                "Only stage misclassification bias","Relative abundance of all parasites",
                "Absolute abundance of all parasites")
focal.best$Bias = factor(focal.best$Bias,
                         levels=levels.list,
                         labels=labels.list)

focal.best$xc = factor(focal.best$xc,levels=c(0.0,0.5),labels=c("Indiscriminate clearance","Late-stage clearance"))

focal.best.0 = filter(focal.best,xc=="Indiscriminate clearance")
focal.best.0.5 = filter(focal.best,xc=="Late-stage clearance")


##### Plotting #####
#Create heatmaps to display estimated values of developmental rate as a function of extent misclassification (epsilon) and sequestration rate (phi) when all biases are applied OR only stage bias is applied.

## All parasites cleared
lambdaA.scale_fill = scale_fill_gradient2(low="#F0E442",mid="#009E73",midpoint=27,high="#CC79A7",name="Cycle \nlength",limits=c(24,30),breaks=c(24,27,30),labels=c("24","27","30"))
kN.scale_fill = scale_fill_gradient2(low="#F0E442",mid="#009E73",midpoint=0.08,high="#CC79A7",name=expression(k[N]),limits=c(0,0.16),breaks=c(0,0.08,0.16),labels=c("0","0.08","0.16"))
kA.scale_fill = scale_fill_gradient2(low="#F0E442",mid="#009E73",midpoint=0.25,high="#CC79A7",name=expression(k[A]),limits=c(0,0.50),breaks=c(0,0.25,0.50),labels=c("0","0.25","0.5"))
c.scale_fill = scale_fill_gradient2(low="#F0E442",mid="#009E73",midpoint=0.45,high="#CC79A7",name=expression(italic(c)),limits=c(0,0.9),breaks=c(0,0.45,0.9),labels=c("0","0.45","0.90"))

theme.nolegend = theme(
  strip.text = element_blank(),
  strip.background = element_blank(),
  axis.text = element_text(size=12),
  axis.title = element_blank(),
  legend.position = "none"
)
theme.legend =   theme(
  strip.text = element_blank(),
  strip.background = element_blank(),
  axis.text = element_text(size=12),
  axis.title = element_blank(),
  legend.text = element_text(size=12),
  legend.title = element_text(size=15)
)

### Facet plot for all biases
plot.lambdaA.all.0 = ggplot(data=filter(focal.best.0,Bias%in%c("All biases")),
                       aes(x=epsilon,y=phi,fill=24/lambdaA.fit))+
  geom_tile()+facet_grid(.~Bias,scales="free")+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+lambdaA.scale_fill+theme.nolegend

plot.kN.all.0 = ggplot(data=filter(focal.best.0,Bias%in%c("All biases")),
                       aes(x=epsilon,y=phi,fill=kN.fit))+
  geom_tile()+facet_grid(.~Bias,scales="free")+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+kN.scale_fill+theme.nolegend

plot.kA.all.0 = ggplot(data=filter(focal.best.0,Bias%in%c("All biases")),
                       aes(x=epsilon,y=phi,fill=kA.fit))+
  geom_tile()+facet_grid(.~Bias,scales="free")+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+kA.scale_fill+theme.nolegend

plot.cN.all.0 = ggplot(data=filter(focal.best.0,Bias%in%c("All biases")),
                       aes(x=epsilon,y=phi,fill=c.fit))+
  geom_tile()+facet_grid(.~Bias,scales="free")+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+c.scale_fill+theme.nolegend

### Facet plot for stage bias only
plot.lambdaA.stages.0 = ggplot(data=filter(focal.best.0,Bias%in%c("Stage misclassification")),
                          aes(x=epsilon,y=phi,fill=24/lambdaA.fit))+
  geom_tile()+facet_grid(.~Bias,scales="free")+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+
  lambdaA.scale_fill+theme.legend

plot.kN.stages.0 = ggplot(data=filter(focal.best.0,Bias%in%c("Stage misclassification")),
                          aes(x=epsilon,y=phi,fill=kN.fit))+
  geom_tile()+facet_grid(.~Bias,scales="free")+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+
  kN.scale_fill+theme.legend

plot.kA.stages.0 = ggplot(data=filter(focal.best.0,Bias%in%c("Stage misclassification")),
                          aes(x=epsilon,y=phi,fill=kA.fit))+
  geom_tile()+facet_grid(.~Bias,scales="free")+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+
  kA.scale_fill+theme.legend

plot.cN.stages.0 = ggplot(data=filter(focal.best.0,Bias%in%c("Stage misclassification")),
                          aes(x=epsilon,y=phi,fill=c.fit))+
  geom_tile()+facet_grid(.~Bias,scales="free")+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+
  c.scale_fill+theme.legend

## Only late-stage parasites cleared
### Facet plot for all biases
plot.lambdaA.all.0.5 = ggplot(data=filter(focal.best.0.5,Bias%in%c("All biases")),aes(x=epsilon,y=phi,fill=24/lambdaA.fit))+
  geom_tile()+facet_grid(.~Bias,scales="free")+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100")) + lambdaA.scale_fill + theme.nolegend

plot.kN.all.0.5 = ggplot(data=filter(focal.best.0.5,Bias%in%c("All biases")),
                         aes(x=epsilon,y=phi,fill=kN.fit))+
  geom_tile()+facet_grid(.~Bias,scales="free")+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100")) + kN.scale_fill + theme.nolegend

plot.kA.all.0.5 = ggplot(data=filter(focal.best.0.5,Bias%in%c("All biases")),
                         aes(x=epsilon,y=phi,fill=kA.fit))+
  geom_tile()+facet_grid(.~Bias,scales="free")+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100")) + kA.scale_fill + theme.nolegend

plot.cN.all.0.5 = ggplot(data=filter(focal.best.0.5,Bias%in%c("All biases")),
                         aes(x=epsilon,y=phi,fill=c.fit))+
  geom_tile()+facet_grid(.~Bias,scales="free")+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100")) + c.scale_fill + theme.nolegend

### Facet plot for stage bias only

plot.lambdaA.stages.0.5 = ggplot(data=filter(focal.best.0.5,Bias%in%c("Stage misclassification")),
                            aes(x=epsilon,y=phi,fill=24/lambdaA.fit))+
  geom_tile()+facet_grid(.~Bias,scales="free")+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100")) + lambdaA.scale_fill + theme.nolegend

plot.kN.stages.0.5 = ggplot(data=filter(focal.best.0.5,Bias%in%c("Stage misclassification")),
                            aes(x=epsilon,y=phi,fill=kN.fit))+
  geom_tile()+facet_grid(.~Bias,scales="free")+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100")) + kN.scale_fill + theme.nolegend

plot.kA.stages.0.5 = ggplot(data=filter(focal.best.0.5,Bias%in%c("Stage misclassification")),
                            aes(x=epsilon,y=phi,fill=kA.fit))+
  geom_tile()+facet_grid(.~Bias,scales="free")+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100")) + kA.scale_fill + theme.nolegend

plot.cN.stages.0.5 = ggplot(data=filter(focal.best.0.5,Bias%in%c("Stage misclassification")),
                            aes(x=epsilon,y=phi,fill=c.fit))+
  geom_tile()+facet_grid(.~Bias,scales="free")+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100")) + c.scale_fill + theme.nolegend

## Plot for main text showing estimating clearance and replication rates as a function of stage misclassification and sequestration rate (6A-C).
theme.nolegend = theme(
  strip.text = element_blank(),
  strip.background = element_blank(),
  axis.text = element_text(size=12),
  legend.position = "none"
)
theme.legend =   theme(
  strip.text = element_blank(),
  strip.background = element_blank(),
  axis.text = element_text(size=10),
  legend.text = element_text(size=12),
  legend.title = element_text(size=15),
  axis.title = element_blank(),
  legend.position="top"
)

plot.kN.main = ggplot(data=filter(focal.best.0.5,Bias%in%c("All biases")),
                      aes(x=epsilon,y=phi,fill=kN.fit))+
  geom_tile()+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+kN.scale_fill+theme.legend

plot.kA.main = ggplot(data=filter(focal.best.0.5,Bias%in%c("All biases")),
                      aes(x=epsilon,y=phi,fill=kA.fit))+
  geom_tile()+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+kA.scale_fill+theme.legend

plot.cN.main = ggplot(data=filter(focal.best.0.5,Bias%in%c("All biases")),
                      aes(x=epsilon,y=phi,fill=c.fit))+
  geom_tile()+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+c.scale_fill+theme.legend


plot_facet = plot_grid(plot.cN.main,plot.kN.main,plot.kA.main,nrow=1)

x.text = ggdraw() + 
  draw_label(
    parse(text="Stage~misclassification~epsilon~('%')"),
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
    expression("Sequestration rate"~phi),
    #fontface = 'bold',
    x = 0.5,
    hjust = 0.5,angle=90,size=15
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7),
    #plot.background = element_rect(fill="#56B4E9",color=NA)
  )

plot_stack = plot_grid(plot_facet,x.text,rel_heights=c(1,0.06),nrow=2)

plot_row = plot_grid(y.text,plot_stack,rel_widths=c(0.06,1),nrow=1)
plot_row

top_row = plot_row + draw_plot_label(label = c("A", "B", "C"), size = 15,
                                     x = c(0.075,0.39,0.7), y = c(0.95,0.95,0.95))

## Create heatmap figure showing estimated developmental rate as a function of stage misclassification and sequestration rate.
lambdaA.scale_fill = scale_fill_gradient2(low="#F0E442",mid="#009E73",midpoint=31,high="#CC79A7",name="Cycle \nlength \n(hours)",limits=c(24,38),breaks=c(24,32,38),labels=c("24","32","38"))

### Facet plot for all biases
epsilon.phi.lambdaA = ggplot()+
  geom_tile(data=filter(focal.best,Bias%in%c("All biases","Stage misclassification")),aes(x=epsilon,y=phi,fill=24/lambdaA.fit))+facet_grid(.~Bias,scales="free")+
  geom_text(data=filter(focal.best,Bias%in%c("All biases","Stage misclassification"),lambdaA.fit!=1),aes(x=epsilon,y=phi,label=round(24/lambdaA.fit,0)),size=5)+
  theme_classic()+facet_grid(xc~Bias)+
  scale_size_continuous(guide=FALSE)+
  xlab(expression("Stage misclassification"~epsilon~"(%)"))+
  ylab(expression("Sequestration rate"~phi))+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+lambdaA.scale_fill+
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size=12),
    axis.title = element_text(size=18),
    legend.position = "none"
  )

## Facet plot for supplemental materials

### Plot for supplemental showing estimating clearance and replication rates as a function of stage misclassification and sequestration rate.
kN.scale_fill = scale_fill_gradient2(low="#F0E442",mid="#009E73",midpoint=0.08,high="#CC79A7",name=expression(k[N]),limits=c(0,0.16),breaks=c(0,0.08,0.16),labels=c("0","0.08","0.16"))
kA.scale_fill = scale_fill_gradient2(low="#F0E442",mid="#009E73",midpoint=0.25,high="#CC79A7",name=expression(k[A]),limits=c(0,0.5),breaks=c(0,0.25,0.5),labels=c("0","0.25","0.5"))
c.scale_fill = scale_fill_gradient2(low="#F0E442",mid="#009E73",midpoint=0.5,high="#CC79A7",name=expression(italic(c)),limits=c(0,1),breaks=c(0,0.5,1),labels=c("0","0.5","1"))

theme.nolegend = theme(
  strip.text = element_blank(),
  strip.background = element_blank(),
  axis.text = element_text(size=12),
  legend.position = "none"
)
theme.legend =   theme(
  strip.text = element_blank(),
  strip.background = element_blank(),
  axis.text = element_text(size=10),
  legend.text = element_text(size=12),
  legend.title = element_text(size=15),
  axis.title = element_blank(),
  legend.position="top"
)

plot.kN.main = ggplot(data=filter(focal.best.0.5,Bias%in%c("Only stage misclassification bias")),
                      aes(x=epsilon,y=phi,fill=kN.fit))+
  geom_tile()+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+kN.scale_fill+theme.legend

plot.kA.main = ggplot(data=filter(focal.best.0.5,Bias%in%c("Only stage misclassification bias")),
                      aes(x=epsilon,y=phi,fill=kA.fit))+
  geom_tile()+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+kA.scale_fill+theme.legend

plot.cN.main = ggplot(data=filter(focal.best.0.5,Bias%in%c("Only stage misclassification bias")),
                      aes(x=epsilon,y=phi,fill=c.fit))+
  geom_tile()+
  theme_classic()+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+c.scale_fill+theme.legend

plot_facet = plot_grid(plot.cN.main,plot.kN.main,plot.kA.main,nrow=1)

x.text = ggdraw() + 
  draw_label(
    parse(text="Stage~misclassification~epsilon~('%')"),
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
    expression("Sequestration rate"~phi),
    #fontface = 'bold',
    x = 0.5,
    hjust = 0.5,angle=90,size=15
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7),
    #plot.background = element_rect(fill="#56B4E9",color=NA)
  )


plot_stack = plot_grid(plot_facet,x.text,rel_heights=c(1,0.06),nrow=2)

plot_row = plot_grid(y.text,plot_stack,rel_widths=c(0.06,1),nrow=1)

top_row_stage = plot_row + draw_plot_label(label = c("A", "B", "C"), size = 15,
                                           x = c(0.075,0.39,0.7), y = c(0.95,0.95,0.95))

#### Statistics ####
#Use step-wise model selection to determine best models (lowest AIC scores) for predicting developmental rate, invasion ratios and clearance rate. The data sets are too small to split into training and testing data, so Beta coefficients and p-values are not calculated.

model.sel.df = data.frame()
model.sel.single = as.data.frame(matrix(nrow=4,ncol=14))
names(model.sel.single) = c("Clearance","Bias","Estimate","Predictor","B","p","SE","Beta","t","model.p","R2","F","df1","df2")

for (v in 1:4){
  
  variable.fit = c("lambdaA.fit","kN.fit","kA.fit","c.fit")[v]
  variable = c("lambdaA","kN","kA","c")[v]
  
  if (v == 1){b.list = c(1,2)} else {b.list = 1:5}
  
  for (b in b.list){
    
    if (b %in% c(1,2)) {
      f.full = as.formula(paste(variable.fit,paste(c("epsilon","phi"), collapse = "*"),sep = " ~ "))
      
      filter.name = unique(focal.best.0$Bias)[b]
      df.name = c("All","Stages")[b]
      
      ##### Standard clearance
      model.sel.single = as.data.frame(matrix(nrow=4,ncol=14))
      model.sel.single[] = matrix(sapply(1:nrow(model.sel.single),function(r){c("Standard",df.name,variable,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)}),nrow=4,byrow=TRUE)
      names(model.sel.single) = c("Clearance","Bias","Estimate","Predictor","B","p","SE","Beta","t","model.p","R2","F","df1","df2")
      model.sel.single[,4] = c("Intercept","epsilon","phi","epsilon:phi")
      
      model = lm(f.full, data=filter(focal.best.0,Bias==filter.name))
      
      model.sel.single[1:4,5] = summary(model)$coefficients[1:4]
      model.sel.single[1:4,6] = summary(model)$coefficients[13:16]
      model.sel.single[1:4,7] = summary(model)$coefficients[5:8]
      
      model.sel.single[1:4,9] = beta(model)$coefficients[9:12]
      model.sel.single[1:4,8] = beta(model)$coefficients[1:4]
      
      model.sel.single[1:4,10] = c(1 -  pf(summary(model)$fstatistic[1],summary(model)$fstatistic[2],summary(model)$fstatistic[3]),NA,NA,NA)
      model.sel.single[1:4,11] = c(summary(model)$adj.r.squared,NA,NA,NA)
      model.sel.single[1:4,12] = c(summary(model)$fstatistic[1],NA,NA,NA)
      model.sel.single[1:4,13] = c(summary(model)$fstatistic[2],NA,NA,NA)
      model.sel.single[1:4,14] = c(summary(model)$fstatistic[3],NA,NA,NA)
      
      model.sel.df = rbind(model.sel.df,model.sel.single)
      
      ##### Alternate clearance
      model.sel.single = as.data.frame(matrix(nrow=4,ncol=14))
      model.sel.single[] = matrix(sapply(1:nrow(model.sel.single),function(r){c("Alternate",df.name,variable,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)}),nrow=4,byrow=TRUE)
      names(model.sel.single) = c("Clearance","Bias","Estimate","Predictor","B","p","SE","Beta","t","model.p","R2","F","df1","df2")
      model.sel.single[,4] = c("Intercept","epsilon","phi","epsilon:phi")
      
      model = lm(f.full, data=filter(focal.best.0.5,Bias==filter.name))
      
      model.sel.single[1:4,5] = summary(model)$coefficients[1:4]
      model.sel.single[1:4,6] = summary(model)$coefficients[13:16]
      model.sel.single[1:4,7] = summary(model)$coefficients[5:8]
      
      model.sel.single[1:4,9] = beta(model)$coefficients[9:12]
      model.sel.single[1:4,8] = beta(model)$coefficients[1:4]
      
      model.sel.single[1:4,10] = c(1 -  pf(summary(model)$fstatistic[1],summary(model)$fstatistic[2],summary(model)$fstatistic[3]),NA,NA,NA)
      model.sel.single[1:4,11] = c(summary(model)$adj.r.squared,NA,NA,NA)
      model.sel.single[1:4,12] = c(summary(model)$fstatistic[1],NA,NA,NA)
      model.sel.single[1:4,13] = c(summary(model)$fstatistic[2],NA,NA,NA)
      model.sel.single[1:4,14] = c(summary(model)$fstatistic[3],NA,NA,NA)
      
      model.sel.df = rbind(model.sel.df,model.sel.single)
      
      model.sel.single = as.data.frame(matrix(nrow=4,ncol=14))
      names(model.sel.single) = c("Clearance","Bias","Estimate","Predictor","B","p","SE","Beta","t","model.p","R2","F","df1","df2")
      
    } else {
      
      f.full = as.formula(paste(variable.fit,paste(c("phi"), collapse = "*"),sep = " ~ "))
      
      filter.name = unique(focal.best.0$Bias)[b]
      df.name = c("All","Stages","Sequestration","Parasitemia","None")[b]
      
      ##### Standard clearance
      model.sel.single = as.data.frame(matrix(nrow=2,ncol=14))
      model.sel.single[] = matrix(sapply(1:nrow(model.sel.single),function(r){c("Standard",df.name,variable,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)}),nrow=2,byrow=TRUE)
      names(model.sel.single) = c("Clearance","Bias","Estimate","Predictor","B","p","SE","Beta","t","model.p","R2","F","df1","df2")
      model.sel.single[,4] = c("Intercept","phi")
      
      model = lm(f.full, data=filter(focal.best.0,Bias==filter.name))
      
      model.sel.single[1:2,5] = summary(model)$coefficients[1:2]
      model.sel.single[1:2,6] = summary(model)$coefficients[7:8]
      model.sel.single[1:2,7] = summary(model)$coefficients[3:4]
      
      model.sel.single[1:2,9] = beta(model)$coefficients[5:6]
      model.sel.single[1:2,8] = beta(model)$coefficients[1:2]
      
      model.sel.single[1:2,10] = c(1 -  pf(summary(model)$fstatistic[1],summary(model)$fstatistic[2],summary(model)$fstatistic[3]),NA)
      model.sel.single[1:2,11] = c(summary(model)$adj.r.squared,NA)
      model.sel.single[1:2,12] = c(summary(model)$fstatistic[1],NA)
      model.sel.single[1:2,13] = c(summary(model)$fstatistic[2],NA)
      model.sel.single[1:2,14] = c(summary(model)$fstatistic[3],NA)
      
      model.sel.df = rbind(model.sel.df,model.sel.single)
      
      ##### Alternate clearance
      model.sel.single = as.data.frame(matrix(nrow=2,ncol=14))
      model.sel.single[] = matrix(sapply(1:nrow(model.sel.single),function(r){c("Alternate",df.name,variable,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)}),nrow=2,byrow=TRUE)
      names(model.sel.single) = c("Clearance","Bias","Estimate","Predictor","B","p","SE","Beta","t","model.p","R2","F","df1","df2")
      model.sel.single[,4] = c("Intercept","phi")
      
      model = lm(f.full, data=filter(focal.best.0.5,Bias==filter.name))
      
      model.sel.single[1:2,5] = summary(model)$coefficients[1:2]
      model.sel.single[1:2,6] = summary(model)$coefficients[7:8]
      model.sel.single[1:2,7] = summary(model)$coefficients[3:4]
      
      model.sel.single[1:2,9] = beta(model)$coefficients[5:6]
      model.sel.single[1:2,8] = beta(model)$coefficients[1:2]
      
      model.sel.single[1:2,10] = c(1 -  pf(summary(model)$fstatistic[1],summary(model)$fstatistic[2],summary(model)$fstatistic[3]),NA)
      model.sel.single[1:2,11] = c(summary(model)$adj.r.squared,NA)
      model.sel.single[1:2,12] = c(summary(model)$fstatistic[1],NA)
      model.sel.single[1:2,13] = c(summary(model)$fstatistic[2],NA)
      model.sel.single[1:2,14] = c(summary(model)$fstatistic[3],NA)
      
      model.sel.df = rbind(model.sel.df,model.sel.single)
      
      names(model.sel.single) = c("Clearance","Bias","Estimate","Predictor","B","p","SE","Beta","t","model.p","R2","F","df1","df2")
    }
    
    
    
  } #end loop over biases
} #end loop over variables
