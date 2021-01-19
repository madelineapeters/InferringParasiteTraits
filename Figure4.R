library(tidyverse)
library(gridExtra)
library(cowplot)

#### Developmental rate figure (epsilon vs. phi, epsilon vs. percent retics) ####
lambdaA.scale_fill = scale_fill_gradient2(low="#F0E442",mid="#009E73",midpoint=37,high="#CC79A7",name="Cycle \nlength",limits=c(24,50),breaks=c(24,37,50),labels=c("24","37","50"))

### Facet plot for all biases
epsilon.phi.lambdaA = ggplot()+
  geom_tile(data=filter(focal.best,Bias%in%c("All biases","Only stage misclassification bias"),xc=="Late-stage clearance"),aes(x=epsilon,y=phi,fill=24/lambdaA.fit))+facet_wrap(Bias~.,scales="free",nrow=2)+
  geom_text(data=filter(focal.best,Bias%in%c("All biases","Only stage misclassification bias",xc=="Late-stage clerance"),lambdaA.fit!=1),aes(x=epsilon,y=phi,label=round(24/lambdaA.fit,0)),size=5)+
  theme_classic()+
  ylab(expression("Sequestration rate"~phi))+
  xlab(expression("Stage misclassification"~epsilon~"(%)"))+
  scale_size_continuous(guide=FALSE)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+lambdaA.scale_fill+
  theme(
    strip.text = element_text(size=12),
    strip.background = element_blank(),
    axis.text = element_text(size=10),
    axis.title = element_text(size=13),
    legend.position = "none"
  )

epsilon.phi.lambdaA

### Facet plot for stage misclassification with misclassified transfused parasites
inoc.best$sens = FALSE
inoc.best$sens[which(inoc.best$epsilon==1&inoc.best$index==7)]=TRUE
inoc.best$sens[which(inoc.best$epsilon==0.5&inoc.best$index==7)]=TRUE
inoc.best$sens[which(inoc.best$epsilon==0.25&inoc.best$index==10)]=TRUE

epsilon.frac.lambdaA = ggplot()+
  geom_tile(data=filter(inoc.best,index<11),aes(x=epsilon,y=index,fill=24/lambdaA.fit))+
  geom_text(data=filter(inoc.best,lambdaA.fit!=1,index<11),aes(x=epsilon,y=index,label=round(24/lambdaA.fit,0)),size=5)+
  theme_classic()+
  ylab("% reticulocytes in \ntransfused parasite cohort")+
  xlab(expression("Stage misclassification"~epsilon~"(%)"))+
  #scale_size_continuous(guide=FALSE)+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10),labels=c("0","3","10","20","30","40","50","60","70","80"))+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+lambdaA.scale_fill+
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size=10),
    legend.text = element_text(size=12),
    axis.title = element_text(size=13),
    legend.title = element_text(size=15)
  )+
  geom_tile(data=filter(inoc.best,index<11), aes(x=epsilon,y=index,colour=factor(sens, c(TRUE, FALSE)), size=factor(sens, c(TRUE, FALSE))), alpha=0) + 
  scale_colour_manual("z", values=c("black", NA)) + 
  scale_size_manual("z", values=c(1, 0))+guides(colour=FALSE,size = FALSE)

epsilon.frac.lambdaA


hlay.main = rbind(c(1,2))

my.plots = list(epsilon.phi.lambdaA,epsilon.frac.lambdaA)
arrangedEstlA = plot_grid(plotlist=my.plots, labels=c("A", "B"),label_size=20,ncol = 2, nrow = 1, rel_widths = c(0.6,0.9))
h = arrangeGrob(arrangedEstlA)

ggsave(file="~/Desktop/Mideo.lab2/Parasite.maturation/TimeOpt/EstimatinglA.pdf", h, units="cm", height=16, width=24)
ggsave(file="~/Desktop/Mideo.lab2/Parasite.maturation/TimeOpt/EstimatinglA.png", h, units="cm", height=16, width=24)


#### Minimum sum of square figure ####
MSS.scale_fill = scale_fill_gradient2(low="#F0E442",mid="#009E73",midpoint=1.75,high="#CC79A7",name=expression(Delta~"MSS"),limits=c(-0.005,3.5),breaks=c(0,1.75,3.5),labels=c("0","1.75","3.5"))

epsilon.frac.fit = ggplot()+
  geom_tile(data=filter(inoc.best,index<11),aes(x=epsilon,y=index,fill=MSS))+
  geom_text(data=filter(inoc.best,index<11),aes(x=epsilon,y=index,label=round(Value,2)),size=4)+
  theme_classic()+
  ylab("% reticulocytes in \ntransfused parasite cohort")+
  xlab(expression("Stage misclassification"~epsilon~"(%)"))+
  #scale_size_continuous(guide=FALSE)+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10),labels=c("0","3","10","20","30","40","50","60","70","80"))+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+MSS.scale_fill+
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size=10),
    legend.text = element_text(size=12),
    axis.title = element_text(size=13),
    legend.title = element_text(size=15)
  )

epsilon.frac.fit

ggsave(file="~/Desktop/Mideo.lab2/Parasite.maturation/TimeOpt/EpsilonFracFit.png", units="cm", height=13, width=12)


#### Proliferation parameter estimate facet plot ####

c.scale_fill = scale_fill_gradient2(low="#F0E442",mid="#009E73",midpoint=0.5,high="#CC79A7",name=expression(italic(c)),limits=c(0,1),breaks=c(0,0.5,1),labels=c("0","0.5","1"))

epsilon.frac.c = ggplot()+
  geom_tile(data=filter(inoc.best,index<11),aes(x=epsilon,y=index,fill=c.fit))+
  #geom_text(data=filter(inoc.best,c.fit>0.01,index<11),aes(x=epsilon,y=index,label=round(c.fit,2)),size=4)+
  theme_classic()+
  ylab("% reticulocytes in \ntransfused parasite cohort")+
  xlab(expression("Stage misclassification"~epsilon~"(%)"))+
  #scale_size_continuous(guide=FALSE)+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10),labels=c("0","3","10","20","30","40","50","60","70","80"))+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+c.scale_fill+
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size=10),
    legend.text = element_text(size=12),
    axis.title = element_blank(),
    legend.title = element_text(size=15),
    legend.position="top"
  )

epsilon.frac.c

kA.scale_fill = scale_fill_gradient2(low="#F0E442",mid="#009E73",midpoint=0.25,high="#CC79A7",name=expression(k[A]),limits=c(0,0.5),breaks=c(0,0.25,0.5),labels=c("0","0.25","0.5"))

epsilon.frac.kA = ggplot()+
  geom_tile(data=filter(inoc.best,index<11),aes(x=epsilon,y=index,fill=kA.fit))+
  #geom_text(data=filter(inoc.best,index<11),aes(x=epsilon,y=index,label=round(kA.fit,2)),size=4)+
  theme_classic()+
  ylab("% reticulocytes in \ntransfused parasite cohort")+
  xlab(expression("Stage misclassification"~epsilon~"(%)"))+
  #scale_size_continuous(guide=FALSE)+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10),labels=c("0","3","10","20","30","40","50","60","70","80"))+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+kA.scale_fill+
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size=10),
    legend.text = element_text(size=12),
    axis.title = element_blank(),
    legend.title = element_text(size=15),
    legend.position="top"
  )

epsilon.frac.kA

kN.scale_fill = scale_fill_gradient2(low="#F0E442",mid="#009E73",midpoint=0.08,high="#CC79A7",name=expression(k[N]),limits=c(0,0.16),breaks=c(0,0.08,0.16),labels=c("0","0.08","0.16"))

epsilon.frac.kN = ggplot()+
  geom_tile(data=filter(inoc.best,index<11),aes(x=epsilon,y=index,fill=kN.fit))+
  theme_classic()+
  ylab("% reticulocytes in \ntransfused parasite cohort")+
  xlab(expression("Stage misclassification"~epsilon~"(%)"))+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10),labels=c("0","3","10","20","30","40","50","60","70","80"))+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","50","100"))+kN.scale_fill+
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size=10),
    legend.text = element_text(size=12),
    axis.title = element_blank(),
    legend.title = element_text(size=15),
    legend.position="top"
  )

epsilon.frac.kN

plot_facet = plot_grid(epsilon.frac.c,epsilon.frac.kN,epsilon.frac.kA,nrow=1)

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
    expression("     % reticulocytes in\ntransfused parasite cohort"),
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

y.text


plot_stack = plot_grid(plot_facet,x.text,rel_heights=c(1,0.06),nrow=2)

plot_row = plot_grid(y.text,plot_stack,rel_widths=c(0.06,1),nrow=1)
plot_row

bottom_row = plot_row + draw_plot_label(label = c("D", "E", "F"), size = 15,
                           x = c(0.09,0.40,0.715), y = c(0.95,0.95,0.95))

plot_grid(top_row,bottom_row,nrow=2)
