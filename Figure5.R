sens.best = bind_rows(sens.optim.full,sens.optim.red)
sens.best$BestModel = factor(sens.best$BestModel,levels=c("full","reduced"),labels=c("Full model \n(delayed maturation)","Reduced \nmodel"))

sens.best$Set = factor(sens.best$Set,levels=c(5,2,3),labels=c("1.0 / 50*'%'","0.5 / 50*'%'","0.25 / 80*'%'"))

#### Developmental rate sensitivity ####
p = ggplot()+
  geom_point(data=filter(sens.best,Parameter=="sigma[1]",Value>0.5,BestModel=="Full model \n(delayed maturation)"),aes(x=Value,y=24/lA.fit),size=1.5)+
  geom_point(data=filter(sens.best,Parameter=="sigma[1]",Value>0.5,BestModel!="Full model \n(delayed maturation)"),aes(x=Value,y=24/lA.fit),col="grey",size=1.5)+
  #geom_smooth(data=filter(sens.best,Parameter=="sigma[1]",Value>0.5),aes(x=Value,y=24/lA.fit,col=factor(Set)),size=1.5,se=FALSE)+
  geom_point(data=filter(sens.best,Parameter!="sigma[1]",BestModel=="Full model \n(delayed maturation)"),aes(x=Value,y=24/lA.fit),size=1.5)+
  geom_point(data=filter(sens.best,Parameter!="sigma[1]",BestModel!="Full model \n(delayed maturation)"),aes(x=Value,y=24/lA.fit),col="grey",size=1.5)+
  #geom_smooth(data=filter(sens.best,Parameter!="sigma[2]",Parameter!="sigma[1]"),aes(x=Value,y=24/lA.fit,col=factor(Set)),se=FALSE,size=1.5)+
  #geom_vline(aes(xintercept=Value,col="Simulated \nvalue"),vline.df,size=2.1)+
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

ggsave(file="~/Desktop/Mideo.lab2/Parasite.maturation/TimeOpt/FacettedSensitivitylA.png", units="cm", height=12, width=28)


#### Replication rate sensitivity ####
q = ggplot()+
  geom_point(data=filter(sens.best,Parameter=="sigma[1]",Value>0.5,BestModel=="Full model \n(delayed maturation)"),aes(x=Value,y=BN.fit),size=1.5)+
  geom_point(data=filter(sens.best,Parameter=="sigma[1]",Value>0.5,BestModel!="Full model \n(delayed maturation)"),aes(x=Value,y=BN.fit),col="grey",size=1.5)+
  #geom_smooth(data=filter(sens.best,Parameter=="sigma[1]",Value>0.5),aes(x=Value,y=24/lA.fit,col=factor(Set)),size=1.5,se=FALSE)+
  geom_point(data=filter(sens.best,Parameter!="sigma[1]",Parameter!="kappa[N]",BestModel=="Full model \n(delayed maturation)"),aes(x=Value,y=BN.fit),size=1.5)+
  geom_point(data=filter(sens.best,Parameter!="sigma[1]",Parameter!="kappa[N]",BestModel!="Full model \n(delayed maturation)"),aes(x=Value,y=BN.fit),col="grey",size=1.5)+
  #geom_smooth(data=filter(sens.best,Parameter!="sigma[2]",Parameter!="sigma[1]"),aes(x=Value,y=24/lA.fit,col=factor(Set)),se=FALSE,size=1.5)+
  #geom_vline(aes(xintercept=Value,col="Simulated \nvalue"),vline.df,size=2.1)+
  facet_grid(Set~Parameter,labeller=label_parsed,scales="free_x")+
  theme_bw()+
  ylab(expression("Estimated replication rate "~k[N]~ ( "day" ^ -1)))+
  xlab("Parameter value")+
  #scale_y_continuous(breaks=c(0.04,0.06,0.08,0.1),limits=c(0.04,0.1))+
  scale_alpha_manual(values=c(1,1))+
  scale_x_continuous(breaks=equal_breaks(n=2,s=0.05),expand=c(0.05,0))+
  theme(
    strip.background=element_blank(),
    strip.text = element_text(size=14),
    #strip.text.x = element_text(color="#009E73"),
    #strip.text.y = element_text(color="#56B4E9"),
    axis.title = element_text(size=14),
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
q

r = ggplot()+
  geom_point(data=filter(sens.best,Parameter=="sigma[1]",Value>0.5,BestModel=="Full model \n(delayed maturation)"),aes(x=Value,y=BA.fit),size=1.5)+
  geom_point(data=filter(sens.best,Parameter=="sigma[1]",Value>0.5,BestModel!="Full model \n(delayed maturation)"),aes(x=Value,y=BA.fit),col="grey",size=1.5)+
  #geom_smooth(data=filter(sens.best,Parameter=="sigma[1]",Value>0.5),aes(x=Value,y=24/lA.fit,col=factor(Set)),size=1.5,se=FALSE)+
  geom_point(data=filter(sens.best,Parameter!="sigma[1]",Parameter!="kappa[A]",BestModel=="Full model \n(delayed maturation)"),aes(x=Value,y=BA.fit),size=1.5)+
  geom_point(data=filter(sens.best,Parameter!="sigma[1]",Parameter!="kappa[A]",BestModel!="Full model \n(delayed maturation)"),aes(x=Value,y=BA.fit),col="grey",size=1.5)+
  #geom_smooth(data=filter(sens.best,Parameter!="sigma[2]",Parameter!="sigma[1]"),aes(x=Value,y=24/lA.fit,col=factor(Set)),se=FALSE,size=1.5)+
  #geom_vline(aes(xintercept=Value,col="Simulated \nvalue"),vline.df,size=2.1)+
  facet_grid(Set~Parameter,labeller=label_parsed,scales="free_x")+
  theme_bw()+
  ylab(expression("Estimated replication rate "~k[A]~ ( "day" ^ -1)))+
  xlab("Parameter value")+
  #scale_y_continuous(breaks=c(0.04,0.06,0.08,0.1),limits=c(0.04,0.1))+
  scale_alpha_manual(values=c(1,1))+
  scale_x_continuous(breaks=equal_breaks(n=2,s=0.05),expand=c(0.05,0))+
  theme(
    strip.background=element_blank(),
    strip.text = element_text(size=14),
    #strip.text.x = element_text(color="#009E73"),
    #strip.text.y = element_text(color="#56B4E9"),
    axis.title = element_text(size=14),
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
r

plot_grid(q,r,nrow=2)



ggsave(file="~/Desktop/Mideo.lab2/Parasite.maturation/TimeOpt/FacettedSensitivityPro.png", units="cm", height=22, width=28)
 