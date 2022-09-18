# point to wherever your data files are
setwd("/Users/uricchio/projects/bootslabOld/anther_smut/manuscript/v2/software/simAnthSmut/") 
dat<- read.csv("data/Transect_1_all.csv", header=T)

library(boot)
library(ggplot2)
library(cowplot)
library(wesanderson)
#quad quality
head(dat,1) 
#  Section.07 Y.orig Xorig Y X X.2 Yplot Xplot X.3 Dp.05 HHp.05 MSp.05 Up.05 Di.05 HHi.05 MSi.05 Ui.05
#  Krum.05 Rock.05 X.4 Dp.07 HHp.07 MSp.07 Up.07 Di.07 HHi.07 MSi.07 Ui.07 Krum.07 Rock.07 Tp.07 X.5
#  X.6 Dp.14 HHp.14 MSp.14 Up.14 Di.14 HHi.14 MSi.14 Ui.14 Krum.14 Tp.14 X.7 Section.15 local.x.15
#  local.y.15 yplot.15 xplot.15 X.8 Krum.15 Dp.15 HHp.15 MSp.15 Up.15 Di.15 Hhi.15 MSi.15 Ui.15
#  Rock.15 X.9 X.10 X.11 Section.18 X.1 Y.1 Dp.18 HHp.18 MSp.18 Up.18 Di.18 HHi.18 MSi.18 Ui.18 X.12
#  Tp.18

 
 # creating a data frame for each sampled year from 2005 to 2018
 # Dp is the diseased plants, total plants is the sum of the Dp, HHp, etc
 # prevalence is just the mean 
 
df0<-data.frame(x=tapply(round(dat$Dp.05,0)+ round(dat$HHp.05,0)+round(dat$MSp.05,0)+round(dat$Up.05,0),round(dat$Dp.05,0)+ round(dat$HHp.05,0)+round(dat$MSp.05,0)+round(dat$Up.05,0),mean),y=tapply(round(dat$Dp.05,0)/(round(dat$Dp.05,0)+ round(dat$HHp.05,0)+round(dat$MSp.05,0)+round(dat$Up.05,0)),round(dat$Dp.05,0)+ round(dat$HHp.05,0)+round(dat$MSp.05,0)+round(dat$Up.05,0),mean),z=tapply(round(dat$Dp.05,0)+ round(dat$HHp.05,0)+round(dat$MSp.05,0)+round(dat$Up.05,0),round(dat$Dp.05,0)+ round(dat$HHp.05,0)+round(dat$MSp.05,0)+round(dat$Up.05,0),length))
df0<-cbind(df0,data.frame(year=rep(2005,length(df0$x))))

df1<-data.frame(x=tapply(dat$Dp.07+ dat$HHp.07+dat$MSp.07+dat$Up.07,dat$Dp.07+ dat$HHp.07+dat$MSp.07+dat$Up.07,mean),y=tapply(dat$Dp.07/(dat$Dp.07+dat$HHp.07+dat$MSp.07+dat$Up.07),dat$Dp.07+ dat$HHp.07+dat$MSp.07+dat$Up.07,mean),z=tapply(dat$Dp.07+ dat$HHp.07+dat$MSp.07+dat$Up.07,dat$Dp.07+ dat$HHp.07+dat$MSp.07+dat$Up.07,length))
df1<-cbind(df1,data.frame(year=rep(2007,length(df1$x))))

df2<-data.frame(x=tapply(dat$Dp.14+ dat$HHp.14+dat$MSp.14+dat$Up.14,dat$Dp.14+ dat$HHp.14+dat$MSp.14+dat$Up.14,mean),y=tapply(dat$Dp.14/(dat$Dp.14+dat$HHp.14+dat$MSp.14+dat$Up.14),dat$Dp.14+ dat$HHp.14+dat$MSp.14+dat$Up.14,mean),z=tapply(dat$Dp.14/(dat$Dp.14+dat$HHp.14+dat$MSp.14+dat$Up.14),dat$Dp.14+ dat$HHp.14+dat$MSp.14+dat$Up.14,length))
df2<-cbind(df2,data.frame(year=rep(2014,length(df2$x))))

df3<-data.frame(x=tapply(dat$Dp.15+ dat$HHp.15+dat$MSp.15+dat$Up.15,dat$Dp.15+ dat$HHp.15+dat$MSp.15+dat$Up.15,mean),y=tapply(dat$Dp.15/(dat$Dp.15+dat$HHp.15+dat$MSp.15+dat$Up.15),dat$Dp.15+ dat$HHp.15+dat$MSp.15+dat$Up.15,mean),z=tapply(dat$Dp.15/(dat$Dp.15+dat$HHp.15+dat$MSp.15+dat$Up.15),dat$Dp.15+ dat$HHp.15+dat$MSp.15+dat$Up.15,length))
df3<-cbind(df3,data.frame(year=rep(2015,length(df3$x))))

df4<-data.frame(x=tapply(dat$Dp.18+ dat$HHp.18+dat$MSp.18+dat$Up.18,dat$Dp.18+ dat$HHp.18+dat$MSp.18+dat$Up.18,mean),y=tapply(dat$Dp.18/(dat$Dp.18+dat$HHp.18+dat$MSp.18+dat$Up.18),dat$Dp.18+ dat$HHp.18+dat$MSp.18+dat$Up.18,mean),z=tapply(dat$Dp.18/(dat$Dp.18+dat$HHp.18+dat$MSp.18+dat$Up.18),dat$Dp.18+ dat$HHp.18+dat$MSp.18+dat$Up.18,length))
df4<-cbind(df4,data.frame(year=rep(2018,length(df4$x))))


# combining all the data
df<-rbind(df0,df1,df2,df3,df4)

# plotting the data 
plA<-ggplot(df,aes(x,y,color=as.factor(year),size=sqrt(z)))+geom_smooth(method='lm',se=FALSE,mapping = aes(weight =z))+geom_point(alpha=0.4)+scale_color_manual(values=wes_palette("Darjeeling2"),name="year")+scale_x_continuous(limits=c(1,50),name="density",breaks=c(1,5,10,15,20,25,30,35,40,45,50))+scale_y_continuous(limits=c(0,1),name="prevalence")+scale_size(guide='none')+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))

plA

summary(lm(df$y[df$year==2005]~df$x[df$year==2005],weights=df$z[df$year==2005]))
summary(lm(df$y[df$year==2007]~df$x[df$year==2007],weights=df$z[df$year==2007]))
summary(lm(df$y[df$year==2014]~df$x[df$year==2014],weights=df$z[df$year==2014]))
summary(lm(df$y[df$year==2015]~df$x[df$year==2015],weights=df$z[df$year==2015]))
summary(lm(df$y[df$year==2018]~df$x[df$year==2018],weights=df$z[df$year==2018]))


new_line<-lm(df$y[df$year==2014]~df$x[df$year==2014],weights=df$z[df$year==2014])
predict(new_line,data.frame(x=seq(1,27)))

ggsave("/Users/uricchio/projects/bootslabOld/anther_smut/manuscript/v2/figures/Fig4.pdf",height=5,width=7)



