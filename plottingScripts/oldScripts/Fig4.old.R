# point to wherever your data files are
setwd("/Users/uricchio/projects/bootslabOld/anther_smut/manuscript/v2/software/simAnthSmut/") 
dat<- read.csv("data/Transect_1_all.csv", header=T)
head(dat)

# this script is just a bunch of manipulation of the associated data table to get the percent krummholz cover and
# total number of plants for each year in the dataset. We don't have % Krummholz data for 2018 so that year is not 
# included.


library(ggplot2)
library(cowplot)
library(wesanderson)

# this is just a list of the names of the variables in the file
# Section.07 Y.orig Xorig Y X X.2 Yplot Xplot X.3 Dp.05 HHp.05 MSp.05 Up.05 Di.05 HHi.05 MSi.05 Ui.05 Krum.05 Rock.05  
#  X.4 Dp.07 HHp.07 MSp.07 Up.07 Di.07 HHi.07 MSi.07 Ui.07 Krum.07 Rock.07 Tp.07 X.5 X.6 Dp.14 HHp.14 MSp.14 Up.14
#  Di.14 HHi.14 MSi.14 Ui.14 Krum.14 Tp.14 X.7 Section.15 local.x.15 local.y.15 yplot.15 xplot.15 X.8 Krum.15 Dp.15
#  HHp.15 MSp.15 Up.15 Di.15 Hhi.15 MSi..15 Ui.15 Rock.15 X.9 X.10 X.11 Section.18 X.1 Y.1 Dp.18 HHp.18 MSp.18 Up.18
#  Di.18 HHi.18 MSi.18 Ui.18 X.12 Tp.18


attach(dat)
#####################################
col05="dark gray";
col07="brown";
col14="blue";
col15="dark green";
col18="red";


# get the 05 data. Tp is total plants
dat$Tp.05=Dp.05 +HHp.05 +MSp.05 +Up.05

krumclass=(trunc(dat$Krum.05*10))*10 # note that data in 05 were scaled to 1, hence the difference in the operation here vs the subsequent years which divide by 10.
krumave=tapply(dat$Tp.05,krumclass,mean)
krumcount=tapply(dat$Tp.05,krumclass,length)

# xpoints are the distances along the linear transect in meters
xpoints=c(0,10,20,30,40,50,60,70,80,90,100)


# model of plant abundance as function of krummholz
mod1=lm(krumave~xpoints+I(xpoints^2),weight=krumcount)
pred1=predict(mod1,type="response")
lines(xpoints,pred1,col=col05,lwd=2)
summary(mod1)

df<-data.frame(x=xpoints,y=krumave,z=sqrt(krumcount/10),w=krumcount,year=rep("2005",11))

# now same for 07
krumclass=(trunc(dat$Krum.07/10))*10
krumave=tapply(dat$Tp.07,krumclass,mean)
krumcount=tapply(dat$Tp.07,krumclass,length)

xpoints=c(0,10,20,30,40, 50, 60,70,80,90,100)

mod1=lm(krumave~xpoints+I(xpoints^2),weight=krumcount)
pred1=predict(mod1,type="response")
lines(xpoints,pred1,col=col07,lwd=2)
summary(mod1)

df<-rbind(df,data.frame(x=xpoints,y=krumave,z=sqrt(krumcount/10),w=krumcount,year=rep("2007",11)))


# now same for 14
dat$Tp.14=dat$Dp.14+dat$HHp.14 +dat$MSp.14 +dat$Up.14

krumclass=(trunc(dat$Krum.14/10))*10
krumave=tapply(dat$Tp.14,krumclass,mean)
krumcount=tapply(dat$Tp.14,krumclass,length)

mod1=lm(krumave~xpoints+I(xpoints^2),weight=krumcount)
pred1=predict(mod1,type="response")
lines(xpoints,pred1,col=col14,lwd=2)
summary(mod1)

df<-rbind(df,data.frame(x=xpoints,y=krumave,z=sqrt(krumcount/10),w=krumcount,year=rep("2014",11)))

# now same for 15
dat$Tp.15=dat$Dp.15+dat$HHp.15+dat$MSp.15+dat$Up.15 
krumclass=(trunc(dat$Krum.15/10))*10
krumave=tapply(dat$Tp.15,krumclass,mean)
krumcount=tapply(dat$Tp.15,krumclass,length)

mod1=lm(krumave~xpoints+I(xpoints^2),weight=krumcount)
pred1=predict(mod1,type="response")
lines(xpoints,pred1,col=col15,lwd=2)
summary(mod1)

df<-rbind(df,data.frame(x=xpoints,y=krumave,z=sqrt(krumcount/10),w=krumcount,year=rep("2015",11)))

 #now plot all the data
ggplot(df,aes(x,y,size=z,col=year,weight=w))+geom_point()+stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 0.6,alpha=0.05)+scale_color_manual(values=wes_palette("Darjeeling2"))+ylim(c(0,17))+xlab("Krummholz % cover")+ylab(expression("plants per " * m^2))+scale_size_continuous(guide=F)+theme(legend.position="NA")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))


ggsave("/Users/uricchio/projects/bootslabOld/anther_smut/manuscript/v2/software/simAnthSmut/figures/Fig3New.pdf",height=4.5,width=8)


