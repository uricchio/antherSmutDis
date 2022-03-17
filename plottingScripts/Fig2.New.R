#install.packages("mgcv") # only need to do this once, if you haven't installed the package
library("mgcv") 
library(ggplot2)
library(cowplot)
library(wesanderson)

# point to wherever your data files are
setwd("/Users/uricchio/projects/bootslabOld/anther_smut/manuscript/v2/software/simAnthSmut/") 
dat<- read.csv("data/Transect_1_all.csv", header=T)

attach(dat)


##DENSITY#####################################
dens.05=round(Dp.05+HHp.05+MSp.05+Up.05,0)
dens.07=Dp.07+HHp.07+MSp.07+Up.07
dens.14=Dp.14+HHp.14+MSp.14+Up.14
dens.15=Dp.15+HHp.15+MSp.15+Up.15 
dens.18=Dp.18+HHp.18+MSp.18+Up.18

dens.data=as.matrix(cbind(Xplot, Yplot, dens.05, dens.07, dens.14, dens.15, dens.18))
dens.data.05=data.frame(cbind(x=Xplot, y=Yplot, z=dens.05))
dens.data.15=data.frame(cbind(x=Xplot, y=Yplot, z=dens.15))
dens.data.14=data.frame(cbind(x=Xplot,y=Yplot, z=dens.14))
dens.data.07=data.frame(cbind(x=Xplot, y=Yplot, z=dens.07))
dens.data.18=data.frame(cbind(x=Xplot,y=Yplot, z=dens.18))
dis.data.05=data.frame(cbind(x=Xplot, y=Yplot, z=Dp.05))
dis.data.07=data.frame(cbind(x=Xplot, y=Yplot, z=Dp.07))
dis.data.18=data.frame(cbind(x=Xplot, y=Yplot, z=Dp.18))
dis.data.14=data.frame(cbind(x=Xplot, y=Yplot, z=Dp.14))
dis.data.15=data.frame(cbind(x=Xplot, y=Yplot, z=Dp.15))

# code block for printing out summary files of Fh,Fd, etc
# not needed to make the figure, I'm retaining this (ugly) code since I used some of these internally
#for (i in seq(1,length(dens.data.05$x))) {
#    write (c(dens.data.05$x[i],dens.data.05$y[i],dens.data.05$z[i]),file="/Users/uricchio/projects/bootslab/anther_smut/simAnthSmut/data/dens.05.txt",append=TRUE)
#}
#
#for (i in seq(1,length(dis.data.05$x))) {
#    write (c(dis.data.05$x[i],dis.data.05$y[i],dis.data.05$z[i]), file="/Users/uricchio/projects/bootslab/anther_smut/simAnthSmut/data/dis.05.txt", append=TRUE)
#}
#
#for (i in seq(1,length(dens.data.07$x))) {
#    write (c(dens.data.07$x[i],dens.data.07$y[i],dens.data.07$z[i]), file="/Users/uricchio/projects/bootslab/anther_smut/simAnthSmut/data/dens.07.txt",append=TRUE)
#}
#for (i in seq(1,length(dis.data.07$x))) {
#    write (c(dis.data.07$x[i],dis.data.07$y[i],dis.data.07$z[i]), file="/Users/uricchio/projects/bootslab/anther_smut/simAnthSmut/data/dis.07.txt",append=TRUE)
#}
#
#for (i in seq(1,length(dens.data.18$x))) {
#    write (c(dens.data.18$x[i],dens.data.18$y[i],dens.data.18$z[i]), file="/Users/uricchio/projects/bootslab/anther_smut/simAnthSmut/data/dens.18.txt",append=TRUE)
#}
#for (i in seq(1,length(dis.data.18$x))) {
#    write (c(dis.data.18$x[i],dis.data.18$y[i],dis.data.18$z[i]),file="/Users/uricchio/projects/bootslab/anther_smut/simAnthSmut/data/dis.18.txt",append=TRUE)
#}
#
#for (i in seq(1,length(dens.data.14$x))) {
#    write (c(dens.data.14$x[i],dens.data.14$y[i],dens.data.14$z[i]), file="/Users/uricchio/projects/bootslab/anther_smut/simAnthSmut/data/dens.14.txt",append=TRUE)
#}
#for (i in seq(1,length(dis.data.14$x))) {
#    write (c(dis.data.14$x[i],dis.data.14$y[i],dis.data.14$z[i]),file="/Users/uricchio/projects/bootslab/anther_smut/simAnthSmut/data/dis.14.txt",append=TRUE)
#}
#for (i in seq(1,length(dens.data.15$x))) {
#    write (c(dens.data.15$x[i],dens.data.15$y[i],dens.data.15$z[i]), file="/Users/uricchio/projects/bootslab/anther_smut/simAnthSmut/data/dens.15.txt",append=TRUE)
#}
#for (i in seq(1,length(dis.data.15$x))) {
#    write (c(dis.data.15$x[i],dis.data.15$y[i],dis.data.15$z[i]),file="/Users/uricchio/projects/bootslab/anther_smut/simAnthSmut/data/dis.15.txt",append=TRUE)
#}
#for (i in seq(1,length(dens.data.15$x))) {
#    print (c(dens.data.15$x[i],dens.data.15$y[i],dens.data.15$z[i]))
#}
#for (i in seq(1,length(dis.data.15$x))) {
#    print (c(dis.data.15$x[i],dis.data.15$y[i],dis.data.15$z[i]))
#}
#for (i in seq(1,length(dens.data.14$x))) {
#    print (c(dens.data.14$x[i],dens.data.14$y[i],dens.data.14$z[i]))
#}
#for (i in seq(1,length(dis.data.14$x))) {
#    print (c(dis.data.14$x[i],dis.data.14$y[i],dis.data.14$z[i]))
#}

#quit()


dens.transect=aggregate(dens.data[,c(3:7)], by=list(Xplot), mean, na.rm=T) 

df<-data.frame(x=c(1:100),y=dens.transect[,2],year=rep("2005",100))
df<-rbind(df,data.frame(x=c(1:100),y=dens.transect[,3],year=rep("2007",100)))
df<-rbind(df,data.frame(x=c(1:100),y=dens.transect[,4],year=rep("2014",100)))
df<-rbind(df,data.frame(x=c(1:100),y=dens.transect[,5],year=rep("2015",100)))
df<-rbind(df,data.frame(x=c(1:100),y=dens.transect[,6],year=rep("2018",100)))

plA<-ggplot(df,aes(x,y,col=year))+geom_point(alpha=0.3)+scale_color_manual(values=wes_palette("Darjeeling2"))

Set.05=as.matrix(cbind(dens.05,Xplot))
Set.05=Set.05[!is.na(Set.05[,1]),]
mod.05 <-gam(Set.05[,1]~s(Set.05[,2]),family=poisson,gamma=0.25)
predgam.05=predict.gam(mod.05);
dfFit<-data.frame(x=unique(Set.05[,2]),y=c(t(t(unique(exp(predgam.05))))),year=rep("2005",length(unique(Set.05[,2]))))

Set.07=as.matrix(cbind(dens.07,Xplot))
Set.07=Set.07[!is.na(Set.07[,1]),]
mod.07 <-gam(Set.07[,1]~s(Set.07[,2]),family=poisson,gamma=0.25)  #unused model parameters:  ,bs="cr",k=20
predgam.07=predict.gam(mod.07); #lines(Set.07[,2], exp(predgam.07),col="black",lwd=2)
dfFit<-rbind(dfFit,data.frame(x=unique(Set.07[,2]),y=c(t(t(unique(exp(predgam.07))))),year=rep("2007",length(unique(Set.07[,2])))))
	
mod.14 <-gam(dens.14~s(Xplot),family=poisson,gamma=0.25) 
predgam.14=predict.gam(mod.14);#lines(Xplot, exp(predgam.14),col="dark green",lwd=2)
dfFit<-rbind(dfFit,data.frame(x=unique(Xplot),y=c(t(t(unique(exp(predgam.14))))),year=rep("2014",100)))

mod.15 <-gam(dens.15~s(Xplot),family=poisson,gamma=0.25)
predgam.15=predict.gam(mod.15);#lines(Xplot, exp(predgam.15),col="brown",lwd=2)
dfFit<-rbind(dfFit,data.frame(x=unique(Xplot),y=c(t(t(unique(exp(predgam.15))))),year=rep("2015",100)))

mod.18 <-gam(dens.18~s(Xplot),family=poisson,gamma=0.25)
predgam.18=predict.gam(mod.18);#lines(Xplot, exp(predgam.18),col="red",lwd=2)
dfFit<-rbind(dfFit,data.frame(x=unique(Xplot),y=c(t(t(unique(exp(predgam.18))))),year=rep("2018",100)))

plA<-plA+geom_line(data=dfFit,aes(x,y,col=year),size=0.8)+scale_y_continuous(limits=c(0,15))+xlab("position (m)")+ylab(expression("Density per " * m^2))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))


#KRUMMHOLZ########################### 
krumr.05=100*(Krum.05+Rock.05)
krumr.05=ifelse(krumr.05>100,100,krumr.05)
Rock.07 =ifelse(is.na(Rock.07),0,Rock.07)
krumr.07=Krum.07+Rock.07
krumr.07=ifelse(krumr.07>100,100,krumr.07)
krumr.14=Krum.14+Rock.07
krumr.14=ifelse(krumr.14>100,100,krumr.14)
krumr.15=Krum.15+Rock.07 
krumr.15=ifelse(krumr.15>100,100,krumr.15)


krum.data=as.matrix(cbind(Xplot, Yplot, krumr.05, krumr.07, krumr.14, krumr.15))
krumr.transect=aggregate(krum.data[,c(3:6)], by=list(Xplot), mean, na.rm=T)


df<-data.frame(x=c(1:100),y=krumr.transect[,2],year=rep("2005",100))
df<-rbind(df,data.frame(x=c(1:100),y=krumr.transect[,3],year=rep("2007",100)))
df<-rbind(df,data.frame(x=c(1:100),y=krumr.transect[,4],year=rep("2014",100)))
df<-rbind(df,data.frame(x=c(1:100),y=krumr.transect[,5],year=rep("2015",100)))
# no data for 2018

plB<-ggplot(df,aes(x,y,col=year))+geom_point(alpha=0.3)+scale_color_manual(values=wes_palette("Darjeeling2"))

knots=30
gampar=0.25

Set.05=as.matrix(cbind(krumr.05,Xplot))
Set.05=Set.05[!is.na(Set.05[,1]),]
mod.05 <-gam(Set.05[,1]~s(Set.05[,2],k=knots/2),family=gaussian, gamma=gampar)
predgam.05=predict.gam(mod.05);
dfFit<-data.frame(x=unique(Set.05[,2]),y=c(t(t(unique(predgam.05)))),year=rep("2005",length(unique(Set.05[,2]))))

Set.07=as.matrix(cbind(krumr.07,Xplot))
Set.07=Set.07[!is.na(Set.07[,1]),]
mod.07 <-gam(Set.07[,1]~s(Set.07[,2],k=knots),family=gaussian,gamma=gampar)  #unused model parameters:  ,bs="cr",k=20
predgam.07=predict.gam(mod.07); #lines(Set.07[,2], predgam.07,col="black",lwd=2)
dfFit<-rbind(dfFit,data.frame(x=unique(Set.07[,2]),y=c(t(t(unique(predgam.07)))),year=rep("2007",length(unique(Set.07[,2])))))


Set.14=as.matrix(cbind(krumr.14,Xplot))
Set.14=Set.14[!is.na(Set.14[,1]),]	
mod.14 <-gam(Set.14[,1]~s(Set.14[,2],k=knots),family=gaussian, gamma=gampar)
predgam.14=predict.gam(mod.14);#lines(Set.14[,2], predgam.14,col="dark green",lwd=2)
dfFit<-rbind(dfFit,data.frame(x=unique(Set.14[,2]),y=c(t(t(unique(predgam.14)))),year=rep("2014",length(unique(Set.14[,2])))))

mod.15 <-gam(krumr.15~s(Xplot,k=knots),family=gaussian, gamma=gampar)
predgam.15=predict.gam(mod.15);#lines(Xplot, predgam.15,col="brown",lwd=2)
dfFit<-rbind(dfFit,data.frame(x=unique(Xplot),y=c(t(t(unique(predgam.15)))),year=rep("2015",length(unique(Xplot)))))

plB<-plB+geom_line(data=dfFit,aes(x,y,col=year),size=0.8)+scale_y_continuous(limits=c(0,105))+xlab("position (m)")+ylab(expression("% cover Krummholz"))+theme(legend.position="NA")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))

#PREVALENCE##########################################################

dis.data=as.matrix(cbind(Xplot, Yplot, round(Dp.05,0), Dp.07, Dp.14, Dp.15, Dp.18))
dis.transect=aggregate(dis.data[,c(3:7)], by=list(Xplot), mean, na.rm=T) 
dens.data=as.matrix(cbind(Xplot, Yplot,dens.05, dens.07, dens.14, dens.15, dens.18))
 
prev.transect=dis.transect/dens.transect



df<-data.frame(x=c(1:100),y=prev.transect[,2],year=rep("2005",100))
df<-rbind(df,data.frame(x=c(1:100),y=prev.transect[,3],year=rep("2007",100)))
df<-rbind(df,data.frame(x=c(1:100),y=prev.transect[,4],year=rep("2014",100)))
df<-rbind(df,data.frame(x=c(1:100),y=prev.transect[,5],year=rep("2015",100)))
df<-rbind(df,data.frame(x=c(1:100),y=prev.transect[,6],year=rep("2018",100)))

plC<-ggplot(df,aes(x,y,col=year))+geom_point(alpha=0.3)+scale_color_manual(values=wes_palette("Darjeeling2"))

dis.data=as.matrix(cbind(Xplot, Yplot, round(Dp.05,0), Dp.07, Dp.14, Dp.15, Dp.18))
dens.data=as.matrix(cbind(Xplot, Yplot,dens.05, dens.07, dens.14, dens.15, dens.18))
prev.data=dis.data/dens.data



Pset.05=cbind(Xplot,prev.data[,3],dens.data[,3])
Pset.05=Pset.05[complete.cases(Pset.05),] 
mod.05 <-gam(Pset.05[,2]~s(Pset.05[,1]),family=binomial(link="logit"),weights=Pset.05[,3],gamma=0.25) 
predgam.05=predict.gam(mod.05);
back=exp(predgam.05)/(1+exp(predgam.05))
dfFit<-data.frame(x=unique(Pset.05[,1]),y=c(t(t(unique(back)))),year=rep("2005",length(unique(Pset.05[,1]))))

Pset.07=cbind(Xplot,prev.data[,4],dens.data[,4])
Pset.07=Pset.07[complete.cases(Pset.07),] 
mod.07 <-gam(Pset.07[,2]~s(Pset.07[,1]),family=binomial(link="logit"),weights=Pset.07[,3],gamma=0.25) 
predgam.07=predict.gam(mod.07);
back=exp(predgam.07)/(1+exp(predgam.07))
dfFit<-rbind(dfFit,data.frame(x=unique(Pset.07[,1]),y=c(t(t(unique(back)))),year=rep("2007",length(unique(Pset.07[,1])))))


Pset.14=cbind(Xplot,prev.data[,5],dens.data[,5])
Pset.14=Pset.14[complete.cases(Pset.14),] 
mod.14 <-gam(Pset.14[,2]~s(Pset.14[,1]),family=binomial(link="logit"),weights=Pset.14[,3],gamma=0.25) 
predgam.14=predict.gam(mod.14);
back=exp(predgam.14)/(1+exp(predgam.14))
dfFit<-rbind(dfFit,data.frame(x=unique(Pset.14[,1]),y=c(t(t(unique(back)))),year=rep("2014",length(unique(Pset.14[,1])))))

Pset.15=cbind(Xplot,prev.data[,6],dens.data[,6])
Pset.15=Pset.15[complete.cases(Pset.15),] 
mod.15 <-gam(Pset.15[,2]~s(Pset.15[,1]),family=binomial(link="logit"),weights=Pset.15[,3],gamma=0.25) 
predgam.15=predict.gam(mod.15);
back=exp(predgam.15)/(1+exp(predgam.15))
dfFit<-rbind(dfFit,data.frame(x=unique(Pset.15[,1]),y=c(t(t(unique(back)))),year=rep("2015",length(unique(Pset.15[,1])))))

Pset.18=cbind(Xplot,prev.data[,7],dens.data[,7])
Pset.18=Pset.18[complete.cases(Pset.18),] 
mod.18 <-gam(Pset.18[,2]~s(Pset.18[,1]),family=binomial(link="logit"),weights=Pset.18[,3],gamma=0.25) 
predgam.18=predict.gam(mod.18);
back=exp(predgam.18)/(1+exp(predgam.18))
dfFit<-rbind(dfFit,data.frame(x=unique(Pset.18[,1]),y=c(t(t(unique(back)))),year=rep("2018",length(unique(Pset.18[,1])))))

plC<-plC+geom_line(data=dfFit,aes(x,y,col=year),size=0.8)+xlab("position (m)")+ylab(expression("disease prevalence"))+theme(legend.position="NA")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))

df <- data.frame()
plBl<-ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)+theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())

row2<-plot_grid(plB,plBl,rel_widths=c(1,0.18),labels=c("B",""))
row3<-plot_grid(plC,plBl,rel_widths=c(1,0.18),labels=c("C",""))

plot_grid(plA,row2,row3,labels=c("A","",""),ncol=1)


ggsave("/Users/uricchio/projects/bootslabOld/anther_smut/manuscript/v2/figures/Fig2.pdf",width=7,height=11)
#ggsave("~/Documents/CV/app/2019/boise/ConceptFigPlotsForApp.pdf",width=7,height=11*2/3)



 
 
