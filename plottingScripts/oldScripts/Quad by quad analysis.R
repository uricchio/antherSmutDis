setwd("~/Dropbox/Transect Simulation/data") # point to wherever you're data files are
dat<- read.csv("Transect 1 all.csv", header=T)
library(boot)
#quad quality
head(dat,1) 
#  Section.07 Y.orig Xorig Y X X.2 Yplot Xplot X.3 Dp.05 HHp.05 MSp.05 Up.05 Di.05 HHi.05 MSi.05 Ui.05
#  Krum.05 Rock.05 X.4 Dp.07 HHp.07 MSp.07 Up.07 Di.07 HHi.07 MSi.07 Ui.07 Krum.07 Rock.07 Tp.07 X.5
#  X.6 Dp.14 HHp.14 MSp.14 Up.14 Di.14 HHi.14 MSi.14 Ui.14 Krum.14 Tp.14 X.7 Section.15 local.x.15
#  local.y.15 yplot.15 xplot.15 X.8 Krum.15 Dp.15 HHp.15 MSp.15 Up.15 Di.15 Hhi.15 MSi.15 Ui.15
#  Rock.15 X.9 X.10 X.11 Section.18 X.1 Y.1 Dp.18 HHp.18 MSp.18 Up.18 Di.18 HHi.18 MSi.18 Ui.18 X.12
#  Tp.18
krum.ave=rowMeans(cbind(dat$Krum.07,dat$Krum.14,dat$Krum.15),na.rm=TRUE)
 rock=ifelse(is.na(dat$Rock.07),0,dat$Rock.07)
 krumrock1=rowSums(cbind(krum.ave,rock)) 
 krumrock=ifelse(krumrock1>100,100,krumrock1)
 qual=(exp(2.153 -5.397e-03*krumrock -3.931e-04*krumrock^2))/exp(2.153)
 
 Hold <- matrix(NA, 15, 100, byrow=TRUE) 
 Hold[cbind( dat$Y.plot, dat$X.plot)] = qual
#Results untransformed data, all significant at P<0.001
#intercept 2005 inv.logit(-1.9680254)= 0.1226011
#intercept 2007 inv.logit(-1.0582503)= 0.257644
#intercept 2014 inv.logit(-1.099743)= 0.2497881
#intercept 2015 inv.logit(-0.810303)= 0.3078259
#intercept 2018 inv.logit(-1.72262 )= 0.151534
##################################################################################### 
#untransformed data
# D=dat$Dp.07;HH=dat$HHp.07; MS=dat$MSp.07; year=2007
#win.graph()
plot(c(0,50),c(0,1),type="n",xlab="density",ylab="prevalence", cex.lab=1.5)
 
DENSITY.DEP.1= function( qual, D=dat$Dp.07, HH=dat$HHp.07, MS=dat$MSp.07, U=dat$Up.07, year)
 {
 if (year==2005){color="dark gray"}
 if (year==2007){color="brown"}	
 if (year==2014){color="blue"}	
 if (year==2015){color="dark green"}	
 if (year==2018){color="red"}	
 	 
 Tot= D+HH+MS+U
 P=D/Tot
 P[!is.finite(P)]=NA

 mT=tapply(Tot,Tot,mean)
 nT=tapply(Tot,Tot,length)
 mP=tapply(P,Tot,mean)
 points(mT,mP,cex=sqrt(nT)/2,col=color)
 frame.mT=as.data.frame(mT)
 mod1=glm(P~Tot+I(Tot^2),weights=Tot,family=binomial) 
 print(summary(mod1))
 pred1=predict(mod1,list(Tot=mT),type="response")
 lines(mT,pred1,col=color,lwd=2)
 }
DENSITY.DEP.1( qual, D=round(dat$Dp.05,0), HH=round(dat$HHp.05,0), 
																					MS=round(dat$MSp.05,0), U=round(dat$Up.05,0),2005)
DENSITY.DEP.1( qual, D=dat$Dp.07, HH=dat$HHp.07, MS=dat$MSp.07, U=dat$Up.07,2007)
DENSITY.DEP.1( qual, D=dat$Dp.14, HH=dat$HHp.14, MS=dat$MSp.14, U=dat$Up.14,2014)
DENSITY.DEP.1( qual, D=dat$Dp.15, HH=dat$HHp.15, MS=dat$MSp.15, U=dat$Up.15,2015)
DENSITY.DEP.1( qual, D=dat$Dp.18, HH=dat$HHp.18, MS=dat$MSp.18, U=dat$Up.18,2018)

#################################################################################
# plot as sqrt on density axis
# D=dat$Dp.07;HH=dat$HHp.07; MS=dat$MSp.07; year=2007
#win.graph()
plot(c(0,7),c(0,1),type="n",xlab="sqrt density",ylab="prevalence",cex.lab=1.5)
 
DENSITY.DEP.2= function(qual, D=dat$Dp.07, HH=dat$HHp.07, MS=dat$MSp.07, U=dat$Up.07,year)
 {
 if (year==2005){color="dark gray"}
 if (year==2007){color="brown"}	
 if (year==2014){color="blue"}	
 if (year==2015){color="dark green"}	
 if (year==2018){color="red"}	
 	 
 Tot= D+HH+MS+U
 H=HH+MS
 P= D/Tot
 P[!is.finite(P)]=NA

 mT=tapply(Tot,Tot,mean)
 nT=tapply(Tot,Tot,length)
 mP=tapply(P,Tot,mean)

 points(sqrt(mT),mP,cex=sqrt(nT)/2,col=color)
 frame.mT=as.data.frame(mT)
 mod1=glm(P~Tot+I(Tot^1.5),weights=Tot,family=binomial) 
 summary(mod1)
 pred1=predict(mod1,list(Tot=mT),type="response")
 lines(sqrt(mT),pred1,col=color,lwd=2)
}
DENSITY.DEP.2( qual, D=round(dat$Dp.05,0), HH=round(dat$HHp.05,0), 
																					MS=round(dat$MSp.05,0), U=round(dat$Up.05,0),2005)
DENSITY.DEP.2( qual, D=dat$Dp.07, HH=dat$HHp.07, MS=dat$MSp.07,U=dat$Up.07,2007)
DENSITY.DEP.2( qual, D=dat$Dp.14, HH=dat$HHp.14, MS=dat$MSp.14,U=dat$Up.14,2014)
DENSITY.DEP.2( qual, D=dat$Dp.15, HH=dat$HHp.15, MS=dat$MSp.15,U=dat$Up.15,2015)
DENSITY.DEP.2( qual, D=dat$Dp.18, HH=dat$HHp.18, MS=dat$MSp.18,U=dat$Up.18,2018)
############################################################################
#log scale
#win.graph()
plot(c(0,4),c(0,1),type="n",xlab="log density",ylab="prevalence",cex.lab=1.5)
 
 DENSITY.DEP.3= function( qual, D=dat$Dp.07, HH=dat$HHp.07, MS=dat$MSp.07, U=dat$Up.07, year)
 {
 if (year==2005){color="dark gray"}
 if (year==2007){color="brown"}	
 if (year==2014){color="blue"}	
 if (year==2015){color="dark green"}	
 if (year==2018){color="red"}	
 	 
 Tot= D+HH+MS+U
 H=HH+MS
 P= D/Tot
 P[!is.finite(P)]=NA

 mT=tapply(Tot,Tot,mean)
 nT=tapply(Tot,Tot,length)
 mP=tapply(P,Tot,mean)
 
 #win.graph()
 #plot(c(0,7),c(0,1),type="n",xlab="sqrt density",ylab="prevalence",main=year)
 points(log(mT),mP,cex=sqrt(nT)/2,col=color)
 frame.mT=as.data.frame(mT)
 mod1=glm(P~Tot+I(Tot^2),weights=Tot,family=binomial) #  +I(tot07^3)+qual+I(qual^2)
 summary(mod1)
 pred1=predict(mod1,list(Tot=mT),type="response")
 lines(log(mT),pred1,col=color,lwd=2)
 #lines(mT,pred1,col=color)
 }
DENSITY.DEP.3( qual, D=round(dat$Dp.05,0), HH=round(dat$HHp.05,0), 
																					MS=round(dat$MSp.05,0), U=round(dat$Up.05,0),2005)
DENSITY.DEP.3( qual, D=dat$Dp.07, HH=dat$HHp.07, MS=dat$MSp.07,U=dat$Up.07,2007)
DENSITY.DEP.3( qual, D=dat$Dp.14, HH=dat$HHp.14, MS=dat$MSp.14,U=dat$Up.14,2014)
DENSITY.DEP.3( qual, D=dat$Dp.15, HH=dat$HHp.15, MS=dat$MSp.15,U=dat$Up.15,2015)
DENSITY.DEP.3( qual, D=dat$Dp.18, HH=dat$HHp.18, MS=dat$MSp.18,U=dat$Up.18,2018)
################################################################################
################################################################################
#plot vs. means at low densities (<10)
#win.graph() 
plot(c(0,3.5),c(0,1),type="n",xlab="sqrt density",ylab="mean prevalence",cex.lab=1.5)
D=dat$Dp.07; HH=dat$HHp.07; MS=dat$MSp.07; year=2007 
 
DENSITY.DEP.4= function(qual, D=dat$Dp.07, HH=dat$HHp.07, MS=dat$MSp.07,U=dat$Up.07, year)
 {
 if (year==2005){color="dark gray"}
 if (year==2007){color="brown"}	
 if (year==2014){color="blue"}	
 if (year==2015){color="dark green"}	
 if (year==2018){color="red"}	
 	 
 Tot= D+HH+MS+U
 H=HH+MS
 P= D/Tot
 P[!is.finite(P)]=NA

 mT=tapply(Tot,Tot,mean)
 nT=tapply(Tot,Tot,length)
 mP=tapply(P,Tot,mean)
 
 mT=mT[1:11]
 nT=nT[1:11]
 mP=mP[1:11]
 
 #win.graph()
 #plot(c(0,7),c(0,1),type="n",xlab="sqrt density",ylab="prevalence",main=year)
 points(sqrt(mT),mP,cex=sqrt(nT/2),col=color)
 #points(mT,mP,cex=sqrt(nT)/2,col=color)
 #frame.mT=as.data.frame(mT)
 sr.mT=sqrt(mT)
 mod1=glm(mP~sr.mT, weights=nT, family="gaussian") #+I(sr.mT^2   
 print( summary(mod1))
 pred1=predict(mod1,list(sr.mT=sr.mT),type="response")
 lines(sr.mT,pred1,col=color,lwd=2)
 }
DENSITY.DEP.4( qual, D=round(dat$Dp.05,0), HH=round(dat$HHp.05,0), 
																					MS=round(dat$MSp.05,0), U=round(dat$Up.05,0),2005)
DENSITY.DEP.4( qual, D=dat$Dp.07, HH=dat$HHp.07, MS=dat$MSp.07,U=dat$Up.07,2007)
DENSITY.DEP.4( qual, D=dat$Dp.14, HH=dat$HHp.14, MS=dat$MSp.14,U=dat$Up.14,2014)
DENSITY.DEP.4( qual, D=dat$Dp.15, HH=dat$HHp.15, MS=dat$MSp.15,U=dat$Up.15,2015)
DENSITY.DEP.4( qual, D=dat$Dp.18, HH=dat$HHp.18, MS=dat$MSp.18,U=dat$Up.18,2018)

#note - non linear term only significant for 2014, P=0.023

#Results data on means
#intercept 2005 = 0.033 P=0.687
#intercept 2007 = 0.064 P=0.591
#intercept 2014 = 0.118 P=0.0516
#intercept 2015 = 0.149 P= 0.052
#intercept 2018 = 0.059 P= 0.347