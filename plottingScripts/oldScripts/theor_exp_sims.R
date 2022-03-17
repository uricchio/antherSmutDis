library(ggplot2)
library(cowplot)
library(wesanderson)


setwd("/Users/uricchio/projects/bootslabOld/anther_smut/manuscript/v2/software/simAnthSmut/")

read.table("TheorExpSims/All.freq.1.dens.1.txt")->tfd      # both freq dependent and dens trans
read.table("TheorExpSims/All.freq.1.dens.0.txt")->tf   # freq dependent
read.table("TheorExpSims/All.freq.0.dens.1.txt")->td   # dens dependent 

# FD

meansFD<-c()
stdFD<-c()
meansFDPrev<-c()
stdFDPrev<-c()

indThresh<-0

for (i in seq(1,100)) {
	meansFD<-c(meansFD,mean(tfd[,i]+tfd[,i+100],na.rm=T))
	stdFD<-c(stdFD,sd(tfd[,i]+tfd[,i+100],na.rm=T))
}
plot(meansFD)


for (i in seq(101,200)) {
	meansFDPrev<-c(meansFDPrev,mean(tfd[,i]/(tfd[,i-100]+tfd[,i]),na.rm=T))
	stdFDPrev<-c(stdFDPrev,sd(tfd[,i]/(tfd[,i-100]+tfd[,i]),na.rm=T))
}
plot(meansFD,meansFDPrev)

# Dens

meansD<-c()
stdD<-c()
meansDPrev<-c()
stdDPrev<-c()

indThresh<-0

for (i in seq(1,100)) {
	meansD<-c(meansD,mean(td[,i]+td[,i+100],na.rm=T))
	stdD<-c(stdD,sd(td[,i]+td[,i+100],na.rm=T))
}
plot(meansD)


for (i in seq(101,200)) {
	meansDPrev<-c(meansDPrev,mean(td[,i]/(td[,i-100]+td[,i]),na.rm=T))
	stdDPrev<-c(stdDPrev,sd(td[,i]/(td[,i-100]+td[,i]),na.rm=T))
}
plot(meansD,meansDPrev)

# Freq

meansF<-c()
stdF<-c()
meansFPrev<-c()
stdFPrev<-c()



for (i in seq(1,100)) {
	meansF<-c(meansF,mean(tf[,i]+tf[,i+100],na.rm=T))
	stdF<-c(stdF,sd(tf[,i]+tf[,i+100],na.rm=T))
}
plot(meansD)


for (i in seq(101,200)) {
	meansFPrev<-c(meansFPrev,mean(tf[,i]/(tf[,i-100]+tf[,i]),na.rm=T))
	stdFPrev<-c(stdFPrev,sd(tf[,i]/(tf[,i-100]+tf[,i]),na.rm=T))
}
plot(meansF,meansFPrev)




mydf<-data.frame(pos=c(seq(1,100),seq(1,100),seq(1,100)),abun=c(meansFD,meansD,meansF),prev=c(meansFDPrev,meansDPrev,meansFPrev),maxAbun<-c(meansFD+stdFD,meansD+stdD,meansF+stdF),minAbun<-c(pmax(meansFD-stdFD,0),pmax(meansD-stdD,0),pmax(meansF-stdF,0)),minPrev<-c(pmax(meansFDPrev-stdFDPrev,0),pmax(meansDPrev-stdDPrev,0),pmax(meansFPrev-stdFPrev,0)),maxPrev<-c(meansFDPrev+stdFDPrev,meansDPrev+stdDPrev,meansFPrev+stdFPrev),type=c(rep("both",100),rep("density-dependent",100),rep("frequency-dependent",100)))

plA<-ggplot()+geom_ribbon(data=mydf,aes(x=pos,ymin=minAbun,ymax=maxAbun,fill=type),alpha=0.3)+scale_color_manual(values=wes_palette("Darjeeling1"))+xlab("position")+ylab("Abundance")+geom_line(data=mydf,aes(pos,abun,col=type),size=1)+theme(legend.position="none")+scale_fill_manual(values=wes_palette("Darjeeling1"))

plB<-ggplot()+geom_ribbon(data=mydf,aes(x=pos,ymin=minPrev,ymax=maxPrev,fill=type),alpha=0.5)+scale_color_manual(values=wes_palette("Darjeeling1"))+xlab("position")+ylab("Prevalence")+geom_line(data=mydf,aes(pos,prev,col=type),size=1)+theme(legend.position="none")+scale_fill_manual(values=wes_palette("Darjeeling1"))+ylim(c(0,0.5))

plC<-ggplot()+geom_ribbon(data=mydf,aes(x=abun,ymin=minPrev,ymax=maxPrev,fill=type),alpha=0.5)+scale_color_manual(values=wes_palette("Darjeeling1"),name="")+xlab("abundance")+ylab("prevalence")+geom_line(data=mydf,aes(abun,prev,col=type),size=1)+scale_fill_manual(values=wes_palette("Darjeeling1"),name="")+xlim(c(0,7))+ylim(c(0,0.5))


row1<-plot_grid(plA,plB,plC,ncol=3,labels=c("A","B","C"),rel_widths=c(1,1,1.6))

# Now with real qual


read.table("TheorExpSims/All.qual.freq.1.dens.1.txt")->tfd      # both freq dependent and dens trans
read.table("TheorExpSims/All.qual.freq.1.dens.0.txt")->tf   # freq dependent
read.table("TheorExpSims/All.qual.freq.0.dens.1.txt")->td   # dens dependent 

# FD

meansFD<-c()
stdFD<-c()
meansFDPrev<-c()
stdFDPrev<-c()

indThresh<-0

for (i in seq(1,100)) {
	meansFD<-c(meansFD,mean(tfd[,i]+tfd[,i+100],na.rm=T))
	stdFD<-c(stdFD,sd(tfd[,i]+tfd[,i+100],na.rm=T))
}
plot(meansFD)


for (i in seq(101,200)) {
	meansFDPrev<-c(meansFDPrev,mean(tfd[,i]/(tfd[,i-100]+tfd[,i]),na.rm=T))
	stdFDPrev<-c(stdFDPrev,sd(tfd[,i]/(tfd[,i-100]+tfd[,i]),na.rm=T))
}
plot(meansFD,meansFDPrev)

# Dens

meansD<-c()
stdD<-c()
meansDPrev<-c()
stdDPrev<-c()

indThresh<-0

for (i in seq(1,100)) {
	meansD<-c(meansD,mean(td[,i]+td[,i+100],na.rm=T))
	stdD<-c(stdD,sd(td[,i]+td[,i+100],na.rm=T))
}
plot(meansD)


for (i in seq(101,200)) {
	meansDPrev<-c(meansDPrev,mean(td[,i]/(td[,i-100]+td[,i]),na.rm=T))
	stdDPrev<-c(stdDPrev,sd(td[,i]/(td[,i-100]+td[,i]),na.rm=T))
}
plot(meansD,meansDPrev)

# Freq

meansF<-c()
stdF<-c()
meansFPrev<-c()
stdFPrev<-c()



for (i in seq(1,100)) {
	meansF<-c(meansF,mean(tf[,i]+tf[,i+100],na.rm=T))
	stdF<-c(stdF,sd(tf[,i]+tf[,i+100],na.rm=T))
}
plot(meansD)


for (i in seq(101,200)) {
	meansFPrev<-c(meansFPrev,mean(tf[,i]/(tf[,i-100]+tf[,i]),na.rm=T))
	stdFPrev<-c(stdFPrev,sd(tf[,i]/(tf[,i-100]+tf[,i]),na.rm=T))
}
plot(meansF,meansFPrev)




mydf<-data.frame(pos=c(seq(1,100),seq(1,100),seq(1,100)),abun=c(meansFD,meansD,meansF),prev=c(meansFDPrev,meansDPrev,meansFPrev),maxAbun<-c(meansFD+stdFD,meansD+stdD,meansF+stdF),minAbun<-c(pmax(meansFD-stdFD,0),pmax(meansD-stdD,0),pmax(meansF-stdF,0)),minPrev<-c(pmax(meansFDPrev-stdFDPrev,0),pmax(meansDPrev-stdDPrev,0),pmax(meansFPrev-stdFPrev,0)),maxPrev<-c(meansFDPrev+stdFDPrev,meansDPrev+stdDPrev,meansFPrev+stdFPrev),type=c(rep("both",100),rep("density-dependent",100),rep("frequency-dependent",100)))

plD<-ggplot()+geom_ribbon(data=mydf,aes(x=pos,ymin=minAbun,ymax=maxAbun,fill=type),alpha=0.3)+scale_color_manual(values=wes_palette("Darjeeling1"))+xlab("position")+ylab("Abundance")+geom_line(data=mydf,aes(pos,abun,col=type),size=1)+theme(legend.position="none")+scale_fill_manual(values=wes_palette("Darjeeling1"))

plE<-ggplot()+geom_ribbon(data=mydf,aes(x=pos,ymin=minPrev,ymax=maxPrev,fill=type),alpha=0.5)+scale_color_manual(values=wes_palette("Darjeeling1"))+xlab("position")+ylab("Prevalence")+geom_line(data=mydf,aes(pos,prev,col=type),size=1)+theme(legend.position="none")+scale_fill_manual(values=wes_palette("Darjeeling1"))+ylim(c(0,0.5))

plF<-ggplot()+geom_ribbon(data=mydf,aes(x=abun,ymin=minPrev,ymax=maxPrev,fill=type),alpha=0.5)+scale_color_manual(values=wes_palette("Darjeeling1"),name="")+xlab("abundance")+ylab("prevalence")+geom_line(data=mydf,aes(abun,prev,col=type),size=1)+scale_fill_manual(values=wes_palette("Darjeeling1"),name="")+xlim(c(0,4.5))+ylim(c(0,0.5))


row2<-plot_grid(plD,plE,plF,ncol=3,labels=c("D","E","F"),rel_widths=c(1,1,1.6))

fin<-plot_grid(row1,row2,ncol=1)

fin

ggsave("/Users/uricchio/projects/bootslabOld/anther_smut/manuscript/v2/figures/Fig5.pdf",fin,width=14,height=7)
