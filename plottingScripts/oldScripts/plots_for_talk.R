library(ggplot2)
library(cowplot)
library(wesanderson)

bear_stand<-runif(1000)
bear_sit<-bear_stand+rnorm(1000,0.,0.1)
palette<-wes_palette("Zissou1", 1000, type = "continuous")
plA<-ggplot(data.frame(x=bear_stand,y=bear_sit,z=as.factor(bear_stand-bear_sit)),aes(x,y,color=z))+geom_point()+xlim(0,1)+ylim(0,1)+xlab(" frequency")+ylab(" frequency")+scale_color_manual(values=palette)+theme(legend.position="NA")


#div<-c(0.)
#for (i in seq(2,10000)) {
#	ne<-div[i-1]+rnorm(1,0,0.1)
	#ne<-max(0,ne)
	#ne<-min(1,ne)
#	div<-c(div,ne)
#}

plB<-ggplot(data.frame(x=seq(1,1000000,100),y=div,z=as.factor(ifelse(abs(div)>5,5,0))),aes(x,y,color=z))+geom_point(size=0.5)+xlab("genomic position")+ylab("Standardized XP-EHH")+geom_hline(yintercept=0,lty=3)+geom_hline(yintercept=-5,lty=2)+geom_hline(yintercept=5,lty=2)+scale_color_manual(values=c("black","red"))+theme(legend.position="NA")


vals<-data.frame(ts=c(rnorm(1000,-1,1),rnorm(1000,0,1),rnorm(1000,1,1)),genotype=c(rep("AA",1000),rep("AT",1000),rep("TT",1000)))

plC<-ggplot(vals,aes(genotype,ts))+geom_boxplot(fill=wes_palette("Darjeeling1")[4])+ylab("Standardized standing time")

plot_grid(plC,plA,plB,labels=c("A","B","C"),ncol=3)
ggsave("~/projects/interviews/cornell/hyp_data.pdf",height=3*1.2,width=14*1.2)


# now data slide


phenoBef<-function(x) {
	 return(exp(-x^2*10))
}

phenoAf<-function(x) {
	 return(exp(-(x-1)^2))
}

phenoBefLv<-function(x) {
	 return(exp(-x^2/2))
}

phenoBefVlv<-function(x) {
	 return(exp(-x^2/4))
}

df<-data.frame(x=seq(-3,3,0.01),y=phenoBef(seq(-3,3,0.01)),type=rep("Before shift",601))
df<-rbind(df,data.frame(x=seq(-3,3,0.01),y=phenoAf(seq(-3,3,0.01)),type=rep("After shift",601)))

dfarrow=data.frame(x=c(0),y=c(1),type=c(""))


#df2<-data.frame()
#for ( i in seq(1,10)) {
#	xinit = runif(1)*0.5
#	xfin=runif(1)*0.5+0.5
#	
 #   df2<-rbind(df2,data.frame(x0=0,xi=xinit,xf=xfin,x1=1,chr=i))	
#}

plA<-ggplot(df2,aes(x0,chr))+geom_segment(aes(xend=x1,yend=chr),size=4.5)+geom_segment(aes(x=xi,y=chr,xend=xf,yend=chr,col="red"),size=4.5)+theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())





plB <- ggplot(df,aes(x,y,fill=type))+geom_area(position="dodge")+ylab("")+scale_y_continuous(breaks=c())+xlab("height")+scale_x_continuous(breaks=c())+geom_segment(data=dfarrow,aes(x=x,xend = x+1,y=y+0.1,yend = y+0.1),arrow = arrow(length = unit(0.1,"cm")))+theme(legend.position="NA")



df<-data.frame(x=seq(-5,5,0.01),y=phenoBef(seq(-5,5,0.01)),type=rep("Total variance",1001))
df<-rbind(df,data.frame(x=seq(-5,5,0.01),y=phenoBefLv(seq(-5,5,0.01)),type=rep("Genetic variance",1001)))
df<-rbind(df,data.frame(x=seq(-5,5,0.01),y=phenoBefVlv(seq(-5,5,0.01)),type=rep("Explained variance",1001)))

plC<- ggplot(df,aes(x,y,fill=type))+geom_area(position="dodge",alpha=0.3)+ylab("")+scale_y_continuous(breaks=c())+xlab("height")+scale_x_continuous(breaks=c())+theme(legend.position="NA")+geom_segment(data=dfarrow,aes(x=-1.2,xend = 1.2,y=0.68,yend =0.68),arrow = arrow(length = unit(0.1,"cm")))+geom_segment(data=dfarrow,aes(x=1.2,xend = -1.2,y=0.68,yend =0.68),arrow = arrow(length = unit(0.1,"cm")))+geom_segment(data=dfarrow,aes(x=-1,xend = 1,y=0.55,yend =0.55),arrow = arrow(length = unit(0.1,"cm")))+geom_segment(data=dfarrow,aes(x=1,xend = -1,y=0.55,yend =0.55),arrow = arrow(length = unit(0.1,"cm")))+geom_segment(data=dfarrow,aes(x=-0.3,xend = 0.3,y=0.43,yend =0.43),arrow = arrow(length = unit(0.1,"cm")))+geom_segment(data=dfarrow,aes(x=0.3,xend = -0.3,y=0.43,yend =0.43),arrow = arrow(length = unit(0.1,"cm"))) 

plC<-plC+geom_text(data=data.frame(a=-3.7,b=0.68,l="V[P]",type=""),aes(a,b,label=l),parse=TRUE,size=6)
plC<-plC+geom_text(data=data.frame(a=-3.7,b=0.55,l="V[G]",type=""),aes(a,b,label=l),parse=TRUE,size=6)
plC<-plC+geom_text(data=data.frame(a=-3.7,b=0.43,l="V[x < 0.05]",type=""),aes(a,b,label=l),parse=TRUE,size=6)

#plC<-plC+geom_label(data=data.frame(a=3.4,b=0.55,l="V[] < V[G]",type=""),aes(a,b,label=l),parse=TRUE,size=6)

plot_grid(plA,plB,plC,ncol=1)
ggsave("~/projects/interviews/cornell/problems.pdf",height=6.5*1.3,width=4.3*1.3)



