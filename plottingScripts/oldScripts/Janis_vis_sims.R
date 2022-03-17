library(ggplot2)
library(cowplot)
library(wesanderson)
library(viridis)

get_data <-function(betaF,betaVJ,index,thresh) {
	file = paste("~/projects/bootslab/anther_smut/simAnthSmut/simsForPlot/sim.betaF.",betaF,".betaVJ.",betaVJ,".txt",sep="")
	read.table(file)->te
	#print (te[,index])
	me = mean(te[,index][te[,10]>thresh])
	se = sd(te[,index][te[,10]>thresh])/sqrt(length(te[,index][te[,10]>thresh]))
	return(list(me=me,su = me+se,sd=me-se))
}

index_list<-c(0.05,0.1,0.2,0.5,1,2,5,10,20)
thresh<-2

# first, disease persistence

VJonly<-c()
VJonlySU<-c()
VJonlySD<-c()
for (i in index_list) {
 VJonly<-c(VJonly,get_data(0,i,10,thresh)$me)
 VJonlySU<-c(VJonlySU,get_data(0,i,10,thresh)$su)
 VJonlySD<-c(VJonlySD,get_data(0,i,10,thresh)$sd)

}

Fonly<-c()
FonlySU<-c()
FonlySD<-c()
for (i in index_list) {
 Fonly<-c(Fonly,get_data(i,0,10,thresh)$me)
 FonlySU<-c(FonlySU,get_data(i,0,10,thresh)$su)
 FonlySD<-c(FonlySD,get_data(i,0,10,thresh)$sd)
}

VJfixandF<-c()
VJfixandFSU<-c()
VJfixandFSD<-c()
for (i in index_list) {
 VJfixandF<-c(VJfixandF,get_data(i,1,10,thresh)$me)
 VJfixandFSU<-c(VJfixandFSU,get_data(i,1,10,thresh)$su)
 VJfixandFSD<-c(VJfixandFSD,get_data(i,1,10,thresh)$sd)
}

FfixandVJ<-c()
FfixandVJSU<-c()
FfixandVJSD<-c()
for (i in index_list) {
 FfixandVJ<-c(FfixandVJ,get_data(1,i,10,thresh)$me)
 FfixandVJSU<-c(FfixandVJSU,get_data(1,i,10,thresh)$su)
 FfixandVJSD<-c(FfixandVJSD,get_data(1,i,10,thresh)$sd)
}
df<-data.frame(x=index_list,y=VJonly,yd=VJonlySD,yu=VJonlySU,type=rep("VJ",9))
df<-rbind(df,data.frame(x=index_list,y=Fonly,yd=FonlySD,yu=FonlySU,type=rep("F",9)))
df<-rbind(df,data.frame(x=index_list,y=FfixandVJ,yd=FfixandVJSD,yu=FfixandVJSU,type=rep("F_VJ",9)))
df<-rbind(df,data.frame(x=index_list,y=VJfixandF,yd=VJfixandFSD,yu=VJfixandFSU,type=rep("VJ_F",9)))


plA<-ggplot(df,aes(x,y,col=type))+geom_point()+geom_line()+scale_x_log10(breaks=c(0.05,0.2,1,5,20))+scale_y_log10(breaks=c(5,10,50,100,500,1000),limits=c(5,1000))+geom_ribbon(aes(ymin=yd,ymax=yu,fill=type,col=type),linetype=0,alpha=0.3)+scale_color_manual(values=wes_palette("Moonrise3"))+scale_fill_manual(values=wes_palette("Moonrise3"))+xlab(expression(beta[rel]))+ylab(expression(italic(T) * " (duration)"))+theme(legend.position="NA")


# next, extinction

VJonly<-c()
VJonlySU<-c()
VJonlySD<-c()
for (i in index_list) {
 VJonly<-c(VJonly,get_data(0,i,11,thresh)$me)
 VJonlySU<-c(VJonlySU,get_data(0,i,11,thresh)$su)
 VJonlySD<-c(VJonlySD,get_data(0,i,11,thresh)$sd)

}

Fonly<-c()
FonlySU<-c()
FonlySD<-c()
for (i in index_list) {
 Fonly<-c(Fonly,get_data(i,0,11,thresh)$me)
 FonlySU<-c(FonlySU,get_data(i,0,11,thresh)$su)
 FonlySD<-c(FonlySD,get_data(i,0,11,thresh)$sd)
}

VJfixandF<-c()
VJfixandFSU<-c()
VJfixandFSD<-c()
for (i in index_list) {
 VJfixandF<-c(VJfixandF,get_data(i,1,11,thresh)$me)
 VJfixandFSU<-c(VJfixandFSU,get_data(i,1,11,thresh)$su)
 VJfixandFSD<-c(VJfixandFSD,get_data(i,1,11,thresh)$sd)
}

FfixandVJ<-c()
FfixandVJSU<-c()
FfixandVJSD<-c()
for (i in index_list) {
 FfixandVJ<-c(FfixandVJ,get_data(1,i,11,thresh)$me)
 FfixandVJSU<-c(FfixandVJSU,get_data(1,i,11,thresh)$su)
 FfixandVJSD<-c(FfixandVJSD,get_data(1,i,11,thresh)$sd)
}

df<-data.frame(x=index_list,y=VJonly,yd=VJonlySD,yu=VJonlySU,type=rep("VJ",9))
df<-rbind(df,data.frame(x=index_list,y=Fonly,yd=FonlySD,yu=FonlySU,type=rep("F",9)))
df<-rbind(df,data.frame(x=index_list,y=FfixandVJ,yd=FfixandVJSD,yu=FfixandVJSU,type=rep("F_VJ",9)))
df<-rbind(df,data.frame(x=index_list,y=VJfixandF,yd=VJfixandFSD,yu=VJfixandFSU,type=rep("VJ_F",9)))



plB<-ggplot(df,aes(x,y,col=type))+geom_point()+geom_line()+scale_x_log10(breaks=c(0.05,0.2,1,5,20))+geom_ribbon(aes(ymin=yd,ymax=yu,fill=type,col=type),linetype=0,alpha=0.3)+scale_color_manual(values=wes_palette("Moonrise3"),name="",
                       labels=c(expression(beta[F] * " = 0, " * beta[VJ] * " = " * beta[rel]), expression(beta[F] * " = " * beta[rel] * ", " * beta[VJ] * " = 0"), expression(beta[F] * " = 1, " * beta[VJ] * " = " * beta[rel]), expression(beta[F] * " = " * beta[rel] * ", " * beta[VJ] * " = 1")))+scale_fill_manual(values=wes_palette("Moonrise3"),name="",
                       labels=c(expression(beta[F] * " = 0, " * beta[VJ] * " = " * beta[rel]), expression(beta[F] * " = " * beta[rel] * ", " * beta[VJ] * " = 0"), expression(beta[F] * " = 1, " * beta[VJ] * " = " * beta[rel] ), expression(beta[F] * " = " * beta[rel] * ", " * beta[VJ] * " = 1")))+xlab(expression(beta[rel]))+ylab(expression(italic(P) * " (persistence)"))


# max change in abundance / year in data is (1409-504)/(18-14) = 226.25
library(ggplot2)
library(cowplot)
library(wesanderson)

plotAndGetMaxSlope<-function(i) {	

  read.table(paste("~/projects/bootslab/anther_smut/simAnthSmut/simsRepGenRandHigh/sim.",i,".txt",sep=""))->t
  scan(paste("~/projects/bootslab/anther_smut/simAnthSmut/simsRepGenRandHigh/sim.",i,".txt",sep=""),nlines=1,what=c("a"))->params
  
  # if epidemic is shorter than actual, skip
  if(length(t$V1) < 10 ) {
  	return (NA)
  }
 
 
  len=min(length(t$V1),1000)
  
  yvals = (t$V1+t$V2)[1:len]
  popSize<-data.frame(x=seq(1,len),y=yvals)

   M<-nls(y~y0/(1+exp(-r*(x-x0)))+ym, start=list(r=0.2,y0=-(max(yvals)),x0=20,ym=max(yvals)),data=popSize,trace=FALSE,lower=c(0.01,-10000,0.01,0.01),algorithm="port"
,control = nls.control(minFactor = 1/10000000,maxiter=1000,warnOnly=TRUE))
  M.df<-data.frame(x=seq(1,len))


  # y = a / (1 + exp(-b *(x-x0))) + c
  M.df$y<-coef(M)[2]/(1+exp(-coef(M)[1]*(M.df$x-coef(M)[3])))+coef(M)[4]

  # Deriv is (a b E^(b (x + x0)))/(E^(b x) + E^(b x0))^2 which is maximized at x0, with a value of a*b/4
  
  myMax = c(t(t(abs(coef(M)[2]*coef(M)[1]/4))))
  
  #if(myMax > 100) {
    #print(myMax)
    #print(coef(M)[3])
    #p<-ggplot(data=popSize,aes(x,y),alpha=0.5)+geom_point()+ylim(c(0,2500))+geom_line(data=M.df,col='red',size=0.7)
    #p+xlim(1,len)
    #print(p)
  #}
  
  tot_final = sum(t[len,][1:5],na.rm=T)
  print(tot_final)
  if(tot_final > 0) {
  	tot_final = 1
  	}
   
  return(data.frame(m=myMax,a=as.numeric(params[2]),b=as.numeric(params[3]),l=len,t=tot_final))
}

data <-data.frame()
emp<-0
for ( i in seq(1,15000)) {
	m<-plotAndGetMaxSlope(i)
	if(is.na(m)) {
		emp<-emp+1
	} else {
	    data<-rbind(data,m)
	}
}

data_plot<-data[data$m<300,]

pl1<-ggplot(data_plot,aes(a,b))+ geom_bin2d(binwidth=c(1,1),aes(fill=cut(m,seq(0,300,15))))+scale_fill_manual(values=wes_palette("Zissou1",21,type="continuous"))


myPdf<-function(df,at,bt,bw) {
	if (length(df[df$a>at & df$b>bt & df$a<at+bw & df$b<bt+bw,]$m) == 0) {
		return(0)
	}
	return(length(df[df$a>at & df$b>bt & df$a<at+bw & df$b<bt+bw & df$m >226.25,]$m)/length(df[df$a>at & df$b>bt & df$a<at+bw & df$b<bt+bw,]$m))
}

myldf<-function(df,at,bt,bw) {
	
	return(mean(df[df$a>at & df$b>bt & df$a<at+bw & df$b<bt+bw ,]$t))
}


pdf<-data.frame()
for (i in seq(0,19,1)) {
	 for (j in seq(0,19,1)) {
	     pdf<-rbind(pdf,data.frame(a=i,b=j,f=myPdf(data_plot,i,j,1),t=myldf(data_plot,i,j,1)))    	
	 } 
}

pl2<-ggplot(pdf,aes(a,b,color=f))+ geom_point(shape=15,size=5)+scale_colour_viridis(option='viridis',name=expression(italic(F)))+xlab(expression(frac(beta[VJ],beta[VJ[0]])))+ylab(expression(frac(beta[F],beta[F[0]])))

pl3<-ggplot(pdf,aes(a,b,color=t))+ geom_point(shape=15,size=5)+scale_colour_viridis(option='magma',name=expression(italic(P)))+xlab(expression(frac(beta[VJ],beta[VJ[0]])))+ylab(expression(frac(beta[F],beta[F[0]])))


plot_grid(pl2,pl3,labels=c("A","B"),rel_widths=c(1,1),ncol=2)



ggsave("~/projects/bootslab/anther_smut/simAnthSmut/figures/figSim.pdf",width=7.5*1.3,height=2.8*1.3)



