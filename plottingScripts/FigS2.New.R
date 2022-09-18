
# max change in abundance / year in data is (1409-504)/(18-14) = 226.25
library(ggplot2)
library(cowplot)
library(wesanderson)
library(viridis)

plotAndGetMaxSlope<-function(i) {	

  read.table(paste("/Users/uricchio/projects/bootslabOld/anther_smut/manuscript/v2/software/simAnthSmut/simsRepGenRandHigh/sim.",i,".txt",sep=""))->t
  scan(paste("/Users/uricchio/projects/bootslabOld/anther_smut/manuscript/v2/software/simAnthSmut/simsRepGenRandHigh/sim.",i,".txt",sep=""),nlines=1,what=c("a"))->params
  
  #read.table(paste("/Users/uricchio/projects/bootslabOld/anther_smut/simAnthSmut/simsRepGenRandHigh/sim.",i,".txt",sep=""))->t
  #scan(paste("/Users/uricchio/projects/bootslabOld/anther_smut/simAnthSmut/simsRepGenRandHigh/sim.",i,".txt",sep=""),nlines=1,what=c("a"))->params

  
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
for ( i in seq(1,5000)) {
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
	
	return(mean(df[df$a>at & df$b>bt & df$a<at+bw & df$b<bt+bw ,]$t,na.rm=T))
}


pdf<-data.frame()
for (i in seq(0,19,1)) {
	 for (j in seq(0,19,1)) {
	     pdf<-rbind(pdf,data.frame(a=i,b=j,f=myPdf(data_plot,i,j,1),t=myldf(data_plot,i,j,1)))    	
	 } 
}

pl2<-ggplot(pdf,aes(a,b,color=f))+ geom_point(shape=15,size=5)+scale_colour_viridis(option='viridis',name=expression(italic(F)))+xlab(expression(frac(beta[d],beta[d[0]])))+ylab(expression(frac(beta[f],beta[f[0]])))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))

pl3<-ggplot(pdf,aes(a,b,color=t))+ geom_point(shape=15,size=5)+scale_colour_viridis(option='magma',name=expression(italic(P)))+xlab(expression(frac(beta[d],beta[d[0]])))+ylab(expression(frac(beta[f],beta[f[0]])))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))


plot_grid(pl2,pl3,labels=c("A","B"),rel_widths=c(1,1),ncol=2)


ggsave("~/projects/bootslabOld/anther_smut/manuscript/v2/figures/FigS1.pdf",width=7.5*1.3,height=2.8*1.3)



