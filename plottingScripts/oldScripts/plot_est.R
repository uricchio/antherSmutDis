library(ggplot2)
library(cowplot)
library(wesanderson)

# get correlations between Fh and unseen variables

read.table("~/projects/bootslab/anther_smut/simAnthSmut/col_corr.txt")->cc
# Fh and Vh
plot(cc$V1,cc$V3)
lm(cc$V3~cc$V1)
# Fh and Vd
plot(cc$V1,cc$V4)
lm(cc$V4~cc$V1)
# Fh and J
plot(cc$V1,cc$V5)
lm(cc$V5~cc$V1)



# get perf of inference assuming we know all relevant variables

read.table("~/projects/bootslab/anther_smut/simAnthSmut/simPerf.BestCase.txt")->sB
read.table("~/projects/bootslab/anther_smut/simAnthSmut/simPerf.Biased.txt")->sF

sF<-sF[sF$V15!="Inf",]

make_df<-function(s) {
	df<-data.frame()
	for(i in seq(1,length(s$V1))) {
		df<-rbind(df,data.frame(x="migSeeds",est=(s$V1[i]-s$V6[i])/s$V6[i]))
	}
	for(i in seq(1,length(s$V1))) {
		df<-rbind(df,data.frame(x="migSpores",est=(s$V2[i]-s$V7[i])/s$V7[i]))
	}
	for(i in seq(1,length(s$V1))) {
		df<-rbind(df,data.frame(x="betaF",est=(s$V3[i]-s$V8[i])/s$V8[i]))
	}
	for(i in seq(1,length(s$V1))) {
		df<-rbind(df,data.frame(x="betaV",est=(s$V4[i]-s$V9[i])/s$V9[i]))
	}
	for(i in seq(1,length(s$V1))) {
		df<-rbind(df,data.frame(x="betaJ",est=(s$V5[i]-s$V10[i])/s$V10[i]))
	}
	return(df)
}

df<-make_df(sB)

plA<-ggplot(df,aes(x,est,fill=as.factor(x)))+geom_boxplot()+scale_fill_manual(values=wes_palette("Darjeeling1"))+geom_hline(yintercept=0,lty=3)+xlab("Parameter")+ylab("Bias")+theme(legend.position="none")+ylim(c(-1,5))

# now assuming we know only Fh and Fd

df<-make_df(sF)

plB<-ggplot(df,aes(x,est,fill=as.factor(x)))+geom_boxplot()+scale_fill_manual(values=wes_palette("Darjeeling1"))+geom_hline(yintercept=0,lty=3)+xlab("Parameter")+ylab("Bias")+theme(legend.position="none")+ylim(c(-1,5))

# now assuming we can't use Vn and Vd and j in the prob calcs

plot_grid(plA,plB,labels=c("A","B"))

# abnundance

read.table("~/projects/bootslab/anther_smut/simAnthSmut/abundance.txt")->ab
ggplot(ab,aes(seq(1,100),V1+V3))+xlab("position")+ylab("abundance")+geom_smooth(span=0.1)+geom_point()+ylim(c(0,200))


# abundance over time

read.table("~/projects/bootslab/anther_smut/simAnthSmut/abundanceDensOnly.txt")->dens
read.table("~/projects/bootslab/anther_smut/simAnthSmut/abundanceDensOnlyNoQual.txt")->dens
read.table("~/projects/bootslab/anther_smut/simAnthSmut/abundanceDisInfParm.txt")->infer
read.table("~/projects/bootslab/anther_smut/simAnthSmut/abundanceFreqOnly.txt")->freq

plA<-ggplot(infer,aes(seq(1,980),V6))+geom_line(col='red')
plB<-ggplot(freq,aes(seq(1,980),V6))+geom_line(col='red')

plot_grid(plA,plB)

