
read.csv('~/projects/bootslab/anther_smut/simAnthSmut/data/Transect census 2019 .csv')->te

mydfH = data.frame(y=te$Y.global+1,x=te$X+1,h=te$N-te$dis)
mydfD = data.frame(y=te$Y.global+1,x=te$X+1,d=te$dis)

write.table(mydfD,'~/projects/bootslab/anther_smut/simAnthSmut/data/dis.19.txt',row.names=FALSE,col.names=FALSE)
write.table(mydfH,'~/projects/bootslab/anther_smut/simAnthSmut/data/dens.19.txt',row.names=FALSE,col.names=FALSE)


pf<-function(p){
	return ((p*(1-p))^0.75)
}

plot(pf(seq(0,1,0.0001)))