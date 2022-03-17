# concept figure for spatial model
library(cowplot)
library(ggplot2)

# replace with your directory path 
setwd("/Users/uricchio/projects/bootslabOld/anther_smut/manuscript/v2/software/simAnthSmut/plottingScripts")

read.table("../data/empDensData.txt")->gr
read.table("../exampleData/simExample.txt")->gs
read.table("../exampleData/qualExample.txt")->qs

p1<-ggplot(gr, aes(V1, V2, fill = V6))+geom_tile(color = "gray")+ylab("")+xlab("")+scale_fill_gradient2(low = "blue", high = "yellow", mid = "white", midpoint = 8,name="Abundance")

p2<-ggplot(gs, aes(V2+1	, V1+1, fill = V3))+geom_tile(color = "white")+ylab(expression("y-coord. (" * italic(m) * ")"))+xlab("")+scale_fill_gradient2(low = "blue", high = "yellow", mid = "white", midpoint = 14,name="Abundance")


p3<-ggplot(qs, aes(V2+1	, V1+1, fill = V3))+geom_tile(color = "white")+xlab(expression("x-coord. (" * italic(m) * ")"))+ylab("")+scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.5,name="Quality       ")


 myGrid<-plot_grid(p1,p2,p3,ncol=1)
 ggsave("~/projects/bootslabOld/anther_smut/simAnthSmut/figures/ConceptFigPlots.pdf",myGrid,width=7,height=6)
