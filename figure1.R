#Figure1
#Wei
#2017.08.22

library(ggplot2)
library(gplots)

#read data
ALL_sign<-read.table("../figure1/significant gene express.txt",sep="\t",header=T)

#calculate fold change
ALL_sign$fold_change<-apply(ALL_sign,1,function(x) mean(as.numeric(x[grepl("NC",colnames(ALL_sign))]))/mean(as.numeric(x[grepl("UC",colnames(ALL_sign))])))

UP<-which(ALL_sign$fold_change<1)[order(ALL_sign[which(ALL_sign$fold_change<1),"upwtest.intersect.up_sig..fc_sig.."])[1:25]]
DOWN<-which(ALL_sign$fold_change>1)[order(ALL_sign[which(ALL_sign$fold_change>1),"upwtest.intersect.up_sig..fc_sig.."])[1:25]]

ALL_sign[UP,"upwtest.intersect.up_sig..fc_sig.."]
ALL_sign[DOWN,"upwtest.intersect.up_sig..fc_sig.."]
#take out the value
GEvalue<-ALL_sign[,grepl("GCRMA",colnames(ALL_sign))]
head(as.matrix(GEvalue))
heatmap.2(as.matrix(GEvalue)[c(UP,DOWN),],scale="row",col = rev(redblue(20)),trace="none")
