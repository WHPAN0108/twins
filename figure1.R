#Figure1
#Wei
#2017.08.22

library(ggplot2)
library(gplots)
library(stringr)

#read data
ALL_sign<-read.table("../figure1/significant gene express.txt",sep="\t",header=T)

#calculate fold change
ALL_sign$fold_change<-apply(ALL_sign,1,function(x) mean(as.numeric(x[grepl("NC",colnames(ALL_sign))]))/mean(as.numeric(x[grepl("UC",colnames(ALL_sign))])))

UP<-which(ALL_sign$fold_change<1)[order(ALL_sign[which(ALL_sign$fold_change<1),"upwtest.intersect.up_sig..fc_sig.."])[1:25]]
DOWN<-which(ALL_sign$fold_change>1)[order(ALL_sign[which(ALL_sign$fold_change>1),"upwtest.intersect.up_sig..fc_sig.."])[1:25]]

ALL_sign[UP,"upwtest.intersect.up_sig..fc_sig.."]
ALL_sign[DOWN,"upwtest.intersect.up_sig..fc_sig.."]

#take out the value
colnames(ALL_sign)[grepl("GCRMA",colnames(ALL_sign))]
label_x<-gsub("(\\d+).*", "\\1", colnames(ALL_sign)[grepl("GCRMA",colnames(ALL_sign))])
label_x<-gsub("NC","H",label_x)
col_side<-(rep(c("darkgreen","orange"),each=10))

GEvalue<-ALL_sign[,grepl("GCRMA",colnames(ALL_sign))]
GEvalue50<-as.matrix(GEvalue)[c(UP,DOWN),]
GEvalue50<-t(apply(GE_value,1,scale))

heatmap.2(GEvalue50,scale="row",col = rev(redblue(11)),trace="none",
          labRow = ALL_sign$Gene.Symbol[c(UP,DOWN)],labCol = label_x,
          ColSideColors = col_side,density.info="none")



