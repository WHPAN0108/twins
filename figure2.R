#figure2
#microbiota analysis
#20170102/20170823
#wei

#library
library(ggplot2)
library(xlsx)
library(fossil)
library(vegan)
library(cowplot)
library(gridExtra)
library(vegan)
library(ape)

#read data
Bact<-read.table("../figure2/bacteria_table_uc.txt",head=T,sep="\t",row.names = 1)

#alpha diversity
Chao_index<-apply(Bact,2,chao1)
Shannon_index<-apply(Bact,2,diversity)
condition<-rep(c("Healthy","UC"),5)
paired<-rep(1:5,each=2)
DF<-data.frame(Chao_index,Shannon_index,condition,paired)

ggplot(DF, aes(y = Chao_index)) +
  geom_boxplot(aes(x = rep(c(-3, 3), 5), group = condition,fill=condition)) +
  geom_point(aes(x = rep(c(-1, 1),5)), size = 5) +xlab("")+ylab("Chao1 index")+
  geom_line(aes(x = rep(c(-1, 1),  5), group = paired))+ggtitle("Chao1 index")

ggplot(DF, aes(y = Shannon_index)) +
  geom_boxplot(aes(x = rep(c(-3, 3), 5), group = condition,fill=condition)) +
  geom_point(aes(x = rep(c(-1, 1),5)), size = 5) +xlab("")+ylab("Shannon_index")+
  geom_line(aes(x = rep(c(-1, 1),  5), group = paired))+ggtitle("Shannon_index")

t.test(DF[c(1,3,5,7,9),1],DF[c(2,4,6,8,10),1],paired=T)
wilcox.test(DF[c(1,3,5,7,9),1],DF[c(2,4,6,8,10),1],paired = T)
mean(DF[c(2,4,6,8,10),2])

#beta diversity
bray=vegdist(t(Bact),"bray")
jaccard=vegdist(t(Bact),"jaccard")
mantel(bray,jaccard)
bray_pcoa<-pcoa(bray)
jacc_pcoa<-pcoa(jaccard)

DF<-data.frame(DF,Bray1=bray_pcoa$vectors[,1],Bray2=bray_pcoa$vectors[,2],
               Bray3=bray_pcoa$vectors[,3],Jacc1=jacc_pcoa$vectors[,1],Jacc2=jacc_pcoa$vectors[,2],Jacc3=jacc_pcoa$vectors[,3])

ggplot(DF,aes(x =Bray1 , y = Bray2,color=condition))+ geom_point(size=5)+geom_line(aes(group=paired),color="steelblue")+
  ggtitle("Pcoa Bray distance")+xlab("Bray1 22.89%")+ylab("Bray2 16.47%")+theme_bw()+ 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=18),axis.title.y = element_text(size=20),axis.title.x = element_text(size=20),
        plot.title = element_text(hjust = 0.5,size=30))

#ggplot(DD,aes(x =Jacc1 , y = Jacc2,color=Week))+ geom_point(size=3)+ggtitle("Jacc")+facet_grid(Organ~.) 
##

