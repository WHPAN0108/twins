#figure4
#20170825
#Wei

library(gplots)
library(ggplot2)
library(xlsx)
library(RColorBrewer)

#read data
COR_matrix<-read.xlsx("../figure4/OTU_gene_table.xlsx",sheetIndex = 1,header=F)
COR_value<-COR_matrix[-c(1,2,3),-1]
COR_value<-t(apply(COR_value,1,as.numeric))
gene_name<-as.character(unlist(COR_matrix[-c(1,2,3),1]))
phyla_name<-as.character(unlist(COR_matrix[1,-1]))
bacteria_name<-as.character(unlist(COR_matrix[2,-1]))
col_side<-(rep(c("darkgreen","orange"),each=10))
table(bacteria_name)

#plot
colours <- colorRampPalette(brewer.pal(7, "Dark2"))(7)[as.numeric(as.factor(bacteria_name))]
heatmap.2(t(COR_value),trace="none",col=colorRampPalette(c("yellow","blue"))(11),labCol = gene_name,
          hclustfun =function(x) hclust(x,method = "complete"),RowSideColors=colours,labRow = "")

#plot2
Bacteria_cor<-apply(t(COR_value),1,function(x) sum(abs(x)>0.5))
colours2 <- colorRampPalette(brewer.pal(7, "Set1"))(5)[as.numeric(as.factor(phyla_name[which(Bacteria_cor>5)]))]
heatmap.2(t(COR_value[,which(Bacteria_cor>5)]),col=colorRampPalette(c("yellow","blue"))(9),trace="none",
          hclustfun =function(x) hclust(x,method = "complete"),labCol = gene_name,
          RowSideColors=colours2,labRow = bacteria_name[which(Bacteria_cor>5)])

