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
row.names(COR_value)<-gene_name
colnames(COR_value)<-bacteria_name
col_side<-(rep(c("darkgreen","orange"),each=10))
table(bacteria_name)

#plot
colours <- colorRampPalette(brewer.pal(7, "Dark2"))(7)[as.numeric(as.factor(bacteria_name))]
heatmap.2(t(COR_value),trace="none",col=colorRampPalette(c("yellow","blue"))(11),labCol = gene_name,
          hclustfun =function(x) hclust(x,method = "complete"),RowSideColors=colours,labRow = "")

#plot2
Bacteria_cor<-apply(t(COR_value),1,function(x) sum(abs(x)>0.5))
Target_bact_gene<-t(COR_value[,which(Bacteria_cor>5)])
row.names(Target_bact_gene)
new_order_row<-c("Otu010_Clostridium_XlVa(95)","Otu052_Dialister(100)","Otu024_Lachnospiracea_incertae_sedis(100)",
                 "Otu042_Blautia(97)","Otu040_Clostridium_XlVa(62)","Otu045_un_Lachnospiraceae(100)",
                 "Otu102_Ruminococcus(100)","Otu133_un_Lachnospiraceae(100)","Otu062_Dorea(100)",
                 "Otu075_un_Ruminococcaceae(100)","Otu026_Barnesiella(100)","Otu050_Alistipes(100)",
                 "Otu004_Bacteroides(100)","Otu113_Bacteroides(100)","Otu088_Alistipes(100)",
                 "Otu029_Bacteroides(100)","Otu207_Desulfovibrio(100)","Otu111_Sutterella(100)",
                 "Otu153_Akkermansia(100)","Otu007_un_Bacteria(100)")
Target_bact_gene[new_order_row,]
colours2 <- colorRampPalette(brewer.pal(7, "Set1"))(5)[as.numeric(as.factor(phyla_name[which(Bacteria_cor>5)]))]
heatmap.2(Target_bact_gene[new_order_row,],col=colorRampPalette(c("yellow","blue"))(9),trace="none",
          hclustfun =function(x) hclust(x,method = "complete"),Rowv = F)

dim(t(COR_value[,which(Bacteria_cor>5)]))

