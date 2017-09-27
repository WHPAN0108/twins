#Figure5
#Wei
#20170925

library(reshape2)

#GENE EXPRESSION

#read all expression
UC_CD<-read.xlsx("../../../validation/validation data/WEI__validation UC_CD_gene list_v3.xlsx", sheetName = "results",startRow=11)
control<-read.xlsx("../../../validation/validation data/WeiHung_BMBF_GeneExpression_2015-09-01.xlsx", sheetName = "results",startRow=11)

#combine two data set and give the NA value
all<-rbind(UC_CD,control)
all[all=="Undetermined"] <- NA

#list the gene and ID name
head(all)
gene.name<-as.character(unique(all$Detector.Name))
ID<-(as.character(unique(all$ID)))

#create the new dataframe for mean value
New_data<-data.frame(matrix(NA,ncol=length(ID),nrow=length(gene.name)),row.names=gene.name)
colnames(New_data)<-ID

#put the value in the dataframe
for(i in gene.name){
  for(j in ID){
    k=which(all$Detector.Name==i & all$ID==j)
    avg<-mean(as.numeric(as.character(all[k,"Ct"])),na.rm=T)
    New_data[i,j]<-avg  
  }  
}

#delta gene expression by normalize by housekeep gene
norm_new_data<-t(apply(New_data[-1,],1,function(x) x-as.numeric(New_data[1,])))
norm_new_data<-norm_new_data[,-which(colnames(norm_new_data)=="water")]
norm_new_data_box<-as.data.frame(t(norm_new_data))
norm_new_data_box$type<-rep(c("UC","CD","control"),each=20)
norm_new_data_box$type <- ordered(norm_new_data_box$type, levels=c("CD", "UC", "control"))

GE_table<-data.frame(CCL11=2^-(norm_new_data_box[c(1:20,41:60),15]),ISG20=2^-(norm_new_data_box[c(1:20,41:60),2]),
                        LYN=2^-(norm_new_data_box[c(1:20,41:60),3]),TNFSF10=2^-(norm_new_data_box[c(1:20,41:60),5]),
                        OAS1=2^-(norm_new_data_box[c(1:20,41:60),6]),AGT=2^-(norm_new_data_box[c(1:20,41:60),9]),
                        CFB=2^-(norm_new_data_box[c(1:20,41:60),11]),S100a8=2^-(norm_new_data_box[c(1:20,41:60),14]),
                        type=norm_new_data_box[c(1:20,41:60),16])
GE_table$type<-gsub("control","Healthy",as.character(GE_table$type))
Mean_GE<-apply(subset(GE_table,type=="UC")[,-9],2,function(x) mean(x,na.rm=T))
GE_table_normalized<-sweep(GE_table[,-9], 2, Mean_GE, `/`)
GE_table_normalized<-cbind(GE_table_normalized,Type=GE_table$type)
GE_table_normalized<-melt(GE_table_normalized)

#plot
library(ggplot2)

All_GE<-ggplot(aes(y=100*value,x=variable,fill=Type),data=GE_table_normalized)+
  ylab("normalized mRNA level")+scale_fill_discrete(guide=FALSE)+xlab("")+
  stat_summary(fun.y=mean, geom="bar",position=position_dodge(0.8),width=0.8) + 
  stat_summary(fun.data = mean_se,geom="errorbar",color="grey40",position=position_dodge(0.8), width=.2)+
  theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   axis.text=element_text(size=18),axis.title.y = element_text(size=20) )

## Methylation and gene expression

#read data
methylation_old<-read.table("../../../validation/first panel/20140616.cpg.txt",sep="\t",header=T)
gene_expression_old<-read.table("../../../validation/first panel/correlated_gene_20140617.txt",sep="\t",header=T)
CFB_cg09583599<-read.table("../../../validation/first panel/CFB_cg09583599.txt",sep="\t",header=T)
CFB_cg01883966<-read.table("../../../validation/first panel/CFB_cg01883966.txt",sep="\t",header=T)
TNFSF10_cg11979312<-read.table("../../../validation/first panel/TNFSF10_cg11979312.txt",sep="\t",header=T)
LYN_cg03973663<-read.table("../../../validation/first panel/LYN_cg03973663.txt",sep="\t",header=T)

#plot all
m_g_combine<-function(methylation,gene_expression,second_panel){
  Condition<-rep("HC",16)
  Condition[seq(from=2,to=16,2)]<-"UC"
  new_data<-cbind(methylation,gene_expression)[,c(1,4)]
  new_data<-data.frame(new_data,Condition)
  colnames(new_data)<-c("methylation","gene_expression","Condition")
  new_data<-new_data[c("Condition","methylation","gene_expression")]
  
  colnames(second_panel)<-c("Condition","methylation","gene_expression")
  all_data<-rbind(new_data,second_panel)
  Group<-c(rep("Twins",16),rep("Validation",40))
  all_data<-data.frame(Group,all_data)  
  
  return(all_data)
}

sort_gene_expression<-function(gene_name){
  gene_value<-gene_expression_old[gene_expression_old$Gene.Symbol==gene_name,7:26]
  sample_name<-do.call(rbind,strsplit(colnames(gene_expression_old)[7:26],"_"))
  new_gene<-data.frame(sample_name[,2],t(gene_value))[-c(4,6,14,16),]
  new_gene<-new_gene[order(new_gene[1]),]
  colnames(new_gene)<-c("sample_number","gene_expression")
  return(new_gene)
}

sort_methylation<-function(cg_number,gene_name){
  new_meth<-methylation_old[which(methylation_old$cgnumber==cg_number & methylation_old$Gene==gene_name),]
  new_meth<-data.frame(t(new_meth[-c(1:3)]))
  new_meth<-data.frame(new_meth[order(rownames(new_meth)),],rownames(new_meth)[order(rownames(new_meth))])
  return(new_meth)
}

CFB_cg09583599_new<-m_g_combine(sort_methylation("cg09583599","CFB"),sort_gene_expression("CFB"),CFB_cg09583599)
CFB_cg09583599_plot<-ggplot(aes(x=methylation,y=log2(gene_expression)),data=subset(CFB_cg09583599_new, methylation>0.6 & Group=="Validation"))+geom_point(size=5,aes(color=Condition))+
  theme_bw()+ylab("log2 gene expression level")+xlab("DNA methylation level")+
  geom_smooth(method=lm, se=FALSE,color="black",linetype="dashed")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=18),axis.title.y = element_text(size=20),axis.title.x = element_text(size=20),
        legend.position="none",plot.title = element_text(hjust = 0.5))

## bacteria and gene expression

#read data
#read data from xlsx file
UC_CD_B_DNA<-read.xlsx("../../../validation/validation data/WP_Richa_UC_CD_BMBF_BaterialAssays.xlsx", sheetName = "UC_CD_DNA_results",startRow=11)
UC_CD_B_DNA_HK<-read.xlsx("../../../validation/validation data/WP_RB_UC_CD_BMBF_BacterialAssays_Part2.xlsx", sheetName = "UC_CD_DNA_results",startRow=11)
UC_CD<-rbind(UC_CD_B_DNA,UC_CD_B_DNA_HK)
control_B<-read.xlsx("../../../validation/validation data/WP_Richa_UC_CD_BMBF_BaterialAssays.xlsx", sheetName = "BMBF_DNA_results",startRow=11)
control_B_HK<-read.xlsx("../../../validation/validation data/WP_RB_UC_CD_BMBF_BacterialAssays_Part2.xlsx", sheetName = "BMBF_results",startRow=11)
control_all<-rbind(control_B,control_B_HK)
control_0<-control_all[grep("A|water", control_all$ID) , ]
all_B<-rbind(UC_CD,control_0)
all_B[all_B=="Undetermined"] <- NA
bacteria.name<-as.character(unique(all_B$Detector.Name))
ID_B<-as.character(unique(all_B$ID))

#Normalized
#put the value in the dataframe
New_data_B<-data.frame(matrix(NA,ncol=length(ID_B),nrow=length(bacteria.name)),row.names=bacteria.name)
colnames(New_data_B)<-ID_B
for(i in bacteria.name){
  for(j in ID_B){
    k=which(all_B$Detector.Name==i & all_B$ID==j)
    avg<-mean(as.numeric(as.character(all_B[k,"Ct"])),na.rm=T)
    New_data_B[i,j]<-avg  
  }  
}

#Normalize with housekeeping gene
N1_data<-t(apply(New_data_B,1,function(x) (x-as.numeric(New_data_B["GAPDH_Hs02758991_g1",]))))
N1_data

#Normalize with all bacteria number
P2_data<-t(apply(N1_data,1,function(x) 2^(x-as.numeric(N1_data["bact_quant_all_bacteria",]))))
colnames(P2_data)<-colnames(New_data_B)
P2_data<-P2_data[-which(grepl("bact_quant_all_bacteria|GAPDH|ACTB",rownames(P2_data))),-which(grepl("water",colnames(P2_data)))]
P2_data["bacteroides",][c(1:20,41:60)]
#plot
DF_CFB<-data.frame(CFB=2^-norm_new_data["CFB",c(1:20,41:60)],bacteroides=log10(P2_data["bacteroides",c(1:20,41:60)]),
                   Condition=rep(c("UC","H"),each=20))

CFB_Bacteroides<-ggplot(aes(y=log2(CFB),x=bacteroides,color=Condition),data=subset(DF_CFB,CFB<0.6 & bacteroides<0))+geom_point(size=5)+
  theme_bw()+ geom_smooth(method=lm, se=FALSE,color="black",linetype="dashed")+
  ylab("log2 gene expression level")+xlab("log10 relative bacteria level")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=18),axis.title.y = element_text(size=20),axis.title.x = element_text(size=20),
        legend.position="none",plot.title = element_text(hjust = 0.5,size=30))


#Combine all three figures
library(ggpubr)
ggarrange(CFB_cg09583599_plot, CFB_Bacteroides,ncol = 2, nrow = 1) 
ggarrange(All_GE, ggarrange(CFB_cg09583599_plot, CFB_Bacteroides, ncol = 2), nrow = 2) 

library("gridExtra")
grid.arrange(All_GE,                             # First row with one plot spaning over 2 columns
             arrangeGrob(CFB_cg09583599_plot, CFB_Bacteroides, ncol = 2), # Second row with 2 plots in 2 different columns
             nrow = 2) 