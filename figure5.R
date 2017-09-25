#Figure5
#Wei
#20170925

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

