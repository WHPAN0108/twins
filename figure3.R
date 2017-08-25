#Figure3
#Wei
#20170823

library(GenomicRanges)
library(Gviz)
library(ggplot2)
library(xlsx)
library(biomaRt)
library(reshape2)
library(gridExtra)
library(rlang)
library(ggpubr)

#read Cpg Info
CpG_grange<-read.xlsx("../figure3/methylation.xlsx",header=T,sheetIndex = 1)
cpgIslands1<-makeGRangesFromDataFrame(CpG_grange[-c(1,2),])
data <- read.table("../figure3/cytoBandIdeo.txt", header=F, sep="\t")
colnames(data) <-c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain')

#plot data

atrack <- AnnotationTrack(cpgIslands1,genome = "hg19",chromosome = "chr1",name="CpG")
plotTracks(atrack)
gtrack <- GenomeAxisTrack()
ideoTrack <- IdeogramTrack(genome="hg19",chromosome = "chr1", bands=data)
plotTracks(list(ideoTrack,gtrack, atrack),showBandId = TRUE)

get.gene.track <- function(chr, start, end, assembly = "mm10") {
  mart <- NULL
  featMap <- Gviz:::.getBMFeatureMap()
  if (assembly == "mm10") {
    mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", 
                    host = "www.ensembl.org")
  } 
  if (assembly == "hg38") {
    mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  } 
  else {
    stop(paste("Assembly", assembly, "not implemented!"))
  }
  
  gene.track <- BiomartGeneRegionTrack(genome = assembly, strand="+",
                                       biomart = mart, featureMap = featMap, chromosome = chr, 
                                       start = start, end = end, name = "Ensembl GENE", 
                                       geneSymbols = TRUE,background.title = "grey50") 
  return(gene.track)
}

gene.track_S100A8<-get.gene.track(chr="chr1",start= 153352508,end=153372508,assembly = "hg38")
plotTracks(gene.track_S100A8)
plotTracks(list(ideoTrack, gtrack,gene.track_S100A8,atrack),showBandId = TRUE)

#correlation plot

#read data
GE<-read.table("../figure3/significant gene express.txt",header=T,sep="\t")
head(GE)
DM<-read.table("../figure3/20140616.cpg.txt",header=T,sep="\t")
head(DM)

GE_S100A8<-subset(GE,Gene.Symbol=="S100A8")[,grepl("GCRMA",colnames(GE))]
GE_S100A8<-melt(GE_S100A8,value.name = "Gene_Expression",variable.name = "sampleID")
GE_samplename_cut<-sapply(strsplit(as.character(GE_S100A8$sampleID), '_'),function(x) x[2])
GE_S100A8$sampleID<-toupper(gsub("^0+([1-9])","\\1",GE_samplename_cut))


DM_cg02813121<-subset(DM,cgnumber=="cg02813121")[,grepl("X",colnames(DM))]
DM_cg02813121<-melt(DM_cg02813121,value.name = "methylation",variable.name = "sampleID")
DM_cg02813121$sampleID<-gsub("X","",DM_cg02813121$sampleID)

GE_DM<-merge(DM_cg02813121,GE_S100A8,by="sampleID",all=T)
GE_DM$condition<-ifelse(grepl("A",GE_DM$sampleID),"H","UC")
COR_plot<-ggplot(aes(y=Gene_Expression,x=methylation,color=condition),data=GE_DM)+
  geom_point(size=3)+ theme(legend.position="bottom")
GE_plot<-ggplot(aes(y=Gene_Expression,x=condition,fill=condition),data=GE_DM)+
  geom_boxplot()+ theme(legend.position="bottom")

ggarrange(COR_plot, GE_plot,ncol = 2, nrow = 1)
