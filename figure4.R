#figure4
#20170825
#Wei

library(gplots)
library(ggplot2)
library(xlsx)

#read data
COR_matrix<-read.xlsx("../figure4/correlation_matrix_OTUvsTranscripts.xlsx",sheetIndex = 3,header=T)
COR_value<-COR_matrix[-1,-c(1,2)]
