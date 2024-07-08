library(CREAM)
library(stringr)
library(dplyr)
setwd("/data3/liangyuan/05_LOCK_and_PMD/")
args = commandArgs(T)
peakFileNameWithPath=args[1]
fileIndex = str_match(peakFileNameWithPath, '.*\\/(.*).peak.bed')[,2]
LOCKPath = 'data_CellLine/04_LOCK/all_LOCKs/'
LongLOCKPath = 'data_CellLine/04_LOCK/long_LOCKs/'
ShortLOCKPath = 'data_CellLine/04_LOCK/short_LOCKs/'

#CREAM
WS = -0.7
LOCK = CREAM(peakFileNameWithPath, WScutoff = WS, MinLength = 1000, peakNumMin = 2)
write.table(LOCK,file = paste0(LOCKPath,fileIndex,'.LOCK_-0.7.bed'),
            sep = '\t',row.names = F,col.names = F,quote = F)
colnames(LOCK) = c('chrom','start','end')
LOCK$d = LOCK$end - LOCK$start
Long_LOCK = filter(LOCK, d > 1e5); Long_LOCK$d = NULL
Short_LOCK = filter(LOCK, d <= 1e5); Short_LOCK$d = NULL
write.table(Long_LOCK, file = paste0(LongLOCKPath,fileIndex,'.Long_LOCK_-0.7.bed'), 
            sep = '\t',row.names = F,col.names = F,quote = F)
write.table(Short_LOCK, file = paste0(ShortLOCKPath,fileIndex,'.Short_LOCK_-0.7.bed'), 
            sep = '\t',row.names = F,col.names = F,quote = F)