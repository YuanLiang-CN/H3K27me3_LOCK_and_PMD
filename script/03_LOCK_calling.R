args = commandArgs(T)
EID = args[1]
HisMark = args[2]

library(CREAM)

peak_dir = '/data3/liangyuan/05_LOCK_and_PMD/data_Roadmap/01_narrowPeak_sorted/'
peak_file = paste0(peak_dir,'/',HisMark,'/',EID,'_',HisMark,'_peak.bed')

LOCK_dir = '/data3/liangyuan/05_LOCK_and_PMD/data_Roadmap/03_LOCK/'
if (!file.exists(paste0(LOCK_dir,'/',HisMark,'/',EID))) {dir.create(paste0(LOCK_dir,'/',HisMark,'/',EID))}
LOCK_path = paste0(LOCK_dir,'/',HisMark,'/',EID,'/')

for (ctoff in seq(-1,1,0.1)) {
  print(c(EID,HisMark,ctoff))
  tryCatch({
    LOCK = CREAM(peak_file, WScutoff = ctoff, MinLength = 1000, peakNumMin = 2)
    write.table(LOCK,file = paste0(LOCK_path,'/',EID,'_',HisMark,'_ctoff',ctoff,'.bed'),sep = '\t',row.names = F,col.names = F,quote = F)
  },error = function(e){
    print('error')
    write.table('error',file = paste0(LOCK_path,'/',EID,'_',HisMark,'_ctoff',ctoff,'_error.txt'),sep = '\t',row.names = F,col.names = F,quote = F)
  })
}