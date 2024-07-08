EID=$1
hisMark=$2
/opt/Anaconda3/bin/computeMatrix scale-regions -p 3 -S /data3/liangyuan/05_LOCK_and_PMD/data_Roadmap/02_ChIPSeq_bw/${hisMark}/${EID}-${hisMark}.pval.signal.bigwig \
-R /data3/liangyuan/05_LOCK_and_PMD/data_Roadmap/01_narrowPeak_sorted/${hisMark}/${EID}_${hisMark}_peak.bed \
-b 0 -a 0 --regionBodyLength 200 \
-o /data3/liangyuan/05_LOCK_and_PMD/data_Roadmap/05_peak_signal/${hisMark}/${EID}_${hisMark}PeakSignal.gz
