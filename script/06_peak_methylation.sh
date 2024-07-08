EID=$1

/opt/Anaconda3/bin/computeMatrix scale-regions -p 3 -S /data3/liangyuan/05_LOCK_and_PMD/data_Roadmap/06_WGBS_bw/rmCGI_WGBS_bw/${EID}_WGBS_subtrCGI.bw \
-R /data3/liangyuan/05_LOCK_and_PMD/data_Roadmap/01_narrowPeak_sorted/H3K27me3/${EID}_H3K27me3_peak.bed \
-b 0 -a 0 --regionBodyLength 200 \
-o /data3/liangyuan/05_LOCK_and_PMD/data_Roadmap/07_peak_DNA_methylation/H3K27me3_peak_met_signal/${EID}_H3K27me3PeakMetSignal.gz
