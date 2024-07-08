EID=$1

cd /data3/liangyuan/05_LOCK_and_PMD

/opt/Anaconda3/bin/computeMatrix scale-regions -p 3 -S data_Roadmap/06_WGBS_bw/rmCGI_WGBS_bw/${EID}_WGBS_subtrCGI.bw \
-R data_Roadmap/03_LOCK/H3K27me3/${EID}/${EID}_H3K27me3_ctoff-0.7.bed \
-b 0 -a 0 --regionBodyLength 10 \
-o data_Roadmap/08_LOCK_DNA_methylation/H3K27me3_LOCK_met_signal/${EID}_-0.7_LOCK_met_signal.gz
