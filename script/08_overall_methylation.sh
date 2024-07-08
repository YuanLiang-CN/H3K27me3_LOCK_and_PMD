EID=$1

cd /data3/liangyuan/05_LOCK_and_PMD

/opt/Anaconda3/bin/computeMatrix scale-regions -p 3 -S data_Roadmap/06_WGBS_bw/rmCGI_WGBS_bw/${EID}_WGBS_subtrCGI.bw \
-R data_Roadmap/08_overall_DNA_methylation/SIL_overall_met/bed/S_PMD.bed data_Roadmap/08_overall_DNA_methylation/SIL_overall_met/bed/I_PMD.bed \
data_Roadmap/08_overall_DNA_methylation/SIL_overall_met/bed/L_PMD.bed data_Roadmap/08_overall_DNA_methylation/SIL_overall_met/bed/HMD.bed\
-b 0 -a 0 --regionBodyLength 10 \
-o data_Roadmap/08_overall_DNA_methylation/SIL_overall_met/gz/${EID}_SIL_met_signal.gz
