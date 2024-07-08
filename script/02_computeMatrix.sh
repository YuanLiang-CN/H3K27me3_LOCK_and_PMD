cd /data3/liangyuan/05_LOCK_and_PMD/data_Roadmap/00_commonPMD_subgroup/
Repli_bw=$1 
computeMatrix scale-regions -p 8 \
-R 02_merged_PMD_HMD/commonPMD_hg19_merge.bed \
-S 01_repli_bw/${Repli_bw} \
-b 0 -a 0 --regionBodyLength 10000000 -bs 1000 -o 03_computeMatrix/${Repli_bw}.gz
