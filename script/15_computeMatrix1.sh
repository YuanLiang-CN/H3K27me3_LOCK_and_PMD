TumorAndNormalLOCK_fileNameWithPath=$1
bw_fileNameWithPath=$2
outputFileNameWithPath=$3
logFileNameWithPath=$4

/opt/Anaconda3/bin//computeMatrix scale-regions -p25 -R ${TumorAndNormalLOCK_fileNameWithPath} \
-S ${bw_fileNameWithPath} \
--beforeRegionStartLength 50000 \
--regionBodyLength 100000 \
--afterRegionStartLength 50000 \
--binSize 10 \
 -o ${outputFileNameWithPath} > ${logFileNameWithPath} 2>&1
