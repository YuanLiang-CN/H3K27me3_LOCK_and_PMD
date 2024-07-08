TumorAndNormalLOCK_fileNameWithPath=$1
bw_fileNameWithPath=$2
outputFileNameWithPath=$3
logFileNameWithPath=$4

/opt/Anaconda3/bin//computeMatrix scale-regions -p25 -R ${TumorAndNormalLOCK_fileNameWithPath} \
-S ${bw_fileNameWithPath} \
-b 0 -a 0 --regionBodyLength 100000 -bs 100 \
-o ${outputFileNameWithPath} > ${logFileNameWithPath} 2>&1
