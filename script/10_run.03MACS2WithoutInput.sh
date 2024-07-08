#!/bin/sh

ipIndex=$1
outputName=${ipIndex}
GSE=`echo ${ipIndex} | sed 's/_.*$//g'`


thread=5
inputFilePath=/data3/liangyuan/01_tumor_and_normal/data/03_ChIPSeq_procedure/${GSE}
outputPath=/data3/liangyuan/01_tumor_and_normal/data/03_ChIPSeq_result/${GSE}

echo -e "${ipIndex}"
macs2 callpeak -t ${inputFilePath}/${ipIndex}.sorted.mkdup.rmblacklist.bam \
      -g 2.7e9 -q 0.01 --keep-dup 1 --nomodel --extsize 146 -n ${outputName} --outdir ${outputPath}

bamCoverage -b ${inputFilePath}/${ipIndex}.sorted.mkdup.rmblacklist.bam --normalizeUsing CPM --binSize 20 \
      --numberOfProcessors ${thread} --ignoreDuplicates --extendReads 146 \
      -o ${outputPath}/${outputName}.CPM.bw
