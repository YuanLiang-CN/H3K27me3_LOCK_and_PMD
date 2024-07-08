#!/bin/sh
ALL_INDEX=$1
ipIndex=`echo $ALL_INDEX | sed 's/__.*$//g'`
controlIndex=`echo $ALL_INDEX | sed 's/^.*__//g'`
outputName=${ipIndex}
GSE=`echo ${ipIndex} | sed 's/_.*$//g'`

thread=10
inputFilePath=/data3/liangyuan/01_tumor_and_normal/data/03_ChIPSeq_procedure/${GSE}
outputPath=/data3/liangyuan/01_tumor_and_normal/data/03_ChIPSeq_result/${GSE}

echo -e "${ipIndex}\t${controlIndex}"
macs2 callpeak -t ${inputFilePath}/${ipIndex}.sorted.mkdup.rmblacklist.bam -c ${inputFilePath}/${controlIndex}.sorted.mkdup.rmblacklist.bam \
      -g 2.7e9 -q 0.01 --keep-dup 1 --nomodel --extsize 146 -n ${outputName} --outdir ${outputPath}

bamCompare -b1 ${inputFilePath}/${ipIndex}.sorted.mkdup.rmblacklist.bam -b2 ${inputFilePath}/${controlIndex}.sorted.mkdup.rmblacklist.bam \
           --operation subtract -o ${outputPath}/${outputName}.SubtractCPM.bw \
           --binSize 20 --numberOfProcessors ${thread} --scaleFactorsMethod None --normalizeUsing CPM \
           --ignoreDuplicates --extendReads 146
