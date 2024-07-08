#!/bin/sh
fileIndex=$1
ends=$2

filePath=/data3/liangyuan/05_LOCK_and_PMD/data_CellLine/02_fastq/WGBS_MethylCSeq/
procedurePath=/data3/liangyuan/05_LOCK_and_PMD/data_CellLine/06_WGBS_run/
resultPath=/data3/liangyuan/05_LOCK_and_PMD/data_CellLine/06_WGBS_run/

mkdir -p ${procedurePath}
mkdir -p ${resultPath}

if [[ ${ends} -eq 1 ]]
then
bismark --genome /data1/reference/hg38/indexed_bismark/ --parallel 10 --prefix ${fileIndex} -o ${procedurePath} ${filePath}/${fileIndex}.fastq.gz
/opt/Anaconda3/bin//deduplicate_bismark --bam --output_dir ${procedurePath} ${procedurePath}/${fileIndex}.${fileIndex}_bismark_bt2.bam 
/opt/Anaconda3/bin//bismark_methylation_extractor --gzip --bedGraph --parallel 10 -o ${resultPath} ${procedurePath}/${fileIndex}.${fileIndex}_bismark_bt2.deduplicated.bam
else
bismark --genome /data1/reference/hg38/indexed_bismark/ --parallel 10 --prefix ${fileIndex} -o ${procedurePath} -1 ${filePath}/${fileIndex}_1.fastq.gz -2 ${filePath}/${fileIndex}_2.fastq.gz
deduplicate_bismark --bam --output_dir ${procedurePath} ${procedurePath}/${fileIndex}.${fileIndex}_1_bismark_bt2_pe.bam 
bismark_methylation_extractor --gzip --bedGraph --parallel 10 -o ${resultPath} ${procedurePath}/${fileIndex}.${fileIndex}_1_bismark_bt2_pe.deduplicated.bam
fi