#!/bin/sh
inputFile=$1
outputPath=/data3/liangyuan/05_LOCK_and_PMD/data_CellLine/03_ChIPSeq_run/00.factqc/

mkdir -p ${outputPath}
fastqc -o ${outputPath} /data3/liangyuan/05_LOCK_and_PMD/data_CellLine/02_fastq/ChIPSeq/${inputFile}
