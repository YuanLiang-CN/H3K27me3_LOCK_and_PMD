#!/bin/sh

inputFileName=$1
echo ${inputFileName}

inputFilePath=/data3/liangyuan/05_LOCK_and_PMD/data_CellLine/03_ChIPSeq_run/00.fastq/
outputPath=/data3/liangyuan/05_LOCK_and_PMD/data_CellLine/03_ChIPSeq_run/01.trim_galore/
logPath=/data3/liangyuan/05_LOCK_and_PMD/data_CellLine/03_ChIPSeq_run/log/

mkdir -p ${outputPath}
mkdir -p ${outputPath}/fastqc/
mkdir -p ${logPath}

trim_galore --gzip --length 30 --fastqc_args "-o ${outputPath}/fastqc/"  \
   --phred33 --illumina ${inputFilePath}/${inputFileName}.fastq.gz \
   -o ${outputPath} 2> ${logPath}/${inputFileName}.trim_galore.log
