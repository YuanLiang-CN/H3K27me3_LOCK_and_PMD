#!/bin/sh

thread=10
inputFileName=$1
echo ${inputFileName}

inputFilePath=/data3/liangyuan/05_LOCK_and_PMD/data_CellLine/03_ChIPSeq_run/01.trim_galore/
outputPath=/data3/liangyuan/05_LOCK_and_PMD/data_CellLine/03_ChIPSeq_run/02.mapping/
mkdir -p ${outputPath}

genomeIndex="/data1/reference/hg38/indexed_BWA/hg38.fa"
blackList="/data1/reference/hg38/hg38-blacklist.v2.bed"

bwa mem -t ${thread} $genomeIndex ${inputFilePath}/${inputFileName}_trimmed.fq.gz > ${outputPath}/${inputFileName}.raw.sam
samtools flagstat ${outputPath}/${inputFileName}.raw.sam > ${outputPath}/${inputFileName}.raw.flagstat.txt
samtools view -q 1 -bS ${outputPath}/${inputFileName}.raw.sam  > ${outputPath}/${inputFileName}.bam
samtools sort -@ ${thread} ${outputPath}/${inputFileName}.bam -o ${outputPath}/${inputFileName}.sorted.bam
sambamba markdup -t ${thread} --tmpdir=/data1/tmp/ ${outputPath}/${inputFileName}.sorted.bam ${outputPath}/${inputFileName}.sorted.mkdup.bam
bedtools intersect -v -a ${outputPath}/${inputFileName}.sorted.mkdup.bam -b ${blackList} > ${outputPath}/${inputFileName}.sorted.mkdup.rmblacklist.bam
samtools index ${outputPath}/${inputFileName}.sorted.mkdup.rmblacklist.bam
rm ${outputPath}/${inputFileName}.raw.sam ${outputPath}/${inputFileName}.bam 
rm ${outputPath}/${inputFileName}.sorted.bam ${outputPath}/${inputFileName}.sorted.mkdup.bam
