#!/bin/sh
fileIndex=$1
ends=$2

echo ${fileIndex}
cd /data3/liangyuan/05_LOCK_and_PMD/data_CellLine/02_fastq/RNASeq/

outputPath_procedure=/data3/liangyuan/05_LOCK_and_PMD/data_CellLine/05_RNASeq_run/
outputPath_result=/data3/liangyuan/05_LOCK_and_PMD/data_CellLine/05_RNASeq/counts/
mkdir -p ${outputPath_procedure}
mkdir -p ${outputPath_result}

if [[ ${ends} -eq 1 ]]
then
echo "single-ended"
hisat2 -p 5 -x /data1/reference/hg38/grch38_snp_tran/genome_snp_tran -U ${fileIndex}.fastq.gz | \
       samtools view -bS - > ${outputPath_procedure}/${fileIndex}.bam 
htseq-count -f bam -s no ${outputPath_procedure}/${fileIndex}.bam /data1/liangyuan/linux_deal/03_LOCK_PMD/01_data/10_gtf/Homo_sapiens.GRCh38.104.gtf > ${outputPath_result}/${fileIndex}.count.txt
else
echo "double-ended"
hisat2 -p 5 -x /data1/reference/hg38/grch38_snp_tran/genome_snp_tran -1 ${fileIndex}_1.fastq.gz -2 ${fileIndex}_2.fastq.gz | \
       samtools view -bS - > ${outputPath_procedure}/${fileIndex}.bam 
htseq-count -f bam -s no ${outputPath_procedure}/${fileIndex}.bam /data1/liangyuan/linux_deal/03_LOCK_PMD/01_data/10_gtf/Homo_sapiens.GRCh38.104.gtf > ${outputPath_result}/${fileIndex}.count.txt
fi


