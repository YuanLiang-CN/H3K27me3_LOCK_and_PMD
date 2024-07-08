WGBS_bwFile=$1
inputPath="/data3/liangyuan/05_LOCK_and_PMD/data_CellLine/06_WGBS_bw/WGBS_bw/"
outputPath_bw="/data3/liangyuan/05_LOCK_and_PMD/data_CellLine/06_WGBS_bw/rmCGI_WGBS_bw/"
outputPath_bdg="/data3/liangyuan/05_LOCK_and_PMD/data_CellLine/06_WGBS_bw/rmCGI_WGBS_bdg/"

/opt/Anaconda3/bin//bigWigToBedGraph ${inputPath}/${WGBS_bwFile} ${outputPath_bdg}/${WGBS_bwFile%%bw}bdg
bedtools intersect -v -a ${outputPath_bdg}/${WGBS_bwFile%%bw}bdg -b /data3/liangyuan/05_LOCK_and_PMD/meta/CGI/CGI_merge_hg38.bed > ${outputPath_bdg}/${WGBS_bwFile%%bw}rmCGI.bdg
/opt/Anaconda3/bin//bedGraphToBigWig ${outputPath_bdg}/${WGBS_bwFile%%bw}rmCGI.bdg /data1/liangyuan/linux_deal/03_LOCK_PMD/01_data/05_chrom_size/hg38.chrom.sizes ${outputPath_bw}/${WGBS_bwFile%%bw}rmCGI.bw