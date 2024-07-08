# KYSE450_shortLOCK_EMT
tissue_type=$1
cd /data3/liangyuan/05_LOCK_and_PMD/data_CellLine/18_scanMotifGenomeWide.pl/
bedtools getfasta -fi /data1/reference/hg38/hg38.fasta \
-bed /data3/liangyuan/05_LOCK_and_PMD/data_CellLine/16_motif/input/NE2_KYSE450.fg.bed > CGI_promoter_KYSE450_loss.fasta
scanMotifGenomeWide.pl /opt/homer/./data/knownTFs/known.motifs CGI_promoter_KYSE450_loss.fasta > CGI_promoter_KYSE450_loss.txt