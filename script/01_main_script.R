suppressMessages(library(valr))
suppressMessages(library(CREAM))
suppressMessages(library(ggpubr))
suppressMessages(library(ggrepel))
suppressMessages(library(cowplot))
suppressMessages(library(circlize))
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ComplexHeatmap))
setwd('/data3/liangyuan/05_LOCK_and_PMD/')

# PART 1: metadata----
# 1. CGI----
list.files('meta')
# [1] "blacklist"           "CGI"                 "chrom_size"          "commonPMD_commonHMD" "commonPMD_subgroup"  "gtf"                
# [7] "liftover"            "LOLACore"            "RefSeq_gene"         "Roadmap_table"       "TSG_OG"

UCSC_hg38 = read.table(file = 'meta/CGI/CGI_UCSC_hg38.bed',sep = '\t',header = F)[,1:3] %>% rename(chrom = 1,start = 2,end = 3)
UCSC_hg19 = read.table(file = 'meta/CGI/CGI_UCSC_hg19.bed',sep = '\t',header = F)[,1:3] %>% rename(chrom = 1,start = 2,end = 3)
Irizarry_hg38 = read.table(file = 'meta/CGI/Irizarry2009-model-based-cpg-islands-hg38.bed',sep = '\t',header = F)[,1:3] %>% rename(chrom = 1,start = 2,end = 3)
Irizarry_hg19 = read.table(file = 'meta/CGI/Irizarry2009-model-based-cpg-islands-hg19.bed',sep = '\t',header = F)[,1:3] %>% rename(chrom = 1,start = 2,end = 3)
Takai_Jones_hg38 = read.table(file = 'meta/CGI/Takai_Jones_from_Fei.hg38.bed',sep = '\t',header = F) %>% rename(chrom = 1,start = 2,end = 3)
Takai_Jones_hg19 = read.table(file = 'meta/CGI/Takai_Jones_from_Fei.hg19.bed',sep = '\t',header = F) %>% rename(chrom = 1,start = 2,end = 3)

CGI_hg38 = rbind(UCSC_hg38,Irizarry_hg38,Takai_Jones_hg38) %>% bed_sort() %>% bed_merge()
CGI_hg19 = rbind(UCSC_hg19,Irizarry_hg19,Takai_Jones_hg19) %>% bed_sort() %>% bed_merge()

write.table(CGI_hg38,file = 'meta/CGI/CGI_merge_hg38.bed',sep = '\t',quote = F,row.names = F,col.names = F)
write.table(CGI_hg19,file = 'meta/CGI/CGI_merge_hg19.bed',sep = '\t',quote = F,row.names = F,col.names = F)

list.files('meta/CGI/')
# [1] "CGI_merge_hg19.bed"                            "CGI_merge_hg38.bed"                            "CGI_UCSC_hg19.bed"                             "CGI_UCSC_hg38.bed"                            
# [5] "Irizarry2009-model-based-cpg-islands-hg19.bed" "Irizarry2009-model-based-cpg-islands-hg38.bed" "Takai_Jones_from_Fei.hg19.bed"                 "Takai_Jones_from_Fei.hg38.bed" 

# 2. Separate common PMDs into S/I/L-PMDs----
# 1) merge common PMD and common HMD
commonPMD_HMD = read.table(file = 'meta/commonPMD_commonHMD/PMD_coordinates_hg19.bed', sep = '\t',header = F)
get_merged_MD = function(MD){
  PMD_hg19 = filter(commonPMD_HMD, V6 == common)
  PMD_hg19 = PMD_hg19[,1:4]
  colnames(PMD_hg19) = c('chrom','start','end','variation')
  PMD_hg19$end = PMD_hg19$end+1
  PMD_merge = bed_merge(PMD_hg19)
  PMD_hg19_merge = bed_intersect(PMD_merge,PMD_hg19,suffix = c('','.y'))
  PMD_hg19_merge = group_by(PMD_hg19_merge,chrom,start,end) %>% 
    summarise(variation = mean(variation.y,na.rm = T)) %>% 
    ungroup()
  return(PMD_hg19_merge)
}
commonPMD_hg19_merge = get_merged_MD('commonPMD')
commonHMD_hg19_merge = get_merged_MD('commonHMD')
write.table(commonPMD_hg19_merge, file = 'meta/commonPMD_subgroup/02_merged_PMD_HMD/commonPMD_hg19_merge.bed', 
            sep = '\t',quote = F, row.names = F, col.names = T)
write.table(commonHMD_hg19_merge, file = 'meta/commonPMD_subgroup/02_merged_PMD_HMD/commonHMD_hg19_merge.bed', 
            sep = '\t',quote = F, row.names = F, col.names = T)

# 2) run computeMatrix
system('ls meta/commonPMD_subgroup/01_repli_bw | xargs -iFile -P5 bash -c "sh script/02_computeMatrix.sh File "')

# 3) get S/I/L PMD
# read computeMatrix result
repli_res = list()
for (i in 1:length(list.files('meta/commonPMD_subgroup/03_computeMatrix', pattern = 'gz'))) {
  computeMatrix_result = read.table(list.files('meta/commonPMD_subgroup/03_computeMatrix', pattern = 'gz', full.names = T)[i],
                                    sep = '\t',header = F,skip = 1)
  repli_res[[i]] = rowMeans(computeMatrix_result[,7:10006],na.rm = T)
}
names(repli_res) = str_match(list.files('meta/commonPMD_subgroup/03_computeMatrix', pattern = 'gz'),'(.*)_.*?')[,2]
repli_cells = str_match(list.files('meta/commonPMD_subgroup/03_computeMatrix', pattern = 'gz'),'(.*)_.*?')[,2]
repli_res = do.call(cbind, repli_res)
hg19_PMD_features = cbind(mutate(commonPMD_hg19_merge, d = end - start, .after = end), repli_res)

#quantile normalization
library(devtools)
library(Biobase)
library(preprocessCore)

edata = hg19_PMD_features[,repli_cells]
norm_edata = normalize.quantiles(as.matrix(edata)) %>% as.data.frame()
colnames(norm_edata) = colnames(edata)
hg19_PMD_features[,repli_cells] = norm_edata

hg19_PMD_features$Replication_timing = rowMeans(hg19_PMD_features[,repli_cells])
hg19_PMD_features$PMD_ID = 1:nrow(hg19_PMD_features)
rownames(hg19_PMD_features) = hg19_PMD_features$PMD_ID
Repli_cluster_plotData = hg19_PMD_features[,c('variation','Replication_timing')]
Repli_cluster_plotData = scale(Repli_cluster_plotData)
Repli_cluster_plotData = as.data.frame(Repli_cluster_plotData)
Repli_cluster_plotData$Replication_timing = -Repli_cluster_plotData$Replication_timing

#clustering
getCluster = function(Repli_cluster_plotData){
  colnames(Repli_cluster_plotData) = c("Methylation\nvariability", "Replication\ntiming")
  
  Data_FigS2D = Repli_cluster_plotData
  colnames(Data_FigS2D) = c("Methylation variability", "Replication timing")
  write.table(Data_FigS2D, file = 'plot/FigS2D/Data_FigS2D.txt', sep = '\t', quote = F, row.names = F, col.names = T)
  
  data = t(Repli_cluster_plotData)
  col_fun = colorRamp2(c(-2.5, 2.5), c("blue", "yellow"))
  col_fun(seq(-3, 3))
  HM = Heatmap(data,
               column_km = 3, column_km_repeats = 100,
               show_column_names = F,
               cluster_rows = F, col = col_fun,
               show_parent_dend_line = FALSE,
               column_title = letters[1:3],
               heatmap_legend_param = list(
                 title = "", at = c(-3, 3), 
                 labels = c("low variability\n/Early", "high variability\n/Late")
               ))
  HM = draw(HM)
  
  c.dend <- column_dend(HM)
  ccl.list <- column_order(HM)
  test = Repli_cluster_plotData
  for (i in 1:length(column_order(HM))){
    if (i == 1) {
      clu <- t(t(row.names(test[column_order(HM)[[i]],])))
      out <- cbind(clu, paste("cluster", i, sep=""))
      colnames(out) <- c("PMD_ID", "Cluster")
    } else {
      clu <- t(t(row.names(test[column_order(HM)[[i]],])))
      clu <- cbind(clu, paste("cluster", i, sep=""))
      out <- rbind(out, clu)
    }
  }
  out = as.data.frame(out)
  cluster1 = filter(out, Cluster == 'cluster1')
  cluster2 = filter(out, Cluster == 'cluster2')
  cluster3 = filter(out, Cluster == 'cluster3')
  
  hg19_PMD_features_insideFun = hg19_PMD_features
  hg19_PMD_features_insideFun$cluster = NA
  hg19_PMD_features_insideFun[cluster1$PMD_ID,]$cluster = 'cluster1'
  hg19_PMD_features_insideFun[cluster2$PMD_ID,]$cluster = 'cluster2'
  hg19_PMD_features_insideFun[cluster3$PMD_ID,]$cluster = 'cluster3'
  scale_Repli = cbind(hg19_PMD_features_insideFun[,'cluster'], Repli_cluster_plotData) %>% dplyr::rename(cluster=1)
  cluster_mean_repli = group_by(scale_Repli,cluster) %>% summarise(Replication_timing = mean(`Replication\ntiming`)) %>% ungroup()
  cluster_mean_repli = arrange(cluster_mean_repli,Replication_timing) %>%
    mutate(PMD_length = c('S','I','L')) %>% 
    arrange(cluster)
  myColumnTitle = cluster_mean_repli$PMD_length
  HM_2 = Heatmap(data,
                 column_km = 3, column_km_repeats = 100,
                 show_column_names = F,
                 cluster_rows = F, col = col_fun,
                 show_parent_dend_line = FALSE, #rm dashed line,
                 column_title = myColumnTitle,
                 heatmap_legend_param = list(
                   title = "", at = c(-3, 3), 
                   labels = c("low variability\n/Early", "high variability\n/Late")
                 ))
  hg19_PMD_features_insideFun = merge(hg19_PMD_features_insideFun,cluster_mean_repli[,c("cluster","PMD_length")],by = 'cluster')
  plotData = hg19_PMD_features_insideFun[,c('d','cluster','PMD_length')]
  plotData$PMD_length = factor(plotData$PMD_length, levels = c('S','I','L'))
  p = ggplot(plotData,aes(x = PMD_length, y = d))+
    geom_boxplot(outlier.size = 0.5,color='black')+
    # geom_jitter(aes(color = PMD_length),size = 0.05, alpha = 0.8)+
    geom_half_point(aes(color = PMD_length),size = 0.05, alpha = 0.8)+
    scale_color_manual(values=c("#92D2C3", "#DDDF98", "#8FBAD9"))+
    theme_classic()+
    guides(color = 'none')+
    xlab('commonPMD class')+
    ylab('Domain size')+
    scale_y_continuous(breaks = c(2500000,7500000),labels = c('2.5Mb','7.5Mb'))
  FigS2E = p
  
  Data_FigS2E = select(plotData, PMD_length, d) %>% 
    dplyr::rename(`commonPMD class` = 1, `Domain size` = 2)
  write.table(Data_FigS2E, file = 'plot/FigS2E/Data_FigS2E.txt', sep = '\t', quote = F, row.names = F, col.names = T)
  
  return(list(hg19_PMD_features_insideFun,HM_2,p))
}
PMD_length_cluster = getCluster(Repli_cluster_plotData)
View(PMD_length_cluster[[1]])

pdf('plot/FigS2D/FigS2D.pdf', width = 8, height = 3)
PMD_length_cluster[[2]]
dev.off()

pdf('plot/FigS2E/FigS2E.pdf', width = 3, height = 4)
PMD_length_cluster[[3]]
dev.off()

table(PMD_length_cluster[[1]]$PMD_length)
#   I   L   S 
# 614 485 759
group_by(PMD_length_cluster[[1]], PMD_length) %>% summarise(whole_size = sum(d))
# # A tibble: 3 × 2
# PMD_length whole_size
# <chr>           <int>
#   1 I           289200000
# 2 L           845500000
# 3 S           178000000

save(PMD_length_cluster,file = 'meta/commonPMD_subgroup/04_PMD_length_Plots/PMD_length_cluster.RData')

#4) save SIL_region_and_annotation
commonPMD_hg19_merge$MD_anno = 'commonPMD'
commonHMD_hg19_merge$MD_anno = 'commonHMD'
commonPMD_or_HMD = rbind(commonPMD_hg19_merge,commonHMD_hg19_merge)

anno_to_intersect = PMD_length_cluster[[1]][,c('chrom','start','end','PMD_length')]
anno_to_intersect$start = anno_to_intersect$start+1 # Otherwise, S and common HMD will have a 1bp overlap, causing common HMD to be annotated as S.
anno_to_intersect$end = anno_to_intersect$end-1
length1 = bed_intersect(commonPMD_or_HMD, anno_to_intersect,suffix = c('','.length')) %>% 
  select(chrom,start,end,MD_anno,PMD_length.length) %>% 
  dplyr::rename(PMD_length_anno = 5)
length2 = bed_intersect(commonPMD_or_HMD, anno_to_intersect, invert = T) %>% 
  select(chrom,start,end,MD_anno) %>% 
  mutate(PMD_length_anno = NA)
SIL_region_and_annotation = rbind(length1,length2) %>% bed_sort() %>% 
  mutate(d_MD = end - start) %>% 
  select(chrom,start,end,d_MD,MD_anno,PMD_length_anno)
write.table(SIL_region_and_annotation, file = 'meta/commonPMD_subgroup/04_PMD_length_Plots/SIL_region_and_annotation.bed',
            sep = '\t', quote = F, row.names = F, col.names = T)

#5) get SIL_hg19 and SIL_hg38 
system("sed '1d' meta/commonPMD_subgroup/04_PMD_length_Plots/SIL_region_and_annotation.bed > meta/commonPMD_commonHMD/SIL_hg19.bed")
system("liftOver meta/commonPMD_commonHMD/SIL_hg19.bed meta/liftover/hg19ToHg38.over.chain.gz meta/commonPMD_commonHMD/SIL_hg38.bed unMapped")


# 3. hg19_RefSeq_gene----
hg19_RefSeq_ucsc = read.table(file = 'meta/RefSeq_gene/hg19_RefSeq.ucsc.txt',
                              sep = '\t',header = T, comment.char = '')
hg19_Refseq_gene = hg19_RefSeq_ucsc[,c('chrom', 'txStart', 'txEnd', 'name2')] %>% 
  dplyr::rename(start = 2, end = 3, gene_name = 4) %>% 
  unique() %>% 
  bed_sort() %>% 
  group_by(gene_name) %>% 
  bed_merge() %>% 
  ungroup() %>% 
  filter(chrom %in% c(paste0('chr', 1:20), 'chrX', 'chrY'))%>% 
  mutate(d = end - start, .after = end) %>% 
  bed_sort()

#add PMD/HMD and S/I/L-PMD annotation
SIL_region_and_annotation = read.table(file = 'meta/commonPMD_subgroup/04_PMD_length_Plots/SIL_region_and_annotation.bed', header = T)
gene_MD = bed_intersect(hg19_Refseq_gene, SIL_region_and_annotation, suffix = c('','.anno')) %>% 
  group_by(chrom, start, end, d, gene_name, MD_anno.anno) %>% 
  summarise(overlap_MD = sum(.overlap), .groups = 'drop') %>% 
  mutate(overlap_MD_ratio = overlap_MD/d) 
gene_MD$MD_anno = as.character(NA)
gene_MD[gene_MD$overlap_MD_ratio > 0.5,]$MD_anno = gene_MD[gene_MD$overlap_MD_ratio > 0.5,]$MD_anno.anno
gene_MD = gene_MD[,c('gene_name','MD_anno')] %>% 
  left_join(hg19_Refseq_gene, ., by = 'gene_name')
gene_MD = arrange(gene_MD, gene_name, MD_anno) %>% 
  distinct(gene_name, .keep_all = T) %>% 
  bed_sort()

gene_SIL = bed_intersect(gene_MD[,1:5], SIL_region_and_annotation, suffix = c('','.anno')) %>% 
  group_by(chrom, start, end, d, gene_name, PMD_length_anno.anno) %>% 
  summarise(overlap_PMD_length = sum(.overlap)) %>% 
  ungroup() %>% 
  mutate(overlap_PMD_length_ratio = overlap_PMD_length/d) 
gene_SIL$PMD_length_anno = as.character(NA)
gene_SIL[gene_SIL$overlap_PMD_length_ratio > 0.5,]$PMD_length_anno = 
  gene_SIL[gene_SIL$overlap_PMD_length_ratio > 0.5,]$PMD_length_anno.anno
gene_SIL = gene_SIL[,c('gene_name','PMD_length_anno')] %>% 
  left_join(gene_MD, ., by = 'gene_name')
gene_SIL = arrange(gene_SIL, gene_name, PMD_length_anno) %>% 
  distinct(gene_name, .keep_all = T) %>% 
  bed_sort()

hg19_Refseq_gene = gene_SIL
save(hg19_Refseq_gene, file = 'meta/RefSeq_gene/hg19_Refseq_gene.RData')

# 4. hg19 promoters (promoter.bed, promoter_with_anno.bed, CGI_promoter.bed, non_CGI_promoter.bed)----
Refseq_transcript = hg19_RefSeq_ucsc[,c('chrom', 'txStart', 'txEnd', 'strand','name2')] %>% 
  dplyr::rename(start = 2, end = 3, gene_name = 5) %>% 
  unique()

transcript = Refseq_transcript
transcript$promoter_start = as.numeric(NA)
transcript$promoter_end = as.numeric(NA)
transcript$TSS = as.numeric(NA)

upstream = 500 # Upstream of the TSS on the sense strand
downstream = 500 # Downstream of the TSS on the sense strand

transcript[transcript$strand=="+",]$promoter_start = transcript[transcript$strand=="+",]$start-upstream
transcript[transcript$strand=="+",]$promoter_end = transcript[transcript$strand=="+",]$start+downstream
transcript[transcript$strand=="-",]$promoter_start = transcript[transcript$strand=="-",]$end-downstream
transcript[transcript$strand=="-",]$promoter_end = transcript[transcript$strand=="-",]$end+upstream
transcript[transcript$strand=='+',]$TSS = transcript[transcript$strand=='+',]$start
transcript[transcript$strand=='-',]$TSS = transcript[transcript$strand=='-',]$end

#Set promoter_start values less than 0 to 0.
transcript[transcript$promoter_start<0,]$promoter_start=0

promoter1 = transcript %>% 
  .[,c('chrom', 'promoter_start', 'promoter_end', 'strand', 'gene_name', 'TSS')] %>% 
  dplyr::rename(start = 2, end = 3) %>% 
  bed_sort()

# Different transcripts of the same gene are merged if they overlap by more than 500 bp.
prom_cluster = bed_cluster(promoter1, max_dist = -500)
promoter_with_anno = group_by(prom_cluster, chrom, gene_name, .id) %>%
  summarise(start = min(start),
            end = max(end), 
            TSSMean = round(mean(TSS), 0),
            TSS = paste(TSS,collapse = ';'))%>% 
  ungroup() %>% 
  select(chrom, start, end, gene_name, TSS, TSSMean, .id) %>% 
  bed_sort() %>% 
  filter(chrom %in% c(paste0('chr',1:20), 'chrX', 'chrY'))
promoter = promoter_with_anno[,1:3]
write.table(promoter, file = 'meta/RefSeq_gene/hg19_promoter.bed', sep = '\t', quote = F, row.names = F, col.names = T)
write.table(promoter_with_anno, file = 'meta/RefSeq_gene/hg19_promoter_with_anno.bed', sep = '\t', quote = F, row.names = F, col.names = T)

# CGI/non-CGI promoters
CGI_hg19 = fread('meta/CGI/CGI_merge_hg19.bed')
colnames(CGI_hg19) = c('chrom', 'start', 'end')
CGI_promoter = bed_intersect(promoter_with_anno, CGI_hg19, suffix = c('', '.y')) %>% 
  select(-start.y, -end.y, -.overlap) %>% 
  distinct()
write.table(CGI_promoter[,1:3], file = 'meta/RefSeq_gene/hg19_CGI_promoter.bed', sep = '\t', quote = F, row.names = F, col.names = F)


# 5. hg38_RefSeq_gene----
hg38_RefSeq_ucsc = read.table(file = 'meta/RefSeq_gene/hg38_RefSeq.ucsc.txt', sep = '\t', header = T, comment.char = '')
hg38_Refseq_gene = hg38_RefSeq_ucsc[,c('chrom', 'txStart', 'txEnd', 'name2')] %>% 
  dplyr::rename(start = 2, end = 3, gene_name = 4) %>% 
  unique() %>% 
  bed_sort() %>% 
  group_by(gene_name) %>% 
  bed_merge() %>% 
  ungroup() %>% 
  filter(chrom %in% c(paste0('chr', 1:22), 'chrX', 'chrY'))
genesWithMultipleSites = names(table(hg38_Refseq_gene$gene_name)[table(hg38_Refseq_gene$gene_name)!=1])
hg38_Refseq_gene = filter(hg38_Refseq_gene, !gene_name %in% genesWithMultipleSites & chrom %in% c(paste0('chr',1:22), 'chrX', 'chrY'))
save(hg38_Refseq_gene, file = 'meta/RefSeq_gene/hg38_Refseq_gene.RData')

setwd('meta/RefSeq_gene')
system('liftOver hg19_promoter_with_anno.bed ../liftover/hg19ToHg38.over.chain.gz hg38_promoter_with_anno.bed')

list.files()
# [1] "hg19_CGI_promoter.bed"       "hg19_promoter_with_anno.bed" "hg19_promoter.bed"          
# [4] "hg19_Refseq_gene.RData"      "hg19_RefSeq.ucsc.txt"        "hg38_promoter_with_anno.bed"
# [7] "hg38_Refseq_gene.RData"      "hg38_RefSeq.ucsc.txt" 

setwd('../..')

# 6. 109 Roadmap samples(Remove ESC/iPSC samples and tumor samples)----
Roadmap_table = read.csv(file = 'meta/Roadmap_table/Roadmap_table.csv',header = T)
tumor_samples = c('E114','E115','E117','E118','E123')
Roadmap_table2 = filter(Roadmap_table,!Group %in% c('ESC','iPSC') & !EID %in% tumor_samples)
Roadmap_table2$summary = paste0(Roadmap_table2$EID,': ',Roadmap_table2$Epigenome.name..from.EDACC.Release.9.directory.) 
Roadmap_table2 = arrange(Roadmap_table2,Order)
write.csv(Roadmap_table2, file = 'meta/Roadmap_LOCK/data/00_Roadmap_table/Roadmap_table2.csv')

# PART 2: Roadmap data preprocessing----
# 1. ChIP-seq data of the Roadmap samples----
# 1.1 narrowPeak----
# For the narrowPeak files from Roadmap: keep the chromosomes of chr1 to chr22 and remove the peaks overlapping with the blacklist(hg19 version)
list.files('data_Roadmap/01_narrowPeak_sorted')
# [1] "H3K27me3" "H3K4me3"  "H3K9me3"
head(list.files('data_Roadmap/01_narrowPeak_sorted/H3K27me3'))
# [1] "E004_H3K27me3_peak.bed" "E005_H3K27me3_peak.bed" "E006_H3K27me3_peak.bed" "E007_H3K27me3_peak.bed" "E009_H3K27me3_peak.bed"
# [6] "E010_H3K27me3_peak.bed"

# 1.2 BigWig----
head(list.files('data_Roadmap/05_peak_signal/H3K27me3/'))
# [1] "E004_H3K27me3PeakSignal.gz" "E005_H3K27me3PeakSignal.gz" "E006_H3K27me3PeakSignal.gz" "E007_H3K27me3PeakSignal.gz" "E009_H3K27me3PeakSignal.gz" "E010_H3K27me3PeakSignal.gz"

# 1.3 Identification of LOCKs----
Roadmap_table2 = read.csv(file = 'meta/Roadmap_table/Roadmap_table2.csv')
system("ls data_Roadmap/01_narrowPeak_sorted/H3K27me3/ | sed 's/_H3K27me3_peak.bed/ H3K27me3/' | xargs -iFile -P20 bash -c \"Rscript script/03_LOCK_calling.R File\"")
system("ls data_Roadmap/01_narrowPeak_sorted/H3K9me3/ | sed 's/_H3K9me3_peak.bed/ H3K9me3/' | xargs -iFile -P20 bash -c \"Rscript script/03_LOCK_calling.R File\"")
head(list.files('data_Roadmap/03_LOCK/H3K27me3/E004'))
# [1] "E004_H3K27me3_ctoff-0.1.bed" "E004_H3K27me3_ctoff-0.2.bed" "E004_H3K27me3_ctoff-0.3.bed" "E004_H3K27me3_ctoff-0.4.bed" "E004_H3K27me3_ctoff-0.5.bed" "E004_H3K27me3_ctoff-0.6.bed"

# 1.4 annotatePeaks----
system("ls data_Roadmap/01_narrowPeak_sorted/H3K27me3/ | sed 's/_H3K27me3_peak.bed/ H3K27me3/' | xargs -iFile -P20 bash -c \"sh script/04_annotatePeaks.pl.sh File\"")
system("ls data_Roadmap/01_narrowPeak_sorted/H3K9me3/ | sed 's/_H3K9me3_peak.bed/ H3K9me3/' | xargs -iFile -P20 bash -c \"sh script/04_annotatePeaks.pl.sh File\"")
head(list.files('data_Roadmap//04_annotatePeaks.pl/H3K27me3'))
# [1] "E004_annotatePeaks.txt" "E005_annotatePeaks.txt" "E006_annotatePeaks.txt" "E007_annotatePeaks.txt" "E009_annotatePeaks.txt" "E010_annotatePeaks.txt"

# 2. WGBS data of the Roadmap samples----
# 2.1 WGBS bw (not used)----
head(list.files('data_Roadmap/06_WGBS_bw/WGBS_bw/'))
# [1] "E006_WGBS.bw" "E050_WGBS.bw" "E058_WGBS.bw" "E065_WGBS.bw" "E066_WGBS.bw" "E079_WGBS.bw"

# 2.2 WGBS bdg (not used)----
WBGS_EIDs = sub('_WGBS.bw','',list.files('data_Roadmap/06_WGBS_bw/WGBS_bw/'))
for (EID in WBGS_EIDs$EID) {
  system(paste0('bigWigToBedGraph data_Roadmap/06_WGBS_bw/WGBS_bw/',EID,'_WGBS.bw data_Roadmap/06_WGBS_bw/WGBS_bdg/',EID,'_WGBS.bdg'))
}

# 2.3 WGBS bdg (remove CGI)----
for (EID in WBGS_EIDs) {
  system(paste0('bedtools subtract -a data_Roadmap/06_WGBS_bw/WGBS_bdg/',EID,'_WGBS.bdg -b meta/CGI/CGI_merge_hg19.bed > data_Roadmap/06_WGBS_bw/rmCGI_WGBS_bdg/',EID,'_WGBS_subtrCGI.bdg'))
  system(paste0('bedtools sort -i data_Roadmap/06_WGBS_bw/rmCGI_WGBS_bdg/',EID,'_WGBS_subtrCGI.bdg > data_Roadmap/06_WGBS_bw/rmCGI_WGBS_bdg/',EID,'_WGBS_subtrCGI.sort.bdg'))
  system(paste0('rm data_Roadmap/06_WGBS_bw/rmCGI_WGBS_bdg/',EID,'_WGBS_subtrCGI.bdg'))
}

# 2.4 WGBS bw (remove CGI)----
# hg19.chrom.sizes is from UCSC database
hg19.chrom.sizes = read.table(file = 'meta/chrom_size/hg19.chrom.sizes')
hg19_ly.chrom.sizes = filter(hg19.chrom.sizes, V1 %in% c(paste0('chr',1:22),'chrX','chrY'))
write.table(hg19_ly.chrom.sizes, file = paste0('meta/hg19_ly.chrom.sizes'), sep = '\t', quote = F, row.names = F, col.names = F)
for (EID in WBGS_EIDs) {
  system(paste0('bedGraphToBigWig data_Roadmap/WGBS_bw/rmCGI_WGBS_bdg/',EID,'_WGBS_subtrCGI.sort.bdg meta/chrom_size/hg19_ly.chrom.sizes data_Roadmap/WGBS_bw/rmCGI_WGBS_bw/',EID,'_WGBS_subtrCGI.bw'))
}

# 3. Calculate the peak intensity----
for (i in 1:nrow(Roadmap_table2)) {
  system(paste0('sh script/05_peak_signal.sh ',Roadmap_table2$EID[i],' H3K27me3'))
}

# 4. Calculate the DNA methylation levels of the peaks and LOCKs----
# 1) peak
WGBS_samples = sub('_.*','',list.files('data_Roadmap/06_WGBS_bw/rmCGI_WGBS_bw/'))
for (i in 1:length(WGBS_samples)) {
  system(paste0('sh script/06_peak_methylation.sh ',WGBS_samples[i]))
}

# 2) LOCK
for (i in 1:length(WGBS_samples)) {
  system(paste0('sh script/07_LOCK_methylation.sh ',WGBS_samples[i]))
}

# 5. RDatas (For ease of subsequent analysis, some data has been pre-processed and saved.)----
# 5.1 samples.RData----
Roadmap_table2 = read.csv(file = 'meta/Roadmap_table/Roadmap_table2.csv')
Roadmap_WGBS_EID = substr(list.files('/data3/liangyuan/03_Roadmap_LOCK/data/06_WGBS_bw/soloWCGW_rmCGI_bw'),1,4)
save(Roadmap_table2, Roadmap_WGBS_EID, file = 'data_Roadmap//RDatas/samples.RData')

# 5.2 peak_109----
# 1) get peak signal data
path = paste0('data_Roadmap/05_peak_signal/H3K27me3//')
peak_109 = list()
for (s109 in 1:length(list.files(path))) {
  print(s109)
  computeMatrix_result = read.table(file = paste0(path, Roadmap_table2$EID[s109],'_H3K27me3PeakSignal.gz'),sep = '\t',header = F,skip = 1)
  peak_109[[s109]] = data.frame(chrom = computeMatrix_result$V1,
                                    start = computeMatrix_result$V2,
                                    end = computeMatrix_result$V3,
                                    d = computeMatrix_result$V3 - computeMatrix_result$V2,
                                    signal = rowMeans(computeMatrix_result[,7:26],na.rm = T))
  
}
names(peak_109) = Roadmap_table2$EID

# 2) add annotatePeaks.pl result
for (s109 in 1:length(peak_109)) {
  print(s109)
  annotatePeaks = data.table::fread(paste0('data_Roadmap/04_annotatePeaks.pl/H3K27me3/',Roadmap_table2[s109,'EID'],'_annotatePeaks.txt'))[,c(2:4,8,15,16)]
  colnames(annotatePeaks) = c('chrom','start','end','Annotation','Nearest_Ensembl','Nearest_Gene_Name')
  annotatePeaks$PeakAnnotation = paste(annotatePeaks$Annotation,annotatePeaks$Nearest_Ensembl,annotatePeaks$Nearest_Gene_Name,sep = '///')
  annotatePeaks = annotatePeaks[,c("chrom","start","end","PeakAnnotation")]
  peak_109[[s109]] = bed_intersect(peak_109[[s109]],annotatePeaks,suffix = c('','.y'))[,-c(6,7,9)] %>% dplyr::rename(PeakAnnotation = 6)
}

# 3) Annotate the methylation domains to which the peak belongs
MD_and_PMD_length_anno = function(df){
  df = mutate(df, ID = paste0(chrom,':',start,'-',end))
  #MD_anno
  df_MD_anno = df %>% 
    bed_intersect(.,SIL_region_and_annotation,suffix = c('','.y')) %>%
    group_by(ID, d, MD_anno.y) %>% 
    summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
    filter(.overlap > d/2) %>% 
    select(ID, MD_anno.y) %>% 
    dplyr::rename(MD_anno = MD_anno.y)
  df = left_join(df, df_MD_anno, by = 'ID')
   
  
  #SIL_anno
  df_MD_anno = df %>% 
    bed_intersect(.,SIL_region_and_annotation,suffix = c('','.y')) %>%
    group_by(ID, d, PMD_length_anno.y) %>% 
    summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
    filter(.overlap > d/2) %>% 
    select(ID, PMD_length_anno.y) %>% 
    dplyr::rename(PMD_length_anno = PMD_length_anno.y)
  df2 = left_join(df, df_MD_anno, by = 'ID')
  return(df2)
}

for (s109 in 1:length(peak_109)) {
  print(s109)
  peak_109[[s109]] = MD_and_PMD_length_anno(peak_109[[s109]])
}

save(peak_109 ,file = 'data_Roadmap/RDatas/peak_109.RData')

# 5.3 LOCK_109 (annotations for all LOCKs in 109 Roadmap samples)----
# 1) prepare for gene size annotation
gene = read.table(file = 'meta/RefSeq_gene/hg19_RefSeq.ucsc.txt', header = T, comment.char = " ")
gene = select(gene, chrom, txStart, txEnd)
colnames(gene) = c('chrom','start','end')
gene_merge = bed_merge(gene[,1:3])

# 2) get LOCK_109
for (s109 in 1:nrow(Roadmap_table2)) {
  print(s109)
  EID = Roadmap_table2$EID[s109]
  LOCK = read.table(file = paste0('data_Roadmap/03_LOCK/H3K27me3/',EID,'/',EID,'_H3K27me3_ctoff-0.7.bed'))
  colnames(LOCK) = c('chrom','start','end')
  LOCK$d = LOCK$end - LOCK$start
  LOCK$ID = paste0(LOCK$chrom,':',LOCK$start,'-',LOCK$end)
  
  # add MD_anno
  df_MD_anno = LOCK %>% 
    bed_intersect(.,SIL_region_and_annotation,suffix = c('','.y')) %>%
    group_by(ID, d, MD_anno.y) %>% 
    summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
    filter(.overlap > d/2) %>% 
    select(ID, MD_anno.y) %>% 
    dplyr::rename(MD_anno = MD_anno.y)
  LOCK = left_join(LOCK, df_MD_anno, by = 'ID')
  
  
  # add SIL_anno
  df_MD_anno = LOCK %>% 
    bed_intersect(.,SIL_region_and_annotation,suffix = c('','.y')) %>%
    group_by(ID, d, PMD_length_anno.y) %>% 
    summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
    filter(.overlap > d/2) %>% 
    select(ID, PMD_length_anno.y) %>% 
    dplyr::rename(PMD_length_anno = PMD_length_anno.y)
  LOCK = left_join(LOCK, df_MD_anno, by = 'ID')
  
  # add peak size
  add_size = bed_intersect(LOCK, peak_109[[s109]], suffix = c('','.y')) %>% 
    group_by(ID) %>% 
    summarise(d_peak = sum(d.y), .groups = 'drop')
  LOCK = left_join(LOCK, add_size, by = 'ID')
  
  # add LOCK_length
  LOCK = dplyr::rename(LOCK, d_LOCK = d)
  LOCK = mutate(LOCK, LOCK_length = ifelse(d_LOCK > 100000, 'Long_LOCK', 'Short_LOCK'))
  
  # add gene_size
  add_geneSize = bed_intersect(LOCK, gene_merge, suffix = c('','.y')) %>% 
    group_by(ID) %>% 
    summarise(geneSize = sum(.overlap), .groups = 'drop')
  LOCK = left_join(LOCK, add_geneSize, by = 'ID')
  LOCK_109[[s109]] = LOCK
}
names(LOCK_109) = Roadmap_table2$EID
save(LOCK_109, file = 'data_Roadmap/RDatas/LOCK_109.RData')

# 5.4 LOCK_peak_109 (annotations for all peaks in 109 Roadmap samples)----
for (s109 in 1:nrow(Roadmap_table2)) {
  print(s109)
  LOCK_peak = bed_intersect(peak_109[[s109]], LOCK_109[[s109]], suffix = c('','.LOCK'))
  LOCK_peak$LOCK_length.peak = NA
  LOCK_peak[LOCK_peak$LOCK_length.LOCK == 'Long_LOCK',]$LOCK_length.peak = 'Peaks in long LOCKs'
  LOCK_peak[LOCK_peak$LOCK_length.LOCK == 'Short_LOCK',]$LOCK_length.peak = 'Peaks in short LOCKs'
  LOCK_peak$.overlap = NULL
  LOCK_peak$.source = NULL
  
  LOCK_peak = left_join(peak_109[[s109]], LOCK_peak, by = colnames(peak_109[[s109]]))
  LOCK_peak[is.na(LOCK_peak$LOCK_length.peak),]$LOCK_length.peak = 'Typical peaks'
  
  LOCK_peak_109[[s109]] = LOCK_peak
}
names(LOCK_peak_109) = Roadmap_table2$EID
save(LOCK_peak_109, file = 'data_Roadmap/RDatas/LOCK_peak_109.RData')

# 5.5 LOCK_109_H3K9me3 and LOCK_peak_109_H3K9me3----
LOCK_109_H3K9me3 = list()
LOCK_peak_109_H3K9me3 = list()
for (s109 in 1:nrow(Roadmap_table2)) {
  print(s109)
  EID = Roadmap_table2$EID[s109]
  peak = read.table(file = paste0('data_Roadmap/01_narrowPeak_sorted/H3K9me3/',EID,'_H3K9me3_peak.bed'))
  colnames(peak) = c('chrom','start','end')
  peak$d = peak$end - peak$start
  peak$ID = paste0(peak$chrom,':',peak$start,'-',peak$end)
  
  LOCK = read.table(file = paste0('data_Roadmap/03_LOCK/H3K9me3/',EID,'/',EID,'_H3K9me3_ctoff-0.7.bed'))
  colnames(LOCK) = c('chrom','start','end')
  LOCK$d = LOCK$end - LOCK$start
  LOCK$LOCK_length = 'Short_LOCK'
  LOCK[LOCK$d > 100000,]$LOCK_length = 'Long_LOCK'
  LOCK = dplyr::rename(LOCK, d_LOCK = d)
  
  LOCK_peak = bed_intersect(peak, LOCK, suffix = c('','.LOCK'))
  LOCK_peak$.overlap = NULL
  LOCK_peak$LOCK_length.peak = NA
  LOCK_peak[LOCK_peak$LOCK_length.LOCK == 'Long_LOCK',]$LOCK_length.peak = 'Peaks in long LOCKs'
  LOCK_peak[LOCK_peak$LOCK_length.LOCK == 'Short_LOCK',]$LOCK_length.peak = 'Peaks in short LOCKs'
  LOCK_peak = left_join(peak, LOCK_peak, by = colnames(peak))
  LOCK_peak[is.na(LOCK_peak$LOCK_length.peak),]$LOCK_length.peak = 'Typical peaks'
  LOCK_peak = bed_sort(LOCK_peak)
  
  LOCK_109_H3K9me3[[s109]] = LOCK
  LOCK_peak_109_H3K9me3[[s109]] = LOCK_peak
}
names(LOCK_109_H3K9me3) = Roadmap_table2$EID
names(LOCK_peak_109_H3K9me3) = Roadmap_table2$EID

save(LOCK_109_H3K9me3, file = 'data_Roadmap/RDatas/LOCK_109_H3K9me3.RData')
save(LOCK_peak_109_H3K9me3, file = 'data_Roadmap/RDatas/LOCK_peak_109_H3K9me3.RData')

# 5.6 LOCK_peak_20 (annotations for all peaks in 20 Roadmap samples with WGBS data)----
load('data_Roadmap/RDatas/samples.RData') #Includes Roadmap_WGBS_EID and Roadmap_table2 objects
LOCK_peak_20 = list()
for (s20 in 1:20) {
  print(s20)
  EID = Roadmap_WGBS_EID[s20]
  x = which(names(LOCK_peak_109) == EID)
  LOCK_peak = LOCK_peak_109[[x]]
  
  computeMatrix_result = read.table(file = paste0(path, EID,'_H3K27me3PeakMetSignal.gz'),sep = '\t',header = F,skip = 1)
  methyl_siganal = data.frame(ID = paste0(computeMatrix_result$V1,':',computeMatrix_result$V2,'-',computeMatrix_result$V3),
                              metSignal = rowMeans(computeMatrix_result[,7:26],na.rm = T))
  LOCK_peak_20[[s20]] = left_join(LOCK_peak, methyl_siganal, by = 'ID')
}
names(LOCK_peak_20) = Roadmap_WGBS_EID
save(LOCK_peak_20, file = 'data_Roadmap/RDatas/LOCK_peak_20.RData')

# 5.7 LOCK_20----
LOCK_20 = list()
for (s20 in 1:length(Roadmap_WGBS_EID)) {
  print(s20)
  EID = Roadmap_WGBS_EID[s20]
  x = which(names(LOCK_109) == EID)
  LOCK = LOCK_109[[x]]
  computeMatrix_LOCKMet = fread(file = paste0('data_Roadmap/08_LOCK_DNA_methylation/H3K27me3_LOCK_met_signal/',
                                              EID,'_-0.7_LOCK_met_signal.gz'),skip = 1)[,c(1:3,7)]
  
  methyl_siganal = data.frame(ID = paste0(computeMatrix_LOCKMet$V1,':',computeMatrix_LOCKMet$V2,'-',computeMatrix_LOCKMet$V3),
                              LOCK_met_signal = computeMatrix_LOCKMet$V7)
  LOCK = left_join(LOCK, methyl_siganal, by = 'ID')
  LOCK_20[[s20]] = LOCK
}
names(LOCK_20) = Roadmap_WGBS_EID
save(LOCK_20, file = 'data_Roadmap/RDatas/LOCK_20.RData')

# 5.9 LOCK_peak_exp_49 (annotations for all peaks in 20 Roadmap samples with RNA-seq data)----
RNASeq = read.table(file = 'data_Roadmap/09_RNASeq/57epigenomes.RPKM.pc.gz',header = F,sep = '\t',skip = 1)[,1:58]
header = read.table(file = 'data_Roadmap/09_RNASeq/57epigenomes.RPKM.pc.gz',header = F,sep = '\t',nrows = 1)
colnames(RNASeq) = header
RNASeq = column_to_rownames(RNASeq,'gene_id')
RNASeq = RNASeq[,colnames(RNASeq) %in% Roadmap_table2$EID]
log_RPKM = log10(RNASeq+0.1)
save(log_RPKM, file = 'data_Roadmap/RDatas//log_RPKM.RData')

LOCK_peak_exp_49 = list()
exp_sample = colnames(log_RPKM)
for (s49 in 1:49) {
  print(s49)
  EID = exp_sample[s49]
  x = which(Roadmap_table2$EID == EID)
  df = LOCK_peak_109[[x]]
  gene_id = str_match(df$PeakAnnotation, '.*///(.*?)///.*$')[,2]
  df = mutate(df, exp = log_RPKM[gene_id,EID])
  LOCK_peak_exp_49[[s49]] = df
}
names(LOCK_peak_exp_49) = exp_sample
save(LOCK_peak_exp_49,file = 'data_Roadmap/RDatas/LOCK_peak_exp_49.RData')

# 5.10 overall DNA methylation of each methylation domain----
# 1) save region
SIL_region_and_annotation = read.table(file = 'meta/commonPMD_subgroup/04_PMD_length_Plots/SIL_region_and_annotation.bed', sep = '\t', header = T)
SIL_region_and_annotation[is.na(SIL_region_and_annotation$PMD_length_anno),]$PMD_length_anno = SIL_region_and_annotation[is.na(SIL_region_and_annotation$PMD_length_anno),]$MD_anno
SIL_region_and_annotation = SIL_region_and_annotation[,c('chrom','start','end','PMD_length_anno')]

S_PMD = filter(SIL_region_and_annotation, PMD_length_anno == 'S')
I_PMD = filter(SIL_region_and_annotation, PMD_length_anno == 'I')
L_PMD = filter(SIL_region_and_annotation, PMD_length_anno == 'L')
HMD = filter(SIL_region_and_annotation, PMD_length_anno == 'commonHMD')

write.table(S_PMD[,1:3], file = 'data_Roadmap/08_overall_DNA_methylation/SIL_overall_met/bed/S_PMD.bed', sep = '\t',quote = F,row.names = F,col.names = F)
write.table(I_PMD[,1:3], file = 'data_Roadmap/08_overall_DNA_methylation/SIL_overall_met/bed/I_PMD.bed', sep = '\t',quote = F,row.names = F,col.names = F)
write.table(L_PMD[,1:3], file = 'data_Roadmap/08_overall_DNA_methylation/SIL_overall_met/bed/L_PMD.bed', sep = '\t',quote = F,row.names = F,col.names = F)
write.table(HMD[,1:3], file = 'data_Roadmap/08_overall_DNA_methylation/SIL_overall_met/bed/HMD.bed', sep = '\t',quote = F,row.names = F,col.names = F)

# 2) run computeMatrix
for (i in 1:length(Roadmap_WGBS_EID)) {
  system(paste0('sh script/08_overall_methylation.sh ',Roadmap_WGBS_EID[i]))
}

# 3) region annotation of computeMatrix result
EID = Roadmap_WGBS_EID[1]
computeMatrix_res = read.table(file = paste0('data_Roadmap/08_overall_DNA_methylation/SIL_overall_met/gz/',EID,'_SIL_met_signal.gz'),sep = '\t',header = F,skip = 1)[,c(1:3,7)]
colnames(computeMatrix_res) = c('chrom','start','end','metSignal')
computeMatrix_FirstLine = read.table(file = paste0('data_Roadmap/08_overall_DNA_methylation/SIL_overall_met/gz/',EID,'_SIL_met_signal.gz'),sep = '\t',header = F,nrows = 1)
group_boundaries = str_match(computeMatrix_FirstLine[1,1],'.*group_boundaries:\\[(.*?)\\].*')[,2]
group_boundaries = as.numeric(strsplit(group_boundaries,split = ',')[[1]])
group_n = group_boundaries[2:5]-group_boundaries[1:4]
SIL_anno_col = c(rep('S',group_n[1]),
                 rep('I',group_n[2]),
                 rep('L',group_n[3]),
                 rep('commonHMD',group_n[4]))

# 4) read and save computeMatrix result
computeMatrix_res = list()
for (s in 1:20) {
  print(s)
  EID = Roadmap_WGBS_EID[s]
  df = read.table(file = paste0('08_overall_DNA_methylation/SIL_overall_met/gz/',EID,'_SIL_met_signal.gz'),sep = '\t',header = F,skip = 1)[,c(1:3,7)] %>% 
    rename(chrom=1,start=2,end=3,metSignal=4) %>% 
    mutate(SIL_anno = SIL_anno_col,sample = EID)
  df_PMD = data.frame(sample = EID,
                      SIL_anno = 'commonPMD',
                      metSignal = median(df[df$SIL_anno != 'commonHMD','metSignal',drop = T], na.rm = T))
  computeMatrix_res[[s]] = group_by(df,sample,SIL_anno) %>% 
    summarise(metSignal = median(metSignal,na.rm = T),.groups = 'drop') %>% 
    rbind(df_PMD,.) 
}
SIL_met_20 = do.call(rbind, computeMatrix_res)
SIL_met_20$SIL_anno = factor(SIL_met_20$SIL_anno, levels = c('S','I','L','commonPMD','commonHMD'))
save(SIL_met_20,file = 'data_Roadmap/RDatas/SIL_met_20.RData')

# 5.11 gene_20----
load('meta/RefSeq_gene/hg19_Refseq_gene.RData')
load('data_Roadmap/RDatas/LOCK_20.RData')
load(file = 'data_Roadmap/RDatas/LOCK_peak_20.RData')

gene_20 = list()
for (s20 in 1:length(Roadmap_WGBS_EID)) {
  print(s20)
  EID = Roadmap_WGBS_EID[s20]
  
  # Annotate if the gene is in LOCK
  LOCK_df = LOCK_20[[s20]][,c('chrom', 'start', 'end', 'LOCK_length')]
  gene_df = mutate(hg19_Refseq_gene, ID = paste0(chrom,':',start,'-',end))
  gene_LOCK = bed_intersect(gene_df, LOCK_df, suffix = c('','.y')) %>% 
    group_by(ID, d, LOCK_length.y) %>% 
    summarise(overlap = sum(.overlap), .groups = 'drop') %>% 
    filter(overlap > d/2) %>% 
    select(ID, LOCK_length.y) %>% 
    dplyr::rename(LOCK_length = LOCK_length.y)
  df = left_join(gene_df, gene_LOCK, by = 'ID')
  
  # Annotate if the gene is a peak nearest gene
  df_peak = LOCK_peak_20[[s20]]
  
  GenePeakAnnoTable = df_peak[,c('PeakAnnotation','LOCK_length.peak')] %>% 
    mutate(nearestGene = str_match(PeakAnnotation,'.*///.*///(.*)')[,2]) %>% 
    select(-PeakAnnotation) %>% 
    unique() 
  whether_one_group = table(GenePeakAnnoTable$nearestGene) == 1
  geneList = names(whether_one_group[whether_one_group == T])
  GenePeakAnnoTable = filter(GenePeakAnnoTable,nearestGene %in% geneList)
  
  gene_20[[s20]] = left_join(df, GenePeakAnnoTable, by = c('gene_name'= 'nearestGene'))
}
names(gene_20) = Roadmap_WGBS_EID
save(gene_20, file = 'data_Roadmap/RDatas/gene_20.RData')

# 5.12. gene_109----
load('meta/RefSeq_gene/hg19_Refseq_gene.RData')
load('data_Roadmap/RDatas/LOCK_109.RData')
load(file = 'data_Roadmap/RDatas/LOCK_peak_109.RData')

gene_109 = list()
for (s109 in 1:109) {
  print(s109)
  EID = Roadmap_WGBS_EID[s109]
  
  # Annotate if the gene is in LOCK
  LOCK_df = LOCK_109[[s109]][,c('chrom', 'start', 'end', 'LOCK_length')]
  gene_df = mutate(hg19_Refseq_gene, ID = paste0(chrom,':',start,'-',end))
  gene_LOCK = bed_intersect(gene_df, LOCK_df, suffix = c('','.y')) %>% 
    group_by(ID, d, LOCK_length.y) %>% 
    summarise(overlap = sum(.overlap), .groups = 'drop') %>% 
    filter(overlap > d/2) %>% 
    select(ID, LOCK_length.y) %>% 
    dplyr::rename(LOCK_length = LOCK_length.y)
  df = left_join(gene_df, gene_LOCK, by = 'ID')
  
  # Annotate if the gene is a peak nearest gene
  df_peak = LOCK_peak_109[[s109]]
  
  GenePeakAnnoTable = df_peak[,c('PeakAnnotation','LOCK_length.peak')] %>% 
    mutate(nearestGene = str_match(PeakAnnotation,'.*///.*///(.*)')[,2]) %>% 
    select(-PeakAnnotation) %>% 
    unique() 
  whether_one_group = table(GenePeakAnnoTable$nearestGene) == 1
  geneList = names(whether_one_group[whether_one_group == T])
  GenePeakAnnoTable = filter(GenePeakAnnoTable,nearestGene %in% geneList)
  
  gene_109[[s109]] = left_join(df, GenePeakAnnoTable, by = c('gene_name'= 'nearestGene'))
}
names(gene_109) = Roadmap_table2$EID
save(gene_109, file = 'data_Roadmap/RDatas/gene_109.RData')

# 5.13. gene_49----
# 1) Convert gene_id to gene_name for log_RPKM
load('data_Roadmap/RDatas/log_RPKM.RData')
log_RPKM = log_RPKM[,colnames(log_RPKM) %in% Roadmap_table2$EID]
exp_sample = colnames(log_RPKM)

gtf_hg19 = read.table('meta/gtf/gencode.v19.annotation.gtf.gz',comment.char = '#', header = F, sep = '\t')
gene_name_anno = filter(gtf_hg19, V3 == 'gene')
gene_name_anno = data.frame(gene_id = str_match(gene_name_anno$V9, pattern = '.*gene_id (.*?)\\..*')[,2],
                            gene_name = str_match(gene_name_anno$V9, pattern = '.*gene_name (.*?);.*')[,2]) %>% unique() %>% 
  column_to_rownames('gene_id') %>% 
  .[rownames(log_RPKM),]

log_RPKM = cbind(gene_name_anno, log_RPKM)
colnames(log_RPKM)[1] = 'gene_name'
rownames(log_RPKM) = NULL
log_RPKM = group_by(log_RPKM, gene_name) %>% 
  summarise_all(mean) %>% 
  ungroup()

# 2) get gene_49
for (s49 in 1:49) {
  print(s49)
  EID = exp_sample[s49]
  x = which(names(gene_109) == EID)
  df = gene_109[[x]]
  exp_df = log_RPKM[,c('gene_name',EID)] %>% 
      dplyr::rename(exp = 2)
  gene_49[[s49]] = left_join(df, exp_df, by = 'gene_name')
}
names(gene_49) = exp_sample
save(gene_49,file = 'data_Roadmap/RDatas/gene_49.RData')

# 6. Metascape----
# Metascape was used for gene function annotation, which was performed within a Docker container.
pathOfDockerContainer = '/data1/liangyuan/linux_deal/software/MSBio'
# 1) input
load('data_Roadmap/RDatas/gene_109.RData')
for (s109 in 1:nrow(Roadmap_table2)) {
  print(s109)
  EID = Roadmap_table2$EID[s109]
  df = gene_109[[s109]]
  
  Long_LOCK = unique(filter(df, LOCK_length == 'Long_LOCK' & LOCK_length.peak == 'Peaks in long LOCKs')$gene_name)
  Short_LOCK = unique(filter(df, LOCK_length == 'Short_LOCK' & LOCK_length.peak == 'Peaks in short LOCKs')$gene_name)
  Typical_peaks = unique(filter(df, LOCK_length.peak == 'Typical peaks')$gene_name)
  
  gene_number = tribble(
    ~sample, ~group, ~gene_number,
    EID,'Long_LOCK',length(Long_LOCK),
    EID,'Short_LOCK',length(Short_LOCK),
    EID,'Typical_peaks',length(Typical_peaks)
  )
  
  if(length(Long_LOCK) > 3000){set.seed(1314); Long_LOCK = sample(Long_LOCK,3000)}
  if(length(Short_LOCK) > 3000){set.seed(1314); Short_LOCK = sample(Short_LOCK,3000)}
  if(length(Typical_peaks) > 3000){set.seed(1314); Typical_peaks = sample(Typical_peaks,3000)}
  
  Metascape_input = data.frame(names = c('Long_LOCK','Short_LOCK','Typical_peaks'),
                               genes = c(paste(Long_LOCK,collapse = ','),
                                         paste(Short_LOCK,collapse = ','),
                                         paste(Typical_peaks,collapse = ',')))
  
  write.table(Metascape_input, file = paste0('data_Roadmap/10_metascape_and_LOLA/Metascape_input/H3K27me3_',EID,'.txt'), 
              sep = '\t', quote = F, row.names = F, col.names = T)
}
# system(paste0('cp -r data_Roadmap/10_metascape_and_LOLA/Metascape_input ',
#               paste0(pathOfDockerContainer,'data/20221124_1/batch_metascape_data')))

# 2) run Metascape
setwd(pathOfDockerContainer)
list.files()
# [1] "bin"                   "data"                  "intersection_gene.log" "license"               "Metascatpe.sh"        
# [6] "run.sh"                "winbin" 
system('cat run.sh')
# #/bin/sh
# for file in `ls data/20221124_1/batch_metascape_data/`
# do
# echo ${file}
# ./Metascatpe.sh ${file}
# done
system('cat Metascatpe.sh')
# #!/bin/sh
# file=$1
# bin/up.sh 1
# bin/ms.sh 1 -o /data/20221124_1/batch_metascape_result/${file%%.txt} /data/20221124_1/batch_metascape_data/${file}
# bin/down.sh 1

setwd('/data3/liangyuan/05_LOCK_and_PMD/data_Roadmap/10_metascape_and_LOLA')
system('ln -s  /data1/liangyuan/linux_deal/software/MSBio/data/20221212_1/batch_metascape_result ./Metascape_output')
setwd('/data3/liangyuan/05_LOCK_and_PMD')

# 7. LOLA----
library(LOLA)
regionDB = loadRegionDB('meta/LOLACore/hg19')
lolaResults_peak_level = list()
for (s109 in 1:nrow(Roadmap_table2)) {
  print(s109)
  df = LOCK_peak_109[[s109]]
  EID = names(LOCK_peak_109)[s109]
  df = df[,c("chrom","start","end","MD_anno","PMD_length_anno","LOCK_length.peak")]
  Peaks_in_long_LOCKs = filter(df, LOCK_length.peak == 'Peaks in long LOCKs')
  Peaks_in_short_LOCKs = filter(df, LOCK_length.peak == 'Peaks in short LOCKs')
  Typical_peaks = filter(df, LOCK_length.peak == 'Typical peaks')
  
  getGRanges <- function(bed){GRanges(seqnames = Rle(bed[,1,drop = T]),ranges = IRanges(bed[,2,drop = T],bed[,3,drop = T]))}
  userSets = GRangesList(
    getGRanges(Peaks_in_long_LOCKs),
    getGRanges(Peaks_in_short_LOCKs),
    getGRanges(Typical_peaks)
  )
  userUniverse = getGRanges(df)
  group_to_merge = data.frame(userSet = 1:3,
                              group = c('Peaks_in_long_LOCKs','Peaks_in_short_LOCKs','Typical_peaks'))
  lolaResults_peak_level[[s109]] = runLOLA(userSets, userUniverse, regionDB, cores=1) %>% 
    mutate(sample = EID) %>% 
    merge(.,group_to_merge,by = 'userSet')
  
}
names(lolaResults_peak_level) = Roadmap_table2$EID
save(lolaResults_peak_level,file = 'data_Roadmap/10_metascape_and_LOLA/lolaResults_peak_level.RData')
# PART 3 Cell line data preprocessing----
# 1. samples----
ESCC_samples = read.csv(file = 'data_CellLine/01_sample_infomation/ESCC_samples.csv', check.names = F)
BRCA_samples = read.csv(file = 'data_CellLine/01_sample_infomation/BRCA_samples.csv', check.names = F)
HNSC_samples = read.csv(file = 'data_CellLine/01_sample_infomation/HNSC_samples.csv', check.names = F)

# 2. pairs----
pairs_for_longLOCKAnalysis = c('NE2_KYSE450','NE2_KYSE510','breastNormalVsHer2Pos','breastNormalVsLuminalA','breastNormalVsTN_Basal','breastNormalVsTN_CL')
pairs_for_shortLOCKAnalysis = c('NE2_KYSE450','NE2_KYSE510','MCF10A_normal_HCC1954_Her2Pos','MCF10A_normal_MCF7_LuminalA','MCF10A_normal_HCC1937_TN_Basal','MCF10A_normal_MB231_TN_CL')
all_pairs = data.frame(cancer = c(rep('ESCC',2), rep('BRCA',8)),
                       pairs = unique(c(pairs_for_longLOCKAnalysis, pairs_for_shortLOCKAnalysis)))
save(pairs_for_longLOCKAnalysis, pairs_for_shortLOCKAnalysis, all_pairs, file = 'data_CellLine/RDatas/Pairs.RData')

# 3. fastq files----
# From .sra files to .fastq files, as shown in the code：'script/09_sra_to_fastq.sh'
list.files('data_CellLine/02_fastq/')
# [1] "ChIPSeq"         "RNASeq"          "WGBS_MethylCSeq"
head(list.files('data_CellLine/02_fastq/ChIPSeq/'))
# [1] "NE2_KYSE450.N.NE2.H3K27ac_1.fastq.gz"  "NE2_KYSE450.N.NE2.H3K27ac_2.fastq.gz"  "NE2_KYSE450.N.NE2.H3K27me3_1.fastq.gz"
# [4] "NE2_KYSE450.N.NE2.H3K27me3_2.fastq.gz" "NE2_KYSE450.N.NE2.H3K4me3_1.fastq.gz"  "NE2_KYSE450.N.NE2.H3K4me3_2.fastq.gz" 

# 4. ChIP-seq----
# 1) fastqc for each .fastq file (fastqc program)
# Example: 
system('script/10_run.00factqc.sh NE2_KYSE450.T.KYSE450.Input_1.fastq.gz') # (for KYSE450 input)
system('script/10_run.00factqc.sh NE2_KYSE450.T.KYSE450.H3K27me3_1.fastq.gz') # (for KYSE450 H3K27me3)

# 2) trim the fastq file and redo fastqc (trim_galore program)
# For the singel-end files, run 'script/10_run.01trimSingle.sh'
# For the paired-end files, run 'script/run.01trimPaired.sh'
# Example: 
system('script/10_run.02mappingPaired.sh NE2_KYSE450.T.KYSE450.H3K27me3')

# 3) mapping (bwa, samtools, sambamba and bedtools)
# For the singel-end files, run 'script/10_run.02mappingSingle.sh'
# For the paired-end files, run 'script/10_run.02mappingPaired.sh'
#Example: 
system('script/10_run.02mappingPaired.sh NE2_KYSE450.T.KYSE450.H3K27me3') 

# 4) mapping (macs2)
# For the sample with IP and input, run 'script/10_run.03MACS2WithInput.sh'
# For the sample with IP but no input, run 'script/10_run.03MACS2WithoutInput.sh'
# Example: 
system('script/10_run.03MACS2WithInput.sh NE2_KYSE450.T.KYSE450.H3K27me3 NE2_KYSE450.T.KYSE450.Input')

# 5) Create symbolic links for the bw files from the macs2 results into the following directory.
list.files('data_CellLine/03_ChIPSeq_Bw')
# [1] "H3K27ac"  "H3K27me3" "H3K4me3"  "H3K9me3" 

# 6) Sort the narrowPeak file, retain only chromosomes 1-22, remove the blacklist regions, and save the processed peaks.
peak_list = list.files('data_CellLine/03_ChIPSeq_run/03.MACS2/', pattern = 'H3K27me3.*narrowPeak', full.names = T)
blacklist_v2_hg38 = read.table('meta/blacklist/hg38-blacklist.v2.bed', sep = '\t', header = F)
for (i in 1:length(peak_list)) {
  print(i)
  df = read.table(peak_list, sep = '\t', header = F)[,1:3]
  colnames(df) = c('chrom', 'start', 'end')
  df = bed_sort(df)
  df = bed_intersect(df, blacklist_v2_hg38, invert = T)
  df = filter(df, chrom %in% c(paste0('chr',1:22)))
  index = sub('_peaks.narrowPeak','.peak.bed',basename(peak_list))
  write.table(df, file = paste0('data_CellLine/03_ChIPSeq_Peak/H3K27me3/',index,'.peak.bed'), 
              sep = '\t', quote = F, col.names = F, row.names = F)
}

peak_list = list.files('data_CellLine/03_ChIPSeq_run/03.MACS2/', pattern = 'H3K9me3.*narrowPeak', full.names = T)
blacklist_v2_hg38 = read.table('meta/blacklist/hg38-blacklist.v2.bed', sep = '\t', header = F)
for (i in 1:length(peak_list)) {
  print(i)
  df = read.table(peak_list, sep = '\t', header = F)[,1:3]
  colnames(df) = c('chrom', 'start', 'end')
  df = bed_sort(df)
  df = bed_intersect(df, blacklist_v2_hg38, invert = T)
  df = filter(df, chrom %in% c(paste0('chr',1:22)))
  index = sub('_peaks.narrowPeak','.peak.bed',basename(peak_list))
  write.table(df, file = paste0('data_CellLine/03_ChIPSeq_Peak/H3K9me3/',index,'.peak.bed'), 
              sep = '\t', quote = F, col.names = F, row.names = F)
}

# For the narrowPeak files identified by macs2, retain only the peaks on chr1-22. 
# Save the results in the following directory.
list.files('data_CellLine/03_ChIPSeq_Peak')
# [1] "H3K27me3" "H3K9me3"

# 5. LOCKs----
# 1) Identification of LOCKs and separate them into long or short LOCKs
peakFiles = c(list.files('data_CellLine/03_ChIPSeq_Peak/H3K27me3/', full.names = T),
              list.files('data_CellLine/03_ChIPSeq_Peak/H3K9me3/', full.names = T))
for (i in 1:length(peakFiles)) {
  print(i)
  system(paste0('Rscript script/03_LOCK_calling_cellLine.R ',peakFiles[i]))
}  

# 2) annotate the regions
PMD_coordinates_hg38 = read.table('meta/commonPMD_commonHMD/PMD_coordinates_hg38.bed', sep = '\t')
colnames(PMD_coordinates_hg38) = c('chrom','start','end','variation','MD_anno','PMD_length_anno')
commonPMD_hg38 = read.table('meta/commonPMD_commonHMD/commonPMD_hg38.bed', sep = '\t') %>% dplyr::rename(chrom=1, start=2, end=3)
commonHMD_hg38 = read.table('meta/commonPMD_commonHMD/commonHMD_hg38.bed', sep = '\t') %>% dplyr::rename(chrom=1, start=2, end=3)
a = c('H3K27me3','H3K9me3')
for (h in 1:2) {
  LOCK_TN[[h]] = list()
  for (s in 1:length(list.files('data_CellLine/04_LOCK/all_LOCKs/'))) {
    #LOCK region
    LOCK_df = read.table(file = list.files('data_CellLine/04_LOCK/all_LOCKs/', full.names = T)[i],
                         sep = '\t', header = F) %>% 
      dplyr::rename(chrom=1, start=2, end=3) %>% 
      bed_sort()
    
    #LOCK length
    LOCK_df = mutate(LOCK_df, 
                     d = end - start,
                     LOCK_ID = paste0(chrom,':',start,'-',end))
    LOCK_df$LOCK_length = 'Long_LOCK'
    LOCK_df[LOCK_df$d <= 1e5,]$LOCK_length = 'Short_LOCK'
    #PMD HMD annotation
    LOCK_df1 = bed_intersect(LOCK_df, PMD_coordinates_hg38, suffix = c('','.y')) %>% 
      select(-start.y, -end.y) %>% 
      dplyr::rename(MD_anno = MD_anno.y) %>% 
      group_by( d, LOCK_ID, MD_anno) %>% 
      summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
      filter(.overlap > d/2) %>% 
      select(LOCK_ID, MD_anno) %>% 
      column_to_rownames('LOCK_ID')
    LOCK_df = cbind(LOCK_df, LOCK_df1[LOCK_df$LOCK_ID,,drop=F])
    rownames(LOCK_df) = NULL
    LOCK_df[is.na(LOCK_df$MD_anno),]$MD_anno = 'Other'
    #PMD length annotation
    LOCK_df2 = bed_intersect(LOCK_df, SIL_region_and_annotation, suffix = c('','.y')) %>% 
      select(-start.y, -end.y) %>% 
      dplyr::rename(PMD_length_anno = PMD_length_anno.y) %>% 
      group_by(LOCK_ID, d, PMD_length_anno) %>% 
      summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
      filter(.overlap > d/2) %>% 
      select(LOCK_ID, PMD_length_anno) %>% 
      column_to_rownames('LOCK_ID')
    LOCK_df = cbind(LOCK_df, LOCK_df2[LOCK_df$LOCK_ID,,drop=F])
    rownames(LOCK_df) = NULL
    write.table(LOCK_df, file = paste0('data_CellLine/04_LOCK/withAnnotations/',list.files('data_CellLine/04_LOCK/all_LOCKs/')[i]),
                sep = '\t', quote = F, row.names = F, col.names = T)
  }
}

# 6. RNA-seq----
# 1) RNA-seq analysis was based on 'script/12_RNASeqMapping.sh'.
# Example:
system('sh script/12_RNASeqMapping.sh NE2_KYSE450.T.KYSE450_rep1 2')

head(list.files('data_CellLine/05_RNASeq/counts/'))
# [1] "breastNormalVsHer2Pos.N.76NF2V_rep1.count.txt" "breastNormalVsHer2Pos.N.76NF2V_rep2.count.txt"
# [3] "breastNormalVsHer2Pos.N.76NF2V_rep3.count.txt" "breastNormalVsHer2Pos.N.76NF2V_rep4.count.txt"
# [5] "breastNormalVsHer2Pos.N.MCF10A_rep1.count.txt" "breastNormalVsHer2Pos.N.MCF10A_rep2.count.txt"

# 2) FPKM_singleSample
exon_length = data.table::fread('meta/gtf/exon_length.GRCh38.84.gene.txt') %>% 
  as.data.frame() %>% 
  dplyr::rename(gene_id = 1, gene_name = 2, exon_length = 3)
countPath = 'data_CellLine/05_RNASeq/counts/'
FPKMPath = 'data_CellLine/05_RNASeq/FPKM_singleSample/'
for (i in 1:length(list.files(countPath))) {
  print(i)
  FPKM_df = read.table(file = list.files(countPath, full.names = T)[i], sep = '\t', header = F) %>% 
    dplyr::rename(gene_id = 1, count =2) %>% 
    left_join(exon_length, by = 'gene_id') %>% 
    select(-gene_id) %>% 
    group_by(gene_name) %>% 
    summarise(count = sum(count, na.rm = T),
              exon_length = sum(exon_length, na.rm = T),
              .groups = 'drop') %>% 
    na.omit() %>% 
    mutate(FPKM = (count/sum(count))/exon_length*1e9) %>% 
    select(gene_name, FPKM) %>% 
    arrange(gene_name)
  write.table(FPKM_df, file = paste0(FPKMPath, sub('count','FPKM',list.files(countPath)[i])))
}

# 3) FPKM_multiSample (Combine the FPKM results according to the tumor-normal pairs)
singleSamplePath = 'data_CellLine/05_RNASeq/FPKM_singleSample/'
multiSamplePath = 'data_CellLine/05_RNASeq/FPKM_multiSample/'
for (i in 1:nrow(all_pairs)) {
  i = 1
  tissue_type = all_pairs$pairs[i]
  files = list.files(singleSamplePath, pattern = tissue_type, full.names = T)
  FPKM = list()
  for (j in 1:length(files)) {
    FPKM[[j]] = read.table(files[j], header = T)
    rownames(FPKM[[j]]) = NULL
    FPKM[[j]] = column_to_rownames(FPKM[[j]], 'gene_name')
  }
  FPKM = do.call(cbind, FPKM)
  colnames(FPKM) = sub('.FPKM.txt','',basename(files))
  forRep0 = sub('(_rep.)', '', colnames(FPKM))
  forRep = table(forRep0)
  FPKM2 = as.data.frame(matrix(nrow = nrow(FPKM), ncol = length(forRep)))
  colnames(FPKM2) = names(forRep)
  rownames(FPKM2) = rownames(FPKM)
  for (j in 1:length(forRep)) {
    FPKM2[,j] = rowMeans(FPKM[,which(forRep0 == names(forRep)[j])])
  }
  write.csv(FPKM2, file = paste0(multiSamplePath,tissue_type,'_FPKM.csv'))
}

# 7. DNA methylation (WGBS, MethylC-seq and EPIC)----
# 1) WGBS or MethylC-seq analysis were based on 'script/13_runWGBS.sh'.
# Example:
system('sh script/13_runWGBS.sh MCF10A_normal_HCC1954_Her2Pos.T.HCC1954 1')

# 2) EPIC data analysis
library(sesame)
sesameDataCache()
setwd('data_CellLine/07_breast_EPIC_GSE171956/')
idat_dir = 'RawData'
betas = openSesame(idat_dir) 
write.table(betas, file = 'procedure/sesame_res_betas.txt', sep = '\t', quote = F, row.names = T, col.names = T)

#The annotation files are from https://zwdzwd.github.io/InfiniumAnnotation.
anno = fread(file = 'mapping_and_masking/EPIC.hg38.manifest.tsv.gz', select = c("CpG_chrm", "CpG_beg", "CpG_end", "probeID")) %>% as.data.frame()
masking = fread(file = 'mapping_and_masking/EPIC.hg38.mask.tsv.gz', select = c("probeID", "MASK_general")) %>% as.data.frame()

anno1 = merge(anno, masking, by = 'probeID') %>% 
  filter(MASK_general == FALSE)
betas1 = betas[anno1$probeID,]

for (i in 1:ncol(betas1)) {
  print(i)
  bdg_df = cbind(anno1[,c("CpG_chrm", "CpG_beg", "CpG_end")],
                 beta = betas1[,i]) %>% 
    na.omit() %>% 
    dplyr::rename(chrom=1, start=2, end=3) %>% 
    bed_sort()
  write.table(bdg_df, file = paste0('bdg/',colnames(betas1)[i],'.bdg'), sep = '\t', quote = F, row.names = F, col.names = F)
  system(paste0('bedtools merge -i ',paste0('bdg/',colnames(betas1)[i],'.bdg'),' -c 4 -o mean -d -1 > ',paste0('bdg/',colnames(betas1)[i],'.merge.bdg')))
  system(paste0('/opt/Anaconda3/bin//bedGraphToBigWig ','bdg/',colnames(betas1)[i],'.merge.bdg ',
                '/data3/liangyuan/05_LOCK_and_PMD/meta/chrom_size/hg38_ly.chrom.sizes ',
                'bw/',colnames(betas1)[i],'.bw'))
}
setwd('/data3/liangyuan/05_LOCK_and_PMD/')
# After completing the previous analysis, we obtained all the bdg and bw files for WGBS, MethylC-seq, and EPIC data. 
# Organize them as symbolic links into the 'data_CellLine/06_WGBS_bw/WGBS_bdg/' and 'data_CellLine/06_WGBS_bw/WGBS_bw/' directories

# 3) remove CpGs within CGIs
bwFiles = list.files('data_CellLine/06_WGBS_bw/WGBS_bw/')
for (i in 1:length(bwFiles)) {
  system(paste0('sh script/14_WGBS_bw_rmCGI.sh ',bwFiles[i]))
}

head(list.files('data_CellLine/06_WGBS_bw/rmCGI_WGBS_bw/'))
# [1] "breastNormalVsHer2Pos.N.HMEC.WGBS.rmCGI.bw"         "breastNormalVsHer2Pos.N.MCF10A_96850.WGBS.rmCGI.bw"
# [3] "breastNormalVsHer2Pos.N.MCF10A_rep3.EPIC.rmCGI.bw"  "breastNormalVsHer2Pos.N.MCF10A_rep4.EPIC.rmCGI.bw" 
# [5] "breastNormalVsHer2Pos.N.MCF10A_rep5.EPIC.rmCGI.bw"  "breastNormalVsHer2Pos.N.MCF10A_rep6.EPIC.rmCGI.bw" 
head(list.files('data_CellLine/06_WGBS_bw/rmCGI_WGBS_bdg/'))
# [1] "breastNormalVsHer2Pos.N.HMEC.WGBS.rmCGI.bdg"         "breastNormalVsHer2Pos.N.MCF10A_96850.WGBS.rmCGI.bdg"
# [3] "breastNormalVsHer2Pos.N.MCF10A_rep3.EPIC.rmCGI.bdg"  "breastNormalVsHer2Pos.N.MCF10A_rep4.EPIC.rmCGI.bdg" 
# [5] "breastNormalVsHer2Pos.N.MCF10A_rep5.EPIC.rmCGI.bdg"  "breastNormalVsHer2Pos.N.MCF10A_rep6.EPIC.rmCGI.bdg" 

# 8. Get Tumor-loss/shared/gain LOCKs----
for (i in 1:length(all_pairs$pairs)) {
  # 1) run bedtools multiinter
  print(i)
  tissue_type = all_tissue_types[i]
  # LOCK_files_long = list.files('data_CellLine/04_LOCK/long_LOCKs/', pattern = tissue_type, full.names = T)
  # LOCK_files_short = list.files('data_CellLine/04_LOCK/short_LOCKs/', pattern = tissue_type, full.names = T)
  # LOCK_files = c(LOCK_files_long, LOCK_files_short)
  # system(paste0('bedtools multiinter -header -i ', paste(LOCK_files, collapse = ' '),
  #               ' > data_CellLine/08_multiinter_LOCKs/',tissue_type,'_multiinter_LOCK_H3K27me3.txt'))
  
  # 2) read multiinter LOCKs
  fileName = paste0('data_CellLine/08_multiinter_LOCKs/',tissue_type,'_multiinter_LOCK_H3K27me3.txt')
  multiinter_LOCK = read.table(file = fileName, sep = '\t', header = T)
  colnames(multiinter_LOCK) = sub('_.0.7.bed','',sub('.H3K.*me3','',sub('X.data3.liangyuan.01_tumor_and_normal.data.03_LOCK.sorted.WS_.0.7.','',colnames(multiinter_LOCK))))
  
  # 3) get Tumor-loss/shared/gain LOCKs
  df = multiinter_LOCK
  get_Freq_anno = function(df, LOCK_length){
    sampleType = str_match(colnames(df), paste0('.*\\.([NT].*',LOCK_length,')'))[,2]
    T_cols = which(grepl('T',sampleType)); T_cols_Length = length(T_cols)
    N_cols = which(grepl('N',sampleType)); N_cols_Length = length(N_cols)
    df$T_Freq = rowSums(df[,T_cols,drop=F])
    df$N_Freq = rowSums(df[,N_cols,drop=F])
    df$T_Freq_ratio = df$T_Freq / T_cols_Length
    df$N_Freq_ratio = df$N_Freq / N_cols_Length
    colnames(df)[which(colnames(df) == 'T_Freq')] = paste0('T_Freq','_',LOCK_length)
    colnames(df)[which(colnames(df) == 'N_Freq')] = paste0('N_Freq','_',LOCK_length)
    colnames(df)[which(colnames(df) == 'T_Freq_ratio')] = paste0('T_Freq_ratio','_',LOCK_length)
    colnames(df)[which(colnames(df) == 'N_Freq_ratio')] = paste0('N_Freq_ratio','_',LOCK_length)
    return(df)
  }
  df = get_Freq_anno(df, 'longLOCK')
  df = get_Freq_anno(df, 'shortLOCK')
  df$group = NA
  ratio_cutoff = 0.667
  ratio_cutoff_complement = 1 - ratio_cutoff
  #long LOCKs
  if(sum(df$T_Freq_ratio_longLOCK>ratio_cutoff &
         df$N_Freq_ratio_longLOCK<ratio_cutoff_complement &
         df$T_Freq_ratio_shortLOCK<ratio_cutoff_complement &
         df$N_Freq_ratio_shortLOCK<ratio_cutoff_complement)>0){
    df[df$T_Freq_ratio_longLOCK>ratio_cutoff &
         df$N_Freq_ratio_longLOCK<ratio_cutoff_complement &
         df$T_Freq_ratio_shortLOCK<ratio_cutoff_complement &
         df$N_Freq_ratio_shortLOCK<ratio_cutoff_complement,]$group = 'Tumor-gain LOCKs'
  }
  if(sum(df$N_Freq_ratio_longLOCK>ratio_cutoff & 
         df$T_Freq_ratio_longLOCK<ratio_cutoff_complement &
         df$T_Freq_ratio_shortLOCK<ratio_cutoff_complement &
         df$N_Freq_ratio_shortLOCK<ratio_cutoff_complement)>0){
    df[df$N_Freq_ratio_longLOCK>ratio_cutoff & 
         df$T_Freq_ratio_longLOCK<ratio_cutoff_complement &
         df$T_Freq_ratio_shortLOCK<ratio_cutoff_complement &
         df$N_Freq_ratio_shortLOCK<ratio_cutoff_complement,]$group = 'Tumor-loss LOCKs'
  }
  if(sum(df$N_Freq_ratio_longLOCK>ratio_cutoff & 
         df$T_Freq_ratio_longLOCK>ratio_cutoff)>0){
    df[df$N_Freq_ratio_longLOCK>ratio_cutoff & 
         df$T_Freq_ratio_longLOCK>ratio_cutoff,]$group = 'Shared LOCKs'
  }
  #short LOCKs
  if(sum(df$T_Freq_ratio_shortLOCK>ratio_cutoff &
         df$N_Freq_ratio_shortLOCK<ratio_cutoff_complement &
         df$T_Freq_ratio_longLOCK<ratio_cutoff_complement &
         df$N_Freq_ratio_longLOCK<ratio_cutoff_complement)>0){
    df[df$T_Freq_ratio_shortLOCK>ratio_cutoff &
         df$N_Freq_ratio_shortLOCK<ratio_cutoff_complement &
         df$T_Freq_ratio_longLOCK<ratio_cutoff_complement &
         df$N_Freq_ratio_longLOCK<ratio_cutoff_complement,]$group = 'Tumor-gain short LOCKs'
  }
  if(sum(df$N_Freq_ratio_shortLOCK>ratio_cutoff & 
         df$T_Freq_ratio_shortLOCK<ratio_cutoff_complement &
         df$T_Freq_ratio_longLOCK<ratio_cutoff_complement &
         df$N_Freq_ratio_longLOCK<ratio_cutoff_complement)>0){
    df[df$N_Freq_ratio_shortLOCK>ratio_cutoff & 
         df$T_Freq_ratio_shortLOCK<ratio_cutoff_complement &
         df$T_Freq_ratio_longLOCK<ratio_cutoff_complement &
         df$N_Freq_ratio_longLOCK<ratio_cutoff_complement,]$group = 'Tumor-loss short LOCKs'
  }
  if(sum(df$N_Freq_ratio_shortLOCK>ratio_cutoff & 
         df$T_Freq_ratio_shortLOCK>ratio_cutoff)>0){
    df[df$N_Freq_ratio_shortLOCK>ratio_cutoff & 
         df$T_Freq_ratio_shortLOCK>ratio_cutoff,]$group = 'Shared short LOCKs'
  }
  df = select(df, chrom, start, end, group)
  TumorAndNormalLOCK_long = filter(df, group %in% c('Tumor-loss LOCKs','Shared LOCKs','Tumor-gain LOCKs'))
  TumorAndNormalLOCK_short = filter(df, group %in% c('Tumor-loss short LOCKs','Shared short LOCKs','Tumor-gain short LOCKs'))
  
  # 4) Annotate for the TumorAndNormalLOCKs
  Annotate_PMD_HMD = function(TumorAndNormalLOCK_input){
    TumorAndNormalLOCK__inFun1 = select(TumorAndNormalLOCK_input, chrom, start, end, group) %>% na.omit() %>% 
      group_by(group) %>% bed_merge() %>% ungroup() %>% 
      mutate(multiLOCK_ID = paste0(chrom,':',start,'-',end),
             d = end - start)
    #MD_anno
    PMD_HMD = rbind(mutate(commonPMD_hg38, MD_anno = 'common PMDs'),
                    mutate(commonHMD_hg38, MD_anno = 'common HMDs'))
    df1 = bed_intersect(TumorAndNormalLOCK__inFun1, PMD_HMD, suffix = c('', '.y')) %>% 
      group_by(multiLOCK_ID, d, MD_anno.y) %>% 
      summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
      filter(.overlap > d/2) %>% 
      select(-d, -.overlap) %>% 
      dplyr::rename(MD_anno = MD_anno.y)
    TumorAndNormalLOCK__inFun1 = left_join(TumorAndNormalLOCK__inFun1, df1, by = 'multiLOCK_ID')
    TumorAndNormalLOCK__inFun1[is.na(TumorAndNormalLOCK__inFun1$MD_anno),]$MD_anno = 'Others'
    return(TumorAndNormalLOCK__inFun1)
  }
  TumorAndNormalLOCK_long = Annotate_PMD_HMD(TumorAndNormalLOCK_long)
  TumorAndNormalLOCK_short = Annotate_PMD_HMD(TumorAndNormalLOCK_short)
  write.table(TumorAndNormalLOCK_long, file = paste0('data_CellLine/09_TumorAndNormalLOCKs/Long_LOCKs/',tissue_type,'.TumorAndNormalLOCK_long.txt'),
              sep = '\t', quote = F, row.names = F, col.names = T)
  write.table(TumorAndNormalLOCK_short, file = paste0('data_CellLine/09_TumorAndNormalLOCKs/Short_LOCKs/',tissue_type,'.TumorAndNormalLOCK_short.txt'),
              sep = '\t', quote = F, row.names = F, col.names = T)
}

# 9. computeMatirx analysis to compute signals for tumor-loss/shared/gain long LOCKs----
# 1) function
Fun_computeMatrix = function(script,bed,bw,output,log){
  system(paste('sh',script,bed,bw,output,log,sep = ' '))
}

# 2) H3K27me3 intensity
all_tissue_types = pairs_for_longLOCKAnalysis
for (t in 1:length(all_tissue_types)) {
  tissue_type = all_tissue_types[t]
  
  ChIPSeq_bw_files = list.files('data_CellLine/03_ChIPSeq_Bw/H3K27me3/', pattern = paste0(tissue_type,'.*H3K27me3'), full.names = T)
  
  Fun_computeMatrix(script = 'script/15_computeMatrix1.sh',
                    bed = paste0('data_CellLine/09_TumorAndNormalLOCKs/Long_LOCKs/',tissue_type,'.TumorAndNormalLOCK_long.txt'),
                    bw = ChIPSeq_bw_files[i],
                    output = paste0('data_CellLine/10_computeMatrix_longLOCKs/H3K27me3/',tissue_type,'.multiinter_LOCK.',fileName,'.RC_0.667.plotProfile.gz'),
                    log = paste0('data_CellLine/10_computeMatrix_longLOCKs/H3K27me3/log/',tissue_type,'.multiinter_LOCK.',fileName,'.RC_0.667.plotProfile.log')
  )
}

# 3) H3K27me3 and H3K9me3 signal (different bash scripts; only for BRCA)
for (t in 3:length(all_tissue_types)) {
  tissue_type = all_tissue_types[t]
  
  ChIPSeq_bw_files1 = list.files('data_CellLine/03_ChIPSeq_Bw/H3K27me3/', pattern = paste0(tissue_type,'.*H3K27me3'), full.names = T)
  ChIPSeq_bw_files2 = list.files('data_CellLine/03_ChIPSeq_Bw/H3K9me3/', pattern = paste0(tissue_type,'.*H3K9me3'), full.names = T)
  ChIPSeq_bw_files = c(ChIPSeq_bw_files1, ChIPSeq_bw_files2)
  Fun_computeMatrix(script = 'script/15_computeMatrix2.sh',
                    bed = paste0('data_CellLine/09_TumorAndNormalLOCKs/Long_LOCKs/',tissue_type,'.TumorAndNormalLOCK_long.txt'),
                    bw = ChIPSeq_bw_files[i],
                    output = paste0('data_CellLine/10_computeMatrix_longLOCKs/H3K27me3_H3K9me3/',tissue_type,'.multiinter_LOCK.',fileName,'.RC_0.667.plotProfile.gz'),
                    log = paste0('data_CellLine/10_computeMatrix_longLOCKs/H3K27me3_H3K9me3/log/',tissue_type,'.multiinter_LOCK.',fileName,'.RC_0.667.plotProfile.log')
  )
}

# 10. DNA methylation levels of the tumor-loss/shared/gain long LOCKs----
blacklist_hg38 = read.table(file = 'meta/blacklist/hg38-blacklist.v2.bed', sep = '\t', header = F)[,c(1:3)]; colnames(blacklist_hg38) = c('chrom', 'start', 'end')
Genome_hg38 = read.table('meta/chrom_size/hg38_ly.chrom.sizes', sep = '\t', header = F) %>%
  dplyr::rename(chrom=1, end=2) %>%
  mutate(start=0, .before = 2) %>%
  filter(!chrom %in% c('chrX','chrY')) %>%
  bed_subtract(.,blacklist_hg38)
Other_hg38 = bed_subtract(Genome_hg38, commonPMD_hg38) %>%
  bed_subtract(., commonHMD_hg38)
MD_anno_Overall = rbind(mutate(commonPMD_hg38, MD_anno = 'common PMDs'),
                        mutate(commonHMD_hg38, MD_anno = 'common HMDs'),
                        mutate(Other_hg38, MD_anno = 'Others'),
                        mutate(Genome_hg38, MD_anno = 'Genome'))

for (t in 1:length(all_tissue_types)) {
  tissue_type = all_tissue_types[t]
  TumorAndNormalLOCK_long = read.table(paste0('data_CellLine/09_TumorAndNormalLOCKs/Long_LOCKs/',
                                              tissue_type,'.TumorAndNormalLOCK_long.txt'))
  met_bdg_path = 'data_CellLine/06_WGBS_bw/rmCGI_WGBS_bdg/'
  if (t %in% 5:6) {
    useEPIC = T
  }
  if (useEPIC) {
    met_bdg_files = list.files(met_bdg_path, pattern = 'EPIC.rmCGI')[str_match(list.files(met_bdg_path, pattern = 'EPIC.rmCGI'), '(^.*?)\\..*?\\..*')[,2] == tissue_type]
  }else{
    met_bdg_files = list.files(met_bdg_path, pattern = 'WGBS.rmCGI')[str_match(list.files(met_bdg_path, pattern = 'WGBS.rmCGI'), '(^.*?)\\..*?\\..*')[,2] == tissue_type]
  }
  plotData_met = list()
  
  
  for (i in 1:length(met_bdg_files)) {
    print(i)
    df = read.table(file = paste0(met_bdg_path,met_bdg_files[i]), sep = '\t', header = F)
    colnames(df) = c('chrom', 'start', 'end', 'metSignal')
    df = filter(df, chrom %in% paste0('chr',1:22))
    df = bed_intersect(df, blacklist_hg38, invert = T)
    if (max(df$metSignal, na.rm = T) <= 1) {
      df$metSignal = df$metSignal*100
    }
    
    df = rbind(rbind(bed_intersect(TumorAndNormalLOCK_long, df, suffix = c('', '.y')) %>%
                       group_by(group, MD_anno) %>% 
                       summarise(metSignal = mean(metSignal.y, na.rm = T),
                                 .groups = 'drop'),
                     bed_intersect(TumorAndNormalLOCK_long, df, suffix = c('', '.y')) %>%
                       group_by(group) %>% 
                       summarise(metSignal = mean(metSignal.y, na.rm = T),
                                 .groups = 'drop') %>% 
                       mutate(MD_anno = 'Genome', .before = 2)),
               bed_intersect(MD_anno_Overall, df, suffix = c('', '.y')) %>% 
                 group_by(MD_anno) %>% 
                 summarise(metSignal = mean(metSignal.y, na.rm = T),
                           .groups = 'drop') %>% 
                 mutate(group = 'Overall', .before = 1)) %>% 
      mutate(gzName = met_bdg_files[i]) %>% 
      select(MD_anno, group, metSignal, gzName)
    df$TorN = str_match(df$gzName, '^.*?\\.([NT])\\..*')[,2]
    df$dataType = str_match(df$gzName, '.*\\.(.*?)\\.rmCGI.bdg')[,2]
    plotData_met[[i]] = df
  }
  plotData_met = do.call(rbind, plotData_met)
  plotData_met$MD_anno = factor(plotData_met$MD_anno, levels = c('Genome', 'common HMDs', 'common PMDs', 'Others'))
  plotData_met$group = factor(plotData_met$group, levels = c('Tumor-loss LOCKs', 'Shared LOCKs', 'Tumor-gain LOCKs', 'Overall'))
  plotData_met$TorN = factor(plotData_met$TorN, levels = c('T', 'N'))
  write.csv(plotData_met, file = paste0('data_CellLine/11_DNA_methylation_PMD_HMD/', tissue_type,'.met.csv'))
}
# 11. Summary of the characteristics of promoters associated with short LOCKs----
for (t in 1:length(pairs_for_shortLOCKAnalysis)) {
  print(t)
  tissue_type = all_tissue_types[t]
  # 1) all promoters (converted from hg19 using liftOver)
  promoter_with_anno_hg38 = read.table(file = 'meta/RefSeq_gene/hg38_promoter_with_anno.bed', sep = '\t', header = T) %>% 
    mutate(promoterID = paste0(chrom,':',start,'-',end), .before = 4)
  #2) add ESC_poised_promoters anno
  promoters = list()
  for (i in 1:length(list.files('meta/LOLACore/hg38/ESC_poised_promoters'))) {
    promoters[[i]] = read.table(file = list.files('meta/LOLACore/hg38/ESC_poised_promoters',
                                                  full.names = T)[i])
    colnames(promoters[[i]]) = c('chrom', 'start', 'end')
  }
  names(promoters) = sub('_hg38.bed','',list.files('meta/LOLACore/hg38/ESC_poised_promoters'))
  
  promoter_1 = bed_intersect(promoter_with_anno_hg38, promoters[[2]], suffix = c('','ESC_poised_promoters')) %>% 
    mutate(ESC_poised_promoters = T) %>% 
    select(promoterID, ESC_poised_promoters) %>% 
    unique() %>% 
    left_join(promoter_with_anno_hg38, ., by = 'promoterID')
  promoter_1[is.na(promoter_1$ESC_poised_promoters),]$ESC_poised_promoters = F
  
  # 3) add short LOCK annotation
  TumorAndNormalShortLOCK = read.table(paste0('data_CellLine/09_TumorAndNormalLOCKs/Short_LOCKs/',tissue_type,'.TumorAndNormalLOCK_short.txt'),
                                       sep = '\t',header = T)
  promoter_2 = bed_intersect(promoter_1, TumorAndNormalShortLOCK[,c('chrom','start','end','group','MD_anno')], suffix = c('','.shortLOCK')) %>% 
    select(-start.shortLOCK, -end.shortLOCK, -.source, -.overlap) %>% 
    unique()
  write.table(promoter_2[,1:3], file = paste0('data_CellLine//14_computeMatrix_shortLOCKs/input_promoters_shortLOCK/',tissue_type,'.bed'), 
              sep = '\t', quote = F, row.names = F, col.names = F)#for computeMatrix
  
  # 4) compute ChIP-Seq intensity
  ChIPSeq_bw_path = 'data_CellLine/03_ChIPSeq_Bw/H3K27me3/'
  ChIPSeq_bw_files = list.files(ChIPSeq_bw_path)[str_match(list.files(ChIPSeq_bw_path), '(^.*?)\\..*?\\..*')[,2] == tissue_type]
  computeMatrixFun = function(bed,bw,output,log){
    system(paste('sh ../script/15_computeMatrix3.sh',bed,bw,output,log,sep = ' '))
  }
  outputPath = 'data_CellLine/14_computeMatrix_shortLOCKs/H3K27me3/'
  logPath = paste0(outputPath,'/log/')
  for (fileName in ChIPSeq_bw_files) {
    computeMatrixFun(bed = paste0('data_CellLine//14_computeMatrix_shortLOCKs/input_promoters_shortLOCK/',tissue_type,'.bed'),
                     bw = paste0(ChIPSeq_bw_path,'/',fileName),
                     output = paste0(outputPath,tissue_type,'.multiinter_shortLOCK.',fileName,'.RC_0.667.gz'),
                     log = paste0(logPath,tissue_type,'.multiinter_shortLOCK.',fileName,'.RC_0.667.log'))
  }
  
  # 5) read computeMatrix result for H3K27me3 intensity
  gzFile = list()
  files = list.files(outputPath, pattern = tissue_type, full.names = T)
  for (i in 1:length(files)) {
    fileName = files[i]
    fileName1 = basename(files[i])
    gzFile[[i]] = read.table(fileName,
                             sep = '\t', header = F, skip = 1)[,c('V4','V7')] %>% 
      column_to_rownames('V4')
    colnames(gzFile[[i]]) = fileName1
  }
  gzFile = do.call(cbind, gzFile)
  colnames(gzFile) = str_match(colnames(gzFile), '.*shortLOCK.(.*?).bw')[,2]
  gzFile = gzFile[promoter_2$promoterID,]
  promoter_3 = cbind(promoter_2, gzFile)
  
  # 6) Compute intensity of other markers
  computeOtherMarkerIntensity = function(marker){
    ChIPSeq_bw_path = paste0('data_CellLine/03_ChIPSeq_Bw/',marker)
    ChIPSeq_bw_files = list.files(ChIPSeq_bw_path)[str_match(list.files(ChIPSeq_bw_path), '(^.*?)\\..*?\\..*')[,2] == tissue_type]
    
    computeMatrixFun = function(bed,bw,output,log){
      system(paste('sh ../script/15_computeMatrix3.sh',bed,bw,output,log,sep = ' '))
    }
    outputPath = paste0('data_CellLine/14_computeMatrix_shortLOCKs/',marker)
    logPath = paste0(outputPath,'/log/')
    for (fileName in ChIPSeq_bw_files) {
      computeMatrixFun(bed = paste0('data_CellLine//14_computeMatrix_shortLOCKs/input_promoters_shortLOCK/',tissue_type,'.bed'),
                       bw = paste0(ChIPSeq_bw_path,'/',fileName),
                       output = paste0(outputPath,tissue_type,'.multiinter_shortLOCK.',fileName,'.RC_',ratio_cutoff,'.gz'),
                       log = paste0(logPath,tissue_type,'.multiinter_shortLOCK.',fileName,'.RC_',ratio_cutoff,'.log'))
    }
  }
  computeOtherMarkerIntensity('H3K4me3')
  computeOtherMarkerIntensity('H3K27me3')
  
  # 7) read computeMatrix result of other markers
  gzFile = list()
  files = c(list.files('data_CellLine/14_computeMatrix_shortLOCKs/H3K4me3/', full.names = T, pattern = tissue_type),
            list.files('data_CellLine/14_computeMatrix_shortLOCKs/H3K27ac/', full.names = T, pattern = tissue_type))
  for (i in 1:length(files)) {
    fileName = files[i]
    fileName1 = basename(files[i])
    gzFile[[i]] = read.table(fileName,
                             sep = '\t', header = F, skip = 1)[,c('V4','V7')] %>% 
      column_to_rownames('V4')
    colnames(gzFile[[i]]) = fileName1
  }
  gzFile = do.call(cbind, gzFile)
  colnames(gzFile) = str_match(colnames(gzFile), '.*shortLOCK.(.*?).bw')[,2]
  promoter_3 = cbind(promoter_3, gzFile[promoter_3$promoterID,])
  
  # 8) DNA methylation (from .bdg)
  met_bdg_path = 'data_CellLine/06_WGBS_bw/WGBS_bdg/'
  met_bdg_files = list.files(met_bdg_path)[str_match(list.files(met_bdg_path), '(^.*?)\\..*?\\..*')[,2] == tissue_type]
  bdg_File = list()
  for (i in 1:length(met_bdg_files)) {
    bdg_File[[i]] = read.table(paste0(met_bdg_path,met_bdg_files[i]),
                               sep = '\t', header = F, skip = 1)
    colnames(bdg_File[[i]]) = c('chrom','start','end','metSignal')
    bdg_File[[i]] = bed_intersect(promoter_3[,1:3], bdg_File[[i]], suffix = c('','.y')) %>% 
      group_by(chrom, start, end) %>% 
      summarise(metSignal = median(metSignal.y, na.rm = T), .groups = 'drop') %>% 
      mutate(promoterID = paste0(chrom,':',start,'-',end)) %>% 
      select(promoterID, metSignal) %>% 
      column_to_rownames('promoterID')
    colnames(bdg_File[[i]]) = sub('bdg','metSignal',met_bdg_files[i])
    bdg_File[[i]] = bdg_File[[i]][promoter_3$promoterID,,drop=F]
  }
  bdg_File = do.call(cbind, bdg_File)
  promoter_3 = cbind(promoter_3, bdg_File)
  
  # 9) gene expression
  #FPKM
  FPKM_matrix_path = paste0('data_CellLine/05_RNASeq/FPKM_multiSample/',tissue_type,'_FPKM.csv')
  FPKM_matrix = read.csv(file = FPKM_matrix_path, row.names = 1)
  FPKM_matrix = FPKM_matrix[rownames(FPKM_matrix) %in% Refseq_gene$gene_name,]
  colnames(FPKM_matrix) = sub('RNASeq','FPKM',colnames(FPKM_matrix))
  
  #log2_ZMAD
  log2_ZMAD_matrix = log2(FPKM_matrix+0.1)
  ZMAD = function(x){0.6745*((x - median(x))/(mad(x)))}
  log2_ZMAD_matrix = as.matrix(log2_ZMAD_matrix)
  for (j in 1:ncol(log2_ZMAD_matrix)) {
    log2_ZMAD_matrix[,j] = as.numeric(ZMAD(log2_ZMAD_matrix[,j]))
  }
  log2_ZMAD_matrix = as.data.frame(log2_ZMAD_matrix)
  colnames(log2_ZMAD_matrix) = sub('FPKM','ZMAD_log2_FPKM',colnames(log2_ZMAD_matrix))
  promoter_3 = cbind(promoter_3,
                     log2_ZMAD_matrix[promoter_3$gene_name,])
  promoter_3 = cbind(promoter_3,
                     FPKM_matrix[promoter_3$gene_name,])
  
  #DESeq2
  countFiles = list.files('data_CellLine/05_RNASeq/counts/', pattern = tissue_type, full.names = T)
  countMatrix = list()
  for (i in 1:length(countFiles)) {
    countMatrix[[i]] = read.table(countFiles[i], sep = '\t', header = F) %>% 
      column_to_rownames('V1')
    colnames(countMatrix[[i]]) = sub('.count.txt','',basename(countFiles[i]))
  }
  countMatrix = do.call(cbind, countMatrix)
  countMatrix = countMatrix[!grepl('__',rownames(countMatrix)),]
  countMatrix = filter(countMatrix, rowMeans(countMatrix) > 0.5)
  geneName_geneID = read.table(file = '/data3/liangyuan/03_Roadmap_LOCK/data/00_gtf/geneName_geneID.txt')
  colnames(geneName_geneID) = c('gene_name', 'gene_id')
  geneName_geneID = column_to_rownames(geneName_geneID, 'gene_id')
  geneName_geneID = geneName_geneID[rownames(countMatrix),]
  countMatrix = cbind(countMatrix, geneName_geneID)
  countMatrix = group_by(countMatrix, geneName_geneID) %>% 
    summarise_all(sum, na.rm = T) %>% 
    ungroup() %>% 
    as.data.frame() %>% 
    na.omit() %>% 
    column_to_rownames('geneName_geneID')
  
  conditionVector = str_match(colnames(countMatrix), '.*\\.([NT])\\..*')[,2]
  coldata = data.frame(sample = colnames(countMatrix),
                       condition = factor(conditionVector,levels = c('N','T')))
  suppressMessages(library(DESeq2))
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = coldata,
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm") %>% 
    as.data.frame() %>% 
    rownames_to_column('gene_name')
  res = res[,c('gene_name', 'log2FoldChange', 'padj')]
  colnames(res)[3] = 'padj_DESeq2'
  promoter_3 = left_join(promoter_3, res, by = 'gene_name')
  promoter_shortLOCK = promoter_3
  
  #tumor FPKM, normal FPKM
  FPKMCol_T = promoter_shortLOCK[,grep('\\.T\\..*\\.FPKM',colnames(promoter_shortLOCK)),drop=F]
  FPKMCol_N = promoter_shortLOCK[,grep('\\.N\\..*\\.FPKM',colnames(promoter_shortLOCK)),drop=F]
  if (ncol(FPKMCol_T) > 1) {
    FPKMCol_T = rowMeans(FPKMCol_T)
    FPKMCol_N = rowMeans(FPKMCol_N)
  }else{
    FPKMCol_T = FPKMCol_T[,1]
    FPKMCol_N = FPKMCol_N[,1]
  }
  promoter_shortLOCK$FPKM_T = FPKMCol_T
  promoter_shortLOCK$FPKM_N = FPKMCol_N
  
  #tumor met, normal met
  metCol_T = promoter_shortLOCK[,grep('\\.T\\..*met',colnames(promoter_shortLOCK)),drop=F]
  metCol_N = promoter_shortLOCK[,grep('\\.N\\..*met',colnames(promoter_shortLOCK)),drop=F]
  if (ncol(metCol_T) > 1) {
    metCol_T = rowMeans(metCol_T)
    metCol_N = rowMeans(metCol_N)
  }else{
    metCol_T = metCol_T[,1]
    metCol_N = metCol_N[,1]
  }
  promoter_shortLOCK$met_T = metCol_T
  promoter_shortLOCK$met_N = metCol_N
  
  # 10) groups of promoters(Fig5E)
  # Upregulated group
  # log2FC >1 & pval < 0.05
  # tumor FPKM > 1
  # normal FPKM < 2
  promoter_shortLOCK$group = NA
  promoter_shortLOCK[promoter_shortLOCK$log2FoldChange > 1 & (!is.na(promoter_shortLOCK$log2FoldChange)) &
                       promoter_shortLOCK$padj_DESeq2 < 0.05 & (!is.na(promoter_shortLOCK$padj_DESeq2)) &
                       promoter_shortLOCK$FPKM_T > 1 & (!is.na(promoter_shortLOCK$FPKM_T)) &
                       promoter_shortLOCK$FPKM_N < 2 & (!is.na(promoter_shortLOCK$FPKM_N)),]$group = 'Upregulated'
  # Hypermethylated grouop
  # normal methylation <0.2
  # tumor_met-normal_met >0.1
  # Expression does not increase
  # normal FPKM<2
  promoter_shortLOCK[(!is.na(promoter_shortLOCK$log2FoldChange)) & (!is.na(promoter_shortLOCK$padj_DESeq2)) &
                       (promoter_shortLOCK$log2FoldChange <= 1 | promoter_shortLOCK$padj_DESeq2 > 0.05) & #不显著上调
                       promoter_shortLOCK$FPKM_N < 2 & (!is.na(promoter_shortLOCK$FPKM_N)) & #normal FPKM<2
                       promoter_shortLOCK$met_N < 0.2 & (!is.na(promoter_shortLOCK$met_N)) & #normal methylation <0.2
                       (promoter_shortLOCK$met_T - promoter_shortLOCK$met_N > 0.1) & (!is.na(promoter_shortLOCK$met_T)), # tumor_met-normal_met >0.1
  ]$group = 'hypermethylated'
  write.table(promoter_shortLOCK, file = paste0('data_CellLine/15_short_LOCK_promoter_feature/',tissue_type,'.bed'),
              sep = '\t', quote = F, row.names = F, col.names = T)
}

# PART 4 Main analysis and figures----
# Figure 1/ Figure S1----
# FigS1A. Taking the E009 sample from the Roadmap project as an example, the number of all LOCK peaks and typical peaks at different WScutoff values----
# 1) load data
library(nortest)
commonPMD = read.table(file = 'meta/commonPMD_commonHMD/commonPMD_hg19.bed',sep = '\t',header = F); colnames(commonPMD) = c('chrom','start','end')
commonHMD = read.table(file = 'meta/commonPMD_commonHMD/commonHMD_hg19.bed',sep = '\t',header = F); colnames(commonHMD) = c('chrom','start','end')
load('data_Roadmap/RDatas/samples.RData')
load('data_Roadmap/RDatas/peak_109.RData')

# 2) read LOCKs identified with all WScutoffs across all samples
ctoffLOCK = list()
for (i in 1:nrow(Roadmap_table2)) {
  ctoffLOCK[[i]] = list()
  for (j in 1:length(seq(-1,1,0.1))) {
    print(c(i,j))
    ctoffLOCK[[i]][[j]] =  read.table(file = paste0('data_Roadmap/03_LOCK/H3K27me3/',Roadmap_table2[i,'EID'],'/',Roadmap_table2[i,'EID'],'_',HisMark,'_ctoff',seq(-1,1,0.1)[j],'.bed'),sep = '\t',header = F) %>% 
      dplyr::rename(chrom = 1,start = 2,end = 3)
  }
  names(ctoffLOCK[[i]]) = paste0('ctoffLOCK',seq(-1,1,0.1)) 
}
names(ctoffLOCK) = Roadmap_table2$EID

# 3) identify LOCK peaks and typical peaks
LOCK_peak = list()
typical_peak = list()
for (i in 1:nrow(Roadmap_table2)) {
  LOCK_peak[[i]] = list()
  typical_peak[[i]] = list()
  for (j in 1:length(seq(-1,1,0.1))) {
    print(c(i,j))
    LOCK_peak[[i]][[j]] = bed_intersect(peak_109[[i]][,1:3],ctoffLOCK[[i]][[j]])[,1:3] %>% dplyr::rename(chrom = 1,start = 2,end = 3) %>% unique()
    typical_peak[[i]][[j]] = bed_intersect(peak_109[[i]][,1:3],ctoffLOCK[[i]][[j]],invert = T)
  }
  names(LOCK_peak[[i]]) = paste0('ctoffLOCK',seq(-1,1,0.1)) 
  names(typical_peak[[i]]) = paste0('ctoffLOCK',seq(-1,1,0.1)) 
}
names(LOCK_peak) = Roadmap_table2$EID
names(typical_peak) = Roadmap_table2$EID

# 4) get numbers of LOCK peaks and typical peaks; get plot
LOCK_peak_df = list()
typical_peak_df = list()
plotdata = list()

for (i in 1:nrow(Roadmap_table2)) {
  LOCK_peak_df[[i]] = data.frame(cutoff = as.numeric(sub('ctoffLOCK','',names(LOCK_peak[[i]]))),peak_number = NA,type = 'LOCK_peak')
  typical_peak_df[[i]] = data.frame(cutoff = as.numeric(sub('ctoffLOCK','',names(typical_peak[[i]]))),peak_number = NA,type = 'typical_peak')
  for (j in 1:length(seq(-1,1,0.1))) {
    print(c(i,j))
    
    LOCK_peak_df[[i]][j,'peak_number'] = nrow(LOCK_peak[[i]][[j]])
    typical_peak_df[[i]][j,'peak_number'] = nrow(typical_peak[[i]][[j]])
  }
  plotdata[[i]] = rbind(LOCK_peak_df[[i]],typical_peak_df[[i]])
}
names(LOCK_peak_df) = Roadmap_table2$EID
names(typical_peak_df) = Roadmap_table2$EID
names(plotdata) = Roadmap_table2$EID

plotdata_CREAM_cutoff = plotdata
plot = list()

i = which(names(plotdata) == 'E009')
plotdata_CREAM_cutoff = plotdata[[i]]
plotdata_CREAM_cutoff$type = sub('LOCK_peak','LOCK peaks',plotdata_CREAM_cutoff$type)
plotdata_CREAM_cutoff$type = sub('typical_peak','Typical peaks',plotdata_CREAM_cutoff$type)
FigS1A = ggplot(plotdata_CREAM_cutoff, aes(x=cutoff, y=peak_number, color=type)) +
  geom_point(alpha=0.8)+theme_light()+xlab("CREAM WS cutoff")+ylab("Peak number")+
  theme_classic()+
  ggtitle(Roadmap_table2[i,'EID'])+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.direction = 'horizontal',
        legend.position = 'bottom')+
  scale_x_continuous(breaks = seq(-1,1,0.2))+
  scale_y_continuous(breaks = seq(0,250000,50000))+
  geom_vline(xintercept = -0.7, linetype = 'dashed', color = 'grey')+
  scale_color_manual(values = c('LOCK peaks' = 'red','Typical peaks' = 'blue'))

pdf('plot/FigS1A/FigS1A.pdf', width = 7, height = 6)
FigS1A
dev.off()

Data_FigS1A = dplyr::rename(plotdata_CREAM_cutoff, `CREAM WS cutoff` = cutoff)
write.table(Data_FigS1A, 'plot/FigS1A/Data_FigS1A.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# FigS1B. The bar chart shows the change in LOCK peak numbers with each 0.1 increase in the WScutoff parameter.----
change = as.data.frame(matrix(nrow = 20,ncol = 0))
for (i in 1:nrow(Roadmap_table2)) {
  change[,i] = LOCK_peak_df[[i]]$peak_number[1:20] - LOCK_peak_df[[i]]$peak_number[2:21]
}
rownames(change) = seq(-1,1,0.1)[1:20]
colnames(change) = Roadmap_table2$EID

Mean = rowMeans(change)
Sd = genefilter::rowSds(as.matrix(change))

change_plotData = data.frame(ctoff = as.numeric(rownames(change)), Mean = Mean, Sd = Sd)
change_plotData$ctoff1 = paste0('[',change_plotData$ctoff,', ',(change_plotData$ctoff+0.1),']')
change_plotData$ctoff1 = factor(change_plotData$ctoff1, levels = change_plotData$ctoff1)
change_plotData$ctoff = factor(change_plotData$ctoff,levels = seq(-1,1,0.1)[1:20])
change_plotData$color = 'normal'
change_plotData[change_plotData$ctoff %in% c(-0.7),]$color = 'special'
FigS1B = ggplot(change_plotData,aes(x = ctoff1,y = Mean))+
  geom_bar(stat = 'identity',aes(color = color, fill = color))+
  geom_errorbar(aes(ymax = Mean+Sd, ymin = Mean-Sd))+
  theme_classic()+
  scale_color_manual(values = c('normal' = 'darkgrey', 'special' = 'red'))+
  scale_fill_manual(values = c('normal' = 'darkgrey', 'special' = 'red'))+
  guides(color = 'none',fill = 'none')+
  ylab('Number of\nchanged LOCK peaks')+
  xlab('CREAM WS cutoff')+
  ggtitle('109 samples')+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1))

pdf('plot/FigS1B/FigS1B.pdf', width = 7, height = 6)
FigS1B
dev.off()

Data_FigS1B = select(change_plotData, ctoff1, Mean, Sd, color) %>% 
  dplyr::rename(`CREAM WS cutoff` = ctoff1, Color = color)
rownames(Data_FigS1B) = NULL
write.table(Data_FigS1B, 'plot/FigS1B/Data_FigS1B.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# FigS1C. The count of different LOCKs in 109 samples.----
# load('data_Roadmap/RDatas/LOCK_109.RData')
LOCK_number = list()
for (i in 1:109) {
  LOCK_number[[i]] = as.data.frame(table(LOCK_109[[i]]$LOCK_length)) %>% 
    mutate(EID = names(LOCK_109)[i])
}
LOCK_number1 = do.call(rbind, LOCK_number)
LOCK_number2 = pivot_wider(LOCK_number1, names_from = 'Var1', values_from = 'Freq')
LOCK_number1$Var1 = sub('_',' ',LOCK_number1$Var1)
LOCK_number1$Var1 = factor(LOCK_number1$Var1, levels = c('Short LOCK','Long LOCK'))
FigS1C = ggplot(LOCK_number1, aes(x = Var1, y = Freq))+
  geom_boxplot(aes(color = Var1))+
  geom_jitter(size = 0.1)+
  theme_classic()+
  xlab('')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  ylab('Number of LOCKs')+
  scale_y_continuous(breaks = seq(0,5200,500))+
  theme(axis.text = element_text(size = 13, color = 'black'),
        axis.title = element_text(size = 13),
        plot.title = element_text(size = 15)
        )+
  coord_flip()+
  scale_color_manual(values = c('Long LOCK' = 'red', 'Short LOCK' = 'orange'))+
  geom_signif(comparisons = list(c('Long LOCK','Short LOCK')),
              map_signif_level = T,
              test = "t.test",test.args = c(var.equal = T,paired = T))+
  guides(color = 'none')

pdf('plot/FigS1C/FigS1C.pdf', width = 7, height = 4)
FigS1C
dev.off()

Data_FigS1C = LOCK_number1 %>% 
  dplyr::rename(LOCK_group = Var1, Number_of_LOCKs = Freq)
write.table(Data_FigS1C, 'plot/FigS1C/Data_FigS1C.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# FigS1D. The count of different peaks (D) in 109 samples.----
# load('data_Roadmap/RDatas/LOCK_peak_109.RData')
peak_number = list()
for (i in 1:109) {
  peak_number[[i]] = as.data.frame(table(LOCK_peak_109[[i]]$LOCK_length.peak)) %>% 
    mutate(EID = names(LOCK_peak_109)[i])
}
peak_number1 = do.call(rbind, peak_number)
peak_number2 = pivot_wider(peak_number1, names_from = 'Var1', values_from = 'Freq')
peak_number1$Var1 = factor(peak_number1$Var1, levels = rev(c('Peaks in long LOCKs', 'Peaks in short LOCKs', 'Typical peaks')))
FigS1D = ggplot(peak_number1, aes(x = Var1, y = Freq))+
  geom_boxplot(aes(color = Var1))+
  geom_jitter(size = 0.1)+
  theme_classic()+
  xlab('')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  ylab('Number of peaks')+
  scale_y_continuous(breaks = seq(0,300000,20000))+
  theme(axis.text = element_text(size = 13, color = 'black'),
        axis.title = element_text(size = 13),
        plot.title = element_text(size = 15))+
  coord_flip()+
  scale_color_manual(values = c('Peaks in long LOCKs' = 'red', 'Peaks in short LOCKs' = 'orange', 'Typical peaks' = 'blue'))+
  geom_signif(comparisons = list(c('Peaks in long LOCKs', 'Peaks in short LOCKs'),
                                 c('Peaks in long LOCKs', 'Typical peaks'),
                                 c('Typical peaks', 'Peaks in short LOCKs')
                                 ),
              map_signif_level = T, step_increase = 0.1, 
              test = "t.test",test.args = c(var.equal = T,paired = T))+
  guides(color = 'none')

pdf('plot/FigS1D/FigS1D.pdf', width = 7, height = 4)
FigS1D
dev.off()

Data_FigS1D = peak_number1 %>% 
  dplyr::rename(Number_of_peaks = Freq, Peak_group = Var1)
write.table(Data_FigS1D, 'plot/FigS1D/Data_FigS1D.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig1B. Comparative analysis of H3K27me3 Intensity, peak size, DNA methylation, and gene expression across different peak types.----
# load(file = 'data_Roadmap/RDatas//samples.RData')
# load(file = 'data_Roadmap/RDatas//LOCK_peak_109.RData')
# 1) peak intensity
getPlotData = function(input_list, var){
  plotData = list()
    for (s in 1:length(input_list)) {
      df = input_list[[s]]
      group = names(table(df$LOCK_length.peak))
      plotData[[s]] = data.frame(sample = names(input_list)[s],group = group, var = NA)
      for (g in 1:length(group)) {
        plotData[[s]][g,'var'] = median(filter(df, LOCK_length.peak == group[g])[,var,drop = T],na.rm = T)
      }
    }
  plotData = do.call(rbind, plotData)
  return(plotData)
}
myPlotFun = function(plotData,  
                     ylab ,
                     mytitle = 'pass',
                     scale_auto = T,
                     y_min = 2,  
                     y_max = 7,
                     step_increase = 0.1){
  group_levels = c('Typical peaks', 'Peaks in long LOCKs', 'Peaks in short LOCKs')
  plotData$group = factor(plotData$group, levels = group_levels)
  p = ggplot(plotData,aes(x = group, y = var))+
    geom_boxplot(aes(color = group))+
    ylab(ylab)+
    xlab('')+
    theme_classic()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(color = 'black'),
          axis.ticks.y = element_line(color = 'black'))+
    scale_color_manual(values = c('Peaks in long LOCKs' = 'red',
                                  'Peaks in short LOCKs' = 'orange',
                                  'Typical peaks' = 'blue'))+
    guides(color = 'none')+
    geom_signif(comparisons = list(c('Typical peaks', 'Peaks in long LOCKs'),
                                   c('Typical peaks', 'Peaks in short LOCKs'),
                                   c('Peaks in long LOCKs', 'Peaks in short LOCKs')
    ),
    test = "t.test",
    map_signif_level = T,
    test.args = c(var.equal = F,paired = T),
    step_increase = step_increase)
  if (!scale_auto) {
    p = p+ylim(y_min,y_max)
  }
  if (mytitle != 'pass') {
    p = p + ggtitle(mytitle)
  }
  return(p)
}
PlotData_intensity_109 = getPlotData(input_list = LOCK_peak_109, var = 'signal')
Fig1B1_intensity = myPlotFun(plotData = PlotData_intensity_109,
                             ylab = 'H3K27me3 intensity',
                             mytitle = 'Peak intensity')
# 2) peak size
# load(file = 'data_Roadmap/RDatas//samples.RData')
# load(file = 'data_Roadmap/RDatas//LOCK_peak_109.RData')
PlotData_peakSize_109 = getPlotData(input_list = LOCK_peak_109, var = 'd')
Fig1B2_peakSize = myPlotFun(plotData = PlotData_peakSize_109,
                            ylab = 'Base pair',
                            mytitle = 'Peak size')
# 3) DNA methylation
# load(file = 'data_Roadmap/RDatas//LOCK_peak_20.RData')
PlotData_metSignal_20 = getPlotData(input_list = LOCK_peak_20, var = 'metSignal')
Fig1B3_DNAMethylation = myPlotFun(plotData = PlotData_metSignal_20,
                                  ylab = 'mean β values',
                                  mytitle = 'DNA methylation')
# 4) Expression of nearest genes
# load(file = 'data_Roadmap/RDatas//LOCK_peak_exp_49.RData')
plotData_exp_49 = list()
for (s49 in 1:49) {
  print(s49)
  df = LOCK_peak_exp_49[[s49]]
  EID = names(LOCK_peak_exp_49)[[s49]]
  myExpTable = df[,c('PeakAnnotation','LOCK_length.peak','exp')] %>% 
    mutate(nearestGene = str_match(PeakAnnotation,'.*///(.*)///.*')[,2]) %>% 
    select(-PeakAnnotation) %>% 
    unique()
  whether_one_group = table(myExpTable$nearestGene) == 1
  geneList = names(whether_one_group[whether_one_group == T])# 一个基因多组
  myExpTable = filter(myExpTable,nearestGene %in% geneList)
  
  plotData_exp_49[[s49]] = group_by(myExpTable,LOCK_length.peak) %>% 
    summarise(var = median(exp,na.rm = T), number = n(), .groups = 'drop') %>% 
    mutate(sample = EID) %>% 
    select(sample,LOCK_length.peak,var,number) %>% dplyr::rename(group = 2)
}
plotData_exp_49 = do.call(rbind,plotData_exp_49)

Fig1B4_exp = myPlotFun(plotData = plotData_exp_49,  
                       ylab = 'log10 (RPKM+0.1)',
                       mytitle = 'Expression of\n nearest genes')
# 5) legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend_group = c('Typical peaks', 'Peaks in long LOCKs', 'Peaks in short LOCKs')
legend_color = c('blue', 'red', 'orange')
legend_df = data.frame(group = factor(legend_group, levels = legend_group),
                       color = legend_color,
                       num = 1:length(legend_group))
p1 = ggplot(legend_df, aes(x = group, y = num, color = group))+
  geom_boxplot()+
  scale_color_manual(values = legend_df$color)+
  theme_bw()+
  guides(color = 'legend')+
  theme(legend.title = element_blank()#,
        # legend.text = element_text(size = 13),
  )
Fig1B5_legend <- g_legend(p1 + guides(fill="legend"))

# ~ combine
Fig1B = ggarrange(Fig1B1_intensity, Fig1B2_peakSize, Fig1B3_DNAMethylation, Fig1B4_exp, Fig1B5_legend, nrow = 1)

pdf('plot/Fig1B/Fig1B.pdf', width = 13, height = 3)
Fig1B
dev.off()

Data_Fig1B1_intensity = PlotData_intensity_109 %>% dplyr::rename(EID = sample, Peak_group = group, Peak_intensity = var)
Data_Fig1B2_peakSize = PlotData_peakSize_109 %>% dplyr::rename(EID = sample, Peak_group = group, Peak_size = var)
Data_Fig1B3_DNAMethylation = PlotData_metSignal_20 %>% dplyr::rename(EID = sample, Peak_group = group, DNA_methylation = var)
Data_Fig1B4_exp = plotData_exp_49 %>% dplyr::rename(EID = sample, Peak_group = group, Expresion_of_nearest_genes = var)
write.table(Data_Fig1B1_intensity, file = 'plot/Fig1B/Data_Fig1B1_intensity.txt', sep = '\t', quote = F, row.names = F, col.names = T)
write.table(Data_Fig1B2_peakSize, file = 'plot/Fig1B/Data_Fig1B2_peakSize.txt', sep = '\t', quote = F, row.names = F, col.names = T)
write.table(Data_Fig1B3_DNAMethylation, file = 'plot/Fig1B/Data_Fig1B3_DNAMethylation.txt', sep = '\t', quote = F, row.names = F, col.names = T)
write.table(Data_Fig1B4_exp, file = 'plot/Fig1B/Data_Fig1B4_exp.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig1C. A heatmap showing the top 20 significantly enriched GO terms for the nearest genes of the peaks.----
# load('data_Roadmap/RDatas//samples.RData')
GO_All = list()
for (ex in 1:length(Roadmap_table2$EID)) {
  EID = Roadmap_table2$EID[ex]
  GO_AllLists = read.csv(file = paste0('data_Roadmap/10_metascape_and_LOLA/Metascape_output/H3K27me3_',
                                       EID,'/Enrichment_GO/GO_AllLists.csv'),
                         header = T)
  GO_All[[ex]] = filter(GO_AllLists, GeneList %in% c('Long_LOCK', 'Short_LOCK', 'Typical_peaks')) %>% 
    mutate(EID = EID) %>% 
    mutate(`-LogP` = -LogP, .before = LogP) %>% 
    mutate(`-Log.q.value.` = -Log.q.value., .before = Log.q.value.) %>% 
    filter(`-Log.q.value.` > -log10(0.05)) #The Metascape webpage shows that both LogP and Log.q.value are base 10.
}

GO_All = do.call(rbind, GO_All)
GO_All$term = paste0(GO_All$GO,': ',GO_All$Description)
GO_All = GO_All[grepl('GO',GO_All$GO),]
topTerms = group_by(GO_All, term) %>% 
  summarise(sumLogq = sum(`-Log.q.value.`)) %>% 
  arrange(desc(sumLogq)) %>% 
  .[1:20,'term',drop = T] 
GO_All = filter(GO_All, term %in% topTerms)
GO_All$GeneList = sub('Long_LOCK','Peaks in long LOCKs',GO_All$GeneList)
GO_All$GeneList = sub('Short_LOCK','Peaks in short LOCKs',GO_All$GeneList)
GO_All$GeneList = sub('Typical_peaks','Typical peaks',GO_All$GeneList)
GO_All$GeneList = factor(GO_All$GeneList, levels = c('Peaks in long LOCKs', 'Peaks in short LOCKs', 'Typical peaks'))

GO_All1 = GO_All
GO_All1$EID_GeneList_term = paste(GO_All1$EID, GO_All1$GeneList, GO_All1$term, sep = '_')
EID_GeneList = c(paste0(Roadmap_table2$EID,'_','Peaks in long LOCKs'), paste0(Roadmap_table2$EID,'_','Peaks in short LOCKs'), paste0(Roadmap_table2$EID,'_','Typical peaks'))
AllDots = list()
for (i in 1:length(EID_GeneList)) {
  AllDots[[i]] = paste0(EID_GeneList[i],'_',topTerms)
}
if (!is.function(c)) {rm(c)}
AllDots = do.call(c, AllDots)

addtable1 = data.frame(EID_GeneList_term = AllDots[!AllDots %in% GO_All1$EID_GeneList_term])
addtable1 = mutate(addtable1, toSeparate = EID_GeneList_term) %>% 
  separate(col = 'toSeparate', into = c('EID','GeneList','term'), sep = '_') %>% 
  select(GeneList, EID, term, EID_GeneList_term)
addtable2 = as.data.frame(matrix(data = NA, nrow = nrow(addtable1), ncol = ncol(GO_All1)-4))
colnames(addtable2) = colnames(GO_All1)[1:(length(colnames(GO_All1))-4)]
addtable = cbind(addtable2, addtable1)
GO_All1 = rbind(GO_All1, addtable)
GO_All_heatmap = GO_All1[,c('term','GeneList','EID','Log.q.value.')] %>% 
  unique() %>% 
  pivot_wider(names_from = 'term', values_from = 'Log.q.value.')
GO_All_heatmap$GeneList = factor(GO_All_heatmap$GeneList, levels = c('Typical peaks', 'Peaks in long LOCKs', 'Peaks in short LOCKs'))
GO_All_heatmap = arrange(GO_All_heatmap, GeneList, EID)

dat = GO_All_heatmap[,-1:-2] %>% t()
dat = -dat
dat[dat<log10(0.05)] = NA
rownames(dat) = sub('GO.*: ','',rownames(dat))
color_min = 1.3
color_max = 30
ht = Heatmap(dat,cluster_rows = F,cluster_columns = F,
             col = colorRamp2(c(color_min,color_max),c('white','red')),
             top_annotation = columnAnnotation(group = GO_All_heatmap$GeneList,
                                               col = list(group = c("Peaks in long LOCKs" = "red",
                                                                    "Peaks in short LOCKs" = "orange",
                                                                    "Typical peaks" = "blue"))),
             heatmap_legend_param = list(title = '-log10(q-Value)',
                                         col_fun = colorRamp2(c(color_min,color_max),c('white','red')),
                                         title = "test", at = c(color_min,color_max)),
             # column_title = "Sample",
             column_title_side = 'bottom',
             row_names_side = "left",
             row_names_max_width = max_text_width(
               rownames(dat),
               gp = gpar(fontsize = 12)
             )
             
)
Fig1C = draw(ht, merge_legend = TRUE)

pdf('plot/Fig1C/Fig1C.pdf', width = 7, height = 7)
Fig1C
dev.off()
Data_Fig1C = dat
write.csv(Data_Fig1C, file = 'plot/Fig1C/Data_Fig1C.txt', sep = '\t', quote = F, row.names = F, col.names = F)

# Fig1D. The left stacked bar chart shows how different peaks are distributed across genomic elements in 109 samples.----
# load('data_Roadmap/RDatas//LOCK_peak_109.RData')
# load('data_Roadmap/RDatas//samples.RData')
annotatePeaks_plotData = list()
for (i in 1:109) {
  print(i)
  EID = names(LOCK_peak_109)[i]
  df = LOCK_peak_109[[i]]
  df$Elements = str_match(df$PeakAnnotation, '(.*?)///')[,2]
  df$Elements2 = str_match(df$Elements, '(.*) \\(')[,2]
  df[is.na(df$Elements2),]$Elements2 = df[is.na(df$Elements2),]$Elements
  
  total_size = group_by(df, LOCK_length.peak) %>% 
    summarise(sum_d = sum(d))
  annotatePeaks_plotData[[i]] = df %>% 
    group_by(LOCK_length.peak, Elements2) %>% 
    summarise(d = sum(d), .groups = 'drop') %>% 
    left_join(., total_size, by = 'LOCK_length.peak') %>% 
    mutate(`Proportion of elements` = d/sum_d*100) %>% 
    mutate(EID = EID, .before = 1) %>% 
    select(EID, LOCK_length.peak, Elements2,`Proportion of elements`) %>% 
    dplyr::rename(Elements = Elements2)
}
annotatePeaks_plotData = do.call(rbind, annotatePeaks_plotData)
annotatePeaks_plotData$LOCK_length.peak = gsub('Peaks in long LOCKs','Peaks in\nlong LOCKs', annotatePeaks_plotData$LOCK_length.peak)
annotatePeaks_plotData$LOCK_length.peak = gsub('Peaks in short LOCKs','Peaks in\nshort LOCKs', annotatePeaks_plotData$LOCK_length.peak)
annotatePeaks_plotData$LOCK_length.peak = factor(annotatePeaks_plotData$LOCK_length.peak, levels = c('Typical peaks',
                                                                                                     'Peaks in\nlong LOCKs',
                                                                                                     'Peaks in\nshort LOCKs'))
colourCount = length(unique(Cellratio$Var1))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
Fig1D_left = ggplot(annotatePeaks_plotData, aes(x = EID, y = `Proportion of elements`, fill = Elements))+
  geom_bar(stat = 'identity', position = 'stack', width=1)+
  facet_grid(.~LOCK_length.peak)+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())+
  xlab('Sample')+
  scale_y_continuous(expand = c(0,0))+
  ylab('Proportion of elements (%)')+
  scale_fill_manual(values = col_vector[-6:-8])

annotatePeaks_plotData1 = annotatePeaks_plotData
annotatePeaks_plotData1$LOCK_length.peak = sub('\n',' ',annotatePeaks_plotData1$LOCK_length.peak)
annotatePeaks_plotData1$LOCK_length.peak = factor(annotatePeaks_plotData1$LOCK_length.peak, levels = c('Typical peaks',
                                                                                                       'Peaks in long LOCKs',
                                                                                                       'Peaks in short LOCKs'))
plotData_intergenic = filter(annotatePeaks_plotData1, Elements == 'Intergenic')
plotData_promoter = filter(annotatePeaks_plotData1, Elements == 'promoter-TSS')
getElementBoxplot = function(plotData){
  ggplot(plotData, aes(x = LOCK_length.peak, y = `Proportion of elements`))+
    geom_boxplot(aes(color = LOCK_length.peak), outlier.size = 0.1)+
    facet_grid(.~Elements)+
    xlab('')+
    theme_classic()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          strip.background = element_blank())+
    # BigWordTheme+StripTheme+
    geom_signif(comparisons = list(c('Peaks in long LOCKs', 'Typical peaks'),
                                   c('Peaks in short LOCKs', 'Typical peaks'),
                                   c('Peaks in long LOCKs', 'Peaks in short LOCKs')),
                map_signif_level = T,
                step_increase = 0.1)+
    scale_color_manual(values = c('Peaks in long LOCKs' = 'red',
                                  'Peaks in short LOCKs' = 'orange',
                                  'Typical peaks' = 'blue'
    ))+
    ylab('Proportion of elements (%)')
}
p_proportionOfElementsBox = getElementBoxplot(annotatePeaks_plotData1)
element_types = unique(annotatePeaks_plotData1$Elements)
p = list()
for (i in 1:length(element_types)) {
  plotData = filter(annotatePeaks_plotData1, Elements == element_types[i])
  p[[i]] = getElementBoxplot(plotData)
}
names(p) = element_types
p_proportionOfElementsBox_intergenic = getElementBoxplot(plotData_intergenic)
p_proportionOfElementsBox_promoter = getElementBoxplot(plotData_promoter)
Fig1D_right = ggarrange(p[[4]]+guides(color = 'none'),
                        p[[7]]+guides(color = 'none'), 
                        nrow = 2)
pdf('plot/Fig1D/Fig1D.pdf', width = 8, height = 7)
ggarrange(Fig1D_left,Fig1D_right, nrow = 1, widths = c(2,1))
dev.off()

Data_Fig1D = annotatePeaks_plotData1 %>% 
  dplyr::rename(Peak_group = LOCK_length.peak)
write.table(Data_Fig1D, file = 'plot/Fig1D/Data_Fig1D.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig1E/FigS1E. Enrichment levels of three types of peaks within common PMD or HMD (H3K27me3 and H3K9me3).----
# load(file = 'data_Roadmap/RDatas//samples.RData')
# load(file = 'data_Roadmap/RDatas//LOCK_peak_109.RData')
# load(file = 'data_Roadmap/RDatas//LOCK_peak_109_H3K9me3.RData')
SIL_region_and_annotation = read.table(file = 'meta/commonPMD_subgroup/04_PMD_length_Plots/SIL_region_and_annotation.bed', sep = '\t', header = T)
SIL_region_and_annotation_select1 = SIL_region_and_annotation[,c('chrom','start','end','MD_anno')] %>% 
  mutate(d = end - start)
MD_total_size = group_by(SIL_region_and_annotation_select1, MD_anno) %>% summarise(MD_size = sum(d)/1e8)
SIL_region_and_annotation_select2 = SIL_region_and_annotation[,c('chrom','start','end','PMD_length_anno')] %>% 
  mutate(d = end - start) %>% na.omit()

getEnrichmentPlot = function(LOCK_peak_File){
  plotData_all = list()
  for (s109 in 1:length(LOCK_peak_File)) {
    print(s109)
    df = LOCK_peak_File[[s109]][,c('chrom','start','end','d','LOCK_length.peak')]
    overlap_size = bed_intersect(df, SIL_region_and_annotation_select1, suffix = c('','.y')) %>% 
      group_by(LOCK_length.peak, MD_anno.y) %>% 
      summarise(.overlap = sum(.overlap), .groups = 'drop')
    peak_size = group_by(df, LOCK_length.peak) %>% 
      summarise(peak_size = sum(d), .groups = 'drop')
    plotData = left_join(overlap_size, peak_size, by = 'LOCK_length.peak')
    plotData = dplyr::rename(plotData, MD_anno = MD_anno.y)
    plotData = left_join(plotData, MD_total_size, by = 'MD_anno')
    plotData$EnrichmentScore = plotData$.overlap/plotData$peak_size/plotData$MD_size
    plotData_all[[s109]] = plotData %>% 
      mutate(EID = names(LOCK_peak_File)[s109]) %>% 
      select(EID, LOCK_length.peak, MD_anno, EnrichmentScore) %>% 
      mutate(LOCK_length.peak = sub('Peaks in short LOCKs','Peaks in\nshort\nLOCKs',
                                    sub('Peaks in long LOCKs','Peaks in\nlong\nLOCKs',
                                        sub('Typical peaks','Typical\npeaks',LOCK_length.peak))),
             MD_anno = sub('MD','MDs',sub('common','common\n',MD_anno))) %>% 
      mutate(LOCK_length.peak = factor(LOCK_length.peak, 
                                       levels = c('Typical\npeaks','Peaks in\nlong\nLOCKs','Peaks in\nshort\nLOCKs')),
             MD_anno = factor(MD_anno, levels = c('common\nPMDs','common\nHMDs')))
  }
  plotData_all = do.call(rbind, plotData_all)
  p = ggplot(plotData_all,aes(x = MD_anno, y = EnrichmentScore))+
    geom_boxplot(aes(color = LOCK_length.peak))+
    facet_grid(.~LOCK_length.peak)+
    ylab(ylab)+
    theme_classic()+
    theme(axis.text = element_text(color = 'black'),
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank())+
    xlab('')+
    scale_color_manual(values = c('Peaks in\nlong\nLOCKs' = 'red',
                                  'Peaks in\nshort\nLOCKs' = 'orange',
                                  'Typical\npeaks' = 'blue'))+
    geom_signif(comparisons = list(c('common\nPMDs','common\nHMDs')),
                map_signif_level = T,
                test = "t.test",test.args = c(var.equal = T,paired = T))
  
  colnames(plotData_all) = c('EID','peak_group','Methylation_domain','Enrichment_score')
  return(list(plotData_all, p))
}
peakEnrichment_K27 = getEnrichmentPlot(LOCK_peak_109)
peakEnrichment_K9 = getEnrichmentPlot(LOCK_peak_109_H3K9me3)
Data_Fig1E = peakEnrichment_K27[[1]]
Data_FigS1E = peakEnrichment_K9[[1]]
Fig1E = peakEnrichment_K27[[2]]
FigS1E = peakEnrichment_K9[[2]]+ggtitle('H3K9me3 (Peak level)')

pdf('plot/Fig1E/Fig1E.pdf', width = 8, height = 6)
Fig1E
dev.off()

pdf('plot/FigS1E/FigS1E.pdf', width = 8, height = 6)
FigS1E
dev.off()

write.table(Data_Fig1E, file = 'plot/Fig1E/Data_Fig1E.txt', sep = '\t', quote = F, row.names = F ,col.names = T)
write.table(Data_FigS1E, file = 'plot/FigS1E/Data_FigS1E.txt', sep = '\t', quote = F, row.names = F ,col.names = T)

# Fig1F. A heatmap shows the significance of peak enrichment in DMV regions or transposon regions----
load('10_metascape_and_LOLA/lolaResults_peak_level.RData') 
all_LOLA = list()
for (i in 1:length(lolaResults_peak_level)) {
  all_LOLA[[i]] = lolaResults_peak_level[[i]][,c('sample','group','support','qValue',
                                                 'c','collection','description','cellType',
                                                 'tissue','antibody','treatment','dataSource',
                                                 'filename','size')] %>% 
    mutate(`Odds ratio` = support/(support+c)) %>% 
    mutate(qValueLog = -log10(qValue + 1e-322)) %>% 
    filter(qValueLog > -log10(0.05))
}
all_LOLA = do.call(rbind, all_LOLA)
all_LOLA = distinct(all_LOLA, sample, group, collection, description, .keep_all = T)

LOLA_TE = filter(all_LOLA, collection == 'UCSC_repeat_class')
LOLA_DMV = filter(all_LOLA, collection == 'DMV')
LOLA_TE_And_DMV = filter(all_LOLA, collection %in% c('UCSC_repeat_class', 'DMV'))
LOLA_TE_And_DMV$description = sub('repeat_','',LOLA_TE_And_DMV$description)
LOLA_TE_And_DMV$description = factor(LOLA_TE_And_DMV$description, levels = c("H1_DMVs", "ME_DMVs",  "NPC_DMVs", "TBL_DMVs", "MSC_DMVs",
                                                                             'LINE','LTR','SINE'))
LOLA_TE_And_DMV = filter(LOLA_TE_And_DMV, !is.na(description))
LOLA_TE_And_DMV = left_join(LOLA_TE_And_DMV, Roadmap_table2[,c('EID', 'summary')], by = c('sample' = 'EID'))
LOLA_TE_And_DMV$group = gsub('_',' ',LOLA_TE_And_DMV$group)


point.size.min = 0.01
point.size.max = 1
LOLA_TE_And_DMV1 = LOLA_TE_And_DMV
LOLA_TE_And_DMV1[LOLA_TE_And_DMV1$qValueLog > 50,]$qValueLog = 50


LOLA_TE_And_DMV2 = LOLA_TE_And_DMV1 %>% 
  select(sample, group, qValueLog, description) %>% 
  arrange(description, group, sample) %>% 
  pivot_wider(names_from = description, values_from = qValueLog)
typical = filter(LOLA_TE_And_DMV2, group == 'Typical peaks') %>% 
  arrange(desc(LTR), desc(LINE), desc(SINE))
long = filter(LOLA_TE_And_DMV2, group == 'Peaks in long LOCKs') %>% 
  arrange(desc(SINE), desc(LINE), desc(LTR))
short = filter(LOLA_TE_And_DMV2, group == 'Peaks in short LOCKs') %>% 
  arrange(desc(H1_DMVs), desc(ME_DMVs), desc(NPC_DMVs), desc(TBL_DMVs), desc(MSC_DMVs))
LOLA_TE_And_DMV3 = rbind(typical, long, short)
LOLA_TE_And_DMV3$group = factor(LOLA_TE_And_DMV3$group, levels = c('Typical peaks', 'Peaks in long LOCKs', 'Peaks in short LOCKs'))
dat = LOLA_TE_And_DMV3[,-1:-2] %>%
  t()
p = Heatmap(dat,cluster_rows = F,cluster_columns = F,
            col = colorRamp2(c(1.32,50),c('white','red')),
            top_annotation = columnAnnotation(group = LOLA_TE_And_DMV3$group,
                                              col = list(group = c("Peaks in long LOCKs" = "red",
                                                                   "Peaks in short LOCKs" = "orange",
                                                                   "Typical peaks" = "blue"))),
            heatmap_legend_param = list(title = '-log10 (q-Value)',
                                        col_fun = colorRamp2(c(1.32,50),c('white','red')),
                                        title = "test", at = c(1.32,50),
                                        labels = c("-log10 (0.05)", "≥50")),
            column_title = "Sample",
            column_title_side = 'bottom',
            row_names_side = 'left'
)
Fig1F = draw(p, merge_legend = TRUE)

pdf('plot/Fig1F/Fig1F.pdf', height = 4,width = 7)
Fig1F
dev.off()

Data_Fig1F = LOLA_TE_And_DMV3
write.table(Data_Fig1F, 'plot/Fig1F/Data_Fig1F.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Combine the plots in Figure 1----
BigWordTheme = theme(plot.title = element_text(hjust = 0.5, size = 15),
                     axis.text = element_text(color = 'black', size = 13),
                     axis.title = element_text(size = 13),
                     legend.title = element_blank(),
                     legend.text = element_text(size = 13))
StripTheme = theme(strip.text = element_text(size = 13))
Empty = Fig1C
Figure_1 = ggarrange(Empty,
                    Fig1B+BigWordTheme+StripTheme,
                    ggarrange(Fig1C, Fig1D_left+BigWordTheme+StripTheme, Fig1D_right+BigWordTheme+StripTheme, 
                              nrow = 1, widths = c(1,1,0.5), labels = c('C','D'), font.label = list(size = 23)),
                    ggarrange(Fig1E+BigWordTheme+StripTheme, Fig1F, nrow = 1, widths = c(1,1.5), labels = c('E','F'), font.label = list(size = 23)),
                    nrow = 4, heights = c(1.5,1,1.5,1), labels = c('A','B'), font.label = list(size = 23))
pdf('20_output_plots/Figure_1.pdf',height = 13, width = 12)
Figure_1
dev.off()

pdf('20_output_plots/Fig1C.pdf',height = 4.3,width = 6)
Fig1C
dev.off()

pdf('20_output_plots/Fig1F.pdf',height = 2.9,width = 5.5)
Fig1F
dev.off()
# Combine the plots in Figure S1----
Figure_S1 = ggarrange(FigS1A+BigWordTheme+StripTheme,
                      FigS1B+BigWordTheme+StripTheme,
                      FigS1C+BigWordTheme+StripTheme,
                      FigS1D+BigWordTheme+StripTheme,
                      FigS1E+BigWordTheme+StripTheme,
                      nrow = 3, ncol = 2, heights = c(1.5,1,1.5),
                      labels = LETTERS[1:6],
                      font.label = list(size = 23))
pdf('20_output_plots/Figure_S1.pdf', width = 13, height = 13)
Figure_S1
dev.off()
# Figure 2/ Figure S2----
# Fig2A. DNA methylation in long LOCKs and background within common PMD or HMD.----
load('data_Roadmap/RDatas/LOCK_20.RData')
# 1) DNA methylation of the overall common PMDs or HMDs
load(file = 'data_Roadmap/RDatas/SIL_met_20.RData')
overall_plotData_LOCKMethylation_MD = filter(SIL_met_20, SIL_anno %in% c('commonPMD', 'commonHMD'))[,1:3] %>% mutate(group = 'Overall',number = NA) %>% select(sample, SIL_anno, group, metSignal, number) %>% dplyr::rename(group_MD = 2, var = 4)
# 2) get the plotData
overall_plotData = overall_plotData_LOCKMethylation_MD
plotData = list()
for (s in 1:length(LOCK_20)) {
  print(s)
  df = LOCK_20[[s]]
  df[is.na(df[,MD_anno,drop = T]),MD_anno] = 'NA'
  group_MD = names(table(df[,MD_anno]))
  group = names(table(df$LOCK_length))
  little_table = list()
  for (g_MD in 1:length(group_MD)) {
    little_table[[g_MD]] = data.frame(sample = names(LOCK_20)[s],group_MD = group_MD[g_MD],group = group, LOCK_met_signal = NA)
    
  }
  little_table = do.call(rbind,little_table)
  
  for(i in 1:nrow(little_table)) {
    little_table[i,'LOCK_met_signal'] = median(df[(df$LOCK_length %in% little_table[i,'group'] & df[,MD_anno,drop=T] %in% little_table[i,'group_MD']),][,LOCK_met_signal,drop = T],na.rm = T)
    little_table[i,'number'] = nrow(df[(df$LOCK_length %in% little_table[i,'group'] & df[,MD_anno,drop=T] %in% little_table[i,'group_MD']),])
  }
  plotData[[s]] = little_table
}
plotData = do.call(rbind,plotData)
plotData = rbind(plotData, overall_plotData_LOCKMethylation_MD)
plotData = filter(plotData, !is.na(LOCK_met_signal))
plotData = filter(plotData, group_MD != 'NA')
plotData = filter(plotData, group != 'Short_LOCK')
plotData$group = sub('Long_LOCK','Long LOCKs',plotData$group)
plotData$group = factor(plotData$group, levels = c('Long LOCKs','Overall'))
plotData$group_MD = sub('MD','MDs',sub('common','common\n',plotData$group_MD))
plotData$group_MD = factor(plotData$group_MD, levels = c('common\nPMDs','common\nHMDs'))

# 3) get the plot
library(ggpubr)
library(ggsci)
library(introdataviz)
Fig2A = ggplot(plotData,aes(x = group_MD,y = LOCK_met_signal,fill = group)) +
  # split violin
  geom_split_violin(alpha = .5, trim = F,color = NA,width = 1) +
  # mean point
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  # errorbar
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,
               position = position_dodge(0.2)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,color = 'black',hjust = 1),
        legend.position = 'bottom',
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  xlab('')+
  ylab('Beta value (%)')+
  scale_fill_manual(values = c('red','purple'))+
  ggtitle('DNA methylation')+
  stat_compare_means(aes(group=group), method = 't.test', paired = T,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "NS")),label = "p.signif",
                     label.y = 1, size = 5)+
  ylim(0,1)

pdf('plot/Fig2A/Fig2A.pdf', width = 4, height = 5)
Fig2A
dev.off()

plotData$number = NULL
Data_Fig2A = dplyr::rename(plotData, DNA_methylation = var)
write.table(Data_Fig2A, file = paste0('plot/Fig2A/Data_Fig2A.txt'), sep = '\t', quote = F, row.names = F, col.names = T)

# Fig2B. Gene expression in long LOCKs and background within common PMD or HMD.----
load('data_Roadmap/RDatas/gene_49.RData')
plotData = list()
for (s in 1:length(gene_49)) {
  print(s)
  df = gene_49[[s]]
  withinLOCK = filter(df, LOCK_length == 'Long_LOCK' & LOCK_length.peak == 'Peaks in long LOCKs') %>% 
    group_by(MD_anno) %>% 
    summarise(exp = median(exp, na.rm = T)) %>% 
    mutate(EID = names(gene_49)[s], group = 'Long_LOCK')
  overall = group_by(df, MD_anno) %>% 
    summarise(exp = median(exp, na.rm = T), .groups = 'drop') %>% 
    mutate(EID = names(gene_49)[s], group = 'Overall')
  plotData[[s]] = rbind(withinLOCK,overall)
}
plotData = do.call(rbind,plotData) %>% na.omit()
plotData$MD_anno = sub('MD','MDs',sub('common','common\n',plotData$MD_anno))
plotData$MD_anno = factor(plotData$MD_anno, levels = c('common\nPMDs','common\nHMDs'))
plotData$group = sub('Long_LOCK','Long LOCKs', plotData$group)
Fig2B = ggplot(plotData,aes(x = MD_anno,y = exp,fill = group)) +
  # split violin
  geom_split_violin(alpha = .5, trim = F,color = NA,width = 1) +
  # mean point
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  # errorbar
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,
               position = position_dodge(0.2)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,color = 'black',hjust = 1),
        legend.position = 'bottom',
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  xlab('')+
  ylab('log10 (FPKM+0.1)')+
  scale_fill_manual(values = c('red','purple'))+
  ggtitle('Gene expression')+
  stat_compare_means(aes(group=group), method = 't.test', paired = T,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "NS")),label = "p.signif",
                     label.y = 1.5, size = 5)+
  ylim(-1.5,1.5)

pdf('plot/Fig2B/Fig2B.pdf', width = 4, height = 5)
Fig2B
dev.off()

Data_Fig2B = plotData %>% dplyr::rename(Expression = exp) %>% select(EID, group, MD_anno, Expression)
write.table(Data_Fig2B, file = 'plot/Fig2B/Data_Fig2B.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig2C/Fig2F. Gene density in long LOCKs and background within common PMD/HMD or S/I/L.----
load('meta/RefSeq_gene/hg19_Refseq_gene.RData')
load('data_Roadmap/RDatas/LOCK_109.RData')
# 1) Overall gene density in methylation domains
SIL_hg19 = read.table('meta/commonPMD_commonHMD/SIL_hg19.bed')
S_PMD = filter(SIL_hg19, V6 == 'S')[,1:3] %>% mutate(MD_anno = 'S-PMDs')
I_PMD = filter(SIL_hg19, V6 == 'I')[,1:3] %>% mutate(MD_anno = 'I-PMDs')
L_PMD = filter(SIL_hg19, V6 == 'L')[,1:3] %>% mutate(MD_anno = 'L-PMDs')
commonPMD_hg19 = read.table('meta/commonPMD_commonHMD/commonPMD_hg19.bed') %>% mutate(MD_anno = 'common PMDs')
commonHMD_hg19 = read.table('meta/commonPMD_commonHMD/commonHMD_hg19.bed') %>% mutate(MD_anno = 'common HMDs')

getGeneDensity = function(targetRegion){
  colnames(targetRegion)[1:3] = c('chrom','start','end')
  gene_number_inside_region = bed_intersect(targetRegion, hg19_Refseq_gene[,1:5], suffix = c('', '.gene')) %>% 
    group_by(gene_name.gene, d.gene, MD_anno) %>% 
    summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
    mutate(gene_number = .overlap / d.gene) %>% 
    group_by(MD_anno) %>% 
    summarise(gene_number = sum(gene_number))
  region_size = sum(targetRegion$end - targetRegion$start)
  geneDensity = data.frame(MD_anno = gene_number_inside_region$MD_anno,
                           gene_density = gene_number_inside_region$gene_number/region_size*1e6)
  return(geneDensity)
}

overall_PMD = getGeneDensity(commonPMD_hg19)
overall_HMD = getGeneDensity(commonHMD_hg19)
overall_S = getGeneDensity(S_PMD)
overall_I = getGeneDensity(I_PMD)
overall_L = getGeneDensity(L_PMD)

overall = list(overall_PMD, overall_HMD, overall_S, overall_I, overall_L)
overall = do.call(rbind, overall)
overall = mutate(overall, EID = 'Overall', group = 'Overall', .before = 1)

# 2) Gene density in long LOCKs
longLOCKGeneDensity = list()
for (s109 in 1:109) {
  print(s109)
  df = LOCK_109[[s109]]
  EID = names(LOCK_109)[s109]
  
  long_LOCK_in_PMD = filter(df, LOCK_length == 'Long_LOCK' & MD_anno == 'commonPMD')[,1:3] %>% mutate(MD_anno = 'common PMDs')
  long_LOCK_in_HMD = filter(df, LOCK_length == 'Long_LOCK' & MD_anno == 'commonHMD')[,1:3] %>% mutate(MD_anno = 'common HMDs')
  long_LOCK_in_S = filter(df, LOCK_length == 'Long_LOCK' & PMD_length_anno == 'S')[,1:3] %>% mutate(MD_anno = 'S-PMDs')
  long_LOCK_in_I = filter(df, LOCK_length == 'Long_LOCK' & PMD_length_anno == 'I')[,1:3] %>% mutate(MD_anno = 'I-PMDs')
  long_LOCK_in_L = filter(df, LOCK_length == 'Long_LOCK' & PMD_length_anno == 'L')[,1:3] %>% mutate(MD_anno = 'L-PMDs')
  
  longLOCKGeneDensity0 = list(getGeneDensity(long_LOCK_in_PMD) %>% mutate(EID = EID, group = 'Long LOCKs', .before = 1),
                         getGeneDensity(long_LOCK_in_HMD) %>% mutate(EID = EID, group = 'Long LOCKs', .before = 1),
                         getGeneDensity(long_LOCK_in_S) %>% mutate(EID = EID, group = 'Long LOCKs', .before = 1),
                         getGeneDensity(long_LOCK_in_I) %>% mutate(EID = EID, group = 'Long LOCKs', .before = 1),
                         getGeneDensity(long_LOCK_in_L) %>% mutate(EID = EID, group = 'Long LOCKs', .before = 1))
  longLOCKGeneDensity[[s109]] = do.call(rbind, longLOCKGeneDensity0)
}
longLOCKGeneDensity = do.call(rbind, longLOCKGeneDensity)
plotData_geneDensity = rbind(overall, longLOCKGeneDensity)
plotData_geneDensity$group = factor(plotData_geneDensity$group, levels = c('Long LOCKs', 'Overall'))
plotData_geneDensity_MD = filter(plotData_geneDensity, MD_anno %in% c('common PMDs','common HMDs')) %>% 
  mutate(MD_anno = factor(MD_anno, levels = c('common PMDs','common HMDs')))
plotData_geneDensity_PMD_length = filter(plotData_geneDensity, MD_anno %in% c('S-PMDs','I-PMDs','L-PMDs')) %>% 
  mutate(MD_anno = factor(MD_anno, levels = c('S-PMDs','I-PMDs','L-PMDs'))) %>% 
  dplyr::rename(PMD_length_anno = MD_anno)

# 3) get plots
Fig2C = ggplot(plotData_geneDensity_MD,aes(x = MD_anno, y = log2(gene_density), fill = group)) +
  # split violin
  geom_split_violin(alpha = .5, trim = F,color = NA,width = 1) +
  # mean point
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  # errorbar
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,
               position = position_dodge(0.2)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,color = 'black',hjust = 1),
        legend.position = 'bottom',
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  xlab('')+
  ylab('log2 (gene number/Mb)')+
  scale_fill_manual(values = c('red','purple'))+
  ggtitle('Gene density')+
  ylim(2,5)

Fig2F = ggplot(plotData_geneDensity_PMD_length,aes(x = PMD_length_anno, y = log2(gene_density), fill = group)) +
  # split violin
  geom_split_violin(alpha = .5, trim = F,color = NA,width = 1) +
  # mean point
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  # errorbar
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,
               position = position_dodge(0.2)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,color = 'black',hjust = 1),
        legend.position = 'bottom',
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  xlab('')+
  ylab('log2 (gene number/Mb)')+
  scale_fill_manual(values = c('red','purple'))+
  ggtitle('Gene density')

pdf('plot/Fig2C/Fig2C.pdf', width = 4, height = 5)
Fig2C
dev.off()

Data_Fig2C = plotData_geneDensity_MD
write.table(Data_Fig2C, file = 'plot/Fig2C/Data_Fig2C.txt', sep = '\t', quote = F, row.names = F, col.names = T)

pdf('plot/Fig2F/Fig2F.pdf', width = 4, height = 5)
Fig2F
dev.off()

Data_Fig2F = plotData_geneDensity_PMD_length
write.table(Data_Fig2F, file = 'plot/Fig2F/Data_Fig2F.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig2D. Relative loop density in the long LOCKs and background within whole genome, common PMD, or common HMDs.----

# 1) Download the Hi-C data in the bedpe format.
# All Hi-C data in ENCODE matched the 109 sample types in Roadmap were downloaded.
# The Hi-C files were renamed based on the EID of their matched Roadmap samples.

head(list.files('data_Roadmap/11_Encode_HiC_hg19/'))
# [1] "E017.ENCSR345VTI.ENCFF267VBT_hg19.bed" "E032.ENCSR847RHU.ENCFF281CTW_hg19.bed" "E046.ENCSR971CJS.ENCFF492SAH_hg19.bed"
# [4] "E065.ENCSR554SLA.ENCFF337GAR_hg19.bed" "E065.ENCSR797MWY.ENCFF161RGE_hg19.bed" "E066.ENCSR753CKM.ENCFF498XDH_hg19.bed"

loops_hg19 = list()
hic_files = list.files('data_Roadmap/11_Encode_HiC_hg19/', full.names = T)
for (i in 1:length(hic_files)) {
  loops_hg19[[i]] = read.table(hic_files[i], sep = '\t', header = F)
  colnames(loops_hg19[[i]]) = c('chrom.x','start.x','end.x','chrom.y','start.y','end.y')
  loops_hg19[[i]]$loopID = 1:nrow(loops_hg19[[i]])
}
names(loops_hg19) = basename(hic_files)
save(loops_hg19, file = 'data_CellLine/RDatas/loops_hg19.RData')

# 2) Compute the number of the loops whose both sides are in long LOCKs or random regions.
load('data_Roadmap/RDatas/samples.RData')
load('data_Roadmap/RDatas/LOCK_109.RData')
load(file = 'data_CellLine/RDatas/loops_hg19.RData')
SIL_region_and_annotation = read.table(file = '00_commonPMD_subgroup/04_PMD_length_Plots/SIL_region_and_annotation.bed', header = T)

MD_options = c('Genome', 'commonPMD', 'commonHMD')
loopNumberPerLOCKNumber = list()
loopNumberPerLOCKSize = list()
for (nMD in 1:length(MD_options)) {
  loopNumberPerLOCKNumber[[nMD]] = as.data.frame(matrix(nrow = length(loops_hg19), ncol = 0))
}
names(loopNumberPerLOCKNumber) = MD_options
genome = read.table(file = '00_chrom_size/hg19_ly.chrom.sizes',sep = '\t',header = F) %>% 
  dplyr::rename(chrom=1, size=2) %>% 
  filter(!chrom %in% c('chrX', 'chrY'))

if(F){
  for (i in 1:length(loops_hg19)) {
    print(i)
    EID = substr(names(loops_hg19)[i], 1, 4)
    x = which(names(LOCK_109) == EID)
    all_LOCK = LOCK_109[[x]]
    LOCK = all_LOCK %>% mutate(LOCK_ID = paste0(chrom,':',start,'-',end)) %>% .[,c("chrom", "start", "end", "d_LOCK", "MD_anno","PMD_length_anno", "LOCK_length", "LOCK_ID")]
    LOCK = filter(LOCK, LOCK_length == 'Long_LOCK')
    LOCK = LOCK[,c('chrom', 'start', 'end', 'LOCK_ID')]
    tryCatch( {
      
      shuffleOfLOCK = list()
      k=1
      cycle = 0
      set.seed(1234)
      sampleNumber = sample(10000000,1000)
      while (cycle <= 1000) {
        tryCatch({
          cycle = cycle+1
          #a. use bed_shuffle to get the random regions.
          shuffleOfLOCK[[k]] = bed_shuffle(LOCK[,1:3], genome, seed = sum(sampleNumber[cycle]+i+k), excl = all_LOCK[,1:3]) %>% bed_sort() %>% mutate(LOCK_ID = paste0(chrom,':',start,'-',end))
          k = k+1
        }, error = function(e){e})
        
        if (k == 51) {
          break
        }
      }
      
      getLoopSta = function(LOCK_in_fun){
        loop_A = loops_hg19[[i]][,c('chrom.x','start.x','end.x','loopID')] %>% dplyr::rename(chrom=1, start=2, end=3)
        loop_B = loops_hg19[[i]][,c('chrom.y','start.y','end.y','loopID')] %>% dplyr::rename(chrom=1, start=2, end=3)
        #b. Add MD_anno and PMD_length_anno to loop_A/B.
        add_anno_for_loop_df = function(df){
          df = mutate(df, d = end - start)
          # add MD_anno
          df1 = bed_intersect(df, SIL_region_and_annotation[,c('chrom', 'start', 'end', 'MD_anno')], suffix = c('','.y')) %>% 
            group_by(chrom, start, end, loopID, d, MD_anno.y) %>% 
            dplyr::rename(MD_anno = MD_anno.y) %>% 
            summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
            filter(.overlap > d/2) %>% 
            select(loopID, MD_anno) %>% 
            left_join(df, ., by = 'loopID')
          # add PMD_length_anno
          df2 = bed_intersect(df1, SIL_region_and_annotation[,c('chrom', 'start', 'end', 'PMD_length_anno')], suffix = c('','.y')) %>% 
            group_by(chrom, start, end, loopID, d, PMD_length_anno.y) %>% 
            dplyr::rename(PMD_length_anno = PMD_length_anno.y) %>% 
            summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
            filter(.overlap > d/2) %>% 
            select(loopID, PMD_length_anno) %>% 
            left_join(df1, ., by = 'loopID')
          return(df2)
        }
        loop_A = add_anno_for_loop_df(loop_A)
        loop_B = add_anno_for_loop_df(loop_B)
        #c. add LOCK_ID to loop_A/B
        add_LOCK_ID_for_loop_df = function(df){
          bed_intersect(df, LOCK_in_fun, suffix = c('', '.y')) %>% 
            filter(.overlap > d/2) %>% 
            dplyr::rename(LOCK_ID = LOCK_ID.y) %>% 
            select(loopID, LOCK_ID) %>% 
            left_join(df, ., by = 'loopID')
        }
        loop_A = add_LOCK_ID_for_loop_df(loop_A)
        loop_B = add_LOCK_ID_for_loop_df(loop_B)
        interact_df = merge(loop_A, loop_B, by = 'loopID', suffixes = c('.A', '.B'))
        return(interact_df)
      }
      LOCK_interact_df = getLoopSta(LOCK)
      Shuffle_interact_df = list()
      for (k in 1:50) {
        Shuffle_interact_df[[k]] = getLoopSta(shuffleOfLOCK[[k]])
      }
      
      for (nMD in 1:length(MD_options)) {
        get_interact_df_MD = function(interact_df){
          one_MD = MD_options[nMD]
          if (one_MD %in% 'Genome') {  interact_df_MD = interact_df  }
          if (one_MD %in% c('commonPMD', 'commonHMD')) {  interact_df_MD = filter(interact_df, MD_anno.A == one_MD & MD_anno.B == one_MD)  }
          if (one_MD %in% 'Others') {  interact_df_MD = filter(interact_df, (!MD_anno.A %in% c('commonPMD', 'commonHMD')) & (!MD_anno.B %in% c('commonPMD', 'commonHMD')))  }
          if (one_MD %in% c('S', 'I', 'L')) {  interact_df_MD = filter(interact_df, PMD_length_anno.A == one_MD & PMD_length_anno.B == one_MD)  }
          return(interact_df_MD)
        }
        LOCK_interact_df__MD = get_interact_df_MD(LOCK_interact_df)
        Shuffle_interact_df__MD = list()
        for (k in 1:50) {
          Shuffle_interact_df__MD[[k]] = get_interact_df_MD(Shuffle_interact_df[[k]])
        }
        
        
        loopNumberPerLOCKNumber[[nMD]][i,'LOCK_res'] = nrow(filter(LOCK_interact_df__MD, LOCK_ID.A == LOCK_ID.B))/nrow(LOCK)
        shuffle_res = list()
        for (k in 1:50) {
          shuffle_res[[k]] = nrow(filter(Shuffle_interact_df__MD[[k]], LOCK_ID.A == LOCK_ID.B))/nrow(shuffleOfLOCK[[k]])
        }
        shuffle_res = do.call(mean,shuffle_res)
        loopNumberPerLOCKNumber[[nMD]][i,'shuffle_res'] = shuffle_res
      }
    },
    error = function(e){message('error')})
  }}
for (i in 1:length(loopNumberPerLOCKNumber)) {
  rownames(loopNumberPerLOCKNumber[[i]]) = names(loops_hg19)
  loopNumberPerLOCKNumber[[i]]$EID = str_match(rownames(loopNumberPerLOCKNumber[[i]]), '(.*?)\\..*')[,2]
}
# To simplify the process, in cases where one EID corresponds to multiple Hi-C data, we use the distinct function to select one.
for (i in 1:length(loopNumberPerLOCKNumber)) {
  loopNumberPerLOCKNumber[[i]] = na.omit(loopNumberPerLOCKNumber[[i]])
  loopNumberPerLOCKNumber[[i]] = loopNumberPerLOCKNumber[[i]] %>% 
    distinct(EID, .keep_all = T) %>%  
    select(-EID)
  print(nrow(loopNumberPerLOCKNumber[[i]]))
}
# save(loopNumberPerLOCKNumber, file = 'data_CellLine/RDatas/loopNumberPerLOCKNumber.RData')

load(file = 'data_CellLine/RDatas/loopNumberPerLOCKNumber.RData')

#The following samples were used finally.
HiC_sample = data.frame(samples = rownames(loopNumberPerLOCKNumber[[1]])) %>% 
  separate(col = 'samples', into = c('EID','Experiment_ID','File_ID'), sep = '\\.')
write.table(HiC_sample, file = 'plot/Fig2D/HiC_samples.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# 3) get the plots
plotData = list()
p = list()
for (nMD in 1:length(MD_options)) {
  loopNumberPerLOCKNumber[[nMD]]$id = rownames(loopNumberPerLOCKNumber[[nMD]])
  loopNumberPerLOCKNumber[[nMD]]$compare = 'LOCK > Shuffle'
  if (sum(loopNumberPerLOCKNumber[[nMD]]$LOCK_res <= loopNumberPerLOCKNumber[[nMD]]$shuffle_res) > 0) {
    loopNumberPerLOCKNumber[[nMD]][loopNumberPerLOCKNumber[[nMD]]$LOCK_res <= loopNumberPerLOCKNumber[[nMD]]$shuffle_res,]$compare = 'LOCK <= Shuffle'
  }
  df = loopNumberPerLOCKNumber[[nMD]] %>% 
    pivot_longer(col = 1:2, names_to = 'group', values_to = 'value')
  df$group = sub('_res','',df$group)
  df$group = sub('shuffle','Shuffle',df$group)
  df$group = factor(df$group, levels = c('LOCK','Shuffle'))
  df$facet = MD_options[nMD]
  #Handle outliers using the Winsorize method.
  df_LOCK = filter(df, group == 'LOCK')
  df_Shuffle = filter(df, group == 'Shuffle')
  max_LOCK = mean(df_LOCK$value)+3*sd(df_LOCK$value)
  min_LOCK = mean(df_LOCK$value)-3*sd(df_LOCK$value)
  max_Shuffle = mean(df_Shuffle$value)+3*sd(df_Shuffle$value)
  min_Shuffle = mean(df_Shuffle$value)-3*sd(df_Shuffle$value)
  df_LOCK = filter(df_LOCK, value > min_LOCK & value < max_LOCK)
  df_Shuffle = filter(df_Shuffle, value > min_Shuffle & value < max_Shuffle)
  df2 = rbind(df_LOCK, df_Shuffle)
  
  plotData[[nMD]] = df2
  p[[nMD]] = ggplot(plotData[[nMD]], aes(x = group, y = value))+
    geom_boxplot()+
    theme_classic()+
    xlab('')+
    ylab('Number of loops in LOCKs/\nnumber of LOCKs')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))+
    ggtitle(MD_options[nMD])+
    geom_signif(comparisons = list(c('LOCK', 'Shuffle')),
                test = 't.test', map_signif_level = T, test.args = c(var.equal = T,paired = T))#+
  # ylim(0,35)
}
plot_grid(p[[1]], p[[2]], p[[3]], nrow = 1)


names(plotData) = sub('common','common ',names(loopNumberPerLOCKNumber))
plots = list()
for (i in 1:length(plotData)) {
  plotData[[i]]$compare = factor(plotData[[i]]$compare, levels = c('LOCK > Shuffle', 'LOCK <= Shuffle'))
  plots[[i]] = ggplot(plotData[[i]], aes(x = group, y = value))+
    facet_grid(.~facet)+
    geom_point()+
    geom_line(aes(group=id, color = compare))+
    theme_classic()+
    xlab('')+
    ylab('Number of loops in LOCKs/\nnumber of LOCKs')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom', legend.title = element_blank(),
          strip.background = element_blank(),
          axis.text = element_text(color = 'black'))+
    geom_signif(comparisons = list(c('LOCK', 'Shuffle')),
                test = 't.test', map_signif_level = T, test.args = c(var.equal = T,paired = T))
}
Fig2D = ggarrange(plots[[1]], plots[[2]]+ylab(''), plots[[3]]+ylab(''), nrow = 1, widths = c(1.2,1,1))

pdf('plot/Fig2D/Fig2D.pdf')
Fig2D
dev.off()

Data_Fig2D_Genome = plotData[[1]] %>% dplyr::rename(RelativeLoopDensity = value)
Data_Fig2D_commonPMD = plotData[[2]] %>% dplyr::rename(RelativeLoopDensity = value)
Data_Fig2D_commonHMD = plotData[[3]] %>% dplyr::rename(RelativeLoopDensity = value)

write.table(Data_Fig2D_Genome, file = 'plot/Fig2D/Data_Fig2D_Genome.txt', sep = '\t', quote = F, row.names = F, col.names = T)
write.table(Data_Fig2D_commonPMD, file = 'plot/Fig2D/Data_Fig2D_commonPMD.txt', sep = '\t', quote = F, row.names = F, col.names = T)
write.table(Data_Fig2D_commonHMD, file = 'plot/Fig2D/Data_Fig2D_commonHMD.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# FigS2B/FigS2C/Fig2E/FigS2F. The enrichment of H3K27me3 or H3K9me3 long LOCKs in S/I/L-PMDs and common HMDs.----
# load(file = 'data_Roadmap/RDatas/LOCK_109.RData')
# load(file = 'data_Roadmap/RDatas/LOCK_109_H3K9me3.RData')
SIL_region_and_annotation = read.table(file = 'meta/commonPMD_subgroup/04_PMD_length_Plots/SIL_region_and_annotation.bed', header = T)
SIL_region_and_annotation = SIL_region_and_annotation[,c('chrom','start','end','PMD_length_anno')]
S_PMD = filter(SIL_region_and_annotation, PMD_length_anno == 'S')
I_PMD = filter(SIL_region_and_annotation, PMD_length_anno == 'I')
L_PMD = filter(SIL_region_and_annotation, PMD_length_anno == 'L')
commonPMD = read.table('meta/commonPMD_commonHMD/commonPMD_hg19.bed') %>% 
  dplyr::rename(chrom=1, start=2, end=3) %>% 
  mutate(PMD_length_anno = 'common\nPMDs')
commonHMD = read.table('meta/commonPMD_commonHMD/commonHMD_hg19.bed') %>% 
  dplyr::rename(chrom=1, start=2, end=3) %>% 
  mutate(PMD_length_anno = 'common\nHMDs')
SIL_domain_size = tribble(
  ~group_MD, ~domain_size,
  'S', sum(S_PMD$end - S_PMD$start)/1e8,
  'I', sum(I_PMD$end - I_PMD$start)/1e8,
  'L', sum(L_PMD$end - L_PMD$start)/1e8,
  'common\nPMDs', sum(commonPMD$end - commonPMD$start)/1e8,
  'common\nHMDs', sum(commonHMD$end - commonHMD$start)/1e8
)

plotDataLOCKEnrichSIL = function(LOCK_file){
  plotData_SIL = list()
  for (s109 in 1:109) {
    print(s109)
    EID = names(LOCK_file)[s109]
    df = LOCK_file[[s109]] %>% 
      mutate(LOCK_length = ifelse(d_LOCK>100000, 'Long_LOCK','Short_LOCK'))
    
    all_LOCK_size = group_by(df,LOCK_length) %>% 
      summarise(d = sum(d_LOCK))
    
    S_overlap = bed_intersect(df[,c( "chrom", "start", "end", "LOCK_length")],S_PMD,suffix = c('','.PMD')) %>% 
      group_by(LOCK_length) %>% 
      summarise(overlap = sum(.overlap)) %>% 
      mutate(group_MD = 'S')
    I_overlap = bed_intersect(df[,c( "chrom", "start", "end", "LOCK_length")],I_PMD,suffix = c('','.PMD')) %>% 
      group_by(LOCK_length) %>% 
      summarise(overlap = sum(.overlap)) %>% 
      mutate(group_MD = 'I')
    L_overlap = bed_intersect(df[,c( "chrom", "start", "end", "LOCK_length")],L_PMD,suffix = c('','.PMD')) %>% 
      group_by(LOCK_length) %>% 
      summarise(overlap = sum(.overlap)) %>% 
      mutate(group_MD = 'L')
    commonPMD_overlap = bed_intersect(df[,c( "chrom", "start", "end", "LOCK_length")],commonPMD,suffix = c('','.PMD')) %>% 
      group_by(LOCK_length) %>% 
      summarise(overlap = sum(.overlap)) %>% 
      mutate(group_MD = 'common\nPMDs')
    commonHMD_overlap = bed_intersect(df[,c( "chrom", "start", "end", "LOCK_length")],commonHMD,suffix = c('','.HMD')) %>% 
      group_by(LOCK_length) %>% 
      summarise(overlap = sum(.overlap)) %>% 
      mutate(group_MD = 'common\nHMDs')
    
    overlap_size = rbind(S_overlap,I_overlap,L_overlap, commonPMD_overlap, commonHMD_overlap)
    
    plotData_SIL[[s109]] = left_join(overlap_size,all_LOCK_size,by = 'LOCK_length') %>% 
      left_join(SIL_domain_size, by = "group_MD") %>% 
      mutate(var = (overlap/d)/domain_size, sample = EID) %>% 
      select(sample, group_MD, LOCK_length, var) %>% 
      dplyr::rename( group = 3)
  }
  plotData_SIL = do.call(rbind,plotData_SIL)
  return(plotData_SIL)
}
LOCKEnrichSIL_K27 = plotDataLOCKEnrichSIL(LOCK_109)
LOCKEnrichSIL_K9 = plotDataLOCKEnrichSIL(LOCK_109_H3K9me3)

plotLOCKEnrichPMDHMD = function(plotData) {
  df = plotData
  df$group = sub('_LOCK',' LOCKs',df$group)
  df$group = factor(df$group, levels = c('Long LOCKs','Short LOCKs'))
  df = filter(df, group_MD %in% c('common\nPMDs','common\nHMDs'))
  df$group_MD = factor(df$group_MD, levels = c('common\nPMDs','common\nHMDs'))
  p = ggplot(df, aes(x = group_MD, y = var))+
    geom_boxplot(aes(color = group))+
    facet_grid(.~group)+
    geom_signif(comparisons = list(c('common\nPMDs','common\nHMDs')), 
                map_signif_level = T, step_increase = 0.05, test = 't.test', test.args = c(paired = T, var.equal = T))+
    theme_classic()+
    theme(axis.text = element_text(color = 'black'),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab('')+
    ylab('Enrichment score')+
    ggtitle(paste0(a[h],' (LOCK level)'))+
    guides(color = 'none')+
    scale_color_manual(values = c('Long LOCKs' = 'red', 'Short LOCKs' = 'orange'))
  return(list(df, p))
}
LOCKEnrichPMDHMD_K27 = plotLOCKEnrichPMDHMD(LOCKEnrichSIL_K27)
FigS2B = LOCKEnrichPMDHMD_K27[[2]]; Data_FigS2B = LOCKEnrichPMDHMD_K27[[1]]
LOCKEnrichPMDHMD_K9 = plotLOCKEnrichPMDHMD(LOCKEnrichSIL_K9)
FigS2C = LOCKEnrichPMDHMD_K9[[2]]; Data_FigS2C = LOCKEnrichPMDHMD_K9[[1]]

pdf('plot/FigS2B/FigS2B.pdf', width = 4, height = 5)
FigS2B
dev.off()

pdf('plot/FigS2C/FigS2C.pdf', width = 4, height = 5)
FigS2C
dev.off()

write.table(Data_FigS2B, file = 'plot/FigS2B/Data_FigS2B.txt', sep = '\t', quote = F, row.names = F, col.names = T)
write.table(Data_FigS2C, file = 'plot/FigS2C/Data_FigS2C.txt', sep = '\t', quote = F, row.names = F, col.names = T)

plotLOCKEnrichSIL = function(plotData) {
  df = plotData
  df$group = sub('_LOCK',' LOCKs',df$group)
  df = filter(df, group == 'Long LOCKs')
  df$group_MD = paste0(df$group_MD,'-PMDs')
  df$group_MD = sub('common\nHMDs-PMDs','common\nHMDs',df$group_MD)
  df = filter(df, group_MD %in% c('S-PMDs','I-PMDs','L-PMDs','common\nHMDs'))
  df$group_MD = factor(df$group_MD, levels = c('S-PMDs','I-PMDs','L-PMDs','common\nHMDs'))
  p = ggplot(df,aes(x = group_MD, y = var))+
    geom_violin(aes(fill = group_MD))+
    geom_boxplot(color = 'black', width=0.1, outlier.size = 0.1)+
    ylab('Enrichment score')+
    xlab('')+
    theme_classic()+
    theme(axis.text = element_text(color = 'black'),
          plot.title = element_text(hjust = 0.5))+
    guides(color = 'none', fill = 'none')+
    geom_signif(comparisons = list( c('common\nHMDs','S-PMDs'),
                                    c('common\nHMDs','I-PMDs'),
                                    c('common\nHMDs','L-PMDs'),
                                    c('S-PMDs','I-PMDs'),
                                    c('I-PMDs','L-PMDs'),
                                    c('S-PMDs','L-PMDs')
    ),
    map_signif_level = T,
    step_increase = 0.05,
    test = "t.test",test.args = c(var.equal = T,paired = T))+
    scale_fill_manual(values = c("#92D2C3", "#DDDF98", "#8FBAD9","#EB9293"))
  return(list(df,p))
}

LOCKEnrichSIL_K27 = plotLOCKEnrichSIL(LOCKEnrichSIL_K27)
Fig2F = LOCKEnrichSIL_K27[[2]]; Data_Fig2F = LOCKEnrichSIL_K27[[1]]
LOCKEnrichSIL_K9 = plotLOCKEnrichSIL(LOCKEnrichSIL_K9)
FigS2F = LOCKEnrichSIL_K9[[2]]; Data_FigS2F = LOCKEnrichSIL_K9[[1]]

pdf('plot/Fig2F/Fig2F.pdf', width = 4, height = 5)
Fig2F
dev.off()

pdf('plot/FigS2F/FigS2F.pdf', width = 4, height = 5)
FigS2F
dev.off()

write.table(Data_Fig2F, file = 'plot/Fig2F/Data_Fig2F.txt', sep = '\t', quote = F, row.names = F, col.names = T)
write.table(Data_FigS2F, file = 'plot/FigS2F/Data_FigS2F.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig2G. The top 20 pathways most enriched in genes within long LOCKs in S/I/L-PMDs.----
load('data_Roadmap/RDatas/samples.RData')
# get top terms in S/I/L respectively
GO_All_SIL = list()
for (ex in 1:length(Roadmap_table2$EID)) {
  EID = Roadmap_table2$EID[ex]
  GO_AllLists = read.csv(file = paste0('data_Roadmap/10_metascape_and_LOLA/Metascape_output//H3K27me3_',
                                       EID,'/Enrichment_GO/GO_AllLists.csv'),
                         header = T)
  GO_All_SIL[[ex]] = filter(GO_AllLists, GeneList %in% c('Long_LOCK__S', 'Long_LOCK__I', 'Long_LOCK__L')) %>% 
    mutate(EID = EID)
}

GO_All_SIL = do.call(rbind, GO_All_SIL)
GO_All_SIL$term = paste0(GO_All_SIL$GO,': ',GO_All_SIL$Description)
GO_All_SIL$LogP = -GO_All_SIL$LogP
GO_All_SIL$Log.q.value. = -GO_All_SIL$Log.q.value.
GO_All_SIL = filter(GO_All_SIL, grepl('GO',GO))
topTerms_SIL = group_by(GO_All_SIL, GeneList, term) %>% 
  summarise(sumLogP = sum(Log.q.value.)) %>% 
  mutate(GeneList = factor(GeneList, levels = c('Long_LOCK__S', 'Long_LOCK__I', 'Long_LOCK__L'))) %>% 
  arrange(GeneList, desc(sumLogP)) %>% 
  top_n(10) %>% 
  ungroup() %>% 
  distinct(term) %>% 
  pull(term)
GO_All_SIL = filter(GO_All_SIL, term %in% topTerms_SIL)
GO_All_SIL$GeneList = sub('Long_LOCK__S','Peaks of LOCK in S',GO_All_SIL$GeneList)
GO_All_SIL$GeneList = sub('Long_LOCK__I','Peaks of LOCK in I',GO_All_SIL$GeneList)
GO_All_SIL$GeneList = sub('Long_LOCK__L','Peaks of LOCK in L',GO_All_SIL$GeneList)
GO_All_SIL$GeneList = factor(GO_All_SIL$GeneList, levels = c('Peaks of LOCK in S', 'Peaks of LOCK in I', 'Peaks of LOCK in L'))
GO_All_SIL1 = GO_All_SIL
GO_All_SIL1$EID_GeneList_term = paste(GO_All_SIL1$EID, GO_All_SIL1$GeneList, GO_All_SIL1$term, sep = '.')
EID_GeneList_SIL = c(paste0(Roadmap_table2$EID,'.','Peaks of LOCK in S'),paste0(Roadmap_table2$EID,'.','Peaks of LOCK in I'), paste0(Roadmap_table2$EID,'.','Peaks of LOCK in L'))
AllDots_SIL = list()
for (i in 1:length(EID_GeneList_SIL)) {
  AllDots_SIL[[i]] = paste0(EID_GeneList_SIL[i],'.',topTerms_SIL)
}
if (!is.function(c)) {rm(c)}
AllDots_SIL = do.call(c, AllDots_SIL)

addtable1_SIL = data.frame(EID_GeneList_term = AllDots_SIL[!AllDots_SIL %in% GO_All_SIL1$EID_GeneList_term])
addtable1_SIL = mutate(addtable1_SIL, toSeparate = EID_GeneList_term) %>% 
  separate(col = 'toSeparate', into = c('EID','GeneList','term'), sep = '\\.') %>% 
  select(GeneList, EID, term, EID_GeneList_term)
addtable2_SIL = as.data.frame(matrix(data = NA, nrow = nrow(addtable1_SIL), ncol = ncol(GO_All_SIL1)-4))
colnames(addtable2_SIL) = colnames(GO_All_SIL1)[1:(length(colnames(GO_All_SIL1))-4)]
addtable_SIL = cbind(addtable2_SIL, addtable1_SIL)
GO_All_SIL1 = rbind(GO_All_SIL1, addtable_SIL)
GO_All_SIL_heatmap = GO_All_SIL1[,c('term','GeneList','EID','Log.q.value.')] %>% 
  mutate(term = factor(term, levels = topTerms_SIL)) %>% 
  arrange(term) %>% 
  pivot_wider(names_from = 'term', values_from = 'Log.q.value.')
GO_All_SIL_heatmap$GeneList = factor(GO_All_SIL_heatmap$GeneList, levels = c('Peaks of LOCK in S', 'Peaks of LOCK in I', 'Peaks of LOCK in L'))
GO_All_SIL_heatmap = arrange(GO_All_SIL_heatmap, GeneList, EID)
# sort by P-value sums.
GO_All_SIL_heatmap = GO_All_SIL_heatmap %>% 
  mutate(.,rowSum = rowSums(.[,-1:-2], na.rm = T)) %>% 
  group_by(GeneList) %>% 
  arrange(GeneList, desc(rowSum)) %>% 
  select(!rowSum)
GO_All_SIL_heatmap$GeneList = sub('Peaks of LOCK in ','Long LOCK peaks in ',GO_All_SIL_heatmap$GeneList)
GO_All_SIL_heatmap$GeneList = paste0(GO_All_SIL_heatmap$GeneList, '-PMDs')
dat_SIL = GO_All_SIL_heatmap[,-1:-2] %>% t()
rownames(dat_SIL) = sub('GO.*: ','',rownames(dat_SIL))
ht_SIL = Heatmap(dat_SIL,cluster_rows = F,cluster_columns = F,
                 col = colorRamp2(c(2,73),c('white','red')),
                 top_annotation = columnAnnotation(group = GO_All_SIL_heatmap$GeneList,
                                                   col = list(group = c("Long LOCK peaks in S-PMDs" = "#92d2c3", 
                                                                        "Long LOCK peaks in I-PMDs" = "#dddf98",
                                                                        "Long LOCK peaks in L-PMDs" = "#8fbad9"))),
                 heatmap_legend_param = list(title = '-log10 (q-Value)',
                                             col_fun = colorRamp2(c(2,73),c('white','red')),
                                             title = "test", at = c(2,73)),
                 column_title = "Sample",
                 row_names_side = "left",
                 row_names_max_width = max_text_width(
                   rownames(dat_SIL),
                   gp = gpar(fontsize = 12)
                 )
                 
)
Fig2G = draw(ht_SIL, merge_legend = TRUE)

pdf('plot/Fig2G/Fig2G.pdf',height = 7, width = 13)
Fig2G
dev.off()

Data_Fig2G = GO_All_SIL_heatmap
write.table(Data_Fig2G, file = 'plot/Fig2G/Data_Fig2G.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig2H.  The beeswarm plot shows the proportion of OGs or TSGs in Long LOCKs among all genes in Long LOCKs, separately for common PMDs/HMDs and S/I/L-PMDs.----
load('data_Roadmap/RDatas/samples.RData')
load('data_Roadmap/RDatas/gene_49.RData')
my_N_top = 1500
gene_sampleList = gene_49
allGene_genome = gene_sampleList[[1]]$gene_name
allGene_commonPMD = gene_sampleList[[1]] %>%
  filter(MD_anno == 'commonPMD') %>%
  pull(gene_name)
allGene_commonHMD = gene_sampleList[[1]] %>%
  filter(MD_anno == 'commonHMD') %>%
  pull(gene_name)
get_OG_TSG_exp_plotdata = function(MD, N_top){
  exp_res = list()
  example_OG_in_S = list()
  for (s49 in 1:49) {
    print(s49)
    gene_df = gene_49[[s49]]
    EID = names(gene_49)[s49]
    if (MD == 'commonPMD' | MD == 'commonHMD') {
      gene_df = filter(gene_df, MD_anno == MD)
    }
    if (MD == 'S' | MD == 'I' | MD == 'L') {
      gene_df = filter(gene_df, PMD_length_anno == MD)
    }
    genestoBeFiltered = allGene_genome
    TSG = read.csv(file = 'meta/TSG_OG/TSG.csv',header = T) %>% filter(Gene %in% genestoBeFiltered) %>% .[1:N_top,] %>% pull(Gene)
    OG = read.csv(file = 'meta/TSG_OG/OG.csv',header = T) %>% filter(Gene %in% genestoBeFiltered) %>% .[1:N_top,] %>% pull(Gene)
    TSG_OG_anno = rbind(data.frame(gene_name = TSG, TSG_OG_anno = 'TSG'),
                        data.frame(gene_name = OG, TSG_OG_anno = 'OG'))
    gene_df1 = left_join(gene_df, TSG_OG_anno, by = "gene_name")
    gene_df1[is.na(gene_df1$LOCK_length),]$LOCK_length = 'outside_of_LOCK'
    gene_df1[is.na(gene_df1$TSG_OG_anno),]$TSG_OG_anno = 'other'
    
    #example_OG_in_S
    if (MD == 'S') {
      example_OG_in_S[[s49]] = filter(gene_df1, TSG_OG_anno == 'OG' & LOCK_length == 'Long_LOCK') %>% 
        mutate(EID = EID)
    }
    
    exp_res[[s49]] = rbind(left_join(gene_df1 %>% 
                                                 group_by(LOCK_length, TSG_OG_anno) %>% 
                                                 summarise(exp = median(exp, na.rm = T), .groups = 'drop',
                                                           gene_number = n()) %>% 
                                                 filter(TSG_OG_anno != 'other'),
                                               gene_df1 %>% 
                                                 filter(!is.na(exp)) %>% 
                                                 group_by(LOCK_length, TSG_OG_anno) %>% 
                                                 summarise(gene_number_with_exp = n(), .groups = 'drop') %>% 
                                                 filter(TSG_OG_anno != 'other'),
                                               by = c('LOCK_length', 'TSG_OG_anno')),
                                     
                                     left_join(gene_df1 %>% 
                                                 group_by(LOCK_length) %>% 
                                                 summarise(exp = median(exp, na.rm = T), 
                                                           gene_number = n(), .groups = 'drop'), 
                                               gene_df1 %>% 
                                                 filter(!is.na(exp)) %>% 
                                                 group_by(LOCK_length) %>% 
                                                 summarise(gene_number_with_exp = n(), .groups = 'drop'),
                                               by = 'LOCK_length')%>% 
                                       mutate(TSG_OG_anno = 'all', .before = 2)) %>% 
      mutate(group = paste0(TSG_OG_anno,'__',LOCK_length)) %>% 
      mutate(group = factor(group, levels = c("OG__Long_LOCK", "OG__Short_LOCK", "OG__outside_of_LOCK", 
                                              "TSG__Long_LOCK", "TSG__Short_LOCK", "TSG__outside_of_LOCK", 
                                              "all__Long_LOCK", "all__Short_LOCK", "all__outside_of_LOCK" ))) %>% 
      mutate(EID = EID, .before = 1)
  }
  exp_res = do.call(rbind, exp_res)
  overallPercentage = data.frame(group = c('OG', 'TSG'), 
                                 value = c(nrow(filter(gene_df1, TSG_OG_anno == 'OG'))/nrow(gene_df1),
                                           nrow(filter(gene_df1, TSG_OG_anno == 'TSG'))/nrow(gene_df1))
  )
  
  #example_OG_in_S
  if (MD == 'S') {
    example_OG_in_S = do.call(rbind,example_OG_in_S)
  }
  
  return(list(exp_res, overallPercentage, example_OG_in_S))
}
exp_res = list()
overallPercentage = list()
all_MD = c('Genome', 'commonHMD', 'commonPMD', 'S', 'I', 'L')
for (i in 1:6) {
  MD = all_MD[i]
  resList = get_OG_TSG_exp_plotdata(MD = MD, N_top = my_N_top)
  exp_res[[i]] = resList[[1]] %>% mutate(MD_anno = MD)
  overallPercentage[[i]] = resList[[2]]
  
  if (MD == 'S') {
    example_OG_in_S = resList[[3]]
  }
}
names(exp_res) = all_MD
names(overallPercentage) = all_MD

for (i in 1:length(overallPercentage)) {
  overallPercentage[[i]] = mutate(overallPercentage[[i]], 
                                  EID = 'Overall',
                                  MD_anno = sub('common','common ',names(overallPercentage)[i]),
                                  .before = 1)
  colnames(overallPercentage[[i]]) = c('EID', 'MD_anno', 'TSG_OG_anno','OgOrTsgToAllGenes')
}
OverallPer = rbind(overallPercentage[[2]], overallPercentage[[3]], overallPercentage[[4]], overallPercentage[[5]], overallPercentage[[6]])
OverallPer$MD_anno = factor(OverallPer$MD_anno, levels = c('common PMD', 'common HMD','S','I','L'))

#plotData
plotData = rbind(exp_res[[2]], exp_res[[3]],exp_res[[4]], exp_res[[5]],exp_res[[6]])
plotData$MD_anno = sub('common','common ',plotData$MD_anno)
plotData$MD_anno = factor(plotData$MD_anno, levels = c('common PMD', 'common HMD', 'S','I','L'))

plotData = filter(plotData, LOCK_length == 'Long_LOCK') %>% 
  select(TSG_OG_anno, gene_number, EID, MD_anno) %>% 
  pivot_wider(names_from = 'TSG_OG_anno', values_from = 'gene_number') %>% 
  mutate(OG = OG/all, TSG = TSG/all) %>% 
  select(EID, OG, TSG, MD_anno) %>% 
  pivot_longer(cols = 2:3, names_to = 'TSG_OG_anno', values_to = 'OgOrTsgToAllGenes')

Fig2H = ggplot(plotData, aes(x = TSG_OG_anno, y = OgOrTsgToAllGenes))+
  facet_grid(.~MD_anno)+
  # geom_boxplot()+
  ggbeeswarm::geom_beeswarm(color = '#eb8100',shape = 18, cex = 3, method = "swarm")+
  geom_point(data = OverallPer, aes(x = TSG_OG_anno, y = OgOrTsgToAllGenes), color = 'purple',shape = 19, size = 3)+
  xlab('')+
  ylab('Ratio of OGs or TSGs to\nall genes within long LOCKs')+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        axis.text = element_text(color = 'black'))+
  geom_signif(comparisons = list( c('OG','TSG')),
              map_signif_level = T,
              step_increase = 0.05,
              test = "t.test",test.args = c(var.equal = T,paired = T))

pdf('plot/Fig2H/Fig2H.pdf', height = 5, width = 9)
Fig2H
dev.off()

Data_Fig2H = rbind(OverallPer, plotData)
write.table(Data_Fig2H, file = 'plot/Fig2H/Data_Fig2H.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig2I. The density plot shows the proportion of OGs among all genes in Long LOCKs, separately for S/I/L-PMDs.----
all_MD2 = sub('common','common ',all_MD)
plotData_percentage2 = list()
for (i in 1:6) {
  df = exp_res[[i]]
  plotData_percentage2[[i]] = list()
  for (j in 1:2) {
    TSG_or_OG = c('OG', 'TSG')[j]
    df2 = filter(df, TSG_OG_anno == TSG_or_OG)
    Numerator = filter(df2, LOCK_length == 'Long_LOCK') %>% select(EID, gene_number)
    Denominator = group_by(df2, EID) %>% 
      summarise(gene_nuber = sum(gene_number, na.rm = T), .groups = 'drop')
    plotData_percentage2[[i]][[j]] = merge(Numerator, Denominator, by = 'EID') %>% 
      dplyr::rename(Numerator = 2, Denominator = 3) %>% 
      mutate(value = Numerator/Denominator) %>% 
      select(EID, value) %>% 
      mutate(MD_anno = all_MD2[i], group = TSG_or_OG)
  }
  plotData_percentage2[[i]] = do.call(rbind, plotData_percentage2[[i]])
}
plotData_percentage2 = do.call(rbind, plotData_percentage2)
plotData_percentage2$MD_anno = factor(plotData_percentage2$MD_anno, levels = all_MD2)
plotData_percentage2 = filter(plotData_percentage2, group == 'OG')
p_FengLuan = ggplot(plotData_percentage2, aes(x = MD_anno, y = value, color = MD_anno))+
  facet_grid(.~group)+
  geom_boxplot()+
  theme_classic()+
  xlab('')+
  ylab('Number of OGs in LOCK/\nnumber of all OGs')+
  theme(strip.background = element_blank(), strip.text = element_text(size=14, face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))
library(ggridges)
library(RColorBrewer)
plotData_percentage2$MD_anno = factor(plotData_percentage2$MD_anno, levels = rev(all_MD2))
plotData_percentage2 = filter(plotData_percentage2, MD_anno %in% c('S','I','L'))
plotData_percentage2$MD_anno = factor(plotData_percentage2$MD_anno, levels = rev(c('S','I','L')))

Fig2I = ggplot(plotData_percentage2, aes(x = value, y = MD_anno, fill = ..density..)) + 
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.00,size = 0.3) + 
  scale_fill_gradientn(colours = c('yellow','red'))+
  theme_classic()+
  ylab('')+
  xlab('# of OGs in long LOCKs / # of all OGs')+
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'bottom')

pdf('plot/Fig2I/Fig2I.pdf', width = 4, height = 5)
Fig2I
dev.off()

Data_Fig2I = plotData_percentage2
write.table(Data_Fig2I, file = 'plot/Fig2I/Data_Fig2I.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig2J. Comparison of the expression of OGs or TSGs within long LOCKs to those outside of long or short LOCKs, separately for common PMDs or HMDs.----
get_OG_TSG_exp_plotData = function(df){
  df = filter(df, LOCK_length != 'Short_LOCK')
  df$LOCK_length = sub('Long_LOCK','LOCK',df$LOCK_length)
  df$LOCK_length = sub('outside_of_LOCK','non-LOCK',df$LOCK_length)
  df$LOCK_length = factor(df$LOCK_length, levels = c('LOCK','non-LOCK'))
  df = filter(df, TSG_OG_anno != 'all')
  df$TSG_OG_anno = factor(df$TSG_OG_anno, levels = c('OG','TSG'))
  for (gene in c('OG','TSG')) {
    x = filter(df, TSG_OG_anno == gene & LOCK_length == 'LOCK')$EID
    y = filter(df, TSG_OG_anno == gene & LOCK_length == 'non-LOCK')$EID
    if (length(x) != length(y)) {
      delList = table(c(x,y)) %>% as.data.frame() %>%
        filter(Freq == 1) %>% pull(Var1) %>%
        as.character()
      df = filter(df, !(TSG_OG_anno == gene & EID %in% delList))
    }
  }
  return(df)
}
get_OG_TSG_exp_plot = function(df,title,ymax = 1.5){
  p = ggplot(df, aes(x = LOCK_length, y = exp))+
    geom_boxplot(aes(color = LOCK_length))+
    facet_grid(.~TSG_OG_anno)+
    xlab('')+
    ylab('log10(RPKM+0.1)')+
    theme_classic()+
    theme(axis.text = element_text(colour = 'black'),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank(), strip.text = element_text(size=14))+
    scale_color_manual(values = c('LOCK' = 'red',
                                  'non-LOCK' = '#f9d126'))+
    geom_signif(comparisons = list(c('LOCK','non-LOCK')),
                map_signif_level = T, test = "t.test",test.args = c(var.equal = T,paired = T))+
    guides(color = 'none')+
    ggtitle(title)+
    ylim(-1,ymax)
  return(p)
}
plotData_exp = list()
p_exp = list()
for (i in 1:3) {
  plotData_exp[[i]] = get_OG_TSG_exp_plotData(exp_res[[i]])
  p_exp[[i]] = get_OG_TSG_exp_plot(plotData_exp[[i]], all_MD2[i])
}
Fig2J = ggarrange(p_exp[[3]]+ggtitle('Within common PMD'),
                  p_exp[[2]]+ggtitle('Within common HMD'),
                  nrow = 1)
pdf('plot/Fig2J/Fig2J.pdf', height = 5, width = 9)
Fig2J
dev.off()

Data_Fig2J_commonPMD = plotData_exp[[3]]
Data_Fig2J_commonHMD = plotData_exp[[2]]

write.table(Data_Fig2J_commonPMD, file = 'plot/Fig2J/Data_Fig2J_commonPMD.txt', sep = '\t', quote = F, row.names = F, col.names = T)
write.table(Data_Fig2J_commonHMD, file = 'plot/Fig2J/Data_Fig2J_commonHMD.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig2K. Consistent with the content depicted in Fig.2J, but for OGs in S/I/L-PMDs.----
df = rbind(exp_res[[4]], exp_res[[5]], exp_res[[6]])
df = filter(df, LOCK_length != 'Short_LOCK')
df$LOCK_length = sub('Long_LOCK','LOCK',df$LOCK_length)
df$LOCK_length = sub('outside_of_LOCK','non-LOCK',df$LOCK_length)
df$LOCK_length = factor(df$LOCK_length, levels = c('LOCK','non-LOCK'))
df = filter(df, TSG_OG_anno == 'OG')
x = filter(df, LOCK_length == 'LOCK')$EID
y = filter(df, LOCK_length == 'non-LOCK')$EID
if (length(x) != length(y)) {
  delList = table(c(x,y)) %>% as.data.frame() %>%
    filter(Freq == 1) %>% pull(Var1) %>%
    as.character()
  df = filter(df, !(EID %in% delList))
}
df$MD_anno = paste0('Within ',df$MD_anno)
df$MD_anno = factor(df$MD_anno, levels = paste0('Within ',c('S','I','L')))
Fig2K = ggplot(df, aes(x = LOCK_length, y = exp))+
  geom_boxplot(aes(color = LOCK_length))+
  facet_grid(.~MD_anno)+
  xlab('')+
  ylab('log10 (RPKM+0.1)')+
  theme_classic()+
  theme(axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(), strip.text = element_text(size=14))+
  scale_color_manual(values = c('LOCK' = 'red',
                                'non-LOCK' = '#f9d126'))+
  geom_signif(comparisons = list(c('LOCK','non-LOCK')),
              map_signif_level = T, test = "t.test",test.args = c(var.equal = T,paired = T))+
  guides(color = 'none')+
  ggtitle('OG')

pdf('plot/Fig2K/Fig2K.pdf', height = 5, width = 5)
Fig2K
dev.off()

Data_Fig2K = select(df, EID, group, TSG_OG_anno, exp, MD_anno)
write.table(Data_Fig2K, file = 'plot/Fig2K/Data_Fig2K.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# FigS2D/FigS2E. ----
load('00_commonPMD_subgroup/04_PMD_length_Plots/PMD_length_cluster.RData')
FigS2D = PMD_length_cluster[[2]]
FigS2E = PMD_length_cluster[[3]]

# Combine the plots in Figure 2----
line1 = ggarrange(Fig2A+BigWordTheme+StripTheme, 
                  Fig2B+BigWordTheme+StripTheme,
                  Fig2C+BigWordTheme+StripTheme, 
                  Fig2E+BigWordTheme+StripTheme, 
                  Fig2F+BigWordTheme+StripTheme, 
                  nrow = 1, labels = c('A','B','C','E','F'), font.label = list(size = 23))
line2 = ggarrange(Fig2D+BigWordTheme+StripTheme, Fig2G, nrow = 1, labels = c('D','G'), font.label = list(size = 23), widths = c(1,1.8))
line3_left = ggarrange(ggarrange(Fig2H+BigWordTheme+StripTheme, Fig2I+BigWordTheme+StripTheme, nrow = 1, labels = c('H','I'), font.label = list(size = 23), widths = c(2,1)),
                       ggarrange(Fig2J+BigWordTheme+StripTheme, Fig2K+BigWordTheme+StripTheme, nrow = 1, labels = c('J','K'), font.label = list(size = 23), widths = c(1.5,1)),
                       ncol = 1)
line3 = ggarrange(line3_left, Empty, nrow = 1, widths = c(2,1), labels = c('','L'), font.label = list(size = 23))
Figure_2 = ggarrange(line1, line2, line3, nrow = 3, heights = c(1.5,1.5,2))

pdf('20_output_plots/Figure_2.pdf',width = 13,height = 15)
Figure_2
dev.off()

pdf('20_output_plots/Fig2G.pdf', width = 8,height = 5)
Fig2G
dev.off()

# Combine the plots in Figure S2----
Figure_S2 = ggarrange(
  ggarrange(FigS2B, FigS2C, nrow = 1, labels = c('B','C'), font.label = list(size = 23)),
  FigS2D,
  ggarrange(FigS2E, FigS2F, nrow = 1, labels = c('E','F'), font.label = list(size = 23)),
  ncol = 1, labels = c('','D',''), font.label = list(size = 23)
)

pdf('20_output_plots/Figure_S2.pdf',width = 13,height = 13)
Figure_S2
dev.off()

pdf('20_output_plots//FigS2D.pdf',width = 13,height = 6.5)
FigS2D
dev.off()

# Figure 3/ Figure S3/ Figure S4----
commonPMD_hg38 = read.table(file = 'meta/commonPMD_commonHMD/commonPMD_hg38.bed', sep = '\t', header = F) %>% mutate(MD_anno = 'commonPMD'); colnames(commonPMD_hg38)[1:3] = c('chrom', 'start', 'end')
commonHMD_hg38 = read.table(file = 'meta/commonPMD_commonHMD/commonHMD_hg38.bed', sep = '\t', header = F) %>% mutate(MD_anno = 'commonHMD'); colnames(commonHMD_hg38)[1:3] = c('chrom', 'start', 'end')
blacklist_hg38 = read.table(file = 'meta/blacklist/hg38-blacklist.v2.bed', sep = '\t', header = F)[,c(1:3)]; colnames(blacklist_hg38) = c('chrom', 'start', 'end')
load('meta/RefSeq_gene/hg38_Refseq_gene.RData')
load('data_CellLine/RDatas/Pairs.RData')
all_tissue_types = pairs_for_longLOCKAnalysis
all_tissue_type_names = c('NE2 vs KYSE450 (ESCC)',
                          'NE2 vs KYSE510 (ESCC)',
                          'Normal vs Her2+ (BRCA)',
                          'Normal vs Luminal A (BRCA)',
                          'Normal vs TN-Basal (BRCA)',
                          'Normal vs TN-CL (BRCA)')


# Fig3A/FigS3A. Changes in the enrichment of H3K27me3 or H3K9me3 long LOCKs in common PMDs or HMDs between normal and tumor cell lines.----
NormalizedEnrichment = T
all_plotData = list()
for(t in 1:6){
  LOCK_TN_path = 'data_CellLine/04_LOCK/withAnnotations/'
  tissue_type = all_tissue_types[t]
  files = list.files(LOCK_TN_path, pattern = paste0(tissue_type))
  enrichment_fun = function(LOCK_df,
                            usingLOCK = 'Short_LOCK',
                            usingMetDomain = 'commonPMD'){
    LOCK_bed = filter(LOCK_df, LOCK_length == usingLOCK)[,1:3]
    if (usingMetDomain == 'commonPMD') {MD_bed = commonPMD_hg38[,1:3]}
    if (usingMetDomain == 'commonHMD') {MD_bed = commonHMD_hg38[,1:3]}
    overlap = sum(bed_intersect(LOCK_bed, MD_bed)$.overlap)
    total_LOCK = sum(LOCK_bed$end - LOCK_bed$start)
    total_MD = sum(MD_bed$end - MD_bed$start)
    
    return(list(overlap/total_MD,
                overlap/total_LOCK,
                overlap/total_LOCK/total_MD*1e8))
  }
  plotData = list()
  for (i in 1:length(files)) {
    LOCK_df = read.table(file = paste0(LOCK_TN_path, files[i]), sep = '\t', header = T)
    whichData = ifelse(NormalizedEnrichment, 3, 2)
    res_PMD = enrichment_fun(LOCK_df = LOCK_df, usingLOCK = 'Long_LOCK', usingMetDomain = 'commonPMD')[[whichData]]
    res_HMD = enrichment_fun(LOCK_df = LOCK_df, usingLOCK = 'Long_LOCK', usingMetDomain = 'commonHMD')[[whichData]]
    plotData[[i]] = data.frame(fileName = files[i],
                               enrichmentScore = c(res_HMD, res_PMD),
                               MD_anno = c('common HMD', 'common PMD'))
    
  }
  plotData = do.call(rbind, plotData)
  plotData[,c('TorN','HisMark')] = str_match(plotData$fileName, '.*?\\.(.*?)\\..*?\\.(.*?)\\..*')[,2:3]
  plotData$TorN = sub('N','Normal',sub('T','Tumor', plotData$TorN))
  all_plotData[[t]] = plotData %>% mutate(tissue_type = all_tissue_type_names[t])
}
all_plotData = do.call(rbind, all_plotData)
all_plotData$tissue_type = factor(all_plotData$tissue_type, levels = rev(all_tissue_type_names))
all_plotData$MD_anno = factor(all_plotData$MD_anno, levels = c('common PMD', 'common HMD'))
all_plotData$TorN = paste0(all_plotData$TorN, '\ncell lines')
all_plotData$TorN = factor(all_plotData$TorN, levels = c('Normal\ncell lines', 'Tumor\ncell lines'))
all_plotData_K27 = filter(all_plotData, HisMark == 'H3K27me3')
all_plotData_K9 = filter(all_plotData, HisMark == 'H3K9me3')

Fig3A = ggplot(all_plotData_K27, aes(x = TorN, y = enrichmentScore))+
  facet_grid(.~MD_anno)+
  geom_boxplot()+
  geom_point(color = 'red', size = 0.5)+
  theme_classic()+
  theme(axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank())+
  xlab('')+
  ylab('Enrichment score of\nH3K27me3 long LOCKs')+
  geom_signif(comparisons = list(c('Normal\ncell lines','Tumor\ncell lines')),
              map_signif_level = T, step_increase = T)

FigS3A = ggplot(all_plotData_K9, aes(x = TorN, y = enrichmentScore))+
  facet_grid(.~MD_anno)+
  geom_boxplot(outlier.size = 0.5)+
  geom_point(color = 'red', size = 0.5)+
  theme_classic()+
  theme(axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  xlab('')+
  ylab('Enrichment score of\nH3K9me3 long LOCKs')+
  geom_signif(comparisons = list(c('Normal\ncell lines','Tumor\ncell lines')),
              map_signif_level = T, step_increase = T)+
  ggtitle('H3K9me3')

pdf('plot/Fig3A/Fig3A.pdf', width = 4, height = 4)
Fig3A
dev.off()

pdf('plot/FigS3A/FigS3A.pdf', width = 4, height = 4)
FigS3A
dev.off()

Data_Fig3A = all_plotData_K27 %>% select(tissue_type, TorN, MD_anno, enrichmentScore, HisMark)
Data_FigS3A = all_plotData_K9 %>% select(tissue_type, TorN, MD_anno, enrichmentScore, HisMark)

write.table(Data_Fig3A, file = 'plot/Fig3A/Data_Fig3A.txt', sep = '\t', quote = F, row.names = F, col.names = T)
write.table(Data_FigS3A, file = 'plot/FigS3A/Data_FigS3A.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig3B. Tumor-loss, shared, and gain long LOCKs based on tumor and normal cell line pairs were identified.----
MD_size = data.frame(MD_anno = c('common PMDs','common HMDs','Others'),
                     MD_size = c(sum(commonPMD_hg38$end - commonPMD_hg38$start),
                                 sum(commonHMD_hg38$end - commonHMD_hg38$start),
                                 2.7e9 - sum(commonPMD_hg38$end - commonPMD_hg38$start) - sum(commonHMD_hg38$end - commonHMD_hg38$start)))

plotData_fraction = list()
for (t in 1:6) {
  tissue_type = all_tissue_types[t]
  TumorAndNormalLOCK_long = read.table(paste0('data_CellLine/09_TumorAndNormalLOCKs/Long_LOCKs/',tissue_type,'.TumorAndNormalLOCK_long.txt'),
                                       sep = '\t', header = T)
  df = TumorAndNormalLOCK_long
  
  FractionMD = group_by(df, group, MD_anno) %>% 
    summarise(d = sum(d), .groups = 'drop') %>% 
    mutate(sample = tissue_type, .before = 1) %>% 
    left_join(., MD_size, by = 'MD_anno') %>% 
    mutate(value = d/MD_size) %>% 
    select(-d, -MD_size) %>% 
    rbind(.,
          group_by(., sample, MD_anno ) %>% 
            summarise(value = sum(value), .groups = 'drop') %>% 
            mutate(value = (1-value)) %>% 
            mutate(group = 'Other', .before = 2))
  FractionMD = filter(FractionMD, group != 'Other')
  FractionMD$group = factor(FractionMD$group, levels = rev(c("Tumor-loss LOCKs", "Shared LOCKs", "Tumor-gain LOCKs" )))
  FractionMD$MD_anno = factor(FractionMD$MD_anno, levels = c('common PMDs', 'common HMDs', 'Others'))
  FractionMD = mutate(FractionMD, tissue_type = all_tissue_type_names[t])
  plotData_fraction[[t]] = FractionMD
}
plotData_fraction = do.call(rbind, plotData_fraction)
plotData_fraction$tissue_type = factor(plotData_fraction$tissue_type, levels = all_tissue_type_names)
Fig3B = ggplot(plotData_fraction, aes(x = MD_anno, y = value*100, fill = group))+
  facet_grid(.~tissue_type)+
  geom_bar(stat = 'identity', position = 'stack')+
  theme_classic()+
  ylab('Proportion of LOCKs in\neach part of the genome (%)')+
  xlab('')+
  scale_fill_manual(values = c('red','yellow','blue'))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        axis.text = element_text(color = 'black'),
        strip.background = element_blank(),
        strip.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'))

pdf('plot/Fig3B/Fig3B.pdf', width = 12, height = 7)
Fig3B
dev.off()

Data_Fig3B = plotData_fraction
write.table(Data_Fig3B, file = 'plot/Fig3B/Data_Fig3B.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig3C/FigS3B. Distribution of tumor-loss LOCKs, shared LOCKs, and tumor-gain LOCKs in each pair across different genomic regions.----
all_pieChart = list()
all_plotData = list()
for (t in 1:6) {
  tissue_type = all_tissue_types[t]
  TumorAndNormalLOCK_long = read.table(paste0('data_CellLine/09_TumorAndNormalLOCKs/Long_LOCKs/',tissue_type,'.TumorAndNormalLOCK_long.txt'),
                                       sep = '\t', header = T)
  df = TumorAndNormalLOCK_long
  blank_theme <- theme_minimal()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank(),
          axis.ticks = element_blank(),
          plot.title=element_text(size=14, face="bold"))
  pieChartData = group_by(df, group, MD_anno) %>% 
    summarise(d = sum(d), .groups = 'drop') %>% 
    mutate(MD_anno = factor(MD_anno, levels = c('common HMDs', 'common PMDs', 'Others')))
  pieChartGroup = names(table(pieChartData$group))
  
  pieCharts = list()
  for (g in 1:length(pieChartGroup)) {
    oneGroup = pieChartGroup[g]
    onePieChartData = filter(pieChartData, group %in% oneGroup)
    totalSize = sum(onePieChartData$d)
    onePieChartData$value = onePieChartData$d/totalSize
    onePieChartData$labelPos=1-(c(0, cumsum(onePieChartData$value))[1:nrow(onePieChartData)]+onePieChartData$value/2)
    pieCharts[[g]] = ggplot(onePieChartData, aes(x="", y=value, fill=MD_anno))+
      geom_bar(width = 1, stat = "identity")+
      coord_polar("y", start=0)+
      scale_fill_manual(values=c(`common HMDs` = "#56B4E9", 
                                 `common PMDs` = "#E69F00", 
                                 Others = "#999999"))+
      theme_bw()+
      blank_theme+
      theme(axis.text.x = element_blank(),
            plot.title = element_text(hjust = 0.5, face = 'plain'))+
      geom_label_repel(data = onePieChartData,
                       aes(y = labelPos, label = MD_anno),
                       size = 4.5, nudge_x = 1, show.legend = FALSE)+
      guides(fill = 'none')+
      ggtitle(oneGroup)
  }
  
  if (length(pieCharts) == 1) {
    p1_pieChart__inFun1 = pieCharts[[1]]
  }
  if (length(pieCharts) == 2) {
    p1_pieChart__inFun1 = cowplot::plot_grid(pieCharts[[1]],pieCharts[[2]], nrow = 1)
  }
  if (length(pieCharts) == 3) {
    p1_pieChart__inFun1 = cowplot::plot_grid(pieCharts[[1]],pieCharts[[2]],pieCharts[[3]], nrow = 1)
  }
  all_pieChart[[t]] = p1_pieChart__inFun1
  all_plotData[[t]] = pieChartData
}
Fig3C = all_pieChart[[1]]
FigS3B = all_pieChart[2:6]

pdf('plot/Fig3C/Fig3C.pdf', width = 10, height = 7)
Fig3C
dev.off()

Data_Fig3C = all_plotData[[1]]
write.table(Data_Fig3C, file = 'plot/Fig3C/Data_Fig3C.txt', sep = '\t', quote = F, row.names = F, col.names = T)

for (t in 2:6) {
  print(t)
  tissue_type = all_tissue_types[t]
  pdf(paste0('plot/FigS3B/FigS3B.',tissue_type,'.pdf'), width = 10, height = 7)
  all_pieChart[[t]]
  dev.off()
  
  write.table(all_plotData[[t]], file = paste0('plot/FigS3B/Data_FigS3B.',tissue_type,'.txt'),
              sep = '\t', quote = F, row.names = F, col.names = T)
}

# Fig3D/FigS3C. Average H3K27me3 intensity in tumor-loss LOCKs, shared LOCKs, and tumor-gain LOCKs and their surrounding regions for each pair.----
plotdata_intensity = list()
p_intensity = list()
for (t in 1:6) {
  tissue_type = all_tissue_types[t]
  fileName = list.files('data_CellLine/10_computeMatrix_longLOCKs/H3K27me3/', 
                        pattern = tissue_type, full.names = T)
  gzFile = list()
  for (i in 1:length(fileName)) {
    oneFileName = fileName[i]
    gzFile[[i]] = read.table(file = oneFileName, skip = 1, header = F, sep = '\t')
    gzFile[[i]][which(!grepl('(.*)_.*',gzFile[[i]]$V4)),]$V4 = paste0(gzFile[[i]][which(!grepl('(.*)_.*',gzFile[[i]]$V4)),]$V4,'_0')
    gzFile[[i]]$V4 = str_match(gzFile[[i]]$V4, '(.*)_.*')[,2]
    oneFileName = basename(oneFileName)
    cellName = str_match(oneFileName,'.*?\\..*?\\..*?\\.(.*?\\..*?)\\..*')[,2]
    gzFile[[i]] = mutate(gzFile[[i]], cellName = cellName, .before = 1)
  }
  gzFile = do.call(rbind, gzFile)
  
  gzFile = select(gzFile,-V1,-V2,-V3,-V5,-V6) %>% group_by(cellName,V4) %>% summarise_if(is_numeric,mean,na.rm = T) %>% as.data.frame()
  rownames(gzFile) = paste0(gzFile$cellName,'_',gzFile$V4)
  gzFile = select(gzFile, -cellName, -V4)  
  gzFile = t(gzFile) %>% as.data.frame()
  
  gzData = gzFile
  gzData = mutate(gzData, site = 1:nrow(gzData))
  gzData$site = gzData$site * 10
  gzData = pivot_longer(gzData, cols = 1:(ncol(gzData)-1),names_to = 'group', values_to = 'signal') %>% as.data.frame()
  gzData = separate(gzData, col = 'group', sep = '_', into = c('cellName','group'))
  # gzData$group = sub('N.spec','Tumor-loss LOCKs',gzData$group)
  # gzData$group = sub('T.spec','Tumor-gain LOCKs',gzData$group)
  # gzData$group = sub('SharedByTAndN','Shared LOCKs',gzData$group)
  gzData$group = factor(gzData$group, levels = c('Tumor-loss LOCKs','Shared LOCKs','Tumor-gain LOCKs'))
  gzData[grepl('N\\..*',gzData$cellName),]$cellName = paste0(str_match(gzData[grepl('N\\..*',gzData$cellName),]$cellName, 'N\\.(.*)')[,2],
                                                             '\n(normal cell line)')
  gzData[grepl('T\\..*',gzData$cellName),]$cellName = paste0(str_match(gzData[grepl('T\\..*',gzData$cellName),]$cellName, 'T\\.(.*)')[,2],
                                                             '\n(tumor cell line)')
  gzData$cellName = factor(gzData$cellName, levels = c(unique(gzData$cellName)[grepl('normal',unique(gzData$cellName))],
                                                       unique(gzData$cellName)[grepl('tumor',unique(gzData$cellName))]))
  p25_TN_intensity = ggplot(gzData, aes(x = site, y = signal, color = cellName))+
    geom_line(alpha = 0.7)+
    facet_grid(.~group)+
    theme_classic()+
    ylab('H3K27me3 intensity')+
    xlab('')+
    scale_x_continuous(breaks = c(1,50000,150000,200000),labels = c('-50kb','Start','End','+50kb'))+
    theme(axis.text = element_text(color = 'black'),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.title = element_blank(),
          strip.background = element_blank())+
    scale_color_manual(values = c('#fec4b8', '#e59572','#2470a1','#515bd4','#4cc2c6'))+
    ggtitle(all_tissue_type_names[t])+
    theme(plot.title = element_text(hjust = 0.5))
  if (tissue_type %in% c('NE2_KYSE450', 'NE2_KYSE510')) {
    p25_TN_intensity = p25_TN_intensity+scale_color_manual(values = c('#2694ab','#eb7072'))
  }
  p_intensity[[t]] = p25_TN_intensity
  plotdata_intensity[[t]] = gzData
}
Fig3D = p_intensity[[1]]
FigS3C = p_intensity[2:6]

pdf('plot/Fig3D/Fig3D.pdf', width = 10, height = 6)
Fig3D
dev.off()

Data_Fig3D = plotdata_intensity[[1]]
write.table(Data_Fig3D, file = 'plot/Fig3D/Data_Fig3D.txt', 
            sep = '\t', quote = F, row.names = F, col.names = T)

for (t in 2:6) {
  tissue_type = all_tissue_types[t]
  
  pdf(paste0('plot/FigS3C/FigS3C.',tissue_type,'.pdf'), width = 10, height = 6)
  p_intensity[[t]]
  dev.off()
  
  write.table(plotdata_intensity[[t]], file = paste0('plot/FigS3C/Data_FigS3C.',tissue_type,'.txt'),
              sep = '\t', quote = F, row.names = F, col.names = T)
}


# Fig3E/FigS3D. Heatmap shows the average DNA methylation levels of the three groups of long LOCK regions, across the whole genome and within common PMDs and HMDs.----
p_methylation = list()
plotData_methylation = list()
for (t in 1:6) {
  print(t)
  tissue_type = all_tissue_types[t]
  plotData_met = read.csv(paste0('data_CellLine/11_DNA_methylation_PMD_HMD/',tissue_type,'.met.csv'))
  met_biological_replicates = ifelse(t %in% 5:6, T, F)
  # for biological_replicates
  if (met_biological_replicates) {
    plotData_met$gzName = sub('_rep.*?\\.','.',plotData_met$gzName)
    plotData_met = group_by(plotData_met, gzName, MD_anno, group) %>% 
      mutate(metSignal = mean(metSignal, na.rm = T),
             ZMADNormMetSignal = mean(ZMADNormMetSignal, na.rm = T),
             ZScoreNormMetSignal = mean(ZScoreNormMetSignal, na.rm = T)) %>% 
      select(MD_anno, group, metSignal, ZMADNormMetSignal, ZScoreNormMetSignal, gzName, TorN, dataType)
  }
  # plotdata
  plotData_met$sample = str_match(plotData_met$gzName, '.*\\.(.*?)\\..*?\\.rmCGI.bdg')[,2]
  plotData_met$sample = sub('MCF10A_96850','MCF10A',plotData_met$sample)#for breast
  plotData_met$TorN = sub('T','tumor',plotData_met$TorN)
  plotData_met$TorN = sub('N','normal',plotData_met$TorN)
  plotData_met$sample_TorN = paste0(plotData_met$sample,'\n(',plotData_met$TorN,')')
  xlabLevels = plotData_met[,c('sample','TorN', 'sample_TorN')] %>% unique() %>% arrange(TorN)
  plotData_met$sample = factor(plotData_met$sample, levels = xlabLevels$sample)
  plotData_met$sample_TorN = factor(plotData_met$sample_TorN, levels = xlabLevels$sample_TorN)
  plotData_met = filter(plotData_met, group != 'Overall') %>% 
    filter(MD_anno != 'Other')
  plotData_met$MD_anno = sub('commonPMD','common PMD',plotData_met$MD_anno)
  plotData_met$MD_anno = sub('commonHMD','common HMD',plotData_met$MD_anno)
  plotData_met$MD_anno = factor(plotData_met$MD_anno,levels = c('Genome','common PMD','common HMD'))
  plotData_met$group = sub('T.spec','Tumor-gain LOCK',plotData_met$group)
  plotData_met$group = sub('N.spec','Tumor-loss LOCK',plotData_met$group)
  plotData_met$group = sub('SharedByTAndN','Shared LOCK',plotData_met$group)
  plotData_met$group  = factor(plotData_met$group, levels = rev(c('Tumor-loss LOCK','Shared LOCK','Tumor-gain LOCK')))
  addNA = as.data.frame(table(plotData_met[,c('MD_anno','group','sample_TorN')])) 
  addNA = addNA[which(addNA$Freq == 0),1:3] %>% mutate(metSignal = NA, ZScoreNormMetSignal = NA)
  plotData_met = select(plotData_met, MD_anno, group, sample_TorN, metSignal, ZScoreNormMetSignal) %>% 
    rbind(., addNA)
  
  ##heatmap
  heatmapMetSignal =
    ggplot(plotData_met, aes(sample_TorN, group, fill = metSignal))+
    geom_tile()+
    facet_grid(.~MD_anno)+
    xlab('')+
    ylab('')+
    theme_classic()+
    geom_text(aes(label=round(metSignal, 1)))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color = 'black'),
          strip.background = element_blank(),
          legend.position = 'bottom', legend.direction = 'horizontal')+
    scale_fill_gradient(
      low = "yellow",high = "red")+
    guides(fill = guide_legend(title = 'β value(%)'))+
    ggtitle(all_tissue_type_names[t])+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  p_methylation[[t]] = heatmapMetSignal
  plotData_methylation[[t]] = plotData_met %>% select(!ZScoreNormMetSignal)
}
Fig3E = p_methylation[[1]]
FigS3D = p_methylation[2:6]

pdf('plot/Fig3E/Fig3E.pdf', width = 5, height = 4)
Fig3E
dev.off()

Data_Fig3E = plotData_methylation[[1]]
write.table(Data_Fig3E, file = 'plot/Fig3E/Data_Fig3E.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)

for (t in 2:6) {
  tissue_type = all_tissue_types[t]
  
  pdf(paste0('plot/FigS3D/FigS3D.',tissue_type,'.pdf'), width = 5, height = 4)
  p_methylation[[t]]
  dev.off()
  
  write.table(plotData_methylation[[t]], file = paste0('plot/FigS3D/Data_FigS3D.',tissue_type,'.txt'),
              sep = '\t', quote = F, row.names = F, col.names = T)
}

# Fig3F/FigS4A. Normalized gene expression levels within the three groups of LOCKs in common PMDs----
load('meta/RefSeq_gene/hg38_Refseq_gene.RData')
winsorize = T
p_GeneExp = list()
plotData_GeneExp = list()
for (t in 1:6) {
  print(t)
  # 1) get gene_TN
  tissue_type = all_tissue_types[t]
  FPKM_matrix = read.csv(list.files('data_CellLine/05_RNASeq/FPKM_multiSample/', 
                                    pattern = tissue_type, full.names = T), 
                         header = T, row.names = 1)
  FPKM_matrix = FPKM_matrix[rownames(FPKM_matrix) %in% hg38_Refseq_gene$gene_name,]
  FPKM_matrix = log2(FPKM_matrix+0.1)
  for (j in 1:ncol(FPKM_matrix)) {
    FPKM_matrix[,j] = as.numeric(scale(FPKM_matrix[,j]))
  }
  TumorAndNormalLOCK_long = read.table(paste0('data_CellLine/09_TumorAndNormalLOCKs/Long_LOCKs/',
                                              tissue_type,'.TumorAndNormalLOCK_long.txt'),
                                       sep = '\t', header = T)
  PMD_HMD = rbind(commonPMD_hg38 %>% mutate(MD_anno = 'common PMDs'), 
                  commonHMD_hg38 %>% mutate(MD_anno = 'common HMDs'))
  gene_TN = list()
  for (j in 1:ncol(FPKM_matrix)) {
    df = FPKM_matrix[,j,drop = F] %>% 
      as.data.frame() %>% 
      dplyr::rename(exp = 1) %>% 
      rownames_to_column('gene_name') %>% 
      merge(., hg38_Refseq_gene, by = 'gene_name') %>% 
      select(chrom, start, end, gene_name, exp) %>% 
      mutate(d = end - start)
    #add LOCK group anno
    df = bed_intersect(df, TumorAndNormalLOCK_long, suffix = c('', '.y')) %>% 
      group_by(chrom, start, end, gene_name, exp, d, group.y) %>% 
      summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
      filter(.overlap > d/2) %>% 
      dplyr::rename(LOCK_group = group.y) %>% 
      select(chrom, start, end, LOCK_group) %>% 
      left_join(df, ., by = c('chrom', 'start', 'end')) %>% 
      bed_sort()
    #add MD anno
    df = bed_intersect(df, PMD_HMD, suffix = c('', '.y')) %>% 
      dplyr::rename(MD_anno = MD_anno.y) %>% 
      group_by(chrom, start, end, gene_name, exp, d, LOCK_group, MD_anno) %>% 
      summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
      filter(.overlap > d/2) %>% 
      left_join(df, ., by = c('chrom', 'start', 'end', 'gene_name', 'exp', 'd', 'LOCK_group'))
    df[is.na(df$LOCK_group),]$LOCK_group = 'Others'
    df[is.na(df$MD_anno),]$MD_anno = 'Others'
    df = select(df, -.overlap)
    gene_TN[[j]] = df
  }
  names(gene_TN) = colnames(FPKM_matrix)
  save(gene_TN, file = paste0('data_CellLine/12_gene_TN/',tissue_type,'.gene_TN.RData'))
  
  # 2) get plotData
  plotData = list()
  for (i in 1:length(gene_TN)) {
    df = gene_TN[[i]] %>% filter(!is.na(exp))
    sample = names(gene_TN)[i]
    plotData_part1 = group_by(df, LOCK_group) %>% 
      summarise(value = median(exp, na.rm = T),
                gene_number = n(),
                .groups = 'drop') %>% 
      na.omit() %>%
      rbind(.,
            data.frame(LOCK_group = 'Overall', value = median(df$exp, na.rm = T), gene_number = nrow(df))) %>% 
      mutate(sample = sample, .before = 1) %>% 
      mutate(MD_anno = 'Genome', .before = value)
    
    plotData_part2 = group_by(df, LOCK_group, MD_anno) %>% 
      summarise(value = median(exp, na.rm = T),
                gene_number = n(),
                .groups = 'drop') %>% 
      na.omit() %>% 
      rbind(.,
            group_by(df, MD_anno) %>% 
              summarise(value = median(exp, na.rm = T),
                        gene_number = n(),
                        .groups = 'drop') %>% 
              mutate(LOCK_group = 'Overall')) %>% 
      mutate(sample = sample, .before = 1)
    
    plotData_all = rbind(plotData_part1, plotData_part2) %>% 
      mutate(to_separate = sample) %>% 
      separate(col = to_separate, into = c('tissue_type', 'TorN', 'name', 'data_type'), sep = '\\.') %>% 
      mutate(index = paste(tissue_type, TorN, sep = '.')) %>% 
      filter(LOCK_group != 'Others')
    plotData[[i]] = plotData_all
  }
  plotData = do.call(rbind, plotData)
  plotData$MD_anno = factor(plotData$MD_anno,
                            levels = c('Genome', 'common HMDs', 'common PMDs', 'Others'))
  
  # A dot represents the average of a gene across all samples
  tumor_sample = grep(paste0(tissue_type, '\\.T\\.'), names(gene_TN))
  normal_sample = grep(paste0(tissue_type, '\\.N\\.'), names(gene_TN))
  
  tumor_anno = gene_TN[[tumor_sample[1]]] %>% select(-exp)
  normal_anno = gene_TN[[normal_sample[1]]] %>% select(-exp)
  
  tumor_matrix = list()
  number = 1
  for (n_one_sample in tumor_sample) {
    tumor_matrix[[number]] = gene_TN[[n_one_sample]]['exp']
    number = number+1
  }
  tumor_matrix = do.call(cbind, tumor_matrix)
  tumor_matrix = mutate(tumor_anno, mean_exp = rowMeans(tumor_matrix, na.rm = T))[,c('gene_name', 'LOCK_group', 'MD_anno', 'mean_exp')] %>% 
    mutate(index = paste0(tissue_type,'.T'))
  
  normal_matrix = list()
  number = 1
  for (n_one_sample in normal_sample) {
    normal_matrix[[number]] = gene_TN[[n_one_sample]]['exp']
    number = number+1
  }
  normal_matrix = do.call(cbind, normal_matrix)
  normal_matrix = mutate(normal_anno, mean_exp = rowMeans(normal_matrix, na.rm = T))[,c('gene_name', 'LOCK_group', 'MD_anno', 'mean_exp')] %>% 
    mutate(index = paste0(tissue_type,'.N'))
  
  plotData_gene = rbind(tumor_matrix, normal_matrix)
  plotData_genome = plotData_gene
  plotData_genome$MD_anno = 'Genome'
  plotData_gene = rbind(plotData_genome, plotData_gene)
  plotData_gene = rbind(plotData_gene, mutate(plotData_gene, LOCK_group = 'Overall')) %>% 
    filter(LOCK_group != 'Others')
  plotData_gene$MD_anno = factor(plotData_gene$MD_anno, 
                                 levels = c('Genome', 'common HMDs', 'common PMDs', 'Others'))
  plotData_gene = na.omit(plotData_gene)
  plotData_gene = unique(plotData_gene)
  # 3) plot
  plotData_gene = filter(plotData_gene, LOCK_group != 'Overall')
  plotData_gene = filter(plotData_gene, !is.na(LOCK_group))
  plotData_gene$LOCK_group = factor(plotData_gene$LOCK_group, 
                                    levels = c('Tumor-loss LOCKs','Shared LOCKs','Tumor-gain LOCKs'))
  df_PMD = plotData_gene %>% filter(MD_anno == 'common PMDs')
  df_PMD$TorN = str_match(df_PMD$index, '.*\\.(.*)')[,2]
  df_PMD$TorN = sub('T','Tumor',sub('N','Normal', df_PMD$TorN))
  
  if (winsorize) {
    df_range = group_by(df_PMD, LOCK_group, TorN) %>% 
      summarise(mean = mean(mean_exp), sd = sd(mean_exp), .groups = 'drop') %>% 
      mutate(max = mean+3*sd, min = mean-3*sd) %>% 
      select(LOCK_group, TorN, max, min)
    df_PMD1 = left_join(df_PMD, df_range, by = c('LOCK_group', 'TorN'))
    df_PMD1$mean_exp_winsorize = df_PMD1$mean_exp
    df_PMD2 = df_PMD1
    
    df_PMD1[df_PMD1$mean_exp > df_PMD1$max,]$mean_exp_winsorize = df_PMD1[df_PMD1$mean_exp > df_PMD1$max,]$max
    df_PMD1[df_PMD1$mean_exp < df_PMD1$min,]$mean_exp_winsorize = df_PMD1[df_PMD1$mean_exp < df_PMD1$min,]$min
    
    df_PMD2 = filter(df_PMD1, mean_exp <= max & mean_exp >= min)
    # df_PMD2 is to remove outliers, df_PMD1 is to handle outliers with Winsorize method
    # df_PMD2 for drawing, df_PMD1 for statistics
  }
  
  if (usePairedLine) {
    # Make color annotations for matching lines
    lineColor = select(df_PMD2, gene_name, mean_exp, LOCK_group, TorN) %>% 
      pivot_wider(names_from = 'TorN', values_from = 'mean_exp')
    lineColor = na.omit(lineColor)
    lineColor$compare = 'No change'
    if (sum(lineColor$Normal > lineColor$Tumor) > 0) {lineColor[lineColor$Normal > lineColor$Tumor,]$compare = 'Down in tumor'}
    if (sum(lineColor$Normal < lineColor$Tumor) > 0) {lineColor[lineColor$Normal < lineColor$Tumor,]$compare = 'Up in tumor'}
    lineColor = left_join(df_PMD2, lineColor[,c('gene_name', 'LOCK_group', 'compare')], by = c('gene_name', 'LOCK_group'))
    lineColor$compare = factor(lineColor$compare, levels = c('Up in tumor','Down in tumor','No change'))
    lineColor$singleValue = paste0(lineColor$gene_name,'|', lineColor$LOCK_group)
    
    # Remove the matching points and leave only one value of Tumor or Normal (the other may be removed as an outlier).
    singleValueToRm = as.data.frame(table(lineColor$singleValue)) %>% filter(Freq == 1) %>% pull(Var1)
    lineColor = filter(lineColor, !singleValue %in% singleValueToRm)
    df_PMD2$singleValue = paste0(df_PMD2$gene_name,'|',df_PMD2$LOCK_group)
    df_PMD2 = filter(df_PMD2, !singleValue %in% singleValueToRm)
  }
  
  if (grepl('NE2',tissue_type)) {
    NormalCell = str_match(tissue_type,'(.*)_.*')[,2]
    TumorCell = str_match(tissue_type,'.*_(.*)')[,2]
    df_PMD1$TorN = sub('Normal',NormalCell,sub('Tumor',TumorCell,df_PMD1$TorN))
    df_PMD1$TorN = factor(df_PMD1$TorN, levels = c(NormalCell,TumorCell))
    df_PMD2$TorN = sub('Normal',NormalCell,sub('Tumor',TumorCell,df_PMD2$TorN))
    df_PMD2$TorN = factor(df_PMD2$TorN, levels = c(NormalCell,TumorCell))
    
    lineColor$TorN = sub('Normal',NormalCell,sub('Tumor',TumorCell,lineColor$TorN))
    lineColor$TorN = factor(lineColor$TorN, levels = c(NormalCell,TumorCell))
    p_GeneExp[[t]] = ggplot(df_PMD2, aes(x = TorN, y = mean_exp_winsorize))+
      # geom_violin(aes(colour = TorN), alpha = 0.7)+
      geom_line(data = lineColor, aes(x = TorN, y = mean_exp, group=gene_name, color = compare), alpha = 0.7, size = 0.3)+
      geom_point(data = df_PMD2, aes(x = TorN, y = mean_exp, group=gene_name), size = 0.5)+
      facet_grid(.~LOCK_group)+
      theme_classic()+
      ggtitle(paste0('Genes in common PMD\n', sub('\n',' ',all_tissue_type_names[t])))+
      xlab('')+
      ylab('Z-MAD normalized\nlog10 (FPKM+0.1)')+
      theme(plot.title = element_text(hjust = 0.5),
            legend.title = element_blank(),
            axis.text = element_text(color = 'black'),
            strip.background = element_blank())+
      geom_signif(data = df_PMD1, aes(x = TorN, y = mean_exp_winsorize),
                  comparisons = list(c(NormalCell,TumorCell)),
                  map_signif_level = T,
                  test = "t.test",test.args = c(var.equal = T,paired = T))
    
    
    if (tissue_type == 'NE2_KYSE450') {p_GeneExp[[t]] = p_GeneExp[[t]]+scale_color_manual(values = c('NE2' = '#2694ab','KYSE450' = '#eb7072', 
                                                                                                     'Up in tumor' = 'purple', 'Down in tumor' = 'orange',
                                                                                                     'No change' = 'grey'))}
    if (tissue_type == 'NE2_KYSE510') {p_GeneExp[[t]] = p_GeneExp[[t]]+scale_color_manual(values = c('NE2' = '#2694ab','KYSE510' = '#eb7072', 
                                                                                                     'Up in tumor' = 'purple', 'Down in tumor' = 'orange',
                                                                                                     'No change' = 'grey'))}
  }else{
    p_GeneExp[[t]] = ggplot(df_PMD2, aes(x = TorN, y = mean_exp_winsorize))+
      geom_line(data = lineColor, aes(x = TorN, y = mean_exp, group=gene_name, color = compare), alpha = 0.7, size = 0.3)+
      geom_point(data = df_PMD2, aes(x = TorN, y = mean_exp, group=gene_name), size = 0.5)+
      facet_grid(.~LOCK_group)+
      theme_classic()+
      ggtitle(paste0('Genes in common PMD\n', sub('\n',' ',all_tissue_type_names[t])))+
      xlab('')+
      ylab('Z-MAD normalized\nlog10 (FPKM+0.1)')+
      theme(plot.title = element_text(hjust = 0.5),
            legend.title = element_blank(),
            axis.text = element_text(color = 'black'),
            strip.background = element_blank())+
      geom_signif(data = df_PMD1, aes(x = TorN, y = mean_exp_winsorize),
                  comparisons = list(c('Normal','Tumor')),
                  map_signif_level = T,
                  test = "t.test",test.args = c(var.equal = T,paired = T))+
      scale_color_manual(values = c('Normal' = '#2694ab',
                                    'Tumor' = '#eb7072',
                                    'Up in tumor' = 'purple', 
                                    'Down in tumor' = 'orange',
                                    'No change' = 'grey'))
    
    
  }
  p_GeneExp[[t]] = p_GeneExp[[t]] + guides(color = 'none') + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  plotData_GeneExp[[t]] = select(df_PMD2, gene_name, LOCK_group, MD_anno, mean_exp, index, TorN, mean_exp_winsorize)
}
Fig3F = p_GeneExp[[1]]
FigS4A = p_GeneExp[2:6]

pdf('plot/Fig3F/Fig3F.pdf', width = 5, height = 5)
Fig3F
dev.off()

Data_Fig3F = plotData_GeneExp[[1]]
write.table(Data_Fig3F, file = 'plot/Fig3F/Data_Fig3F.txt', sep = '\t', quote = F, row.names = F, col.names = T)

for (t in 2:6) {
  tissue_type = all_tissue_types[t]
  
  pdf(paste0('plot/FigS4A/FigS4A.',tissue_type,'.pdf'), width = 5, height = 5)
  p_GeneExp[[t]]
  dev.off()
  
  write.table(plotData_GeneExp[[t]], file = paste0('plot/FigS4A/FigS4A.',tissue_type,'.txt'), sep = '\t', quote = F, row.names = F, col.names = T)
}

# Fig3G/FigS4B. Consistent with Fig.3F, but for oncogenes.----
winsorize = T
N_top = 1500
p_OG_exp = list()
plotData_OG_exp = list()
for (t in 1:6) {
  print(t)
  tissue_type = all_tissue_types[t]
  load(paste0('data_CellLine/12_gene_TN/',tissue_type, '.gene_TN.RData'))
  gene_TN_inFun = gene_TN
  get_OG_TSG_exp = function(df, geneList, OG, TSG){
    to_return = c(median(filter(df, gene_name %in% intersect(geneList, OG))$exp, na.rm = T),
                  median(filter(df, gene_name %in% intersect(geneList, TSG))$exp, na.rm = T))
    names(to_return) = c('OG','TSG')
    return(to_return)
  } # if this function is used, the value returned by 'get_OG_TSG_plotData' will be the 'expression'.
  get_OG_TSG_plotData = function(OG_TSG_function, gene_TN_inFun, N_top){
    output = list()
    for (i in 1:length(gene_TN_inFun)) {
      df = gene_TN_inFun[[i]]
      index = names(gene_TN_inFun)[i]
      genestoBeFiltered = df$gene_name
      TSG = read.csv(file = '/data1/liangyuan/linux_deal/03_LOCK_PMD/01_data/15_TSG_OG/TSG.csv',header = T) %>% filter(Gene %in% genestoBeFiltered) %>% .[1:N_top,] %>% pull(Gene)
      OG = read.csv(file = '/data1/liangyuan/linux_deal/03_LOCK_PMD/01_data/15_TSG_OG/OG.csv',header = T) %>% filter(Gene %in% genestoBeFiltered) %>% .[1:N_top,] %>% pull(Gene)
      #overall
      overall_plotData = cbind(data.frame(Genome = OG_TSG_function(df, df$gene_name, OG, TSG)),
                               data.frame(commonPMD = OG_TSG_function(df, df[df$MD_anno == 'commonPMD',]$gene_name, OG, TSG)),
                               data.frame(commonHMD = OG_TSG_function(df, df[df$MD_anno == 'commonHMD',]$gene_name, OG, TSG)),
                               data.frame(Other = OG_TSG_function(df, df[df$MD_anno == 'Other',]$gene_name, OG, TSG))) %>% 
        t() %>% as.data.frame() %>% rownames_to_column('MD_anno') %>% 
        mutate(index = index, .before = 1) %>% 
        mutate(LOCK_group = 'Overall')
      
      #LOCK group, genome
      LOCK_group_plotData_genome = list()
      for (x in 1:4) {
        LOCK_group = c('N.spec', 'SharedByTAndN', 'T.spec', 'Other')[x]
        res = OG_TSG_function(df, df[df$LOCK_group == LOCK_group,]$gene_name, OG, TSG)
        LOCK_group_plotData_genome[[x]] = data.frame(index = index, MD_anno = 'Genome', OG = res[1], TSG = res[2], LOCK_group = LOCK_group)
        rownames(LOCK_group_plotData_genome[[x]]) = NULL
      }
      
      LOCK_group_plotData_genome = do.call(rbind, LOCK_group_plotData_genome)
      
      #LOCK group, MD_anno
      LOCK_group_plotData_MDanno = list()
      for (x in 1:3) {
        MD_group = c('commonPMD', 'commonHMD', 'Other')[x]
        LOCK_group_plotData_MDanno[[x]] = list()
        for (y in 1:4) {
          LOCK_group = c('N.spec', 'SharedByTAndN', 'T.spec', 'Other')[y]
          res = OG_TSG_function(df, df[df$LOCK_group == LOCK_group & df$MD_anno == MD_group,]$gene_name, OG, TSG)
          LOCK_group_plotData_MDanno[[x]][[y]] = data.frame(index = index, MD_anno = MD_group, OG = res[1], TSG = res[2], LOCK_group = LOCK_group)
          rownames(LOCK_group_plotData_MDanno[[x]][[y]]) = NULL
        }
        LOCK_group_plotData_MDanno[[x]] = do.call(rbind, LOCK_group_plotData_MDanno[[x]])
      }
      LOCK_group_plotData_MDanno = do.call(rbind, LOCK_group_plotData_MDanno)
      
      output[[i]] = rbind(overall_plotData, LOCK_group_plotData_genome, LOCK_group_plotData_MDanno) %>% 
        mutate(MD_anno = factor(MD_anno, levels = c('Genome', 'commonHMD', 'commonPMD', 'Other'))) 
    }
    output = do.call(rbind, output) %>% 
      unique() %>% 
      pivot_longer(cols = c('OG', 'TSG'), names_to = 'OG_TSG_anno', values_to = 'value') %>% 
      separate(col = 'index', into = c('tissue_type', 'TorN', 'name', 'data_type'), sep = '\\.') %>% 
      mutate(index = paste(tissue_type, TorN, data_type, sep = '.'))
    return(output)
  }
  plotData_TsgOgExp = get_OG_TSG_plotData(OG_TSG_function = get_OG_TSG_exp, gene_TN_inFun = gene_TN_inFun, N_top = N_top) # expression
  
  # plot of 'expression'
  plotData_TsgOgExp$LOCK_group = factor(plotData_TsgOgExp$LOCK_group, levels = c('Overall', 'N.spec', 'SharedByTAndN', 'T.spec', 'Other'))
  plotData_TsgOgExp$index1 = str_match(plotData_TsgOgExp$index, '(.*\\..*)\\..*')[,2]
  
  # MeanExpOfEachGene
  MeanExpOfEachGene = as.data.frame(matrix(nrow = nrow(gene_TN_inFun[[1]]), ncol = 0))
  for (i in 1:length(gene_TN_inFun)) {
    MeanExpOfEachGene[,i] = gene_TN_inFun[[i]]$exp
  }
  colnames(MeanExpOfEachGene) = names(gene_TN_inFun)
  #get mean in tumor or normal
  sampleClass = str_match(colnames(MeanExpOfEachGene), '.*\\.([NT])\\..*')[,2]
  MeanExpOfEachGene_tumor = MeanExpOfEachGene[,which(sampleClass == 'T')]
  if (class(MeanExpOfEachGene_tumor) != 'numeric') {
    MeanExpOfEachGene_tumor = rowMeans(MeanExpOfEachGene_tumor)
  }
  MeanExpOfEachGene_tumor = gene_TN_inFun[[1]] %>% 
    select(-exp) %>% 
    mutate(MeanExp = MeanExpOfEachGene_tumor, .before = 5) %>% 
    mutate(index = paste0(tissue_type, '.T'))
  MeanExpOfEachGene_normal = MeanExpOfEachGene[,which(sampleClass == 'N')]
  if (class(MeanExpOfEachGene_normal) != 'numeric') {
    MeanExpOfEachGene_normal = rowMeans(MeanExpOfEachGene_normal)
  }
  MeanExpOfEachGene_normal = gene_TN_inFun[[1]] %>% 
    select(-exp) %>% 
    mutate(MeanExp = MeanExpOfEachGene_normal, .before = 5) %>% 
    mutate(index = paste0(tissue_type, '.N'))
  MeanExpOfEachGene = rbind(MeanExpOfEachGene_tumor, MeanExpOfEachGene_normal)
  #OG/TSG annotation
  TSG = read.csv(file = '/data1/liangyuan/linux_deal/03_LOCK_PMD/01_data/15_TSG_OG/TSG.csv',header = T) %>% filter(Gene %in% gene_TN_inFun[[1]]$gene_name) %>% .[1:N_top,] %>% pull(Gene)
  OG = read.csv(file = '/data1/liangyuan/linux_deal/03_LOCK_PMD/01_data/15_TSG_OG/OG.csv',header = T) %>% filter(Gene %in% gene_TN_inFun[[1]]$gene_name) %>% .[1:N_top,] %>% pull(Gene)
  OG_TSG_list = rbind(data.frame(gene_name = OG, OG_TSG_anno = 'OG'),
                      data.frame(gene_name = TSG, OG_TSG_anno = 'TSG'))
  MeanExpOfEachGene = merge(MeanExpOfEachGene, OG_TSG_list, by = "gene_name")
  #LOCK group: add overall
  MeanExpOfEachGene1 = MeanExpOfEachGene
  MeanExpOfEachGene1$LOCK_group = 'Overall'
  MeanExpOfEachGene = rbind(MeanExpOfEachGene, MeanExpOfEachGene1)
  #MD_anno: add Genome
  MeanExpOfEachGene2 = MeanExpOfEachGene
  MeanExpOfEachGene2$MD_anno = 'Genome'
  MeanExpOfEachGene = rbind(MeanExpOfEachGene, MeanExpOfEachGene2)
  using_MeanExpOfEachGene = MeanExpOfEachGene %>% filter(MD_anno == 'common PMDs')
  plotData = filter(using_MeanExpOfEachGene,MD_anno == 'common PMDs' & OG_TSG_anno == 'OG')
  plotData = filter(plotData, LOCK_group %in% c('Tumor-loss LOCKs','Shared LOCKs','Tumor-gain LOCKs'))
  plotData$LOCK_group = factor(plotData$LOCK_group, levels = c('Tumor-loss LOCKs','Shared LOCKs','Tumor-gain LOCKs'))
  plotData$TorN = str_match(plotData$index, '.*\\.(.*)')[,2]
  plotData$TorN = sub('T','Tumor',sub('N','Normal',plotData$TorN))
  plotData$TorN = factor(plotData$TorN, levels = c('Normal','Tumor'))
  OG_plotData = plotData
  if (winsorize) {
    #winsorize
    OG_range = group_by(OG_plotData, LOCK_group, TorN) %>% 
      summarise(mean = mean(MeanExp), sd = sd(MeanExp), .groups = 'drop') %>% 
      mutate(max = mean+3*sd, min = mean-3*sd) %>% 
      select(LOCK_group, TorN, max, min)
    OG_plotData = left_join(OG_plotData, OG_range, by = c('LOCK_group', 'TorN'))
    OG_plotData1 = OG_plotData
    
    OG_plotData[OG_plotData$MeanExp > OG_plotData$max,]$MeanExp = OG_plotData[OG_plotData$MeanExp > OG_plotData$max,]$max
    OG_plotData[OG_plotData$MeanExp < OG_plotData$min,]$MeanExp = OG_plotData[OG_plotData$MeanExp < OG_plotData$min,]$min
    
    OG_plotData1 = filter(OG_plotData, MeanExp <= max & MeanExp >= min)
    #OG_plotData1 removes outliers, OG_plotData handles outliers with the winsorize method
  }
  plotData_OG_exp[[t]] = OG_plotData
  
  # get plot
  #Make color annotations for matching lines
  lineColor = select(OG_plotData1, gene_name, MeanExp, LOCK_group, TorN) %>% 
    pivot_wider(names_from = 'TorN', values_from = 'MeanExp')
  lineColor = na.omit(lineColor)
  lineColor$compare = 'No change'
  if (sum(lineColor$Normal > lineColor$Tumor) > 0) {lineColor[lineColor$Normal > lineColor$Tumor,]$compare = 'Down in tumor'}
  if (sum(lineColor$Normal < lineColor$Tumor) > 0) {lineColor[lineColor$Normal < lineColor$Tumor,]$compare = 'Up in tumor'}
  lineColor = left_join(OG_plotData1, lineColor[,c('gene_name', 'LOCK_group', 'compare')], by = c('gene_name', 'LOCK_group'))
  lineColor$compare = factor(lineColor$compare, levels = c('Up in tumor','Down in tumor','No change'))
  lineColor$singleValue = paste0(lineColor$gene_name,'|', lineColor$LOCK_group)
  #Remove the matching points and leave only one value of Tumor or Normal (the other may be removed as an outlier).
  singleValueToRm = as.data.frame(table(lineColor$singleValue)) %>% filter(Freq == 1) %>% pull(Var1)
  lineColor = filter(lineColor, !singleValue %in% singleValueToRm)
  OG_plotData1$singleValue = paste0(OG_plotData1$gene_name,'|',OG_plotData1$LOCK_group)
  OG_plotData1 = filter(OG_plotData1, !singleValue %in% singleValueToRm)
  
  if (grepl('NE2',tissue_type)) {
    NormalCell = str_match(tissue_type,'(.*)_.*')[,2]
    TumorCell = str_match(tissue_type,'.*_(.*)')[,2]
    OG_plotData$TorN = sub('Normal',NormalCell,sub('Tumor',TumorCell,OG_plotData$TorN))
    OG_plotData$TorN = factor(OG_plotData$TorN, levels = c(NormalCell,TumorCell))
    OG_plotData1$TorN = sub('Normal',NormalCell,sub('Tumor',TumorCell,OG_plotData1$TorN))
    OG_plotData1$TorN = factor(OG_plotData1$TorN, levels = c(NormalCell,TumorCell))
    lineColor$TorN = sub('Normal',NormalCell,sub('Tumor',TumorCell,lineColor$TorN))
    lineColor$TorN = factor(lineColor$TorN, levels = c(NormalCell,TumorCell))
    p27_OG_Tumorloss_PMD = ggplot(OG_plotData1, aes(x = TorN, y = MeanExp))+
      # geom_boxplot(aes(color = TorN), outlier.size = 0.1)+
      geom_line(data = lineColor, aes(x = TorN, y = MeanExp, group=gene_name, color = compare), alpha = 0.7, size = 0.3)+
      geom_point(data = OG_plotData1, aes(x = TorN, y = MeanExp, group=gene_name), size = 0.5)+
      theme_classic()+
      xlab('')+
      facet_grid(.~LOCK_group)+
      ylab('Z-MAD normalized\nlog10 (FPKM+0.1)')+
      theme(plot.title = element_text(hjust = 0.5),
            axis.text = element_text(color = 'black'),
            legend.title = element_blank(),
            strip.background = element_blank())+
      geom_signif(data = OG_plotData, aes(x = TorN, y = MeanExp),
                  comparisons = list(c(NormalCell,TumorCell)),#
                  map_signif_level = T, test = "t.test",test.args = c(var.equal = T,paired = T))+
      ggtitle(paste0('Oncogenes in common PMD\n', sub('\n',' ',all_tissue_type_names[t])))
    if (tissue_type == 'NE2_KYSE450') {p27_OG_Tumorloss_PMD = p27_OG_Tumorloss_PMD+scale_color_manual(values = c('NE2' = '#2694ab','KYSE450' = '#eb7072', 
                                                                                                                 'Up in tumor' = 'purple', 'Down in tumor' = 'orange',
                                                                                                                 'No change' = 'grey'))}
    if (tissue_type == 'NE2_KYSE510') {p27_OG_Tumorloss_PMD = p27_OG_Tumorloss_PMD+scale_color_manual(values = c('NE2' = '#2694ab','KYSE510' = '#eb7072', 
                                                                                                                 'Up in tumor' = 'purple', 'Down in tumor' = 'orange',
                                                                                                                 'No change' = 'grey'))}
  }else{
    p27_OG_Tumorloss_PMD = ggplot(OG_plotData1, aes(x = TorN, y = MeanExp))+
      # geom_boxplot(aes(color = TorN), outlier.size = 0.1)+
      geom_line(data = lineColor, aes(x = TorN, y = MeanExp, group=gene_name, color = compare), alpha = 0.7, size = 0.3)+
      geom_point(data = OG_plotData1, aes(x = TorN, y = MeanExp, group=gene_name), size = 0.5)+
      theme_classic()+
      xlab('')+
      facet_grid(.~LOCK_group)+
      ylab('Z-MAD normalized\nlog10 (FPKM+0.1)')+
      theme(plot.title = element_text(hjust = 0.5),
            axis.text = element_text(color = 'black'),
            legend.title = element_blank(),
            strip.background = element_blank())+
      geom_signif(data = OG_plotData, aes(x = TorN, y = MeanExp),
                  comparisons = list(c('Normal','Tumor')),#
                  map_signif_level = T, test = "t.test",test.args = c(var.equal = T,paired = T))+
      ggtitle(paste0('Oncogenes in common PMD\n', sub('\n',' ',all_tissue_type_names[t])))+
      scale_color_manual(values = c('Normal' = '#2694ab',
                                    'Tumor' = '#eb7072',
                                    'Up in tumor' = 'purple', 
                                    'Down in tumor' = 'orange',
                                    'No change' = 'grey'))
  }
  p_OG_exp[[t]] = p27_OG_Tumorloss_PMD + guides(color = 'none') + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  plotData_OG_exp[[t]] = OG_plotData1
}

Fig3G = p_OG_exp[[1]]
FigS4B = p_OG_exp[2:6]

pdf('plot/Fig3G/Fig3G.pdf', width = 5, height = 5)
Fig3G
dev.off()

Data_Fig3G = plotData_OG_exp[[1]]
write.table(Data_Fig3G, file = 'plot/Fig3G/Data_Fig3G.txt', sep = '\t', quote = F, row.names = F, col.names = T)

for (t in 2:6) {
  tissue_type = all_tissue_types[t]
  
  pdf(paste0('plot/FigS4B/FigS4B.',tissue_type,'.pdf'), width = 5, height = 5)
  p_OG_exp[[t]]
  dev.off()
  
  write.table(plotData_OG_exp[[t]], file = paste0('plot/FigS4B/FigS4B.',tissue_type,'.txt'), sep = '\t', quote = F, row.names = F, col.names = T)
}


# FigS4C. Normalized expression levels of OGs within long LOCKs versus non-LOCK regions, categorized by S/I/L-PMDs.----
SIL = read.table(file = 'meta/commonPMD_commonHMD/SIL_hg38.bed')
colnames(SIL) = c('chrom','start','end','d','MD_anno','PMD_length')
SIL = select(SIL, chrom, start, end, PMD_length)
SIL = filter(SIL, PMD_length != 'N')
# FPKM
j = 0
FPKM = list()
for (t in 1:6) {
  tissue_type = all_tissue_types[t]
  FPKM_file = read.csv(list.files('data_CellLine/05_RNASeq/FPKM_multiSample/', pattern = tissue_type, full.names = T)) %>% 
    column_to_rownames('X')
  for (i in 1:ncol(FPKM_file)) {
    j=j+1
    FPKM[[j]] = FPKM_file[,i,drop=F]
    names(FPKM)[j] = colnames(FPKM_file)[i]
    colnames(FPKM[[j]]) = 'FPKM'
  }
}
cellName = str_match(names(FPKM),'^.*?\\.[NT]\\.(.*?)\\..*')[,2]
cellName = data.frame(cellName = cellName, order = 1:length(cellName))
cellName = distinct(cellName, cellName, .keep_all = T)
FPKM = FPKM[cellName$order]

# ZMAD
ZMAD = list()
for (t in 1:6) {
  tissue_type = all_tissue_types[t]
  load(paste0('data_CellLine/12_gene_TN/',tissue_type,'.gene_TN.RData'))
  ZMAD[[t]] = gene_TN
}
ZMAD = do.call(c, ZMAD)
cellName_ZMAD = str_match(names(ZMAD),'^.*?\\.[NT]\\.(.*?)\\..*')[,2]
cellName_ZMAD = data.frame(cellName_ZMAD = cellName_ZMAD, order = 1:length(cellName_ZMAD))
cellName_ZMAD = distinct(cellName_ZMAD, cellName_ZMAD, .keep_all = T)
ZMAD = ZMAD[cellName_ZMAD$order]

# LOCK
LOCK = list()
j=0
for (t in 1:6) {
  tissue_type = all_tissue_types[t]
  files = list.files('data_CellLine/04_LOCK/withAnnotations/', 
                     pattern = paste0(tissue_type,'.*H3K27me3.LOCK_-0.7.bed'), full.names = T)
  for (i in 1:length(files)) {
    j=j+1
    LOCK[[j]] = read.table(files[i], header = T)
    names(LOCK)[j] = basename(files[i])
  }
}
names(LOCK) = sub('MCF7','MCF_7',names(LOCK))
cellName_LOCK = names(LOCK)
cellName_LOCK = str_match(cellName_LOCK, '.*?\\.[NT]\\.(.*?)\\..*')[,2]
cellName_LOCK = data.frame(cellName_LOCK = cellName_LOCK, order = 1:length(cellName_LOCK))
cellName_LOCK = distinct(cellName_LOCK, cellName_LOCK, .keep_all = T)
cellName_LOCK$cellName_LOCK = factor(cellName_LOCK$cellName_LOCK, levels = cellName$cellName)
cellName_LOCK = arrange(cellName_LOCK, cellName_LOCK)
LOCK = LOCK[cellName_LOCK$order]

#For each cell line, define groups for the genes（LOCK/non-LOCK, MD_anno, PMD_length_anno, OG/TSG）
plotData_list = list()
for (i in 1:length(LOCK)) {
  print(i)
  FPKM_df = FPKM[[i]] %>% rownames_to_column('gene_name')
  LOCK_df = LOCK[[i]][,c('chrom','start','end','LOCK_length')]
  ZMAD_df = ZMAD[[i]]
  #gene region, ZMAD
  plotData_df = ZMAD_df %>% 
    select(chrom, start, end, gene_name, exp, d, MD_anno)
  #add OG and TSG annotation
  genestoBeFiltered = plotData_df$gene_name
  TSG = read.csv(file = '/data1/liangyuan/linux_deal/03_LOCK_PMD/01_data/15_TSG_OG/TSG.csv',header = T) %>% filter(Gene %in% genestoBeFiltered) %>% .[1:N_top,] %>% pull(Gene)
  OG = read.csv(file = '/data1/liangyuan/linux_deal/03_LOCK_PMD/01_data/15_TSG_OG/OG.csv',header = T) %>% filter(Gene %in% genestoBeFiltered) %>% .[1:N_top,] %>% pull(Gene)
  OG_TSG_df = rbind(data.frame(gene_name = OG, OG_TSG_anno = 'OG'),
                    data.frame(gene_name = TSG, OG_TSG_anno = 'TSG'))
  OG_TSG_df = filter(OG_TSG_df, !gene_name %in% names(table(OG_TSG_df$gene_name)[table(OG_TSG_df$gene_name) == 2]))#去掉既是OG又是TSG的
  plotData_df = left_join(plotData_df, OG_TSG_df, by = 'gene_name')
  #add LOCK annotation
  inter1 = bed_intersect(plotData_df, LOCK_df, suffix = c('','.LOCK')) %>% 
    group_by(chrom, start, end, gene_name, exp, d, MD_anno, OG_TSG_anno, LOCK_length.LOCK) %>% 
    summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
    dplyr::rename(LOCK_length = LOCK_length.LOCK) %>% 
    filter(.overlap > d/2) %>% #
    select(!.overlap) 
  inter2 = filter(plotData_df, !gene_name %in% inter1$gene_name) %>% 
    mutate(LOCK_length = 'Outside_of_LOCK')
  plotData_df = rbind(inter1, inter2)
  #add PMD_length_anno
  inter_SIL1 = bed_intersect(plotData_df, SIL, suffix = c('','.SIL')) %>% 
    group_by(chrom, start, end, gene_name, exp, d, MD_anno, OG_TSG_anno, LOCK_length, PMD_length.SIL) %>% 
    summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
    dplyr::rename(PMD_length_anno = PMD_length.SIL) %>% 
    filter(.overlap > d/2) %>% 
    select(!.overlap)
  inter_SIL2 = filter(plotData_df, !gene_name %in% inter_SIL1$gene_name) %>% 
    mutate(PMD_length_anno = NA)
  plotData_df = rbind(inter_SIL1, inter_SIL2)
  #add FPKM
  plotData_df = left_join(plotData_df, FPKM_df, by = 'gene_name')
  plotData_df$sample = names(FPKM)[i]
  plotData_list[[i]] = plotData_df
}

# plotData
plotData = do.call(rbind, plotData_list)
plotData[,c('TorN','cellName')] = str_match(plotData$sample, '^.*?\\.([NT])\\.(.*?)\\..*')[,2:3]
plotData = filter(plotData, !is.na(OG_TSG_anno) & LOCK_length != 'Short_LOCK')
plotData$OG_TSG_anno = factor(plotData$OG_TSG_anno, levels = c('OG','TSG'))
plotData$LOCK_length = factor(plotData$LOCK_length, levels = c('Long_LOCK','Outside_of_LOCK'))
plotData_PMD = filter(plotData, MD_anno == 'common PMDs') %>% 
  group_by(LOCK_length, OG_TSG_anno, sample, TorN) %>% summarize(ZMAD = median(exp, na.rm = T), FPKM = median(FPKM, na.rm = T), .groups = 'drop')
plotData_HMD = filter(plotData, MD_anno == 'common HMDs') %>% 
  group_by(LOCK_length, OG_TSG_anno, sample, TorN) %>% summarize(ZMAD = median(exp, na.rm = T), FPKM = median(FPKM, na.rm = T), .groups = 'drop')
plotData_S = filter(plotData, PMD_length_anno == 'S') %>% 
  group_by(LOCK_length, OG_TSG_anno, sample, TorN) %>% summarize(ZMAD = median(exp, na.rm = T), FPKM = median(FPKM, na.rm = T), .groups = 'drop')
plotData_I = filter(plotData, PMD_length_anno == 'I') %>% 
  group_by(LOCK_length, OG_TSG_anno, sample, TorN) %>% summarize(ZMAD = median(exp, na.rm = T), FPKM = median(FPKM, na.rm = T), .groups = 'drop')
plotData_L = filter(plotData, PMD_length_anno == 'L') %>% 
  group_by(LOCK_length, OG_TSG_anno, sample, TorN) %>% summarize(ZMAD = median(exp, na.rm = T), FPKM = median(FPKM, na.rm = T), .groups = 'drop')

# plots
plotData_SIL_OG = list(plotData_S %>% mutate(PMD_length_anno = 'Within S'),
                       plotData_I %>% mutate(PMD_length_anno = 'Within I'),
                       plotData_L %>% mutate(PMD_length_anno = 'Within L'))
plotData_SIL_OG = do.call(rbind, plotData_SIL_OG)
plotData_SIL_OG = filter(plotData_SIL_OG, OG_TSG_anno == 'OG')
plotData_SIL_OG$PMD_length_anno = factor(plotData_SIL_OG$PMD_length_anno, levels = c('Within S','Within I','Within L'))
plotData_SIL_OG$LOCK_length = sub('Long_LOCK','Within long LOCKs',plotData_SIL_OG$LOCK_length)
plotData_SIL_OG$LOCK_length = sub('Outside_of_LOCK','Outside of long or short LOCKs',plotData_SIL_OG$LOCK_length)
plotData_SIL_OG$LOCK_length = factor(plotData_SIL_OG$LOCK_length, levels = c('Within long LOCKs','Outside of long or short LOCKs'))
PlotData_N = plotData_SIL_OG %>% filter(TorN == 'N')
PlotData_T = plotData_SIL_OG %>% filter(TorN == 'T')

plot_N_ZMAD = ggplot(PlotData_N, aes(x = LOCK_length, y = ZMAD))+ #
  facet_grid(.~PMD_length_anno)+
  geom_boxplot(aes(color = LOCK_length))+
  geom_point(size = 0.5)+
  xlab('')+
  ylab('ZMAD normalized\nlog10 (FPKM+0.1)')+
  theme_classic()+
  theme(strip.background = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank()
  )+
  scale_color_manual(values = c('Within long LOCKs' = '#ff0000', 
                                'Outside of long or short LOCKs' = '#f9d22a'))+
  geom_signif(comparisons = list(c('Within long LOCKs',
                                   'Outside of long or short LOCKs')),
              map_signif_level = T)+
  ggtitle('OG expression in normal cell lines')+
  guides(color = 'none')

plot_T_ZMAD = ggplot(PlotData_T, aes(x = LOCK_length, y = ZMAD))+ #
  facet_grid(.~PMD_length_anno)+
  geom_boxplot(aes(color = LOCK_length))+
  geom_point(size = 0.5)+
  xlab('')+
  ylab('ZMAD normalized\nlog10 (FPKM+0.1)')+
  theme_classic()+
  theme(strip.background = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank()
  )+
  scale_color_manual(values = c('Within long LOCKs' = '#ff0000', 'Outside of long or short LOCKs' = '#f9d22a'))+
  geom_signif(comparisons = list(c('Within long LOCKs','Outside of long or short LOCKs')),
              map_signif_level = T)+
  ggtitle('OG expression in tumor cell lines')+
  guides(color = 'none')

FigS4C = plot_grid(plot_N_ZMAD, plot_T_ZMAD, nrow = 1)

pdf('plot/FigS4C/FigS4C.pdf', width = 10, height = 6)
FigS4C
dev.off()

Data_FigS4C.Normal = PlotData_N %>% select(!ZMAD)
Data_FigS4C.Tumor = PlotData_T %>% select(!ZMAD)
write.table(Data_FigS4C.Normal, file = 'plot/FigS4C/Data_FigS4C.Normal.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)
write.table(Data_FigS4C.Tumor, file = 'plot/FigS4C/Data_FigS4C.Tumor.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)

# Combine the plots in Figure 3----
Figure_3 = ggarrange(
  ggarrange(Fig3A, Fig3B, nrow = 1, labels = c('A','B'), widths = c(0.7,2)),
  ggarrange(Fig3C, Fig3D, labels = c('C', 'D'), widths = c(1,1)),
  ggarrange(Fig3E, Fig3F, Fig3G, nrow = 1, labels = c('E','F','G'), widths = c(1.5,1,1)),
  ncol = 1, heights = c(1.5,1,1.5)
)
pdf('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/Figure_3.pdf', width = 13, height = 14)
Figure_3
dev.off()

# Combine the plots in Figure S3----
Figure_S3 = ggarrange(ggarrange(
  ggarrange(FigS3A, ggarrange(FigS3C[[1]], FigS3C[[2]], FigS3C[[3]], FigS3C[[4]], FigS3C[[5]], ncol = 1), 
            ncol = 1, labels = c('A','C'), font.label = list(size = 23), heights = c(2,5)),
  ggarrange(FigS3B[[1]], FigS3B[[2]], FigS3B[[3]], FigS3B[[4]], FigS3B[[5]], 
            ncol = 1, labels = 'B', font.label = list(size = 23)),
  nrow = 1, widths = c(1,1.2)),
  ggarrange(ggarrange(FigS3D[[2]], FigS3D[[3]], nrow = 1), 
            ggarrange(FigS3D[[1]], FigS3D[[4]], FigS3D[[5]], nrow = 1), 
            ncol = 1, labels = 'D', font.label = list(size = 23)),
  ncol = 1,
  heights = c(2,1)
)
pdf('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/Figure_S3.pdf', width = 13, height = 23)
Figure_S3
dev.off()

# Combine the plots in Figure S4----
Figure_S4 = ggarrange(
  ggarrange(ggarrange(FigS4A[[1]], FigS4A[[2]], FigS4A[[3]], FigS4A[[4]], FigS4A[[5]], nrow = 1),
            ggarrange(FigS4B[[1]], FigS4B[[2]], FigS4B[[3]], FigS4B[[4]], FigS4B[[5]], nrow = 1),
            ncol = 1, labels = c('A','B'), font.label = list(size = 23)),
  ggarrange(FigS4C, ggplot(), nrow = 1, widths = c(4,1)), 
  ncol = 1, labels = c('','C'), font.label = list(size = 23), heights = c(3,1))

pdf('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/Figure_S4.pdf', width = 13, height = 13)
Figure_S4
dev.off()


# Figure 4/ Figure S5----
# Fig4A/FigS5A. The enrichment of H3K27me3 or H3K9me3 long LOCKs in S/I/L-PMDs.----
overall_SIL = mutate(SIL, d = end - start) %>% 
  group_by(PMD_length) %>% 
  summarise(overall_length = sum(d), .groups = 'drop')
files = list.files('data_CellLine/04_LOCK//withAnnotations', full.names = T)
files_name = list.files('data_CellLine/04_LOCK//withAnnotations', full.names = F)
files_name = sub('.LOCK_-0.7.bed','',files_name)
plotData_enrichment = list()
for (i in 1:length(files)) {
  Long_LOCK = read.table(file = files[i], sep = '\t', header = T) %>% 
    filter(LOCK_length == 'Long_LOCK')
  Long_LOCK_size = sum(Long_LOCK$d)
  plotData_enrichment[[i]] = Long_LOCK %>% 
    bed_intersect(., SIL, suffix = c('.x','')) %>% 
    group_by(PMD_length) %>% 
    summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
    left_join(., overall_SIL, by = 'PMD_length') %>% 
    mutate(percent = .overlap/overall_length*100) %>% 
    select(PMD_length, percent) %>% 
    mutate(enrichment_score = percent/Long_LOCK_size) %>% 
    select(PMD_length, enrichment_score) %>% 
    mutate(sample = files_name[i])
}
plotData_enrichment = do.call(rbind, plotData_enrichment)
plotData_enrichment[,c('tissue_type','TorN','cell_name','hisMark')] = str_match(plotData_enrichment$sample, '(.*?)\\.(.*?)\\.(.*?)\\.(.*)')[,2:5]
plotData_enrichment$PMD_length = factor(plotData_enrichment$PMD_length, levels = c('S','I','L'))
plotData_enrichment = filter(plotData_enrichment, tissue_type %in% all_tissue_types)
tissue_type_df = data.frame(all_tissue_types = all_tissue_types, all_tissue_type_names = all_tissue_type_names)
plotData_enrichment = left_join(plotData_enrichment, tissue_type_df, by = c('tissue_type' = 'all_tissue_types'))
plotData_enrichment$all_tissue_type_names = factor(plotData_enrichment$all_tissue_type_names, levels = rev(all_tissue_type_names))

plotData_enrichment_cellLevel = plotData_enrichment
H3K27me3_enrich1 = filter(plotData_enrichment_cellLevel, hisMark == 'H3K27me3')
H3K9me3_enrich1 = filter(plotData_enrichment_cellLevel, hisMark == 'H3K9me3')

Fig4A = ggplot(H3K27me3_enrich1, aes(x = PMD_length, y = enrichment_score))+
  facet_grid(.~TorN)+
  geom_boxplot(outlier.size = 0.5)+
  geom_point(color = 'red', size = 0.5)+
  theme_classic()+
  xlab('')+
  ylab('Enrichment score')+
  theme(axis.text = element_text(color = 'black'),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('H3K27me3')+
  geom_signif(comparisons = list(c('S','I'), c('I','L'), c('S','L')),
              map_signif_level = T, step_increase = T)
FigS5A = ggplot(H3K9me3_enrich1, aes(x = PMD_length, y = enrichment_score))+
  facet_grid(.~TorN)+
  geom_boxplot(outlier.size = 0.5)+
  geom_point(color = 'red', size = 0.5)+
  theme_classic()+
  xlab('')+
  ylab('Enrichment score')+
  theme(axis.text = element_text(color = 'black'),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('H3K9me3')+
  geom_signif(comparisons = list(c('S','I'), c('I','L'), c('S','L')),
              map_signif_level = T, step_increase = T)

pdf('plot/Fig4A/Fig4A.pdf', width = 5, height = 7)
Fig4A
dev.off()

pdf('plot/FigS5A/FigS5A.pdf', width = 5, height = 7)
FigS5A
dev.off()

Data_Fig4A = H3K27me3_enrich1
Data_FigS5A = H3K9me3_enrich1
write.table(Data_Fig4A, file = 'plot/Fig4A/Data_Fig4A.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)
write.table(Data_FigS5A, file = 'plot/FigS5A/Data_FigS5A.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)

# Fig4B. The proportions of the three groups of long LOCKs covering S/I/L-PMDs.----
plotData_SIL_ratio = list()
for (t in 1:6) {
  tissue_type = all_tissue_types[t]
  TumorAndNormalLOCK = read.table(paste0('data_CellLine/09_TumorAndNormalLOCKs/Long_LOCKs/',tissue_type,'.TumorAndNormalLOCK_long.txt'), sep = '\t', header = T)
  inter_SIL = bed_intersect(TumorAndNormalLOCK, SIL, suffix = c('','.SIL')) %>% 
    group_by(group, PMD_length.SIL) %>% 
    summarise(.overlap = sum(.overlap), .groups = 'drop')
  overall_SIL = mutate(SIL, d = end - start) %>% 
    group_by(PMD_length) %>% 
    summarise(overall_length = sum(d), .groups = 'drop')
  inter_SIL1 = left_join(inter_SIL, overall_SIL, by = c('PMD_length.SIL' = 'PMD_length')) %>% 
    mutate(percent = .overlap/overall_length*100,
           tissue_type = all_tissue_type_names[t]) %>% 
    select(group, PMD_length.SIL, percent, tissue_type)
  plotData_SIL_ratio[[t]] = inter_SIL1
}
plotData_SIL_ratio = do.call(rbind, plotData_SIL_ratio)
plotData_SIL_ratio$PMD_length.SIL = factor(plotData_SIL_ratio$PMD_length.SIL, levels = c('S','I','L'))
plotData_SIL_ratio$group = factor(plotData_SIL_ratio$group, 
                                  levels = rev(c('Tumor-loss LOCKs','Shared LOCKs','Tumor-gain LOCKs')))

Fig4B = ggplot(plotData_SIL_ratio, aes(x = PMD_length.SIL, y = percent, fill = group))+
  facet_grid(.~tissue_type)+
  geom_bar(stat = 'identity', position = 'stack')+
  theme_classic()+
  xlab('')+
  ylab('% of LOCKs in each common PMD class')+
  theme(axis.text = element_text(color = 'black'),
        strip.background = element_blank(),
        legend.title = element_blank())+
  scale_fill_manual(values = c('Tumor-gain LOCKs' = 'red', 'Shared LOCKs' = 'yellow', 'Tumor-loss LOCKs' = 'blue'))

pdf('plot/Fig4B/Fig4B.pdf', width = 14, height = 6)
Fig4B
dev.off()

Data_Fig4B = plotData_SIL_ratio
write.table(Data_Fig4B, file = 'plot/Fig4B/Data_Fig4B.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig4C/FigS5B. The heatmap shows average DNA methylation levels of the three groups of long LOCKs within S/I/L-PMDs in normal cell lines and tumor cell lines.----
if (F) {
  methylation = list()
  for (t in 1:6) {
    print(t)
    tissue_type = all_tissue_types[t]
    R_file = paste0('data_CellLine/09_TumorAndNormalLOCKs/Long_LOCKs/',tissue_type,'.TumorAndNormalLOCK_long.txt')
    TumorAndNormalLOCK = read.table(R_file, sep = '\t', header = T)
    colnames(TumorAndNormalLOCK) = c('chrom','start','end','LOCK_group','ID','d','MD_anno')
    
    TumorAndNormalLOCK_SIL = filter(TumorAndNormalLOCK, MD_anno == 'commonPMD') %>% 
      bed_intersect(., SIL, suffix = c('','.SIL')) %>% 
      group_by(chrom, start, end, LOCK_group, ID, d, PMD_length.SIL) %>% 
      summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
      filter(.overlap > d/2) %>% 
      select(ID, PMD_length.SIL)
    
    TumorAndNormalLOCK = left_join(TumorAndNormalLOCK, TumorAndNormalLOCK_SIL, by = 'ID')
    
    dataType = ifelse(tissue_type %in% c('breastNormalVsTN_Basal','breastNormalVsTN_CL'), 'EPIC', 'WGBS')
    bdgFile = list.files('05_WGBS_rmCGI_bdg/', pattern = paste0(tissue_type,'.*',dataType), full.names = T)
    bdg_df = list()
    AvMet_MD_anno = list()
    AvMet_PMD_length.SIL = list()
    for (i in 1:length(bdgFile)) {
      print(i)
      bdg_df[[i]] = read.table(file = bdgFile[i])
      colnames(bdg_df[[i]]) = c('chrom','start','end','signal')
      if (!(max(bdg_df[[i]]$signal, na.rm = T) <= 1)) {
        bdg_df[[i]]$signal = bdg_df[[i]]$signal/100
      }
      bdg_df[[i]] = filter(bdg_df[[i]], chrom %in% paste0('chr',1:22))
      bdg_df[[i]] = bed_intersect(bdg_df[[i]], blacklist_hg38, invert = T)
      
      #Average of all CpG sites
      AvMet_MD_anno[[i]] = bed_intersect(bdg_df[[i]], TumorAndNormalLOCK, suffix = c('.met','')) %>%  #少数LOCK没有CpG位点，不管它们了反正也是NA
        group_by(LOCK_group, MD_anno) %>% # PMD_length.SIL
        summarise(signal = mean(signal.met, na.rm = T), .groups = 'drop') %>% 
        mutate(sample = basename(bdgFile)[i])
      AvMet_PMD_length.SIL[[i]] = bed_intersect(bdg_df[[i]], TumorAndNormalLOCK, suffix = c('.met','')) %>%  #少数LOCK没有CpG位点，不管它们了反正也是NA
        group_by(LOCK_group, PMD_length.SIL) %>% # PMD_length.SIL
        summarise(signal = mean(signal.met, na.rm = T), .groups = 'drop') %>% 
        mutate(sample = basename(bdgFile)[i])
    }
    bdg_df = do.call(rbind, bdg_df)
    AvMet_MD_anno = do.call(rbind, AvMet_MD_anno)
    AvMet_PMD_length.SIL = do.call(rbind, AvMet_PMD_length.SIL)
    
    methylation[[t]] = list(AvMet_MD_anno, AvMet_PMD_length.SIL)
    names(methylation[[t]]) = c('AvMet_MD_anno', 'AvMet_PMD_length.SIL')
  }
  names(methylation) = all_tissue_types
  save(methylation, file = 'data_CellLine/13_DNA_methylation_SIL/TumorAndNormalLOCK_methylation.RData')
}

load(file = 'data_CellLine/13_DNA_methylation_SIL/TumorAndNormalLOCK_methylation.RData')
p_SIL_met = list()
plotData_SIL_met = list()
for (t in 1:6) {
  tissue_type = all_tissue_types[t]
  
  R_file = paste0('data_CellLine/09_TumorAndNormalLOCKs/Long_LOCKs/',tissue_type,'.TumorAndNormalLOCK_long.txt')
  TumorAndNormalLOCK = read.table(R_file, sep = '\t', header = T)
  colnames(TumorAndNormalLOCK) = c('chrom','start','end','LOCK_group','ID','d','MD_anno')
  
  TumorAndNormalLOCK_SIL = filter(TumorAndNormalLOCK, MD_anno == 'common PMDs') %>% 
    bed_intersect(., SIL, suffix = c('','.SIL')) %>% 
    group_by(chrom, start, end, LOCK_group, ID, d, PMD_length.SIL) %>% 
    summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
    filter(.overlap > d/2) %>% 
    select(ID, PMD_length.SIL)
  
  TumorAndNormalLOCK = left_join(TumorAndNormalLOCK, TumorAndNormalLOCK_SIL, by = 'ID')
  
  #DNA methylation
  AvMet_PMD_length.SIL = methylation[[t]][[2]]
  
  getMetPlot = function(met_df = AvMet_MD_anno, t) {
    met_df$LOCK_group = sub('SharedByTAndN','Shared LOCK',sub('T.spec','Tumor-gain LOCK',sub('N.spec','Tumor-loss LOCK',met_df$LOCK_group)))
    met_df$LOCK_group = factor(met_df$LOCK_group, levels = rev(c('Tumor-loss LOCK','Shared LOCK','Tumor-gain LOCK')))
    if ('MD_anno' %in%  colnames(met_df)) {
      met_df$MD_anno = sub('Other','Others',sub('common','common ',met_df$MD_anno))
      met_df$MD_anno = factor(met_df$MD_anno, levels = c('common PMD', 'common HMD', 'Others'))
    }
    if ('PMD_length.SIL' %in%  colnames(met_df)) {
      met_df$PMD_length.SIL = factor(met_df$PMD_length.SIL, levels = c('S','I','L'))
      met_df[,c('TorN','cellName')] = str_match(met_df$sample, '^.*?\\.([NT])\\.(.*?)\\..*')[,2:3]
    }
    met_df[,c('TorN','cellName')] = str_match(met_df$sample, '^.*?\\.([NT])\\.(.*?)\\..*')[,2:3]
    
    if (t %in% 5:6) {
      met_df$cellName = sub('_rep.','',met_df$cellName)
      if ('MD_anno' %in%  colnames(met_df)) {
        met_df = group_by(met_df, LOCK_group, MD_anno, TorN, cellName) %>% 
          summarise(signal = mean(signal), .groups = 'drop')
      }
      if ('PMD_length.SIL' %in%  colnames(met_df)) {
        met_df = group_by(met_df, LOCK_group, PMD_length.SIL, TorN, cellName) %>% 
          summarise(signal = mean(signal), .groups = 'drop')
      }
    }
    met_df$TorN1 = 'normal'
    met_df[met_df$TorN == 'T',]$TorN1 = 'tumor'
    met_df$sample_TorN = paste0(met_df$cellName,'\n(',met_df$TorN1,')')
    met_df$TorN = factor(met_df$TorN, levels = c('N','T'))
    xlabLevels = arrange(unique(met_df[,c('TorN','sample_TorN')]), TorN)$sample_TorN
    met_df$sample_TorN = factor(met_df$sample_TorN, levels = xlabLevels)
    
    if ('MD_anno' %in%  colnames(met_df)) {
      p = ggplot(met_df, aes(sample_TorN, LOCK_group, fill = signal*100))+
        geom_tile()+
        facet_grid(.~MD_anno)+
        xlab('')+
        ylab('')+
        theme_classic()+
        geom_text(aes(label=round(signal*100, 1)))+
        theme(plot.title = element_text(hjust = 0.5),
              axis.text = element_text(color = 'black'),
              strip.background = element_blank())+
        scale_fill_gradient(
          low = "yellow",high = "red")+
        guides(fill = guide_legend(title = 'β value'))+
        ggtitle(all_tissue_type_names[t])+
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust = 1))
    }
    if ('PMD_length.SIL' %in%  colnames(met_df)) {
      met_df = filter(met_df, !is.na(PMD_length.SIL))
      p = ggplot(met_df, aes(sample_TorN, LOCK_group, fill = signal*100))+
        geom_tile()+
        facet_grid(.~PMD_length.SIL)+
        xlab('')+
        ylab('')+
        theme_classic()+
        geom_text(aes(label=round(signal*100, 1)))+
        theme(plot.title = element_text(hjust = 0.5),
              axis.text = element_text(color = 'black'),
              strip.background = element_blank())+
        scale_fill_gradient(
          low = "yellow",high = "red")+
        guides(fill = guide_legend(title = 'β value'))+
        ggtitle(all_tissue_type_names[t])+
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust = 1))
    }
    return(list(p,met_df))
  }
  
  res = getMetPlot(met_df = AvMet_PMD_length.SIL, t = t)
  p_SIL_met[[t]] = res[[1]]
  plotData_SIL_met[[t]] = res[[2]]
  
}
p_SIL_met.1 = plot_grid(p_SIL_met[[1]], p_SIL_met[[2]], p_SIL_met[[3]], p_SIL_met[[4]], p_SIL_met[[5]], p_SIL_met[[6]], nrow = 3)

Fig4C = p_SIL_met[[3]]
FigS5B = p_SIL_met[c(1,2,4:6)]

pdf('plot/Fig4C/Fig4C.pdf', height = 5)
Fig4C
dev.off()

Data_Fig4C = plotData_SIL_met[[1]]
write.table(Data_Fig4C, file = 'plot/Fig4C/Data_Fig4C.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)

for (t in 2:6) {
  tissue_type = all_tissue_types[t]
  pdf(paste0('plot/FigS5B/FigS5B.',tissue_type,'.pdf'), height = 5)
  p_SIL_met[[t]]
  dev.off()
  
  write.table(plotData_SIL_met[[t]],paste0('plot/FigS5B/Data_FigS5B.',tissue_type,'.txt'),
              sep = '\t', quote = F, row.names = F, col.names = T)
}

# Fig4D/FigS5C. The heatmap shows H3K27me3 and H3K9me3 signals of three groups of long LOCKs within S/I/L-PMDs in normal breast cell lines and Her2+ BRCA cell lines.----
p_ht = list()
plotData_ht = list()
heatmap_data = list()
percent_to_label = list()
num = 0
for (t in 3:6) {
  tissue_type = all_tissue_types[t]
  
  #input region
  R_file = paste0('data_CellLine/09_TumorAndNormalLOCKs/Long_LOCKs/',tissue_type,
                  '.TumorAndNormalLOCK_long.txt')
  TumorAndNormalLOCK = read.table(R_file,
                                  sep = '\t', header = T)
  colnames(TumorAndNormalLOCK) = c('chrom','start','end','LOCK_group','ID','d','MD_anno')
  
  TumorAndNormalLOCK_SIL = filter(TumorAndNormalLOCK, MD_anno == 'common PMDs') %>% 
    bed_intersect(., SIL, suffix = c('','.SIL')) %>% 
    group_by(chrom, start, end, LOCK_group, ID, d, PMD_length.SIL) %>% 
    summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
    filter(.overlap > d/2) %>% 
    select(ID, PMD_length.SIL)
  
  TumorAndNormalLOCK = left_join(TumorAndNormalLOCK, TumorAndNormalLOCK_SIL, by = 'ID')
  #only for SIL
  TumorAndNormalLOCK = filter(TumorAndNormalLOCK, !is.na(PMD_length.SIL))
  #Regions smaller than 100bp have no signal because the bin size in computeMatrix is set to 100bp.
  TumorAndNormalLOCK = filter(TumorAndNormalLOCK, d >= 100)
  
  #computeMatrix result
  gzFileName = c(list.files('data_CellLine/10_computeMatrix_longLOCKs/H3K27me3_H3K9me3/', full.names = T, pattern = paste0(tissue_type,'.*H3K27me3')),
                 list.files('data_CellLine/10_computeMatrix_longLOCKs/H3K27me3_H3K9me3/', full.names = T, pattern = paste0(tissue_type,'.*H3K9me3')))
  gzFile = list()
  regionMean = list()
  for (i in 1:length(gzFileName)) {
    gzFile[[i]] = read.table(file = gzFileName[i], sep = '\t', skip = 1) %>% 
      mutate(ID = paste0(V1,':',V2,'-',V3)) %>% 
      column_to_rownames('ID') %>% 
      .[,-1:-6] %>% 
      .[TumorAndNormalLOCK$ID,]
    regionMean[[i]] = as.data.frame(rowMeans(gzFile[[i]], na.rm = T))
    colnames(regionMean[[i]]) = sub('.gz','',basename(gzFileName[i]))
  }
  names(gzFile) = basename(gzFileName)
  regionMean = do.call(cbind, regionMean)
  
  regionMean_K27 = regionMean[,str_match(colnames(regionMean), '.*(H3K.*me3).*')[,2] == 'H3K27me3']
  regionMean_K9 = regionMean[,str_match(colnames(regionMean), '.*(H3K.*me3).*')[,2] == 'H3K9me3']
  getBarDataChange = function(hisMark, regionMean_hisMark){
    TorN = str_match(colnames(regionMean_hisMark), '^.*?\\.([NT])\\..*')[,2]
    regionMean_T = regionMean_hisMark[,which(TorN == 'T')]
    regionMean_N = regionMean_hisMark[,which(TorN == 'N')]
    change_of_signal = as.data.frame(matrix(nrow = nrow(TumorAndNormalLOCK), ncol = 0))
    for (i in 1:nrow(regionMean_hisMark)) {
      x = as.numeric(regionMean_T[i,]) 
      y = as.numeric(regionMean_N[i,])
      if (sum(is.na(c(x,y))) > 0) {
        change_of_signal[i,'t.test_p'] = NA
      }else{
        change_of_signal[i,'t.test_p'] = t.test(x,y)$p.value
      }
      change_of_signal[i,'log2FC'] = log2(mean(x, na.rm = T) / mean(y, na.rm = T))
      change_of_signal[i,'minus'] = mean(x, na.rm = T) - mean(y, na.rm = T)
    }
    change_of_signal$`Change in tumor` = 'No change'
    change_of_signal[change_of_signal$minus > 0.1,]$`Change in tumor` = 'Up'
    change_of_signal[change_of_signal$minus < -0.1,]$`Change in tumor` = 'Down'
    change_of_signal$`Change in tumor` = factor(change_of_signal$`Change in tumor`, levels = c('Down','No change','Up'))
    change_of_signal$t.test_signif = 'non-significant'
    if (sum(change_of_signal$t.test_p < 0.05 & !is.na(change_of_signal$t.test_p)) > 0) {
      change_of_signal[change_of_signal$t.test_p < 0.05 & !is.na(change_of_signal$t.test_p),]$t.test_signif = 'significant'
    }
    change_of_signal = select(change_of_signal, minus, `Change in tumor`, t.test_signif)
    colnames(change_of_signal) = paste0(colnames(change_of_signal),'_',hisMark)
    return(change_of_signal)
  }
  TumorAndNormalLOCK_withRowMean = cbind(TumorAndNormalLOCK, cbind(getBarDataChange('H3K27me3', regionMean_K27), getBarDataChange('H3K9me3', regionMean_K9)))
  
  TumorAndNormalLOCK_withRowMean$PMD_length.SIL = factor(TumorAndNormalLOCK_withRowMean$PMD_length.SIL, levels = c('S','I','L'))
  TumorAndNormalLOCK_withRowMean$LOCK_group = factor(TumorAndNormalLOCK_withRowMean$LOCK_group, 
                                                     levels = c('Tumor-loss LOCKs','Shared LOCKs', 'Tumor-gain LOCKs'))
  TumorAndNormalLOCK_withRowMean = arrange(TumorAndNormalLOCK_withRowMean, PMD_length.SIL, LOCK_group, minus_H3K9me3)
  heatmap_data[[t]] = TumorAndNormalLOCK_withRowMean
  for (i in 1:length(gzFile)) {
    gzFile[[i]] = gzFile[[i]][TumorAndNormalLOCK_withRowMean$ID,]
  }
  plot = list()
  for (i in 1:length(gzFile)) {
    my_col_title = paste(str_match(names(gzFile)[i], '^.*?\\..*?\\.(.*?)\\.(.*?)\\..*')[,3:2], collapse = '\n') 
    htData = gzFile[[i]] %>% as.matrix()
    colnames(htData) = NULL
    rownames(htData) = NULL
    plot[[i]] = Heatmap(htData, cluster_rows = F, cluster_columns = F,
                        col = colorRamp2(c(-0.2, 0.5), c("yellow", "red")),
                        heatmap_legend_param = list(title = 'H3K27me3\nH3K9me3'), show_heatmap_legend = F,
                        column_title = my_col_title
    )
    if (i == 1) {
      plot[[i]] = Heatmap(htData, cluster_rows = F, cluster_columns = F,
                          col = colorRamp2(c(-0.2, 0.5), c("yellow", "red")),
                          heatmap_legend_param = list(title = 'H3K27me3 or\nH3K9me3 intensity'),
                          column_title = my_col_title,
                          row_title = all_tissue_type_names[t],
                          left_annotation = rowAnnotation(`common PMD class` = TumorAndNormalLOCK_withRowMean$PMD_length.SIL, 
                                                          `LOCK group` = TumorAndNormalLOCK_withRowMean$LOCK_group,
                                                          `H3K27me3 change in tumor` = TumorAndNormalLOCK_withRowMean$`Change in tumor_H3K27me3`,
                                                          `H3K9me3 change in tumor` = TumorAndNormalLOCK_withRowMean$`Change in tumor_H3K9me3`,
                                                          col = list(`LOCK group` = c("Tumor-loss LOCKs" = "#702fa7",
                                                                                      "Shared LOCKs" = "#f3747d",
                                                                                      "Tumor-gain LOCKs" = "#f2bd6d"),
                                                                     `common PMD class` = c('S' = '#92d2c3', 'I' = '#dddf98', 'L' = '#8fbad9'),
                                                                     `H3K27me3 change in tumor` = c('Down' = '#2195ce', 'No change' = 'grey', 'Up' = '#f75339'),
                                                                     `H3K9me3 change in tumor` = c('Down' = '#2195ce', 'No change' = 'grey', 'Up' = '#f75339')
                                                          ))
      )
    }
  }
  num = num+1
  if (length(plot) == 8) {
    p_ht[[num]] = plot[[1]]+plot[[2]]+plot[[3]]+plot[[4]]+plot[[5]]+plot[[6]]+plot[[7]]+plot[[8]]
  }
  if (length(plot) == 10) {
    p_ht[[num]] = plot[[1]]+plot[[2]]+plot[[3]]+plot[[4]]+plot[[5]]+plot[[6]]+plot[[7]]+plot[[8]]+plot[[9]]+plot[[10]]
  }
  plotData_ht[[num]] = gzFile
  #Annotate ratio
  region_number1 = group_by(TumorAndNormalLOCK_withRowMean, PMD_length.SIL, LOCK_group) %>% summarise(region_number = n())
  region_number2 = group_by(TumorAndNormalLOCK_withRowMean, PMD_length.SIL, LOCK_group, `Change in tumor_H3K9me3`) %>% summarise(region_number = n())
  region_number = left_join(region_number2, region_number1, by = c('PMD_length.SIL', 'LOCK_group')) %>% 
    mutate(percent = region_number.x / region_number.y * 100)
  percent_to_label[[t]] = region_number %>% 
    filter(PMD_length.SIL %in% c('S','I','L') & LOCK_group == 'Tumor-gain LOCKs')
}
names(percent_to_label)[3:6] = all_tissue_type_names[3:6]
names(plotData_ht) = all_tissue_types[3:6]
Fig4D = p_ht[[1]]
FigS5C = p_ht[2:4]

png('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/Fig4D.png',  width = 2000, height = 3000, res = 300)
Fig4D
dev.off()

png('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/FigS5C_luminalA.png',  width = 2000, height = 3000, res = 300)
FigS5C[[1]]
dev.off()

png('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/FigS5C_TN_Basal.png',  width = 2000, height = 3000, res = 300)
FigS5C[[2]]
dev.off()

png('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/FigS5C_TN_CL.png',  width = 2000, height = 3000, res = 300)
FigS5C[[3]]
dev.off()

Data_Fig4D = plotData_ht[[1]]
save(Data_Fig4D, file = 'plot/Fig4D/Data_Fig4D.RData')

Data_FigS5C.breastNormalVsLuminalA = plotData_ht[[2]]
Data_FigS5C.breastNormalVsTN_Basal = plotData_ht[[3]]
Data_FigS5C.breastNormalVsTN_CL = plotData_ht[[4]]

save(Data_FigS5C.breastNormalVsLuminalA, file = 'plot/FigS5C/Data_FigS5C.breastNormalVsLuminalA.RData')
save(Data_FigS5C.breastNormalVsTN_Basal, file = 'plot/FigS5C/Data_FigS5C.breastNormalVsTN_Basal.RData')
save(Data_FigS5C.breastNormalVsTN_CL, file = 'plot/FigS5C/Data_FigS5C.breastNormalVsTN_CL.RData')

Data_Fig4D_percent = percent_to_label[[3]]
Data_FigS5C_percent.breastNormalVsLuminalA = percent_to_label[[4]]
Data_FigS5C_percent.breastNormalVsTN_Basal = percent_to_label[[5]]
Data_FigS5C_percent.breastNormalVsTN_CL = percent_to_label[[6]]

write.table(Data_Fig4D_percent, file = 'plot/Fig4D/Data_Fig4D_percent.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)
write.table(Data_FigS5C_percent.breastNormalVsLuminalA, file = 'plot/FigS5C/Data_FigS5C_percent.breastNormalVsLuminalA.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)
write.table(Data_FigS5C_percent.breastNormalVsTN_Basal, file = 'plot/FigS5C/Data_FigS5C_percent.breastNormalVsTN_Basal.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)
write.table(Data_FigS5C_percent.breastNormalVsTN_CL, file = 'plot/FigS5C/Data_FigS5C_percent.breastNormalVsTN_CL.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)

# FigS5D. The impact on gene expression of the replacement of H3K9me3 in normal cells with H3K27me3 in tumor cells in I-PMDs and L-PMDs.----
PMD_class = c('I','L')
plotData2 = list()
for (t in 3:6) {
  print(t)
  # 1) Obtain the regions in L with Tumor-gain LOCK where H3K9me3 is decreased.
  target_region = heatmap_data[[t]]
  target_region = filter(target_region, PMD_length.SIL %in% PMD_class & LOCK_group == 'Tumor-gain LOCKs' & `Change in tumor_H3K9me3` == 'Down')
  target_region = target_region[,c('chrom','start','end')]
  # 2) Obtain the genes in the regions.
  Refseq_gene = mutate(hg38_Refseq_gene, d = end - start)
  getGenes = bed_intersect(Refseq_gene, target_region, suffix = c('','.y')) %>% 
    group_by(chrom, start, end, gene_name, d) %>% 
    summarise(.overlap = sum(.overlap), .groups = 'drop') %>% 
    filter(.overlap > d/2)
  target_genes = getGenes$gene_name
  # 3) Get expression of these genes in normal and tumor cells, separated by each cell.
  tissue_type = all_tissue_types[[t]]
  #FPKM
  FPKM = read.csv(paste0('data_CellLine/05_RNASeq/FPKM_multiSample/',tissue_type,'_FPKM.csv'))
  colnames(FPKM)[1] = 'gene_name' 
  FPKM = pivot_longer(FPKM, cols = 2:ncol(FPKM), names_to = 'sample_ID', values_to = 'exp')
  exp_data = FPKM; ylab_exp = 'FPKM'
  exp_data = filter(exp_data, gene_name %in% target_genes)
  plotData = separate(exp_data, col = 'sample_ID', sep = '\\.', into = c('tissue_type', 'TorN', 'cellName', 'dataType'))
  plotData$TorN = sub('T','Tumor cell lines',sub('N','Normal cell lines', plotData$TorN))
  plotData$TorN = factor(plotData$TorN, levels = c('Normal cell lines','Tumor cell lines'))
  plotData = arrange(plotData, TorN, cellName)
  plotData$cellName = factor(plotData$cellName, levels = unique(plotData$cellName))
  # 4) get mean expression in tumor or normal for each gene
  mean_data2 = group_by(plotData, tissue_type, gene_name, TorN) %>% 
    summarise(exp = mean(exp), .groups = 'drop') %>% 
    mutate(tissue_type_name = all_tissue_type_names[t])
  plotData2[[t]] = mean_data2
}
plotData2 = do.call(rbind, plotData2)
plotData2$tissue_type_name = factor(plotData2$tissue_type_name, levels = all_tissue_type_names[3:6])
gene_number_df = unique(plotData2[,c('gene_name','tissue_type_name')])
gene_number = as.data.frame(table(gene_number_df$tissue_type_name))
colnames(gene_number) = c('tissue_type_name','gene_number')
plotData2 = left_join(plotData2, gene_number, by = 'tissue_type_name')
plotData2$tissue_type_name2 = paste0(plotData2$tissue_type_name,'\n(n = ',plotData2$gene_number,')')
FigS5D = ggplot(plotData2, aes(x = TorN, y = log10(exp+0.1)))+
  facet_grid(.~tissue_type_name2)+
  geom_boxplot(outlier.size = 0.5)+
  geom_signif(comparisons = list(c('Normal cell lines','Tumor cell lines')),
              map_signif_level = T, step_increase = T)+
  xlab('')+
  ylab(paste0('log10 (',ylab_exp,' + 0.1)'))+
  theme_classic()+
  theme(axis.text = element_text(color = 'black'),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        strip.background = element_blank())

pdf('plot/FigS5D/FigS5D.pdf', width = 9, height = 9)
FigS5D
dev.off()

Data_FigS5D = plotData2
write.table(Data_FigS5D, file = 'plot/FigS5D/Data_FigS5D.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)

# Combine the plots in Figure 4----
Fig4AB = ggarrange(Fig4A, Fig4B, nrow = 1, widths = c(1,2), labels = c('A','B'), font.label = list(size = 23))

pdf('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/Fig4AB_high.pdf', width = 13, height = 8)
Fig4AB
dev.off()

pdf('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/Fig4AB_low.pdf', width = 13, height = 4)
Fig4AB
dev.off()

pdf('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/Fig4C.pdf',  width = 4, height = 4)
Fig4C
dev.off()

pdf('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/Fig4D.pdf',  width = 6.2, height = 12)
Fig4D
dev.off()

png('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/Fig4D.png',  width = 2000, height = 3000, res = 300)
Fig4D
dev.off()

# Combine the plots in Figure S5----
pdf('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/FigS5_line1_low.pdf', width = 13, height = 6)
ggarrange(FigS5A, FigS5B[[1]],FigS5B[[2]], widths = c(0.8,1,1), labels = c('A','B'), nrow = 1)
dev.off()

pdf('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/FigS5_line1_high.pdf', width = 13, height = 10)
ggarrange(FigS5A, FigS5B[[1]],FigS5B[[2]], widths = c(0.8,1,1), nrow = 1, labels = c('A','B'), font.label = list(size = 23))
dev.off()

pdf('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/FigS5_line2.pdf', width = 13, height = 6)
ggarrange(FigS5B[[3]],FigS5B[[4]],FigS5B[[5]], widths = c(1,1,1), nrow = 1)
dev.off()

png('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/FigS5C_luminalA.png',  width = 2000, height = 3000, res = 300)
FigS5C[[1]]
dev.off()

png('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/FigS5C_TN_Basal.png',  width = 2000, height = 3000, res = 300)
FigS5C[[2]]
dev.off()

png('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/FigS5C_TN_CL.png',  width = 2000, height = 3000, res = 300)
FigS5C[[3]]
dev.off()

pdf('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/FigS5D.pdf', width = 5, height = 6.5)
FigS5D
dev.off()

save(percent_to_label, file = '/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/FigS5C_percent_to_label.RData')

# Figure 5/ Figure S6----
# Fig5A. Peak density within short or long LOCKs, indicating the percentage of the total peak length relative to the total LOCK length, in 109 Roadmap samples.----
load('data_Roadmap/RDatas/LOCK_peak_109.RData')
load('data_Roadmap/RDatas/LOCK_109.RData')
peak_density = list()
for (i in 1:109) {
  peak_density[[i]] = LOCK_109[[i]] %>% 
    mutate(density = d_peak/d_LOCK) %>% 
    group_by(LOCK_length) %>% 
    summarise(density = median(density), .groups = 'drop') %>% 
    mutate(EID = names(LOCK_109)[i])
}
peak_density = do.call(rbind, peak_density)
peak_density$LOCK_length = sub('_LOCK',' LOCKs',peak_density$LOCK_length)
peak_density$LOCK_length = factor(peak_density$LOCK_length, levels = c('Long LOCKs', 'Short LOCKs'))
Fig5A = ggplot(peak_density, aes(x = LOCK_length, y = density*100))+
  geom_boxplot(aes(color = LOCK_length))+
  theme_classic()+
  xlab('')+
  ylab('Peak density (%)')+
  geom_signif(comparisons = list(c('Long LOCK','Short LOCK')),
              map_signif_level = T)+
  theme(axis.text = element_text(color = 'black'),
        plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c('Long LOCKs' = 'red', 'Short LOCKs' = 'orange'))+
  guides(color = 'none')+
  ggtitle('Roadmap samples\n(n = 109)')+
  coord_flip()

pdf('plot/Fig5A/Fig5A.pdf', width = 7, height = 3)
Fig5A
dev.off()

Data_Fig5A = peak_density
write.table(Data_Fig5A, file = 'plot/Fig5A/Data_Fig5A.txt', 
            sep = '\t', quote = F, row.names = F, col.names = T)

# Fig5B. The heatmap shows the top 20 LOLA terms most significantly enriched in peaks within short LOCKs.----
load('data_Roadmap/10_metascape_and_LOLA/lolaResults_peak_level.RData')
plotData_list0 = do.call(rbind, lolaResults_peak_level)
# Remove terms with special treatment.
names(table(plotData_list0$treatment))
plotData_list0 =  plotData_list0%>% 
  filter(treatment %in% c("No treatment", "No treatment(Antibody RAR _A704)", "No treatment(Antibody RXR _sc774)", "None", "regular medium") | is.na(treatment)) %>% 
  filter(!collection %in% c('UCSC_repeat_class', 'DMV'))
# calculate qValueLog
plotData_list0$qValueLog = -log10(plotData_list0$qValue+1e-322)
plotData_list0$oddsRatio[plotData_list0$oddsRatio == Inf] = range(plotData_list0$oddsRatio, finite = T)[2]
plotData_list0 = filter(plotData_list0, group %in% c('Peaks_in_long_LOCKs', 'Peaks_in_short_LOCKs', 'Typical_peaks'))
plotData_list0$group = gsub('_',' ',plotData_list0$group)
plotData_list0$group = factor(plotData_list0$group, levels = c('Typical peaks', 'Peaks in long LOCKs', 'Peaks in short LOCKs'))
plotData_list0 = mutate(plotData_list0, supportRatio = support/(support+c)) %>% 
  select(pValueLog, qValue, qValueLog, oddsRatio, supportRatio, collection,
         description, cellType, tissue, antibody, treatment, sample, group) %>% 
  filter(qValue < 0.05) %>% 
  mutate(all_description = paste(collection, description, cellType, tissue, antibody, treatment, sep = '||'))
plotData_list0 = as.data.frame(plotData_list0)
# only for peaks in short LOCKs
plotData = filter(plotData_list0, group == 'Peaks in short LOCKs')
# Order of descriptions
theOrders = group_by(plotData, all_description) %>% 
  summarise(oddsRatio = sum(oddsRatio, na.rm = T),
            supportRatio = sum(supportRatio, na.rm = T),
            qValueLog = sum(qValueLog, na.rm = T),
            .groups = 'drop')
theOrders_oddsRatio = arrange(theOrders, desc(oddsRatio)) %>% pull(all_description) %>% .[1:20]
theOrders_supportRatio = arrange(theOrders, desc(supportRatio)) %>% pull(all_description) %>% .[1:20]
theOrders_qValueLog = arrange(theOrders, desc(qValueLog)) %>% pull(all_description) %>% .[1:20]
# top20
theOrders__ = theOrders_qValueLog
title_anno = 'qValueLog'
plotData1 = filter(plotData_list0, all_description %in% theOrders__)
plotData1$all_description = factor(plotData1$all_description, levels = rev(theOrders__))
plotData1 = as.data.frame(plotData1)
# descriptions
top_description = unique(select(plotData1, all_description, collection, description, cellType, tissue, antibody, treatment))
top_description$description2 = top_description$description
top_description[8,'description2'] = 'ChIP K562 KDM4A'
top_description[13,'description2'] = 'ChIP K562 RBBP5'
top_description[17,'description2'] = 'Day 5 of in vitro differentiation FOXA2'
top_description = select(top_description, all_description, description2)
plotData1 = left_join(plotData1, top_description, by = 'all_description')
# get heatmap
plotData1 = plotData1[,c('group', 'sample', 'description2', title_anno)] %>% unique() %>% 
  distinct(group, sample, description2, .keep_all = T)
plotData2 = pivot_wider(plotData1, names_from = 'sample', values_from = title_anno)
plotData2$description2 = factor(plotData2$description2, levels = top_description$description2)
plotData2 = arrange(plotData2, group, description2)
sample_order = group_by(plotData1, group, sample) %>% 
  summarise(sum_qValue = sum(.data[[title_anno]]), .groups = 'drop') %>% 
  arrange(sample, desc(sum_qValue)) %>% 
  distinct(sample, .keep_all = T) %>% 
  arrange(group, desc(sum_qValue)) %>% 
  pull(sample)
plotData2 = plotData2[,c('group','description2',sample_order)]
plotData3 = split(plotData2, plotData2$group)
for (i in 1:length(plotData3)) {
  plotData3[[i]] = column_to_rownames(plotData3[[i]], 'description2') %>% 
    select(-group)
}
plotData3 = do.call(cbind, plotData3)
top_group_anno = str_match(colnames(plotData3), '(.*)\\..*')[,2]

dat = as.matrix(plotData3)
dat[dat > 50] = 50

p = Heatmap(dat,cluster_rows = F,cluster_columns = F,
            show_column_names = F,
            col = colorRamp2(c(1.31,50),c('white','red')),
            top_annotation = columnAnnotation(group = top_group_anno,
                                              col = list(group = c("Peaks in long LOCKs" = "red",
                                                                   "Peaks in short LOCKs" = "orange",
                                                                   "Typical peaks" = "blue"))),
            heatmap_legend_param = list(title = '-log10(qValue)',
                                        col_fun = colorRamp2(c(1.32,50),c('white','red')),
                                        title = "test", at = c(1.32,50),
                                        labels = c("-log10(0.05)", "≥50")),
            column_title = "Sample",
            row_names_side = "left",
            column_title_side = 'bottom',
            row_names_max_width = max_text_width(
              rownames(dat), 
              gp = gpar(fontsize = 12)
            )
)
p_LOLA = draw(p, merge = T)
Fig5B = p_LOLA

pdf('plot/Fig5B/Fig5B.pdf', width = 12)
Fig5B
dev.off()

Data_Fig5B = plotData3
write.table(Data_Fig5B, file = 'plot/Fig5B/Data_Fig5B.txt',
            sep = '\t', quote = F, row.names = T, col.names = T)
# Fig5C. % of CGI promoters or poised promoters----
Roadmap_table2 = read.csv(file = 'meta/Roadmap_table/Roadmap_table2.csv')
CGI_promoters = read.table(file = 'meta/RefSeq_gene/hg19_CGI_promoter.bed')
colnames(CGI_promoters) = c('chrom', 'start', 'end')
ESC_poised_promoters = read.table('meta/LOLACore/hg19/ESC_poised_promoters/regions/ESC_poised_promoters.bed')
colnames(ESC_poised_promoters) = c('chrom', 'start', 'end')
ESC_poised_promoters = bed_intersect(CGI_promoters, ESC_poised_promoters, suffix = c('','.y')) %>% 
  select(chrom, start, end) %>% unique()
load('data_Roadmap/RDatas/peak_109.RData')
H3K4me3_peak_path = 'data_Roadmap/01_narrowPeak_sorted/H3K4me3/'
plotData_diff_withinShortLOCK = list()
plotData_toShortLOCK_diff = list()
plotData_toShortLOCK_ESC = list()
plotData_toK27CGI_diff = list()
plotData_diff_all = list()
for (i in 1:length(peak_109)) {
  print(i)
  EID = names(peak_109)[i]
  summary = Roadmap_table2[Roadmap_table2$EID == EID,'summary']
  #Identification of poised promoters
  fileName = paste0(H3K4me3_peak_path,EID,'_H3K4me3_peak.bed')
  H3K4me3 = read.table(fileName)
  colnames(H3K4me3) = c('chrom', 'start', 'end')
  H3K27me3 = peak_109[[i]]
  diff_poised_promoters = bed_intersect(CGI_promoters, H3K27me3[,1:3], suffix = c('','.y')) %>% 
    select(chrom, start, end) %>% 
    unique() %>% 
    bed_intersect(., H3K4me3, suffix = c('','.y')) %>% 
    select(chrom, start, end) %>% 
    unique()
  
  #LOCK, short LOCK
  LOCK = read.table(file = paste0('data_Roadmap/03_LOCK/H3K27me3/',EID,'/',EID,'_H3K27me3_ctoff-0.7.bed')) %>% 
    dplyr::rename(chrom = 1, start = 2, end = 3) %>% 
    mutate(d = end - start)
  short_LOCK = filter(LOCK, d < 100000)
  LOCK$LOCK_length_anno = 'Long_LOCK'
  LOCK[LOCK$d < 100000,]$LOCK_length_anno = 'Short_LOCK'
  LOCK$d = NULL
  
  #peak annotation(long,short,typical)
  LOCK_anno = bed_intersect(H3K27me3[,1:3], LOCK, suffix = c('','.LOCK')) %>% 
    select(chrom, start, end, LOCK_length_anno.LOCK) %>% 
    dplyr::rename(LOCK_length_anno = LOCK_length_anno.LOCK)
  H3K27me3 = left_join(H3K27me3, LOCK_anno, by = c('chrom','start','end'))
  H3K27me3$peakClass = 'Typical peaks'
  H3K27me3[H3K27me3$LOCK_length_anno == 'Long_LOCK' & (!is.na(H3K27me3$LOCK_length_anno)),]$peakClass = 'Peaks in long LOCKs'
  H3K27me3[H3K27me3$LOCK_length_anno == 'Short_LOCK' & (!is.na(H3K27me3$LOCK_length_anno)),]$peakClass = 'Peaks in short LOCKs'
  
  #Poised promoters have H3K27me3, but which type of H3K27me3 are they: long, short, or typical?
  diff_poised_promoters_peakClass = bed_intersect(diff_poised_promoters, H3K27me3[,c('chrom','start','end','peakClass')], suffix = c('','.y')) %>% 
    select(chrom, start, end, peakClass.y) %>% 
    dplyr::rename(peakClass = peakClass.y) %>% 
    mutate(value = 1) %>% 
    unique() %>% 
    pivot_wider(id_cols = everything(), names_from = 'peakClass', values_from = 'value', values_fill = 0)
  diff_poised_promoters_peakClass$`More than one peak class` = rowSums(diff_poised_promoters_peakClass[,4:6])
  diff_poised_promoters_peakClass$`Peak class` = NA
  if (sum(diff_poised_promoters_peakClass$`More than one peak class` > 1) > 0) {
    diff_poised_promoters_peakClass[diff_poised_promoters_peakClass$`More than one peak class` > 1,]$`Peak class` = 'More than one peak class'
  }
  if (sum(diff_poised_promoters_peakClass$`More than one peak class` <= 1 & 
          diff_poised_promoters_peakClass$`Typical peaks` == 1) > 0) {
    diff_poised_promoters_peakClass[diff_poised_promoters_peakClass$`More than one peak class` <= 1 & 
                                      diff_poised_promoters_peakClass$`Typical peaks` == 1,]$`Peak class` = 'Typical peaks'
  }
  if (sum(diff_poised_promoters_peakClass$`More than one peak class` <= 1 & 
          diff_poised_promoters_peakClass$`Peaks in short LOCKs` == 1) > 0) {
    diff_poised_promoters_peakClass[diff_poised_promoters_peakClass$`More than one peak class` <= 1 & 
                                      diff_poised_promoters_peakClass$`Peaks in short LOCKs` == 1,]$`Peak class` = 'Peaks in short LOCKs'
  }
  if (sum(diff_poised_promoters_peakClass$`More than one peak class` <= 1 & 
          diff_poised_promoters_peakClass$`Peaks in long LOCKs` == 1) > 0) {
    diff_poised_promoters_peakClass[diff_poised_promoters_peakClass$`More than one peak class` <= 1 & 
                                      diff_poised_promoters_peakClass$`Peaks in long LOCKs` == 1,]$`Peak class` = 'Peaks in long LOCKs'
  }
  diff_poised_promoters_peakClass = select(diff_poised_promoters_peakClass, chrom, start, end, `Peak class`)
  
  
  #poised promoters within short LOCKs
  diff_poised_promoters1 = bed_intersect(diff_poised_promoters, short_LOCK, suffix = c('','.LOCK')) %>% 
    mutate(d = end - start) %>% 
    group_by(chrom, start, end, d) %>% 
    summarise(overlap = sum(.overlap), .groups = 'drop') %>% 
    filter(overlap > d/2) %>% 
    select(chrom, start, end)
  ESC_poised_promoters1 =  bed_intersect(ESC_poised_promoters, short_LOCK, suffix = c('','.LOCK')) %>% 
    mutate(d = end - start) %>% 
    group_by(chrom, start, end, d) %>% 
    summarise(overlap = sum(.overlap), .groups = 'drop') %>% 
    filter(overlap > d/2) %>% 
    select(chrom, start, end)
  
  #CGI promoters within short LOCKs(with H3K27me3 peaks)
  short_LOCK_CGI_promoter = bed_intersect(CGI_promoters, short_LOCK, suffix = c('','.LOCK')) %>% 
    mutate(d = end - start) %>% 
    group_by(chrom, start, end, d) %>% 
    summarise(overlap = sum(.overlap), .groups = 'drop') %>% 
    filter(overlap > d/2) %>% 
    select(chrom, start, end) %>% 
    bed_intersect(H3K27me3[,1:3], suffix = c('','.y')) %>% 
    select(-start.y, -end.y, -.source, -.overlap) %>% 
    unique()
  #CGI promoters with H3K27me3 peak binding
  K27_CGI_promoters = bed_intersect(CGI_promoters, H3K27me3[,1:3], suffix = c('','.y')) %>% 
    select(chrom, start, end) %>% 
    unique()
  #ESC poised promoter within K27_CGI_promoters(counts if overlap) 
  #(not for diff_poised_promoters, because diff_poised_promoters are all within K27_CGI_promoters)
  ESC_poised_promoters2 = bed_intersect(ESC_poised_promoters, K27_CGI_promoters, suffix = c('','.y')) %>% 
    select(chrom,start,end) %>% 
    unique()
  #The relationship between each sample's poised promoter and ESC poised promoter, regardless of whether they are in the short LOCK.
  inter = bed_intersect(diff_poised_promoters, ESC_poised_promoters)
  diff_overlap = select(inter, chrom, start.x, end.x) %>% unique()%>% nrow(); diff_notOverlap = nrow(diff_poised_promoters) - diff_overlap
  ESC_overlap = select(inter, chrom, start.y, end.y) %>% unique()%>% nrow(); ESC_notOverlap = nrow(ESC_poised_promoters) - ESC_overlap
  
  plotData_diff_all[[i]] = rbind(data.frame(group = 'Overlapping with ESC poised promoters', number = diff_overlap, EID = EID),
                                 data.frame(group = 'Not overlapping with ESC poised promoters', number = diff_notOverlap, EID = EID)) %>% 
    mutate(percent = number/nrow(diff_poised_promoters)*100)
  
  #The relationship between each sample's poised promoter and ESC poised promoter within the short LOCK.
  inter = bed_intersect(diff_poised_promoters1, ESC_poised_promoters1)
  diff_overlap = select(inter, chrom, start.x, end.x) %>% unique()%>% nrow(); diff_notOverlap = nrow(diff_poised_promoters1) - diff_overlap
  ESC_overlap = select(inter, chrom, start.y, end.y) %>% unique()%>% nrow(); ESC_notOverlap = nrow(ESC_poised_promoters1) - ESC_overlap
  
  plotData_diff_withinShortLOCK[[i]] = rbind(data.frame(group = 'Overlapping with ESC poised promoters', number = diff_overlap, EID = EID),
                                             data.frame(group = 'Not overlapping with ESC poised promoters', number = diff_notOverlap, EID = EID)) %>% 
    mutate(percent = number/nrow(diff_poised_promoters1)*100)
  
  
  #The proportion of poised promoters among the CGI promoters with peaks within the short LOCK in each sample .
  plotData_toShortLOCK_diff[[i]] = rbind(data.frame(percent = (nrow(diff_poised_promoters1)/nrow(short_LOCK_CGI_promoter))*100, EID = EID, group = 'poised promoters'),
                                         data.frame(percent = 100-(nrow(diff_poised_promoters1)/nrow(short_LOCK_CGI_promoter)*100), EID = EID, group = 'Other CGI promoters'))
  #The proportion of poised promoters in each sample among all K27 CGI promoters.
  plotData_toK27CGI_diff[[i]] = rbind(data.frame(percent = (nrow(diff_poised_promoters)/nrow(K27_CGI_promoters))*100, EID = EID, group = 'poised promoters'),
                                      data.frame(percent = 100-(nrow(diff_poised_promoters)/nrow(K27_CGI_promoters)*100), EID = EID, group = 'Other CGI promoters'))
  
}
plotData_diff_all  = do.call(rbind, plotData_diff_all)
plotData_diff_withinShortLOCK  = do.call(rbind, plotData_diff_withinShortLOCK)
plotData_toShortLOCK_diff  = do.call(rbind, plotData_toShortLOCK_diff)
plotData_toK27CGI_diff  = do.call(rbind, plotData_toK27CGI_diff)

plotData_percent_of_poised_to_CGI_wtih_peak = rbind(plotData_toK27CGI_diff %>% mutate(facet = 'Genome-wide'),
                                                    plotData_toShortLOCK_diff %>% mutate(facet = 'Within short LOCKs')) %>% 
  filter(group == 'poised promoters') %>% 
  mutate(facet = factor(facet, levels = c('Genome-wide','Within short LOCKs')))

plotData_poised_inter_ESC_poised = rbind(plotData_diff_all %>% mutate(facet = 'Genome-wide'), 
                                         plotData_diff_withinShortLOCK %>% mutate(facet = 'Within short LOCKs')) %>% 
  filter(group == 'Overlapping with ESC poised promoters') %>% 
  mutate(facet = factor(facet, levels = c('Genome-wide', 'Within short LOCKs')))

plot_percent_of_poised_to_CGI_with_peak = ggplot(plotData_percent_of_poised_to_CGI_wtih_peak, aes(x = facet, y = percent))+
  geom_boxplot(aes(fill = facet))+
  theme_classic()+
  xlab('')+
  ylab('% of poised promoters to\nCGI promoters with\nH3K27me3 peaks')+
  geom_signif(comparisons = list(c('Genome-wide','Within short LOCKs')),
              test = 't.test', test.args = c(var.equal = T, paired = T), map_signif_level = T)+
  scale_fill_manual(values = c('Genome-wide' = '#c57ea5', 'Within short LOCKs' = '#cabad4'))+
  guides(fill = 'none')+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = 'black'))

plot_poised_inter_ESC_poised = ggplot(plotData_poised_inter_ESC_poised, aes(x = facet, y = percent))+
  geom_boxplot(aes(fill = facet))+
  theme_classic()+
  xlab('')+
  ylab('% of poised promoters\noverlapping with ESC\npoised promoters')+
  geom_signif(comparisons = list(c('Genome-wide','Within short LOCKs')),
              test = 't.test', test.args = c(var.equal = T, paired = T), map_signif_level = T)+
  scale_fill_manual(values = c('Genome-wide' = '#c57ea5', 'Within short LOCKs' = '#cabad4'))+
  guides(fill = 'none')+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = 'black'))
Fig5C = plot_grid(plot_percent_of_poised_to_CGI_with_peak, plot_poised_inter_ESC_poised, nrow = 2)

pdf('plot/Fig5C/Fig5C.pdf', height = 6, width = 2)
Fig5C
dev.off()

Data_Fig5C_up = plotData_percent_of_poised_to_CGI_wtih_peak
Data_Fig5C_down = plotData_poised_inter_ESC_poised
write.table(Data_Fig5C_up, file = 'plot/Fig5C/Data_Fig5C_up.txt', sep = '\t', quote = F, row.names = F, col.names = T)
write.table(Data_Fig5C_down, file = 'plot/Fig5C/Data_Fig5C_down.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig5D. The number of CGI promoters within Tumor-loss short LOCKs, shared short LOCKs, and tumor-gain short LOCKs identified across different normal-tumor pairs.----
all_tissue_types = pairs_for_shortLOCKAnalysis
all_tissue_type_names = c('NE2 vs KYSE450\n(ESCC)',
                          'NE2 vs KYSE510\n(ESCC)',
                          'MCF10A vs HCC1954\n(BRCA Her2+)',
                          'MCF10A vs MCF7\n(BRCA Luminal A)',
                          'MCF10A vs HCC1937\n(BRCA TN-Basal)',
                          'MCF10A vs MB231\n(BRCA TN-CL)')
plotData_number_of_CGI_promoters = list()
for (t in 1:6) {
  tissue_type = all_tissue_types[t]
  promoter_shortLOCK = read.table(paste0('data_CellLine/15_short_LOCK_promoter_feature/',tissue_type,'.bed'), 
                                  sep = '\t', header = T)
  plotData_number_of_CGI_promoters[[t]] = filter(promoter_shortLOCK, CGI == T) %>% 
    group_by(group.shortLOCK) %>% 
    summarise(prom_number = n()) %>% 
    mutate(tissue_type = all_tissue_type_names[t])
}
plotData_number_of_CGI_promoters = do.call(rbind, plotData_number_of_CGI_promoters)
plotData_number_of_CGI_promoters$group.shortLOCK = sub('LOCK','short LOCK',plotData_number_of_CGI_promoters$group.shortLOCK)
plotData_number_of_CGI_promoters$group.shortLOCK = factor(plotData_number_of_CGI_promoters$group.shortLOCK, 
                                                          levels = rev(c('Tumor-loss short LOCK',
                                                                         'Shared short LOCK',
                                                                         'Tumor-gain short LOCK')))
plotData_number_of_CGI_promoters$tissue_type = factor(plotData_number_of_CGI_promoters$tissue_type, 
                                                      levels = rev(all_tissue_type_names))
Fig5D = ggplot(plotData_number_of_CGI_promoters, aes(y = tissue_type, x = prom_number, fill = group.shortLOCK))+
  geom_bar(stat = 'identity', position = 'stack')+
  scale_fill_manual(values = c('Tumor-gain short LOCK' = '#f9b317',
                               'Shared short LOCK' = '#584439',
                               'Tumor-loss short LOCK' = '#a0302f'))+
  theme_classic()+
  xlab('# of CGI promoters')+
  ylab('')+
  theme(axis.text = element_text(colour = 'black'),
        legend.title = element_blank())

pdf('plot/Fig5D/Fig5D.pdf', width = 7, height = 4)
Fig5D
dev.off()

Data_Fig5D = plotData_number_of_CGI_promoters
write.table(Data_Fig5D,'plot/Fig5D/Data_Fig5D.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig5F. Number of CGI promoters within upregulated or hypermethylated group.----
plotData_prom_number_two_groups = list()
for (t in 1:6) {
  tissue_type = all_tissue_types[t]
  promoter_shortLOCK = read.table(paste0('data_CellLine/15_short_LOCK_promoter_feature/',tissue_type,'.bed'), sep = '\t', header = T)
  plotData_prom_number_two_groups[[t]] = filter(promoter_shortLOCK, CGI == T & group.shortLOCK == 'Tumor-loss LOCK' & !is.na(group)) %>% 
    group_by(group) %>% 
    summarise(prom_number = n()) %>% 
    mutate(tissue_type = all_tissue_type_names[t])
}
plotData_prom_number_two_groups = do.call(rbind, plotData_prom_number_two_groups)
plotData_prom_number_two_groups$tissue_type = factor(plotData_prom_number_two_groups$tissue_type, levels = rev(all_tissue_type_names))
plotData_prom_number_two_groups$group = factor(plotData_prom_number_two_groups$group, levels = c('Upregulated','Hypermethylated'))

Fig5F = ggplot(plotData_prom_number_two_groups, aes(y = tissue_type, x = prom_number, fill = group))+
  geom_bar(stat = 'identity', position = 'dodge', color = 'black')+
  theme_classic()+
  xlab('# of CGI promoters')+
  ylab('')+
  scale_fill_manual(values = c('Upregulated' = '#f2d6a7', 'Hypermethylated' = '#babadb'))+
  geom_text(aes(label = prom_number, group = group), hjust = -0.5, color = "black", size = 4, position = position_dodge(0.9))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = 'black'),
        legend.title = element_blank())

pdf('plot/Fig5F/Fig5F.pdf', width = 7, height = 4)
Fig5F
dev.off()

Data_Fig5F = plotData_prom_number_two_groups
write.table(Data_Fig5F, 'plot/Fig5F/Data_Fig5F.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig5G/Figure_S6. In each pair, the H3K27me3 intensity, DNA methylation, gene expression, H3K27ac intensity, and H3K4me3 intensity of upregulated or hypermethylated CGI promoters. ----
p_featuresTwoGroup = list()
plotData_featuresTwoGroup = list()
for (t in 1:6) {
  tissue_type = all_tissue_types[t]
  plotData_featuresTwoGroup[[t]] = list()
  promoter_shortLOCK = read.table(paste0('data_CellLine/15_short_LOCK_promoter_feature/',tissue_type,'.bed'), sep = '\t', header = T)
  promoter_shortLOCK_CGI_loss = filter(promoter_shortLOCK, CGI == T & group.shortLOCK == 'Tumor-loss LOCK' & !is.na(group))
  promoter_shortLOCK_CGI_loss$group = factor(promoter_shortLOCK_CGI_loss$group, levels = c('Upregulated','Hypermethylated'))
  if (t == 6) {
    colnames(promoter_shortLOCK_CGI_loss) = sub('MDA_MB_','MB',colnames(promoter_shortLOCK_CGI_loss))
  }
  # b)Intensity----
  plotData = promoter_shortLOCK_CGI_loss[,c(colnames(promoter_shortLOCK_CGI_loss)[grep('H3K27me3', colnames(promoter_shortLOCK_CGI_loss))], 'group.shortLOCK', 'group')] %>% 
    pivot_longer(., cols = 1:(ncol(.)-2), names_to = 'sample', values_to = 'intensity')
  plotData[,c('TorN','xlab')] = str_match(plotData$sample, '.*\\.(.*?)\\.(.*?)\\..*?$')[,2:3]
  plotData$group.shortLOCK = factor(plotData$group.shortLOCK, levels = c('Tumor-loss LOCK','Shared LOCK','Tumor-gain LOCK'))
  xlabFactor = plotData[,c('TorN','xlab')] %>% unique() %>% mutate(TorN = factor(TorN, levels = c('N','T'))) %>% arrange(TorN) %>% pull(xlab)
  plotData$xlab = factor(plotData$xlab, levels = xlabFactor)
  plotData$TorN = factor(plotData$TorN, levels = c('N','T'))
  plotData_featuresTwoGroup[[t]][[1]] = plotData
  p_intensity = ggplot(plotData, aes(x = xlab, y = intensity))+
    geom_boxplot(aes(color = TorN), outlier.size = 0.1)+
    facet_grid(.~group)+
    xlab('')+
    ylab('H3K27me3 intensity')+
    theme_classic()+
    theme(strip.background = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_color_manual(values = c('N' = '#2694ab',
                                  'T' = '#eb7072'))+
    geom_signif(comparisons = list(xlabFactor),
                map_signif_level = T)+
    ggtitle('H3K27me3')+
    guides(color = 'none')
  # c)Other markers----
  OtherHis = c('H3K27ac','H3K4me3')
  p_otherHis = list()
  m = 4
  for (o in 1:length(OtherHis)) {
    oneOtherHis = OtherHis[o]
    cols_name = colnames(promoter_shortLOCK_CGI_loss)[grepl(oneOtherHis, colnames(promoter_shortLOCK_CGI_loss))]
    plotData = promoter_shortLOCK_CGI_loss[,c(cols_name, 'group.shortLOCK', 'group')] %>% 
      pivot_longer(., cols = 1:(ncol(.)-2), names_to = 'sample', values_to = 'intensity')
    plotData[,c('TorN','xlab')] = str_match(plotData$sample, '.*\\.(.*?)\\.(.*?)\\..*?$')[,2:3]
    xlabFactor = plotData[,c('TorN','xlab')] %>% unique() %>% mutate(TorN = factor(TorN, levels = c('N','T'))) %>% arrange(TorN) %>% pull(xlab)
    plotData$xlab = factor(plotData$xlab, levels = xlabFactor)
    plotData$TorN = factor(plotData$TorN, levels = c('N','T'))
    plotData_featuresTwoGroup[[t]][[m]] = plotData
    p_otherHis[[o]] = ggplot(plotData, aes(x = xlab, y = intensity))+
      geom_boxplot(aes(color = TorN), outlier.size = 0.1)+
      facet_grid(.~group)+
      xlab('')+
      ylab(paste0(oneOtherHis, ' intensity'))+
      theme_classic()+
      theme(strip.background = element_blank(),
            legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1))+
      scale_color_manual(values = c('N' = '#2694ab',
                                    'T' = '#eb7072'))+
      geom_signif(comparisons = list(xlabFactor),
                  map_signif_level = T)+
      ggtitle(oneOtherHis)+
      guides(color = 'none')
    m = m+1
  }
  p_H3K27ac = p_otherHis[[1]]
  p_H3K4me3 = p_otherHis[[2]]
  # d)DNA methylation----
  pattern = ifelse(grepl('TN_Basal|TN_CL',tissue_type), 'EPIC.metSignal', 'WGBS.metSignal')
  plotData = promoter_shortLOCK_CGI_loss[,c(colnames(promoter_shortLOCK_CGI_loss)[grep(pattern, colnames(promoter_shortLOCK_CGI_loss))], 'group.shortLOCK', 'group')] %>% 
    pivot_longer(., cols = 1:(ncol(.)-2), names_to = 'sample', values_to = 'DNA_methylation')
  plotData[,c('TorN','xlab')] = str_match(plotData$sample, '.*\\.(.*?)\\.(.*?)\\..*?\\.metSignal')[,2:3]
  plotData$group.shortLOCK = factor(plotData$group.shortLOCK, levels = c('Tumor-loss LOCK','Shared LOCK','Tumor-gain LOCK'))
  xlabFactor = plotData[,c('TorN','xlab')] %>% unique() %>% mutate(TorN = factor(TorN, levels = c('N','T'))) %>% arrange(TorN) %>% pull(xlab)
  plotData$xlab = factor(plotData$xlab, levels = xlabFactor)
  plotData$TorN = factor(plotData$TorN, levels = c('N','T'))
  dataType = ifelse(grepl('Basal',tissue_type), '(EPIC)', '(WGBS)')
  plotData_featuresTwoGroup[[t]][[2]] = plotData
  p_methylation = ggplot(plotData, aes(x = xlab, y = DNA_methylation))+
    geom_boxplot(aes(color = TorN), outlier.size = 0.1)+
    facet_grid(.~group)+
    xlab('')+
    ylab('β values')+
    theme_classic()+
    theme(strip.background = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_color_manual(values = c('N' = '#2694ab',
                                  'T' = '#eb7072'))+
    ggtitle(paste0('DNA methylation ',dataType))+
    guides(color = 'none')
  if (!tissue_type %in% c('MCF10A_normal_MCF7_LuminalA', 'MCF10A_normal_HCC1937_TN_Basal', 'MCF10A_normal_MB231_TN_CL')) {
    p_methylation = p_methylation + geom_signif(comparisons = list(xlabFactor), map_signif_level = T)  
  }
  # e)gene expression----
  plotData = promoter_shortLOCK_CGI_loss[,c('FPKM_T', 'FPKM_N', 'group.shortLOCK', 'group')] %>% 
    pivot_longer(., cols = 1:(ncol(.)-2), names_to = 'sample', values_to = 'Gene_expression')
  plotData$TorN = str_match(plotData$sample, 'FPKM_(.*)')[,2]
  df = as.data.frame(str_match(grep('\\.FPKM',colnames(promoter_shortLOCK_CGI_loss), value=T), '\\.*.([NT])\\.(.*?)\\.FPKM')[,2:3]) %>% 
    dplyr::rename(TorN=1, xlab=2) %>% 
    mutate(TorN = factor(TorN, levels = c('N','T'))) %>% 
    arrange(TorN) %>% 
    mutate(xlab = factor(xlab, levels = xlab))
  xlabFactor = levels(df$xlab)
  plotData = left_join(plotData, df, by = 'TorN')
  plotData$group.shortLOCK = factor(plotData$group.shortLOCK, levels = c('Tumor-loss LOCK','Shared LOCK','Tumor-gain LOCK'))
  plotData_featuresTwoGroup[[t]][[3]] = plotData
  p_expression = ggplot(plotData, aes(x = xlab, y = log10(Gene_expression+0.1)))+
    geom_boxplot(aes(color = TorN), outlier.size = 0.1)+
    facet_grid(.~group)+
    xlab('')+
    ylab('log10(FPKM+0.1)')+
    theme_classic()+
    theme(strip.background = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_color_manual(values = c('N' = '#2694ab',
                                  'T' = '#eb7072'))+
    geom_signif(comparisons = list(xlabFactor),
                map_signif_level = T)+
    ggtitle('Gene expression')+
    guides(color = 'none')
  
  p_featuresTwoGroup[[t]] = plot_grid(p_intensity, p_methylation, p_expression, p_H3K27ac, p_H3K4me3, nrow = 1)
  names(plotData_featuresTwoGroup[[t]]) = c('H3K27me3','DNA_methylation','Gene_expression',
                                            'H3K27ac','H3K4me3')
}
Fig5G = p_featuresTwoGroup[[1]]
Figure_S6 = plot_grid(p_featuresTwoGroup[[2]],
                      p_featuresTwoGroup[[3]],
                      p_featuresTwoGroup[[4]],
                      p_featuresTwoGroup[[5]],
                      p_featuresTwoGroup[[6]],
                      ncol = 1)
pdf('plot/Fig5G/Fig5G.pdf', width = 13, height = 4)
Fig5G
dev.off()

pdf('plot/FigS6/FigS6.pdf', width = 13, height = 18)
Figure_S6
dev.off()

Data_Fig5G = plotData_featuresTwoGroup[[1]]
Data_FigS6.NE2_KYSE510 = plotData_featuresTwoGroup[[2]]
Data_FigS6.MCF10A_normal_HCC1954_Her2Pos = plotData_featuresTwoGroup[[3]]
Data_FigS6.MCF10A_normal_MCF7_LuminalA = plotData_featuresTwoGroup[[4]]
Data_FigS6.MCF10A_normal_HCC1937_TN_Basal = plotData_featuresTwoGroup[[5]]
Data_FigS6.MCF10A_normal_MB231_TN_CL = plotData_featuresTwoGroup[[6]]

save(Data_Fig5G, file = 'plot/Fig5G/Data_Fig5G.RData')
save(Data_FigS6.NE2_KYSE510, file = 'plot/FigS6/Data_FigS6.NE2_KYSE510.RData')
save(Data_FigS6.MCF10A_normal_HCC1954_Her2Pos, file = 'plot/FigS6/Data_FigS6.MCF10A_normal_HCC1954_Her2Pos.RData')
save(Data_FigS6.MCF10A_normal_MCF7_LuminalA, file = 'plot/FigS6/Data_FigS6.MCF10A_normal_MCF7_LuminalA.RData')
save(Data_FigS6.MCF10A_normal_HCC1937_TN_Basal, file = 'plot/FigS6/Data_FigS6.MCF10A_normal_HCC1937_TN_Basal.RData')
save(Data_FigS6.MCF10A_normal_MB231_TN_CL, file = 'plot/FigS6/Data_FigS6.MCF10A_normal_MB231_TN_CL.RData')

# Fig5H. Percentage of poised promoters among upregulated or hypermethylated CGI promoters.----
# 1) get MCF10A_poised and NE2_poised
#promoters
CGI_promoters_hg38 = read.table(file = 'meta/LOLACore/hg38/ESC_poised_promoters/CGI_promoters_hg38.bed') %>% dplyr::rename(chrom=1,start=2,end=3)
ESC_poised_promoters_hg38 = read.table(file = 'meta/LOLACore/hg38/ESC_poised_promoters/ESC_poised_promoters_hg38.bed') %>% dplyr::rename(chrom=1,start=2,end=3)#是CGI的子集
#H3K4me3 peaks 
MCF10A_H3K4me3 = read.table(file = 'data_CellLine/03_ChIPSeq_Peak/H3K4me3/MCF10A_H3K4me3_peaks.narrowPeak', skip = 1)[,1:3] %>% dplyr::rename(chrom=1,start=2,end=3)
NE2_H3K4me3 = read.table(file = 'data_CellLine/03_ChIPSeq_Peak/H3K4me3/NE2_H3K4me3_peaks.narrowPeak', skip = 1)[,1:3] %>% dplyr::rename(chrom=1,start=2,end=3)
#H3K27me3 peaks
MCF10A_H3K27me3 = read.table(file = 'data_CellLine/03_ChIPSeq_Peak/H3K27me3/MCF10A_normal_HCC1954_Her2Pos.N.MCF10A.H3K27me3.peak.bed') %>% dplyr::rename(chrom=1,start=2,end=3)
NE2_H3K27me3 = read.table(file = 'data_CellLine/03_ChIPSeq_Peak/H3K27me3//NE2_KYSE450.N.NE2.H3K27me3.peak.bed') %>% dplyr::rename(chrom=1,start=2,end=3)
#poised promoters
MCF10A_poised = bed_intersect(CGI_promoters_hg38, MCF10A_H3K4me3, suffix = c('','.y')) %>% 
  .[,1:3] %>% 
  unique() %>% 
  bed_intersect(., MCF10A_H3K27me3, suffix = c('','.y')) %>% 
  .[,1:3] %>% 
  unique()
NE2_poised = bed_intersect(CGI_promoters_hg38, NE2_H3K4me3, suffix = c('','.y')) %>% 
  .[,1:3] %>% 
  unique() %>% 
  bed_intersect(., NE2_H3K27me3, suffix = c('','.y')) %>% 
  .[,1:3] %>% 
  unique()

# 2) get plot
all_tissue_types = pairs_for_shortLOCKAnalysis
all_tissue_type_names = c('NE2 vs KYSE450\n(ESCC)',
                          'NE2 vs KYSE510\n(ESCC)',
                          'MCF10A vs HCC1954\n(BRCA Her2+)',
                          'MCF10A vs MCF7\n(BRCA Luminal A)',
                          'MCF10A vs HCC1937\n(BRCA TN-Basal)',
                          'MCF10A vs MB231\n(BRCA TN-CL)')
plotData = list()
for (i in 1:length(all_tissue_types)) {
  tissue_type = all_tissue_types[i]
  tissue_type_name = all_tissue_type_names[i]
  file = read.table(paste0('data_CellLine/15_short_LOCK_promoter_feature/',tissue_type,'.bed'), sep = '\t', header = T) %>% 
    filter(CGI == T & group.shortLOCK == 'Tumor-loss LOCK' & !is.na(group)) %>% 
    select(chrom, start, end, group)
  if (i %in% 1:2) {
    cell_poised = NE2_poised
  }else{
    cell_poised = MCF10A_poised
  }
  poised = bed_intersect(file, cell_poised, suffix = c('','.y')) %>% 
    select(chrom, start, end, group) %>% 
    unique()
  df = left_join(as.data.frame(table(file$group)), 
                 as.data.frame(table(poised$group)),
                 by = 'Var1',
                 suffix = c('.CGI_promoter','.poised_promoter')) %>% 
    mutate(poised_ratio = Freq.poised_promoter/Freq.CGI_promoter) %>% 
    dplyr::rename(group = Var1)
  df$group = factor(df$group, levels = c('Upregulated','Hypermethylated','No change'))
  df$tissue_type = tissue_type_name
  plotData[[i]] = df
}
plotData = do.call(rbind, plotData)
plotData$tissue_type = factor(plotData$tissue_type, levels = rev(all_tissue_type_names))
plotData[is.na(plotData)] = 0
Fig5H = ggplot(plotData, aes(x = tissue_type, y = poised_ratio*100, fill = group))+
  geom_bar(stat = 'identity', position = 'dodge', color = 'black')+
  theme_classic()+
  theme(plot.title = element_text(size = 8, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = 'black'),
        legend.title = element_blank())+
  xlab('')+
  ylab('% of poised promoters')+
  scale_fill_manual(values = c('Upregulated' = '#f9e2ae',
                               'Hypermethylated' = '#c8c5e4',
                               'No change' = '#5ca8b8'))+
  coord_flip()

pdf('plot/Fig5H/Fig5H.pdf', width = 7, height = 4)
Fig5H
dev.off()

Data_Fig5H = plotData
write.table(Data_Fig5H, file = 'plot/Fig5H/Data_Fig5H.txt', 
            sep = '\t', quote = F, row.names = F, col.names = T)

# Combine the plots in Figure 5----
Figure_5 = ggarrange(ggarrange(
  ggarrange(ggarrange(Fig5A, Fig5B, ncol = 1, heights = c(1,3), labels = c('A','B'), font.label = list(size = 23)),
            Fig5C,
            nrow = 1, widths = c(2,1)),
  ggplot(), nrow = 1, widths = c(6,5)),
  ggarrange(Fig5D+guides(fill = 'none'), Fig5F+guides(fill = 'none'), Fig5H+guides(fill = 'none'), nrow = 1, labels = c('D','F','H'), font.label = list(size = 23)),
  Fig5G,
  ncol = 1, heights = c(2,1,1), labels = c('','','G'), font.label = list(size = 23))
pdf('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/Figure_5.pdf', width = 13, height = 15)
Figure_5
dev.off()
# Combine the plots in Figure S6----
pdf('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/Figure_S6.pdf', width = 13, height = 18)
Figure_S6
dev.off()

# Figure 6/ Figure S7. ----

# Fig6A. Bubble plot shows motifs identified using Homer software with upregulated CGI-promoters as foreground and hypermethylated CGI-promoters as background.----
# 1) run motif analysis
runMotif = F
for (t in 1:6) {
  tissue_type = all_tissue_types[t]
  motifPath = 'data_CellLine/16_motif/'
  dir.create(paste0(motifPath,'/input'), showWarnings = F)
  dir.create(paste0(motifPath,'/output'), showWarnings = F)
  dir.create(paste0(motifPath,'/log'), showWarnings = F)
  promoter_shortLOCK = read.table(paste0('data_CellLine/15_short_LOCK_promoter_feature/',tissue_type,'.bed'),
                                  sep = '\t', header = T)
  promoter_shortLOCK_CGI_loss = filter(promoter_shortLOCK, CGI == T & group.shortLOCK == 'Tumor-loss LOCK' & !is.na(group))
  #fg
  fg = promoter_shortLOCK_CGI_loss %>% 
    filter(group == 'Upregulated') %>% 
    select(chrom, start, end)
  bg = promoter_shortLOCK_CGI_loss %>% 
    filter(group == 'Hypermethylated') %>% 
    select(chrom, start, end)
  write.table(fg, 
              file = paste0(motifPath,'input/',tissue_type,'.fg.bed'),
              sep = '\t', quote = F, row.names = F, col.names = F)
  write.table(bg, 
              file = paste0(motifPath,'input/',tissue_type,'.bg.bed'),
              sep = '\t', quote = F, row.names = F, col.names = F)
  
  if (runMotif) {
    outputPath=paste0(motifPath,'/output/')
    logPath=paste0(motifPath,'/log/')
    
    path_env <- Sys.getenv("PATH")
    if (!grepl('homer',path_env)) {
      path_to_add = '/opt/homer/.//bin///'
      Sys.setenv(PATH = paste(path_to_add, path_env, sep = ":"))
    }
    runMotif = function(fgFile, bgFile, outputPath, logFile){
      system(paste0('/opt/homer/.//bin///findMotifsGenome.pl ',fgFile,' hg38 ',outputPath,' -p 5 -bg ',bgFile,' -keepOverlappingBg',' > ',logFile))
      
    }
    runMotif(fgFile = paste0(motifPath,'input/',tissue_type,'.fg.bed'),
             bgFile = paste0(motifPath,'input/',tissue_type,'.bg.bed'),
             outputPath = paste0(outputPath,'/',tissue_type),
             logFile = paste0(logPath,'/',tissue_type,'.log'))
  }
}

# 2) read results
knownResList = list()
FPKMList = list()
for (i in 1:length(all_tissue_types)) {
  tissue_type = all_tissue_types[i]
  #motif result
  fileName = paste0('data_CellLine/16_motif/output/',tissue_type,'/knownResults.txt')
  knownRes = read.table(file = fileName, sep = '\t', header = F, skip = 1)[,c(1,5)] %>% 
    dplyr::rename(motif_source = 1, FDR = 2) %>% 
    filter(FDR < 0.05)
  if (nrow(knownRes) == 0) {next}
  knownRes$motif = str_match(knownRes$motif_source, '(.*?)\\/.*')[,2]
  knownRes$motifName = str_match(knownRes$motif, '(.*?)\\(.*')[,2]
  knownRes$geneName = toupper(knownRes$motifName)
  knownResList[[i]] = knownRes %>% mutate(tissue_type = tissue_type)
  #FPKM list
  FPKMList[[i]] = read.csv(file = paste0('data_CellLine/05_RNASeq/FPKM_multiSample/',tissue_type,'_FPKM.csv'))
}
knownResList = do.call(rbind, knownResList)
allGenes = unique(do.call(c, sapply(FPKMList, function(df){rownames(df)})))
knownResList$correctName = knownResList$geneName %in% allGenes
unique(knownResList[knownResList$correctName == F,]$geneName)
# [1] "CRE"         "E2A"         "AP-2ALPHA"   "SLUG"        "SNAIL1"      "AR-HALFSITE" "RXR"         "BMAL1"       "N-MYC"       "TATA-BOX"    "HEB"         "ZFP809" 

# 3) Edit the gene names
knownResList[knownResList$geneName == 'CRE',]$geneName = 'CREB1'
knownResList[knownResList$geneName == 'E2A',]$geneName = 'TCF3'
knownResList[knownResList$geneName == 'AP-2ALPHA',]$geneName = 'TFAP2A'
knownResList[knownResList$geneName == 'SLUG',]$geneName = 'SNAI2'
knownResList[knownResList$geneName == 'AR-HALFSITE',]$geneName = 'AR'
#There are three retinoic X receptors (RXR): RXR-alpha, RXR-beta, and RXR-gamma, encoded by the RXRA, RXRB, RXRG genes, respectively. 
knownResList[knownResList$geneName == 'RXR',]$geneName = 'RXRA'
#HEB factors are the products of the TCF12 gene locus, which encodes both the canonical HEB protein (HEBCan) and a shorter variant (HEBAlt)
knownResList[knownResList$geneName == 'HEB',]$geneName = 'TCF12'
#BMAL1 and ZFP809 are gene names themselves.

# 4) Indicate whether the transcription factors (TFs) are expressed.
knownResList = split(knownResList, knownResList$tissue_type)
for (i in 1:length(knownResList)) {
  #FPKM
  tissue_type = names(knownResList)[i]
  knownResList[[i]]$FPKM.T = read.csv(file = paste0('data_CellLine/05_RNASeq/FPKM_multiSample/',tissue_type,'_FPKM.csv'), row.names = 1)[,2,drop=F] %>% 
    .[knownResList[[i]]$geneName,]
  knownResList[[i]] = filter(knownResList[[i]], FPKM.T > 1)
}
knownResList = do.call(rbind, knownResList)
knownResList$correctName = NULL
rownames(knownResList) = NULL
# 5) get the plot
knownResList$`log10(FPKM)\nin tumor` = log10(knownResList$FPKM.T)
tissue_type_df = data.frame(all_tissue_types = all_tissue_types, 
                            all_tissue_type_names = all_tissue_type_names)
knownResList = left_join(knownResList, tissue_type_df, by = c('tissue_type' = 'all_tissue_types'))
knownResList$all_tissue_type_names = factor(knownResList$all_tissue_type_names, levels = rev(all_tissue_type_names))
knownResList = arrange(knownResList, all_tissue_type_names, desc(FDR), `log10(FPKM)\nin tumor`)
knownResList = dplyr::rename(knownResList, TF = motif)
knownResList$`-log10(FDR)` = -log10(knownResList$FDR + 1e-5)
# xlabLevels
sharedTFs_3 = table(knownResList$TF) %>% as.data.frame() %>% filter(Freq == 3) %>% pull(Var1)
sharedTFs_order_3 = filter(knownResList, TF %in% sharedTFs_3) %>% 
  group_by(TF) %>% 
  summarise(FDR = sum(FDR)) %>% 
  arrange(FDR) %>% 
  pull(TF)
sharedTFs_2 = table(knownResList$TF) %>% as.data.frame() %>% filter(Freq == 2) %>% pull(Var1)
sharedTFs_order_2 = filter(knownResList, TF %in% sharedTFs_2) %>% 
  group_by(TF) %>% 
  summarise(FDR = sum(FDR)) %>% 
  arrange(FDR) %>% 
  pull(TF)
sharedTFs_order = c(sharedTFs_order_3, sharedTFs_order_2)
sharedTFs_order = c('E2A(bHLH),near_PU.1',
                    'Slug(Zf)',
                    'ETV4(ETS)',
                    'AP-2alpha(AP2)',
                    'TEAD3(TEA)',
                    'ETS1(ETS)',
                    'E2A(bHLH)')
UniqTFs = rev(unique(knownResList$TF[!knownResList$TF %in% sharedTFs_order]))
ylabLevels = c(sharedTFs_order, UniqTFs)
knownResList$TF = factor(knownResList$TF, levels = (ylabLevels))
non_EMT_TF = c('AP-2alpha(AP2)','THRb(NR)','Elf4(ETS)',"RXR(NR),DR1","Tgif1(Homeobox)","KLF5(Zf)","KLF10(Zf)","GABPA(ETS)","SPDEF(ETS)","EHF(ETS)","ELF1(ETS)")
EMT_TF = knownResList$TF[!knownResList$TF %in% non_EMT_TF]
highlight = rep('red', length(EMT_TF))
names(highlight) = EMT_TF
Other_TF = knownResList$TF[!knownResList$TF %in% EMT_TF]
others = rep('black', length(Other_TF))
names(others) = Other_TF
y_cols = c(highlight,others)
y_cols = y_cols[levels(knownResList$TF)]
coord_knownResList = knownResList
coord_knownResList$TF = factor(coord_knownResList$TF, levels = rev(levels(coord_knownResList$TF)))
coord_y_cols = rev(y_cols)
coord_knownResList$all_tissue_type_names = factor(coord_knownResList$all_tissue_type_names, levels = rev(levels(coord_knownResList$all_tissue_type_names)))
p1_TF_bubble = ggplot(coord_knownResList, aes(y = all_tissue_type_names, x = TF))+
  geom_point(aes(color = FDR, size = `log10(FPKM)\nin tumor`))+
  scale_color_gradient(low ="red",high ="blue")+
  theme_bw()+
  xlab('')+
  ylab('')+
  theme(axis.text = element_text(color = 'black'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_flip()
if (if_highlight) {
  p1_TF_bubble = p1_TF_bubble + theme(axis.text.y = element_text(colour=coord_y_cols))
}
Fig6A = p1_TF_bubble

pdf('plot/Fig6A/Fig6A.pdf', width = 7, height = 10)
Fig6A
dev.off()

Data_Fig6A = coord_knownResList
write.table(Data_Fig6A, file = 'plot/Fig6A/Data_Fig6A.txt', 
            sep = '\t', quote = F, row.names = F, col.names = T)

# Fig6B. Counting motif subfamily occurrences in all normal-tumor pairs.----
knownResList$family = str_match(knownResList$TF, '.*\\((.*)\\)')[,2]
knownResList_unique = knownResList %>% distinct(TF, .keep_all = T)
family_df = as.data.frame(table(knownResList_unique$family)) %>% arrange(desc(Freq))
colnames(family_df) = c('Motif family','Frequency')
family_df$`Motif family` = factor(family_df$`Motif family`, levels = family_df$`Motif family`)
Fig6B = ggplot(family_df, aes(x = `Motif family`, y = Frequency))+
  geom_bar(stat = 'identity', fill = '#f5a79a')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = 'black'))+
  xlab('')+
  ggtitle('Motif subfamilies')

pdf('plot/Fig6B/Fig6B.pdf', width = 7, height = 5)
Fig6B
dev.off()

Data_Fig6B = family_df
write.table(Data_Fig6B, file = 'plot/Fig6B/Data_Fig6B.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)

# FigS7A. The genomic environment around ETS1 target genes in HNSCC.---- 

# The ChIP-seq analysis processes are the same as in other cell lines. However, if it is not a histone modification ChIP-seq, 
# do not use the --nomodel --extsize 146 parameters in the macs2 callpeak step.

# 1) ChIP-seq analysis (ETS1, in SCC25)
list.files('data_CellLine/17_HNSC/ETS1_ChIPSeq_SCC25/')
# [1] "00.fastq"       "01.trim_galore" "02.mapping"     "03.MACS2"       "log" 

# 2) ChIP-seq analysis (H3K27me3, H3K27ac, H3K4me3, in SCC25)
list.files('data_CellLine/17_HNSC/SCC25_histone_ChIPSeq/')
# [1] "00.fastq"       "01.trim_galore" "02.mapping"     "03.MACS2"       "GSE103554"      "log" 

# FigS7B. The heatmap illustrates changes in the expression of ETS1 target genes after ETS1 knockdown.----
# 1) RNASeq analysis (Knock down ETS1 in SCC25) #The RNA-seq analysis process is the same as in other cell lines. 
list.files('data_CellLine/17_HNSC/ETS1_knockDown_SCC25')
# [1] "countMatrix" "diffExp" 

# 2) get ETS1 expression before and after ETS1 knock down in SCC25
sh_FPKM = read.csv('data_CellLine/17_HNSC/ETS1_knockDown_SCC25/diffExp/shETS1_SCC25_DESeq2Result.csv')
head(sh_FPKM)
# gene_name    baseMean log2FoldChange     lfcSE       pvalue         padj SCC25_shCON_Rep1 SCC25_shCON_Rep2 SCC25_shETS1_Rep1 SCC25_shETS1_Rep2
# 1     IGHV3-74    8.907212       5.980379 2.6455701 9.325687e-05 3.308857e-04       0.00000000       0.00000000          1.605240         0.7220436
# 2          NTS  507.700537       5.945502 0.2842028 2.045829e-97 5.642225e-95       0.40472731       0.46029546         24.212260        25.3733863
# 3 RP11-439C8.2   16.888355       5.333397 1.4127764 1.025622e-04 3.608076e-04       0.08053115       0.00000000          2.988107         1.9228460
# 4        POF1B   64.501557       5.280960 0.6618973 6.493898e-16 9.415796e-15       0.03124105       0.02960865          1.088020         1.3293104
# 5         MMP7 2573.457066       5.236719 0.1288368 0.000000e+00 0.000000e+00       2.47441763       2.34512468         86.346037        80.0307666
# 6       IGSF23   22.062079       5.007410 1.1016269 2.238211e-06 1.046237e-05       0.00000000       0.05904145          1.186171         1.2872230

plotData_ETS1_FPKM = t(sh_FPKM[sh_FPKM$gene_name == 'ETS1',c("SCC25_shCON_Rep1", "SCC25_shETS1_Rep1", "SCC25_shCON_Rep2", "SCC25_shETS1_Rep2")]) %>% 
  as.data.frame()%>% dplyr::rename(ETS1 = 1) %>% rownames_to_column('sample')
plotData_ETS1_FPKM$sample = sub('SCC25_','',plotData_ETS1_FPKM$sample)
plotData_ETS1_FPKM$sample = factor(plotData_ETS1_FPKM$sample, levels = c("shCON_Rep1", "shCON_Rep2", "shETS1_Rep1", "shETS1_Rep2"))
plotData_ETS1_FPKM$fill = str_match(plotData_ETS1_FPKM$sample, '(.*)_.*')[,2]

# 3) get expression of ETS1 target genes
Gene = 'ETS1'
promoter_KYSE450 = read.table('data_CellLine/15_short_LOCK_promoter_feature/NE2_KYSE450.bed', sep = '\t', header = T)
if(F){system('sh script/16_scanMotifGenomeWide.pl.sh')}
scanMotifRes = read.table('data_CellLine/18_scanMotifGenomeWide.pl/CGI_promoter_KYSE450_loss.txt', sep = '\t') %>% 
  select(V1, V2) %>% 
  dplyr::rename(motifInfor=1, promoterID=2) %>% 
  merge(., promoter_KYSE450, by = 'promoterID')
targetGenes = unique(filter(scanMotifRes, grepl(Gene,motifInfor))$gene_name)
geneList2 = c('SYK','EPCAM','C1GALT1C1L','CNTN1','COCH','KCNS3')
targetGenes = targetGenes[targetGenes %in% geneList2]

# Plot the expression of ETS1 and its target genes.
sh_FPKM1 = filter(sh_FPKM, gene_name %in% targetGenes) %>% arrange(log2FoldChange)
ylabLevels = sh_FPKM1$gene_name
sh_FPKM2 = select(sh_FPKM1, gene_name, SCC25_shCON_Rep1, SCC25_shCON_Rep2, SCC25_shETS1_Rep1, SCC25_shETS1_Rep2)
plotData_TargetExp = pivot_longer(sh_FPKM2, cols = 2:5, names_to = 'sample', values_to = 'FPKM')
plotData_TargetExp$gene_name = factor(plotData_TargetExp$gene_name, levels = ylabLevels)
plotData_TargetExp$sample = sub('SCC25_','',plotData_TargetExp$sample)
plotData_TargetExp$sample = factor(plotData_TargetExp$sample, levels = c("shCON_Rep1", "shCON_Rep2", "shETS1_Rep1", "shETS1_Rep2"))
plotData_TargetExp$fill = str_match(plotData_TargetExp$sample, '(.*)_.*')[,2]
heatData1 = plotData_ETS1_FPKM %>% 
  select(!fill) %>% 
  pivot_wider(names_from = 'sample', values_from = 'ETS1') %>% 
  mutate(gene_name = 'ETS1', .before = 1)
heatData2 = plotData_TargetExp %>% 
  select(!fill) %>% 
  pivot_wider(names_from = 'sample', values_from = 'FPKM')
heatData = rbind(heatData1, heatData2) %>% 
  column_to_rownames('gene_name') %>% 
  select("shCON_Rep1", "shCON_Rep2", "shETS1_Rep1", "shETS1_Rep2")
pheatmap(as.matrix(heatData), scale = 'row', cluster_rows = F, cluster_cols = F)
heatData_t = t(heatData) %>% as.data.frame()
heatData_t = heatData_t[,c('ETS1',geneList2[geneList2 %in% colnames(heatData_t)])]
FigS7B = pheatmap(as.matrix(heatData_t), scale = 'column', cluster_rows = F, cluster_cols = F)

pdf('plot/FigS7B/FigS7B.pdf', height = 3, width = 7)
FigS7B
dev.off()

Data_FigS7B = heatData_t
write.table(Data_FigS7B, file = 'plot/FigS7B/Data_FigS7B.txt', sep = '\t', quote = F, row.names = T, col.names = T)

# Fig6C. Normalized expression level of the ETS1 gene in different ESCC cell lines from the CCLE database.----
CCLE = read.csv('data_CellLine/19_DepMap/CCLE_expression.21Q2.csv')
CCLE = CCLE[,c(1,grep('ETS1',colnames(CCLE)), drop = F)]
colnames(CCLE) = c('DepMap_ID','ETS1')
CCLE_sampleinfor_NEW = read.csv('data_CellLine/19_DepMap/CCLE_sampleinfor_NEW.csv')
esophagus = filter(CCLE_sampleinfor_NEW, lineage == 'esophagus')
CCLE_ETS1 = merge(esophagus, CCLE, by = 'DepMap_ID')
ESCC = filter(CCLE_ETS1, tcga == 'ESCC') %>% 
  arrange(desc(ETS1)) %>% 
  mutate(cell_line_name = factor(cell_line_name, levels = (cell_line_name)))
ESCC$text_color = 'black'
ESCC[ESCC$cell_line_name == 'KYSE-150',]$text_color = 'red'
ESCC[ESCC$cell_line_name == 'T.T',]$text_color = 'blue'
text_color = c('blue',rep('black',(length(ESCC$cell_line_name)-2)),'red')
Fig6C = ggplot(ESCC, aes(x = cell_line_name, y = ETS1))+
  geom_bar(fill = '#e2c3c8', stat = 'identity')+
  theme_classic()+
  guides(fill = 'none')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   colour=rev(text_color)),
        plot.title = element_text(hjust = 0.5))+
  xlab('')+
  ylab('TPM of ETS1')+
  ggtitle('ESCC cell lines in DepMap')+
  theme(axis.text = element_text(color = 'black'))

pdf('plot/Fig6C/Fig6C.pdf', width = 7, height = 5)
Fig6C
dev.off()

Data_Fig6C = ESCC
write.table(Data_Fig6C, file = 'plot/Fig6C/Data_Fig6C.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig6D/FigS7C. IGV snapshots display the genomic environment of representative ETS1 target genes----

# The ChIP-seq analysis process is the same as in other cell lines. However, if it is not a histone modification ChIP-seq, 
# do not use the --nomodel --extsize 146 parameters in the macs2 callpeak step.

# 1) ChIP-seq analysis (ETS1 in KYSE150)
list.files('data_CellLine/17_HNSC/ETS1_ChIPSeq_KYSE150/')
# [1] "00.fastq"       "01.trim_galore" "02.mapping"     "03.MACS2"

# 2) ChIP-seq analysis (H3K27ac in KYSE150)
list.files('data_CellLine/17_HNSC/H3K27ac_ChIPSeq_KYSE150/')
# [1] "00.fastq"       "01.trim_galore" "02.mapping"     "03.MACS2"
# Fig6E. Relative expression of ETS1 target genes detected by qRT-PCR after ETS1 knockdown in KYSE150----
sh_data = read.table('data_CellLine/20_qPCR/plotData_shETS1_shScramble.txt', sep = '\t', header = T)
pointplotData = group_by(sh_data, Treatment, Gene, Bio_rep) %>% summarise(RelativeExp = mean(RelativeExp), .groups = 'drop')
barplotData = group_by(pointplotData, Gene, Treatment) %>% 
  summarise(mean_exp = mean(RelativeExp), 
            sd_exp = sd(RelativeExp),
            .groups = 'drop')
gene_levels = c("ETS1", "MNX1", "EPCAM", "BARX2", "OCLN", "SYK", "OTULINL")
pointplotData$Gene = factor(pointplotData$Gene, levels = gene_levels)
barplotData$Gene = factor(barplotData$Gene, levels = gene_levels)
pointplotData$Treatment = factor(pointplotData$Treatment, levels = c('shScramble','shETS1'))
barplotData$Treatment = factor(barplotData$Treatment, levels = c('shScramble','shETS1'))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
Fig6E = ggplot(barplotData, aes(x = Gene, y = mean_exp, group = Treatment))+
  geom_bar(stat = 'identity', position = 'dodge', mapping = aes(fill = Treatment))+
  geom_errorbar(aes(ymin = mean_exp, ymax = mean_exp + sd_exp), position = position_dodge(0.9), width = 0.2)+
  scale_fill_manual(values = col_vector)+
  theme_classic()+
  geom_jitter(data = pointplotData, mapping = aes(x = Gene, y = RelativeExp, group = Treatment, fill = Treatment),
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), shape = 18, size = 2, color = '#695860'
  )+
  xlab('')+
  ylab('Relative expression')+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text = element_text(color = 'black'),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('Knockdown of ETS1 in KYSE150')+
  stat_compare_means(mapping = aes(x = Gene, y = RelativeExp, group = Treatment), data = pointplotData,
                     label = "p.signif", method = 't.test', method.args = list(var.equal=T), vjust = -0.5)

pdf('plot/Fig6E/Fig6E.pdf', width = 7, height = 5)
Fig6E
dev.off()

Data_Fig6E = sh_data
write.table(Data_Fig6E, file = 'plot/Fig6E/Data_Fig6E.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Fig6F. Relative expression of ETS1 target genes detected by qRT-PCR after ETS1 overexpression in T.T----
OE_data = read.table('data_CellLine/20_qPCR/plotData_OEETS1_vector.txt', sep = '\t', header = T)
pointplotData = group_by(OE_data, Treatment, Gene, Bio_rep) %>% summarise(RelativeExp = mean(RelativeExp), .groups = 'drop')
barplotData = group_by(pointplotData, Gene, Treatment) %>% 
  summarise(mean_exp = mean(RelativeExp), 
            sd_exp = sd(RelativeExp),
            .groups = 'drop')
gene_levels = c("ETS1", "MNX1", "EPCAM", "BARX2", "OCLN", "SYK", "OTULINL")
pointplotData$Gene = factor(pointplotData$Gene, levels = gene_levels)
barplotData$Gene = factor(barplotData$Gene, levels = gene_levels)
pointplotData$Treatment = sub('ETS1-OE', 'ETS1\noverexpression', sub('Control', 'Control\nvector', pointplotData$Treatment))
barplotData$Treatment = sub('ETS1-OE', 'ETS1\noverexpression', sub('Control', 'Control\nvector', barplotData$Treatment))
pointplotData$Treatment = factor(pointplotData$Treatment, levels = c('Control\nvector','ETS1\noverexpression'))
barplotData$Treatment = factor(barplotData$Treatment, levels = c('Control\nvector','ETS1\noverexpression'))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
Fig6F = ggplot(barplotData, aes(x = Gene, y = mean_exp, group = Treatment))+
  geom_bar(stat = 'identity', position = 'dodge', mapping = aes(fill = Treatment))+
  geom_errorbar(aes(ymin = mean_exp, ymax = mean_exp + sd_exp), position = position_dodge(0.9), width = 0.2)+
  scale_fill_manual(values = col_vector)+
  theme_classic()+
  geom_jitter(data = pointplotData, mapping = aes(x = Gene, y = RelativeExp, group = Treatment, fill = Treatment),
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), shape = 18, size = 2, color = '#695860'
  )+
  xlab('')+
  ylab('Relative expression')+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text = element_text(color = 'black'),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('Knockdown of ETS1 in KYSE150')+
  stat_compare_means(mapping = aes(x = Gene, y = RelativeExp, group = Treatment), data = pointplotData,
                     label = "p.signif", method = 't.test', method.args = list(var.equal=T), vjust = -0.5)

pdf('plot/Fig6F/Fig6F.pdf', width = 7, height = 5)
Fig6F
dev.off()

Data_Fig6F = OE_data
write.table(Data_Fig6F, file = 'plot/Fig6F/Data_Fig6F.txt', sep = '\t', quote = F, row.names = F, col.names = T)

# Combine the plots in Figure 6----
Figure_6 = ggarrange(Fig6A, ggarrange(Fig6B,Fig6C, nrow = 2, labels = c('B','C'), font.label = list(size = 23)), 
                     nrow = 1, widths = c(1,1),labels = c('A'), font.label = list(size = 23))
pdf('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/Figure_6.pdf', width = 13, height = 8)
Figure_6
dev.off()

# FigS8. Overlap between H3K27me3 and H3K9me3 in Roadmap samples----
load('data_Roadmap/RDatas/LOCK_peak_109.RData')
load('data_Roadmap/RDatas/LOCK_peak_109_H3K9me3.RData')
plotData = list()
for (i in 1:109) {
  K27_peak = LOCK_peak_109[[i]][,1:4]
  K9_peak = LOCK_peak_109_H3K9me3[[i]][,1:4]
  both = sum(bed_intersect(K27_peak, K9_peak)$.overlap)
  K27_all = sum(K27_peak$d)
  K9_all = sum(K9_peak$d)
  K27_only = K27_all - both
  K9_only = K9_all - both
  none = 2.7e9 - (K27_only + K9_only + both)
  df = data.frame(group = c('H3K27me3 only', 'both', 'H3K9me3 only', 'none'),
                  size_bp = c(K27_only, both, K9_only, none))
  df$percent = df$size_bp/2.7e9*100
  plotData[[i]] = df %>% mutate(EID = names(peak_109)[i])
}
plotData = do.call(rbind, plotData)
plotData = filter(plotData, group != 'none')
plotData$group = factor(plotData$group, levels = rev(c('H3K27me3 only', 'both', 'H3K9me3 only')))
Figure_S8 = ggplot(plotData, aes(x = EID, y = percent, fill = group))+
  geom_bar(stat = 'identity', position = 'stack')+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.text = element_text(colour = 'black'),
        plot.title = element_text(hjust = 0.5))+
  xlab('sample')+
  ylab('% of domain size to genome size')+
  scale_fill_manual(values = c('H3K27me3 only' = '#ff9900',
                               'both' = 'red',
                               'H3K9me3 only' = '#0099ff'))+
  ggtitle('Roadmap samples (n = 109)')

pdf('plot/FigS8/FigS8.pdf', width = 13, height = 9)
Figure_S8
dev.off()

Data_FigS8 = plotData
write.table(Data_FigS8, file = 'plot/FigS8/Data_FigS8.txt', sep = '\t', quote = F, row.names = F, col.names = T)

pdf('/data3/liangyuan/03_Roadmap_LOCK/data/20_output_plots/Figure_S8.pdf', width = 13, height = 9)
Figure_S8
dev.off()
