
# ----------------------------------------------------------------------------------------------------------------
# CTCF chip-seq and footprint analysis 
# ----------------------------------------------------------------------------------------------------------------
setwd("/Users/fyu/Documents/GitHub/")
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(data.table)

##
## mast0035_down
##
motif_data <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_down/CTCF_fimo.bed")
colnames(motif_data) <- c("chr", "start", "end", "misc", "score", "strand")
gene_data <- read.table("mecom_var_sankaran/data/down0035_deg_info/01222021-mast0035_down_interaction_gene.bed")
colnames(gene_data) <- c("chr", "start", "end", "type", "peak_rank", "gene_name", "peak_gene_name", "interaction_strength", "interaction_p", "interaction_type")
re_data <- read.table("mecom_var_sankaran/data/down0035_deg_info/01222021-mast0035_down_interaction_re.bed")
colnames(re_data) <- c("chr", "start", "end", "type", "peak_rank", "gene_name", "peak_gene_name", "interaction_strength", "interaction_p", "interaction_type")
factor_cp_hsc <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_down/hsc/hsc_CTCF_fimo.pdf.postpr.txt")[, 1]
gene_peak_obj_g <- makeGRangesFromDataFrame(gene_data, keep.extra.columns = TRUE)
gene_peak_obj_re <- makeGRangesFromDataFrame(re_data, keep.extra.columns = TRUE)
motif_obj_g <- makeGRangesFromDataFrame(motif_data, keep.extra.columns = TRUE)
footprint_obj_g <- makeGRangesFromDataFrame(motif_data[factor_cp_hsc>=0.95, ], keep.extra.columns = TRUE) # 1892/5053
ol_motif_g <- findOverlaps(gene_peak_obj_g, motif_obj_g) # length(ol) 
ol_motif_re <- findOverlaps(gene_peak_obj_re, motif_obj_g) # length(ol) 
ol_fp_g <- findOverlaps(gene_peak_obj_g, footprint_obj_g) # length(ol) 
ol_fp_re <- findOverlaps(gene_peak_obj_re, footprint_obj_g) # length(ol) 

length(table(queryHits(ol_motif_g))) # 4122
length(table(queryHits(ol_motif_re))) # 4001
length(table(queryHits(ol_fp_g))) # 1366
length(table(queryHits(ol_fp_re))) # 2040

length(unique(subjectHits(ol_fp_re)))
# [1] 2292
length(unique(subjectHits(ol_fp_g)))
# [1] 57
motif_data2 <- data.frame(motif_data, factor_cp_hsc)
fp_overlap_re <- data.frame(motif_data2[motif_data2$factor_cp_hsc>=0.95, ][subjectHits(ol_fp_re), ], re_data[queryHits(ol_fp_re), ])
fp_overlap_re_sort <- fp_overlap_re[order(fp_overlap_re$factor_cp_hsc, fp_overlap_re$score, decreasing=T), ]
setwd("mecom_var_sankaran/data/CTCF_enrich")
write.table(fp_overlap_re_sort, "mecom_var_sankaran/data/CTCF_enrich/ctcf_fp_re_sort_bed_info.txt", row.names=F, col.names=F, sep="\t", quote=F)


#### NOT in R
#### the bw files are too large to upload, download them to the bwdir first by yourself
cd mecom_var_sankaran/data/CTCF_enrich
bwdir=mecom_var_sankaran/data/CTCF_enrich
outdir=mecom_var_sankaran/data/CTCF_enrich/deeptools_ana
mkdir -p $outdir
# make sure you have deepTools installed
computeMatrix reference-point \
    -S $bwdir/CD34_D0_CTCF.bw  \
       $bwdir/Ery_cell.bw  \
       $bwdir/T_cell.bw  \
       $bwdir/B_cell.bw  \
       $bwdir/mono_cell.bw \
    -R mecom_var_sankaran/data/CTCF_enrich/ctcf_fp_re_sort_bed_info.txt \
    --referencePoint center \
    -a 1000 -b 1000 --outFileName $outdir/fp_signal5_ctcf_mast0035_down.mat.gz

# RdBu
# YlGnBu
# gist_heat
# Reds
plotHeatmap \
    -m $outdir/fp_signal5_ctcf_mast0035_down.mat.gz \
    -out $outdir/fp_signal5_ctcf_mast0035_down_OrRd.pdf \
    --outFileNameMatrix $outdir/heatmap_value.mat \
    --heatmapHeight 15  \
    --refPointLabel footprint.center \
    --regionsLabel signal \
    --plotTitle 'CTCF signal' \
    --sortRegions no \
    --boxAroundHeatmaps no \
    --colorMap OrRd \
    --xAxisLabel "" \
    --refPointLabel "footprint"


