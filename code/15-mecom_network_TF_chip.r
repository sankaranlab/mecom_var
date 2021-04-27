# ----------------------------------------------------------------------------------------------------------------
# ChIP-seq peak overlapping analysis
# ----------------------------------------------------------------------------------------------------------------
library(regioneR)
library(GenomicRanges)
library(dplyr)
library(data.table)

setwd("/Users/fyu/Documents/GitHub/")

### 
### mecom all, FLI1, GATA2, RUNX1
### 
query_re_all <- toGRanges(data.frame(fread("mecom_var_sankaran/data/down0035_deg_info/all_mast0035_interactions_peak4.bed")))
# TF
chip_peak_FLI1=toGRanges(data.frame(fread("mecom_var_sankaran/data/CTCF_enrich/FLI1_peaks.narrowPeak")))
chip_peak_GATA2=toGRanges(data.frame(fread("mecom_var_sankaran/data/CTCF_enrich/GATA2_peaks.narrowPeak")))
chip_peak_RUNX1=toGRanges(data.frame(fread("mecom_var_sankaran/data/CTCF_enrich/RUNX1_peaks.narrowPeak")))


FLI1_peakOverlap <- permTest(A=query_re_all, B=chip_peak_FLI1, verbose=F, randomize.function=randomizeRegions, evaluate.function=numOverlaps, ntimes=1000, alternative="greater", genome="hg19", force.parallel=TRUE)
GATA2_peakOverlap <- permTest(A=query_re_all, B=chip_peak_GATA2, verbose=F, randomize.function=randomizeRegions, evaluate.function=numOverlaps, ntimes=1000, alternative="greater", genome="hg19", force.parallel=TRUE)
RUNX1_peakOverlap <- permTest(A=query_re_all, B=chip_peak_RUNX1, verbose=F, randomize.function=randomizeRegions, evaluate.function=numOverlaps, ntimes=1000, alternative="greater", genome="hg19", force.parallel=TRUE)
(FLI1_peakOverlap$numOverlaps$observed)/median(FLI1_peakOverlap$numOverlaps$permuted)
# [1] 35.23077
(GATA2_peakOverlap$numOverlaps$observed)/median(GATA2_peakOverlap$numOverlaps$permuted)
# [1] 29.61111
(RUNX1_peakOverlap$numOverlaps$observed)/median(RUNX1_peakOverlap$numOverlaps$permuted)
# [1] 41.408

# # try to sample the re across the entire hematopoietic catalogs
# # validate the num
# uni_re <- toGRanges(data.frame(fread("mecom_var_sankaran/data/hememap/10142020-atac_peak_pro.txt")))
# FLI1_peakOverlap <- permTest(A=query_re_all, B=chip_peak_FLI1, verbose=F, randomize.function=resampleRegions, universe=uni_re, evaluate.function=numOverlaps, ntimes=1000, alternative="greater", genome="hg19", force.parallel=TRUE)
# (FLI1_peakOverlap$numOverlaps$observed)/median(FLI1_peakOverlap$numOverlaps$permuted)
# # [1] 4.287051
# # validate
# sampled_query_re_all <- resampleRegions(query_re_all, universe=uni_re)
# numOverlaps(query_re_all, chip_peak_FLI1)
# numOverlaps(sampled_query_re_all, chip_peak_FLI1)
# GATA2_peakOverlap <- permTest(A=query_re_all, B=chip_peak_GATA2, verbose=F, randomize.function=resampleRegions, universe=uni_re, evaluate.function=numOverlaps, ntimes=1000, alternative="greater", genome="hg19", force.parallel=TRUE)
# RUNX1_peakOverlap <- permTest(A=query_re_all, B=chip_peak_RUNX1, verbose=F, randomize.function=resampleRegions, universe=uni_re, evaluate.function=numOverlaps, ntimes=1000, alternative="greater", genome="hg19", force.parallel=TRUE)
# (GATA2_peakOverlap$numOverlaps$observed)/median(GATA2_peakOverlap$numOverlaps$permuted)
# # 3.5
# (RUNX1_peakOverlap$numOverlaps$observed)/median(RUNX1_peakOverlap$numOverlaps$permuted)
# # 4.7

#
# all
#
TF=rep(c("FLI1", "GATA2", "RUNX1"), times=2)
count_sd=c(rep(0, 3), c(sd(FLI1_peakOverlap$numOverlaps$permuted), sd(GATA2_peakOverlap$numOverlaps$permuted), sd(RUNX1_peakOverlap$numOverlaps$permuted)))
group=rep(c("Observed", "Expected"), each=3)
overlap=c(
          c(FLI1_peakOverlap$numOverlaps$observed, GATA2_peakOverlap$numOverlaps$observed, RUNX1_peakOverlap$numOverlaps$observed), 
          c(median(FLI1_peakOverlap$numOverlaps$permuted), median(GATA2_peakOverlap$numOverlaps$permuted), median(RUNX1_peakOverlap$numOverlaps$permuted))
          )
overlap
df2 <- data.frame(TF, group, overlap, count_sd)
df2$count_sd[df2$overlap == 0] <- 0
library(ggpubr)
p<- ggplot(df2, aes(x=reorder(TF, overlap), y=overlap, fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=overlap-count_sd, ymax=overlap+count_sd), width=.2,
                 position=position_dodge(.9))+coord_flip() +
                 xlab("") +
                 theme_pubr()+
                 scale_fill_manual("", values = c("Expected" = "gray", "Observed" = "steelblue"))
p


### 
### mecom down CTCF
### 
query_re_down <- toGRanges(data.frame(fread("mecom_var_sankaran/data/down0035_deg_info/down_mast0035_interactions_peak4.bed")))
# CTCF from blood paper
chip_peak_CTCF_hsc=toGRanges(data.frame(fread("mecom_var_sankaran/data/CTCF_enrich/CD34_D0.bed")))
CTCF_hsc_peakOverlap <- permTest(A=query_re_down, B=chip_peak_CTCF_hsc, verbose=F, randomize.function=randomizeRegions, evaluate.function=numOverlaps, ntimes=1000, alternative="greater", genome="hg19", force.parallel=TRUE)
# CTCF_hsc_peakOverlap <- permTest(A=query_re_down, B=chip_peak_CTCF_hsc, verbose=F, randomize.function=resampleRegions, universe=uni_re, evaluate.function=numOverlaps, ntimes=1000, alternative="greater", genome="hg19", force.parallel=TRUE)
# (CTCF_hsc_peakOverlap$numOverlaps$observed)/median(CTCF_hsc_peakOverlap$numOverlaps$permuted)

# 

#
# down
#
TF=rep(c("CTCF_hsc"), times=2)
count_sd=c(rep(0, 2), c(sd(ERG_peakOverlap$numOverlaps$permuted), sd(FLI1_peakOverlap$numOverlaps$permuted), sd(SCL_peakOverlap$numOverlaps$permuted), sd(LYL1_peakOverlap$numOverlaps$permuted), sd(GATA2_peakOverlap$numOverlaps$permuted), sd(RUNX1_peakOverlap$numOverlaps$permuted), sd(LMO2_peakOverlap$numOverlaps$permuted), sd(CTCF_hsc_peakOverlap$numOverlaps$permuted), sd(CTCF_ery_peakOverlap$numOverlaps$permuted), sd(CTCF_b_peakOverlap$numOverlaps$permuted), sd(CTCF_t_peakOverlap$numOverlaps$permuted), sd(CTCF_mono_peakOverlap$numOverlaps$permuted)))
group=rep(c("Observed", "Expected"), each=1)
overlap=c(
          c(CTCF_hsc_peakOverlap$numOverlaps$observed), 
          c(median(CTCF_hsc_peakOverlap$numOverlaps$permuted))
          )
overlap
df2 <- data.frame(TF, group, overlap, count_sd)
df2$count_sd[df2$overlap == 0] <- 0
library(ggpubr)
p <- ggplot(df2, aes(x=reorder(TF, overlap), y=overlap, fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=overlap-count_sd, ymax=overlap+count_sd), width=.2,
                 position=position_dodge(.9))+coord_flip() +
                 xlab("") +
                 theme_pubr()+
                 scale_fill_manual("", values = c("Expected" = "gray", "Observed" = "steelblue"))
p




# ----------------------------------------------------------------------------------------------------------------
# chip-seq peak overlap analysis for three tfs
# ----------------------------------------------------------------------------------------------------------------
setwd("mecom_var_sankaran/data/CTCF_enrich")
a_file <- read.table("FLI1_peaks.narrowPeak")[, 1:3]
b_file <- read.table("GATA2_peaks.narrowPeak")[, 1:3]
c_file <- read.table("RUNX1_peaks.narrowPeak")[, 1:3]
colnames(a_file) <- colnames(b_file) <- colnames(c_file) <- c("chr", "start", "end")
a_file_g <- makeGRangesFromDataFrame(a_file)
b_file_g <- makeGRangesFromDataFrame(b_file)
c_file_g <- makeGRangesFromDataFrame(c_file)
###
### all deg
### 
setwd("/Users/fyu/Documents/GitHub/")
mast0035_all_file <- read.table("mecom_var_sankaran/data/down0035_deg_info/all_mast0035_interactions_peak4.bed")[, 1:3]
colnames(mast0035_all_file) <- c("chr", "start", "end")
mast0035_all_file_g <- makeGRangesFromDataFrame(mast0035_all_file)
ol_1 <- findOverlaps(mast0035_all_file_g, a_file_g, select="first") 
ol_2 <- findOverlaps(mast0035_all_file_g, b_file_g, select="first") 
ol_3 <- findOverlaps(mast0035_all_file_g, c_file_g, select="first") 
ol_mat <- data.frame(ol_1, ol_2, ol_3)
ol_mat[is.na(ol_mat)] <- 0
ol_mat[ol_mat>0] <- 1
a_file_num <- (1:nrow(mast0035_all_file))[as.logical(ol_mat[, 1])]
b_file_num <- (1:nrow(mast0035_all_file))[as.logical(ol_mat[, 2])]
c_file_num <- (1:nrow(mast0035_all_file))[as.logical(ol_mat[, 3])]

write.table(data.frame(a_file_num), "mecom_var_sankaran/data/chip_overlapping/TF3_overlap_FLI1_mast0035_all.txt", row.names=F, col.names=F, sep="\t", quote=F)
write.table(data.frame(b_file_num), "mecom_var_sankaran/data/chip_overlapping/TF3_overlap_GATA2_mast0035_all.txt", row.names=F, col.names=F, sep="\t", quote=F)
write.table(data.frame(c_file_num), "mecom_var_sankaran/data/chip_overlapping/TF3_overlap_RUNX1_mast0035_all.txt", row.names=F, col.names=F, sep="\t", quote=F)
