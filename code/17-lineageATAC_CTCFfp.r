# ----------------------------------------------------------------------------------------------------------------
# ATAC signal on the CTCF foorprint (mecom down)
# ----------------------------------------------------------------------------------------------------------------
library(annotables)
library(tidyr)
library(GenomicRanges)
library(stringr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)
setwd("/Users/fyu/Documents/GitHub/")

#### 
#### get REs contain and did not contain CTCF fp 
#### 
cisre <- read.table("mecom_var_sankaran/data/down0035_deg_info/01222021-mast0035_down_interaction_re.bed")
temp_occur_mast0035_down_binary_fp <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_down/temp_occur_mast0035_down_binary_fp.txt", header=T, sep="\t")
cisre <- data.frame(cisre, temp_occur_mast0035_down_binary_fp$CTCF.re_fp)
cisre_ctcf <- cisre[cisre[, 11]==1, ]
cisre_noctcf <- cisre[cisre[, 11]==0, ]
peak_sample_mat <- as.data.frame(data.table::fread("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/09022020-MECOM_Richard/data/02262021-downMast_re-analysis/supplement/ctcf_atac_analysis//10142020-peak_sample_mat_pro.txt", header=T))
peak_sample_mat <- edgeR::cpm(peak_sample_mat) %>% as.data.frame
ctcf_idx <- str_split(cisre_ctcf$V5, "_", simplify = T)[, 2] %>% as.numeric
noctcf_idx <- str_split(cisre_noctcf$V5, "_", simplify = T)[, 2] %>% as.numeric

ctcf_mat <- peak_sample_mat[ctcf_idx, ]
noctcf_mat <- peak_sample_mat[noctcf_idx, ]

ctcf_mat_melt <- melt(ctcf_mat)
noctcf_mat_melt <- melt(noctcf_mat)
ctcf_mat_melt <- data.frame(ctcf_mat_melt, group="CTCF")
noctcf_mat_melt <- data.frame(noctcf_mat_melt, group="no CTCF")
mydata_combine <- rbind.data.frame(ctcf_mat_melt, noctcf_mat_melt)

# ggplot(data = mydata_combine, aes(x = variable, y = value, fill = group))+
#      scale_fill_viridis_d()+
#      # geom_violin(alpha=0.4, position = position_dodge(width = .75), size=1, color="black") +
#     #  geom_point( shape = 16, size=0.2, position = position_jitterdodge(), color="gray", alpha=0.2)+
#           geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=.2, alpha = 0.8)+ ylim(0, 30)+
#         pretty_plot(fontsize = 6) + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
# ery lineage
# boxplot(ctcf_mat[, c("HSC", "MPP", "CMP", "MEP", "Ery")] %>% log1p, outline=FALSE, col=as.character(jdb_palette("brewer_spectra")), main="Ery lineage")
boxplot(ctcf_mat[, c("HSC", "MPP", "CMP", "MEP", "Ery")], outline=FALSE, col=as.character(jdb_palette("brewer_spectra")), main="Erythroid lineage", ylab="Chromatin accessiblity")

# Lymph lineage
# boxplot(ctcf_mat[, c("HSC", "MPP", "CLP", "CD8", "CD4", "B")] %>% log1p, outline=FALSE, col=as.character(jdb_palette("brewer_spectra")), main="Lymph lineage")
boxplot(ctcf_mat[, c("HSC", "MPP", "CLP", "CD8", "CD4", "B")], outline=FALSE, col=as.character(jdb_palette("brewer_spectra")), main="Lymphoid lineage", ylab="Chromatin accessiblity")

# Mye lineage
# boxplot(ctcf_mat[, c("HSC", "MPP", "CMP", "GMP.B", "Mono")] %>% log1p, outline=FALSE, col=as.character(jdb_palette("brewer_spectra")), main="Myeloid lineage")
boxplot(ctcf_mat[, c("HSC", "MPP", "CMP", "GMP.B", "Mono")], outline=FALSE, col=as.character(jdb_palette("brewer_spectra")), main="Myeloid lineage", ylab="Chromatin accessiblity")

