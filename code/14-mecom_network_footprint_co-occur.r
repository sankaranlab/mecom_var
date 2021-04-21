
# --------------------------------------------------------------------------------------------------------------
# footprinting coocurrence & chromatin accessibility on CTCF fp across different lineages
# ---------------------------------------------------------------------------------------------------------
library(annotables)
library(tidyr)
library(GenomicRanges)
library(stringr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)
setwd("/Users/fyu/Documents/GitHub/")

# define functions
# return the occurrence of fp and motif in cisREs, also output the differences and graph of HemeMap scores between fp-containing cisREs and those do not
get_motif_and_fp <- function(motif_data, gene_data, re_data, factor_cp_hsc){
    motif_data <- read.table(motif_data)
    colnames(motif_data) <- c("chr", "start", "end", "misc", "score", "strand")
    gene_data <- read.table(gene_data)
    colnames(gene_data) <- c("chr", "start", "end", "type", "peak_rank", "gene_name", "peak_gene_name", "interaction_strength", "interaction_p", "interaction_type")
    re_data <- read.table(re_data)
    colnames(re_data) <- c("chr", "start", "end", "type", "peak_rank", "gene_name", "peak_gene_name", "interaction_strength", "interaction_p", "interaction_type")
    factor_cp_hsc <- read.table(factor_cp_hsc)[, 1]
    
    gene_peak_obj_g <- makeGRangesFromDataFrame(gene_data, keep.extra.columns = TRUE)
    gene_peak_obj_re <- makeGRangesFromDataFrame(re_data, keep.extra.columns = TRUE)
    motif_obj_g <- makeGRangesFromDataFrame(motif_data, keep.extra.columns = TRUE)
    footprint_obj_g <- makeGRangesFromDataFrame(motif_data[factor_cp_hsc>=0.95, ], keep.extra.columns = TRUE)
    ol_motif_g <- findOverlaps(gene_peak_obj_g, motif_obj_g) # length(ol) 
    ol_motif_re <- findOverlaps(gene_peak_obj_re, motif_obj_g) # length(ol) 
    ol_fp_g <- findOverlaps(gene_peak_obj_g, footprint_obj_g) # length(ol) 
    ol_fp_re <- findOverlaps(gene_peak_obj_re, footprint_obj_g) # length(ol) 

    length(table(queryHits(ol_motif_g))) # 4228
    length(table(queryHits(ol_motif_re))) # 3786
    length(table(queryHits(ol_fp_g))) # 1322
    length(table(queryHits(ol_fp_re))) # 1741

    ol_re_motifnum <- rep(0, length(gene_peak_obj_re))
    ol_re_motifnum[as.numeric(names(table(queryHits(ol_motif_re))))] <- table(queryHits(ol_motif_re))

    ol_re_fpnum <- rep(0, length(gene_peak_obj_re))
    ol_re_fpnum[as.numeric(names(table(queryHits(ol_fp_re))))] <- table(queryHits(ol_fp_re))
    boxplot(log1p(re_data$interaction_strength[ol_re_fpnum!=0]), log1p(re_data$interaction_strength[ol_re_fpnum==0]), col=c("indianred", "lightblue"), 
            main=paste0("Wilcox.test P-value=", wilcox.test(re_data$interaction_strength[ol_re_fpnum!=0], re_data$interaction_strength[ol_re_fpnum==0])$p.value), 
            names=c("contain fp", "not contain fp"), 
            outline=F,
            ylab="log2(interaction strength)",
            boxwex=0.7)
    print(t.test(re_data$interaction_strength[ol_re_motifnum==0], re_data$interaction_strength[ol_re_motifnum!=0]))
    print(t.test(re_data$interaction_strength[ol_re_fpnum==0], re_data$interaction_strength[ol_re_fpnum!=0]))
    print(wilcox.test(re_data$interaction_strength[ol_re_fpnum==0], re_data$interaction_strength[ol_re_fpnum!=0]))
    return(data.frame(ol_re_motifnum, ol_re_fpnum))
}


### 
### mecom all
### 
setwd("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_all")
# TF_names="ETS CTCF RUNX JUN KLF GATA"
# ETS
ETS_motif_fp_occur <- get_motif_and_fp("./ETS_fimo.bed",
                    "../../../down0035_deg_info/01222021-mast0035_all_interaction_gene.bed",
                    "../../../down0035_deg_info/01222021-mast0035_all_interaction_re.bed",
                    "./hsc/hsc_ETS_fimo.pdf.postpr.txt")
colnames(ETS_motif_fp_occur) <- c("ETS-re_motif", "ETS-re_fp")
# ETS_fp_strength_mast0035_all
# CTCF
CTCF_motif_fp_occur <- get_motif_and_fp("./CTCF_fimo.bed",
                    "../../../down0035_deg_info/01222021-mast0035_all_interaction_gene.bed",
                    "../../../down0035_deg_info/01222021-mast0035_all_interaction_re.bed",
                    "./hsc/hsc_CTCF_fimo.pdf.postpr.txt")
colnames(CTCF_motif_fp_occur) <- c("CTCF-re_motif", "CTCF-re_fp")
# CTCF_fp_strength_mast0035_all

# RUNX
RUNX_motif_fp_occur <- get_motif_and_fp("./RUNX_fimo.bed",
                    "../../../down0035_deg_info/01222021-mast0035_all_interaction_gene.bed",
                    "../../../down0035_deg_info/01222021-mast0035_all_interaction_re.bed",
                    "./hsc/hsc_RUNX_fimo.pdf.postpr.txt")
colnames(RUNX_motif_fp_occur) <- c("RUNX-re_motif", "RUNX-re_fp")
# JUN
JUN_motif_fp_occur <- get_motif_and_fp("./JUN_fimo.bed",
                    "../../../down0035_deg_info/01222021-mast0035_all_interaction_gene.bed",
                    "../../../down0035_deg_info/01222021-mast0035_all_interaction_re.bed",
                    "./hsc/hsc_JUN_fimo.pdf.postpr.txt")
colnames(JUN_motif_fp_occur) <- c("JUN-re_motif", "JUN-re_fp")
# KLF
KLF_motif_fp_occur <- get_motif_and_fp("./KLF_fimo.bed",
                    "../../../down0035_deg_info/01222021-mast0035_all_interaction_gene.bed",
                    "../../../down0035_deg_info/01222021-mast0035_all_interaction_re.bed",
                    "./hsc/hsc_KLF_fimo.pdf.postpr.txt")
colnames(KLF_motif_fp_occur) <- c("KLF-re_motif", "KLF-re_fp")
# GATA
GATA_motif_fp_occur <- get_motif_and_fp("./GATA_fimo.bed",
                    "../../../down0035_deg_info/01222021-mast0035_all_interaction_gene.bed",
                    "../../../down0035_deg_info/01222021-mast0035_all_interaction_re.bed",
                    "./hsc/hsc_GATA_fimo.pdf.postpr.txt")
colnames(GATA_motif_fp_occur) <- c("GATA-re_motif", "GATA-re_fp")
temp_occur_mast0035_all <- data.frame(ETS_motif_fp_occur, CTCF_motif_fp_occur, RUNX_motif_fp_occur, JUN_motif_fp_occur, KLF_motif_fp_occur, GATA_motif_fp_occur)
temp_occur_mast0035_all_binary <- temp_occur_mast0035_all; temp_occur_mast0035_all_binary[temp_occur_mast0035_all_binary>1] <- 1
temp_occur_mast0035_all_binary_motif <- temp_occur_mast0035_all_binary[, c(1, 3, 5, 7, 9, 11)]
temp_occur_mast0035_all_binary_fp <- temp_occur_mast0035_all_binary[, c(2, 4, 6, 8, 10, 12)]

out_cor_p_mat <- matrix(0, nrow=ncol(temp_occur_mast0035_all_binary_fp), ncol=ncol(temp_occur_mast0035_all_binary_fp))
for (i in 1:ncol(temp_occur_mast0035_all_binary_fp)){
    for (j in  1:ncol(temp_occur_mast0035_all_binary_fp)){
            if(i ==j){
                out_cor_p_mat[i, j] <- 1
            } else {
                p_value <- phyper(sum(temp_occur_mast0035_all_binary_fp[, i]== 1 & temp_occur_mast0035_all_binary_fp[, j]==1)-1,
                        colSums(temp_occur_mast0035_all_binary_fp)[i], 
                        nrow(temp_occur_mast0035_all_binary_fp)-colSums(temp_occur_mast0035_all_binary_fp)[i], 
                        colSums(temp_occur_mast0035_all_binary_fp)[j], lower.tail = F)

                out_cor_p_mat[i, j] <- as.numeric(formatC(p_value, format = "e", digits = 2))
            }
            
    }
}

TF_names=c("ETS", "CTCF", "RUNX", "JUN", "KLF", "GATA")
colnames(out_cor_p_mat) <- rownames(out_cor_p_mat) <- TF_names
### fp_coocurrence_upper_heatmap
library(corrplot)
library(RColorBrewer)
corrplot(-log10(out_cor_p_mat), is.corr=F, type = "upper", col = brewer.pal(n = 10, name = "RdYlBu"),  order = "FPC")



### 
### mecom down
### 
# focus on these 6 TFs
setwd("/Users/fyu/Documents/GitHub/")
setwd("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_down")
# TF_names="ETS CTCF RUNX JUN KLF"
ETS_motif_fp_occur <- get_motif_and_fp("./ETS_fimo.bed",
                    "../../../down0035_deg_info/01222021-mast0035_down_interaction_gene.bed",
                    "../../../down0035_deg_info/01222021-mast0035_down_interaction_re.bed",
                    "./hsc/hsc_ETS_fimo.pdf.postpr.txt")
colnames(ETS_motif_fp_occur) <- c("ETS-re_motif", "ETS-re_fp")
# ETS_fp_strength_mast0035_down
# CTCF
CTCF_motif_fp_occur <- get_motif_and_fp("./CTCF_fimo.bed",
                    "../../../down0035_deg_info/01222021-mast0035_down_interaction_gene.bed",
                    "../../../down0035_deg_info/01222021-mast0035_down_interaction_re.bed",
                    "./hsc/hsc_CTCF_fimo.pdf.postpr.txt")
colnames(CTCF_motif_fp_occur) <- c("CTCF-re_motif", "CTCF-re_fp")
# CTCF_fp_strength_mast0035_down
# RUNX
RUNX_motif_fp_occur <- get_motif_and_fp("./RUNX_fimo.bed",
                    "../../../down0035_deg_info/01222021-mast0035_down_interaction_gene.bed",
                    "../../../down0035_deg_info/01222021-mast0035_down_interaction_re.bed",
                    "./hsc/hsc_RUNX_fimo.pdf.postpr.txt")
colnames(RUNX_motif_fp_occur) <- c("RUNX-re_motif", "RUNX-re_fp")
# JUN
JUN_motif_fp_occur <- get_motif_and_fp("./JUN_fimo.bed",
                    "../../../down0035_deg_info/01222021-mast0035_down_interaction_gene.bed",
                    "../../../down0035_deg_info/01222021-mast0035_down_interaction_re.bed",
                    "./hsc/hsc_JUN_fimo.pdf.postpr.txt")
colnames(JUN_motif_fp_occur) <- c("JUN-re_motif", "JUN-re_fp")
# KLF
KLF_motif_fp_occur <- get_motif_and_fp("./KLF_fimo.bed",
                    "../../../down0035_deg_info/01222021-mast0035_down_interaction_gene.bed",
                    "../../../down0035_deg_info/01222021-mast0035_down_interaction_re.bed",
                    "./hsc/hsc_KLF_fimo.pdf.postpr.txt")
colnames(KLF_motif_fp_occur) <- c("KLF-re_motif", "KLF-re_fp")

temp_occur_mast0035_down <- data.frame(ETS_motif_fp_occur, CTCF_motif_fp_occur, RUNX_motif_fp_occur, JUN_motif_fp_occur, KLF_motif_fp_occur)
temp_occur_mast0035_down_binary <- temp_occur_mast0035_down; temp_occur_mast0035_down_binary[temp_occur_mast0035_down_binary>1] <- 1
temp_occur_mast0035_down_binary_motif <- temp_occur_mast0035_down_binary[, c(1, 3, 5, 7, 9)]
temp_occur_mast0035_down_binary_fp <- temp_occur_mast0035_down_binary[, c(2, 4, 6, 8, 10)]
write.table(temp_occur_mast0035_down_binary_fp, "temp_occur_mast0035_down_binary_fp.txt", row.names=F, col.names=T, quote=F, sep="\t")

out_cor_p_mat <- matrix(0, nrow=ncol(temp_occur_mast0035_down_binary_fp), ncol=ncol(temp_occur_mast0035_down_binary_fp))
for (i in 1:ncol(temp_occur_mast0035_down_binary_fp)){
    for (j in  1:ncol(temp_occur_mast0035_down_binary_fp)){
            if(i ==j){
                out_cor_p_mat[i, j] <- 1
            } else {
                p_value <- phyper(sum(temp_occur_mast0035_down_binary_fp[, i]== 1 & temp_occur_mast0035_down_binary_fp[, j]==1)-1,
                        colSums(temp_occur_mast0035_down_binary_fp)[i], 
                        nrow(temp_occur_mast0035_down_binary_fp)-colSums(temp_occur_mast0035_down_binary_fp)[i], 
                        colSums(temp_occur_mast0035_down_binary_fp)[j], lower.tail = F)

                out_cor_p_mat[i, j] <- as.numeric(formatC(p_value, format = "e", digits = 2))
            }
            
    }
}

TF_names=c("ETS", "CTCF", "RUNX", "JUN", "KLF")
colnames(out_cor_p_mat) <- rownames(out_cor_p_mat) <- TF_names
### fp_coocurrence_upper_heatmap
library(corrplot)
library(RColorBrewer)
corrplot(-log10(out_cor_p_mat), is.corr=F, type = "upper", col = brewer.pal(n = 10, name = "RdYlBu"),  order = "FPC")



### 
### mecom up
### 
# focus on these 6 TFs
setwd("/Users/fyu/Documents/GitHub/")
setwd("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_up")
# TF_names="ETS CTCF RUNX JUN KLF"
# ETS
ETS_motif_fp_occur <- get_motif_and_fp("./ETS_fimo.bed",
                    "../../../down0035_deg_info/01222021-mast0035_up_interaction_gene.bed",
                    "../../../down0035_deg_info/01222021-mast0035_up_interaction_re.bed",
                    "./hsc/hsc_ETS_fimo.pdf.postpr.txt")
colnames(ETS_motif_fp_occur) <- c("ETS-re_motif", "ETS-re_fp")
# ETS_fp_strength_mast0035_up
# CTCF
CTCF_motif_fp_occur <- get_motif_and_fp("./CTCF_fimo.bed",
                    "../../../down0035_deg_info/01222021-mast0035_up_interaction_gene.bed",
                    "../../../down0035_deg_info/01222021-mast0035_up_interaction_re.bed",
                    "./hsc/hsc_CTCF_fimo.pdf.postpr.txt")
colnames(CTCF_motif_fp_occur) <- c("CTCF-re_motif", "CTCF-re_fp")
# CTCF_fp_strength_mast0035_up
# RUNX
RUNX_motif_fp_occur <- get_motif_and_fp("./RUNX_fimo.bed",
                    "../../../down0035_deg_info/01222021-mast0035_up_interaction_gene.bed",
                    "../../../down0035_deg_info/01222021-mast0035_up_interaction_re.bed",
                    "./hsc/hsc_RUNX_fimo.pdf.postpr.txt")
colnames(RUNX_motif_fp_occur) <- c("RUNX-re_motif", "RUNX-re_fp")
# JUN
JUN_motif_fp_occur <- get_motif_and_fp("./JUN_fimo.bed",
                    "../../../down0035_deg_info/01222021-mast0035_up_interaction_gene.bed",
                    "../../../down0035_deg_info/01222021-mast0035_up_interaction_re.bed",
                    "./hsc/hsc_JUN_fimo.pdf.postpr.txt")
colnames(JUN_motif_fp_occur) <- c("JUN-re_motif", "JUN-re_fp")
# KLF
KLF_motif_fp_occur <- get_motif_and_fp("./KLF_fimo.bed",
                    "../../../down0035_deg_info/01222021-mast0035_up_interaction_gene.bed",
                    "../../../down0035_deg_info/01222021-mast0035_up_interaction_re.bed",
                    "./hsc/hsc_KLF_fimo.pdf.postpr.txt")
colnames(KLF_motif_fp_occur) <- c("KLF-re_motif", "KLF-re_fp")

temp_occur_mast0035_up <- data.frame(ETS_motif_fp_occur, CTCF_motif_fp_occur, RUNX_motif_fp_occur, JUN_motif_fp_occur, KLF_motif_fp_occur)
temp_occur_mast0035_up_binary <- temp_occur_mast0035_up; temp_occur_mast0035_up_binary[temp_occur_mast0035_up_binary>1] <- 1
temp_occur_mast0035_up_binary_motif <- temp_occur_mast0035_up_binary[, c(1, 3, 5, 7, 9)]
temp_occur_mast0035_up_binary_fp <- temp_occur_mast0035_up_binary[, c(2, 4, 6, 8, 10)]

out_cor_p_mat <- matrix(0, nrow=ncol(temp_occur_mast0035_up_binary_fp), ncol=ncol(temp_occur_mast0035_up_binary_fp))
for (i in 1:ncol(temp_occur_mast0035_up_binary_fp)){
    for (j in  1:ncol(temp_occur_mast0035_up_binary_fp)){
            if(i ==j){
                out_cor_p_mat[i, j] <- 1
            } else {
                p_value <- phyper(sum(temp_occur_mast0035_up_binary_fp[, i]== 1 & temp_occur_mast0035_up_binary_fp[, j]==1)-1,
                        colSums(temp_occur_mast0035_up_binary_fp)[i], 
                        nrow(temp_occur_mast0035_up_binary_fp)-colSums(temp_occur_mast0035_up_binary_fp)[i], 
                        colSums(temp_occur_mast0035_up_binary_fp)[j], lower.tail = F)

                out_cor_p_mat[i, j] <- as.numeric(formatC(p_value, format = "e", digits = 2))
            }
            
    }
}

TF_names=c("ETS", "CTCF", "RUNX", "JUN", "KLF")
colnames(out_cor_p_mat) <- rownames(out_cor_p_mat) <- TF_names
### fp_coocurrence_upper_heatmap
library(corrplot)
library(RColorBrewer)
corrplot(-log10(out_cor_p_mat), is.corr=F, type = "upper", col = brewer.pal(n = 10, name = "RdYlBu"),  order = "FPC")

