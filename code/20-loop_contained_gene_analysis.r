
# --------------------------------------------------------------------------------------------------------------
# Loop associated genes analysis with respect to CTCF knockdown
# ---------------------------------------------------------------------------------------------------------
setwd("/Users/fyu/Documents/GitHub/")
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(data.table)

# --------------------------------------------------------------------------------------------------------------
# look into the CTCF knockdown deg enrichment on the loop genes
peak_obj <- as.data.frame(data.table::fread("mecom_var_sankaran/data/LowC_ana/peak_obj_genename.tsv")) #450920
bedpe_file <- read.table("mecom_var_sankaran/data/LowC_ana/mast0035_down/ctcf_fp_eitherloop.bedpe")


if(sum(bedpe_file[, 2] <bedpe_file[, 5]) == nrow(bedpe_file)){
    temp_bed <- bedpe_file[, c(1, 2, 6)]
    colnames(temp_bed) <- c("chr", "start", "end")
    temp_bed_g <- makeGRangesFromDataFrame(temp_bed)
} else {
    stop()
}

peak_obj_g <- peak_obj[, 1:3]
colnames(peak_obj_g) <- c("chr", "start", "end")
peak_obj_g <- makeGRangesFromDataFrame(peak_obj_g)

ol <- findOverlaps(peak_obj_g, temp_bed_g, maxgap=-1L)
overlap_peak <- peak_obj[queryHits(ol), ] # 82180
table(overlap_peak$peak_type) # down genes
#     G    RE 
#  3264 78916 


# determine how many genes in total
peak_gene_num <- peak_obj$peak_gene_name[peak_obj$peak_gene_name != "0"] %>% unique # 18830
lt_hsc_exp <- data.frame(fread("mecom_var_sankaran/data/LowC_ana/GSE125063_shCTCF_RNA_Seq_FPKM.tsv")) #2512
lt_hsc_exp <- lt_hsc_exp[, 1:7] # 24331
# idx2 <- rowSums(lt_hsc_exp[, 2:7]) != 0 # remove duplicates with the wrong with 0 expression all the samples
# lt_hsc_exp <- lt_hsc_exp[idx2, ] # 17206
bdg_gene_name <- lt_hsc_exp$NAME  %>% unique 
loop_genes <- overlap_peak$peak_gene_name[overlap_peak$peak_gene_name %in% bdg_gene_name] %>% unique # 2564 unique 2200

ctcf_deg <- read.table("mecom_var_sankaran/data/LowC_ana/ctcf_knock_deg.txt", sep="\t", header=T)
lt_all_df <- (ctcf_deg$sample_1=="HSC_LacZ" & ctcf_deg$sample_2=="HSC_shCTCF") %>% ctcf_deg[., ] # 861
lt_up_df <- (ctcf_deg$sample_1=="HSC_LacZ" & ctcf_deg$sample_2=="HSC_shCTCF") %>% ctcf_deg[., ] %>% subset(., log2_fold_change>0) # 296
lt_down_df <- (ctcf_deg$sample_1=="HSC_LacZ" & ctcf_deg$sample_2=="HSC_shCTCF") %>% ctcf_deg[., ] %>% subset(., log2_fold_change<0) # 565

# ----------------------------------- for ctcf_up and MECOM down
# intersect(overlap_peak$peak_gene_name, lt_down_df$Gene.Name) %>% length
loop_genes <- overlap_peak$peak_gene_name[overlap_peak$peak_gene_name %in% bdg_gene_name] %>% unique # 2564 unique 1364
mast0035_down <- data.frame(fread("mecom_var_sankaran/data/down0035_deg_info/deg_mast0035_down.txt"))[, 1]
ctcfdeg_loop_num <- intersect(mast0035_down, lt_up_df$Gene.Name) %>% length
notctcfdeg_loop_num <- length(mast0035_down)-ctcfdeg_loop_num
ctcfdeg_notloop_num <- nrow(lt_up_df)-ctcfdeg_loop_num
notctcfdeg_notloop_num <- length(bdg_gene_name)-length(mast0035_down)-ctcfdeg_notloop_num

observed_table <- matrix(c(ctcfdeg_loop_num, notctcfdeg_loop_num, ctcfdeg_notloop_num, notctcfdeg_notloop_num),
                          nrow = 2, ncol = 2, byrow=T)
rownames(observed_table) <- c('mecom down', 'Not mecom down')
colnames(observed_table) <- c('CTCF up', 'Not CTCF up')
observed_table
#                CTCF up Not CTCF up
# mecom down          16         305
# Not mecom down     280       22376
chisq.test(observed_table)
# Pearson's Chi-squared test with Yates' continuity correction

# data:  observed_table
# X-squared = 32.089, df = 1, p-value = 1.473e-08