#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#   HemeMap construction
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
library(annotables)
library(qvalue)
library(tidyselect)
library(tidyr)
library(GenomicRanges)
library(stringr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)
library(doParallel)
library(parallel)
setwd("/Users/fyu/Documents/GitHub/")

#####################################################################################
#   Find the indirect target genes
#####################################################################################
# -----------------------------------
# cisRE <-> cisRE <-> gene
# -----------------------------------

# -----------------------------------
# ATAC-ATAC correlation
# 
search_radius <- 500000 # search RE for each gene within sournding 500kb
# read gene filtered and promoter integrated ATAC peaks
peak_obj <- as.data.frame(data.table::fread("mecom_var_sankaran/data/hememap/10142020-atac_peak_pro.txt")) #450920
peak_obj <- data.frame(peak_obj, paste("peak", seq(1:nrow(peak_obj)), sep="_"))
colnames(peak_obj) <- c("chr", "start", "end", "peak_type", "peak_rank") # no M
peak_sample_mat <- as.data.frame(data.table::fread("mecom_var_sankaran/data/hememap/10142020-peak_sample_mat_pro.txt", header=T))
# peak_sample_mat <- subset(peak_sample_mat, select = -c(mDC, Mega)) # remove two cell types
rownames(peak_sample_mat) <- peak_obj$peak_rank

# normalize count tables for atac - 18 cell types
peak_sample_mat_cpm <- edgeR::cpm(peak_sample_mat)

# find the RE-RE pairs using promoter and atac peaks
# convert to GRange obj
peak_obj_g <- makeGRangesFromDataFrame(peak_obj, keep.extra.columns = TRUE)
ol <- findOverlaps(peak_obj_g, peak_obj_g, maxgap=search_radius) # length(ol) 4114621
rho <- simplify2array(
            parallel::mclapply(1:length(ol), mc.cores = 4, function(x){
                rho <- cor(log2(peak_sample_mat_cpm[queryHits(ol)[x], ]+1), log2(peak_sample_mat_cpm[subjectHits(ol)[x], ]+1))
            }  
        )
    )
# sum(is.na((rho))) # 71335
# Define a function to get pvalue from correlation given a sample size
rhoToP <- function(rho, n = 18){
  t <- abs(rho / sqrt((1-rho)^2/(n-2)))
  2*pt(t, n-1, lower=FALSE)
}
p_value <- rhoToP(rho)
output_df <- data.frame(data.frame(ol), rho, p_value) # 95875260

# rm self correlation and na
output_df2 <- subset(output_df, (queryHits != subjectHits) & (!is.na(rho))) # 95353336

######## Remove duplicates (B-A & A-B ---> A-B); sort A < B
# df is a dataframe
rm_dup <- function(df) {
    cols <- c(1, 2)
    newdf <- df[, cols]
    newdf <- matrix(c(apply(newdf, 1, sort)), ncol=2, byrow=T)
    idx <- !duplicated(newdf)
    return(idx)
}
dup_idx <- rm_dup(output_df2)
output_df2_nodupSort <- output_df2[dup_idx, ]
# already sorted
# inverse_idx <- output_df2_nodupSort$queryHits > output_df2_nodupSort$subjectHits
# sum(inverse_idx)
# [1] 0
candidate_pair1 <- peak_obj[output_df2_nodupSort$queryHits, ]
candidate_pair2 <- peak_obj[output_df2_nodupSort$subjectHits, ]
ATACcor_unique_name <- paste(candidate_pair1$chr, candidate_pair1$start, candidate_pair1$end, candidate_pair2$chr, candidate_pair2$start, candidate_pair2$end, sep="_")
output_df2_nodupSort_allinfo <- data.frame(output_df2_nodupSort, ATACcor_unique_name, peak1_type=candidate_pair1$peak_type, peak2_type=candidate_pair2$peak_type, peak1_rank=candidate_pair1$peak_rank, peak2_rank=candidate_pair2$peak_rank)
gene_idx <- output_df2_nodupSort_allinfo$peak1_type == "G" | output_df2_nodupSort_allinfo$peak2_type == "G"
output_df2_nodupSort_gene <- output_df2_nodupSort_allinfo[gene_idx, ] # 4,176,928
output_df2_nodupSort_nogene <- output_df2_nodupSort_allinfo[!gene_idx, ] # 43,499,740
sum(output_df2_nodupSort_nogene[, 1] > output_df2_nodupSort_nogene[, 2]) # no need to sort again; will be used for merge to build network
# [1] 0
output_df2_nodupSort_nogene_cutoffp005 <- output_df2_nodupSort_nogene[output_df2_nodupSort_nogene$p_value <= 0.05, ] # 15061815 # [1] 0.3453158 0.9940528 # do not need to filter out the negative rho value
quantile(table(output_df2_nodupSort_nogene_cutoffp005$queryHits)) # look at how many elements link to an element
# 0%  25%  50%  75% 100%
# 1   16   28   48  251
length(unique(output_df2_nodupSort_nogene_cutoffp005$queryHits)) # 429276/450920 # # look at how many elements were covered by this analysis

# processing the direct interactions
gene_idx <- output_df2_nodupSort_allinfo$peak1_type == "G" | output_df2_nodupSort_allinfo$peak2_type == "G"
output_df2_nodupSort_gene <- output_df2_nodupSort_allinfo[gene_idx, ] # 4,176,928
output_df2_nodupSort_nogene <- output_df2_nodupSort_allinfo[!gene_idx, ] # 43,499,740
sum(output_df2_nodupSort_nogene[, 1] > output_df2_nodupSort_nogene[, 2]) # no need to sort again; will be used for merge to build network
# [1] 0
# ------------------
output_df2_nodupSort_nogene_cutoffp005 <- output_df2_nodupSort_nogene[output_df2_nodupSort_nogene$p_value <= 0.05, ] # 15,061,815 # [1] 0.3453158 0.9940528 # do not need to filter out the negative rho value
output_df2_nodupSort_nogene_cutoffp005$rho <- round(output_df2_nodupSort_nogene_cutoffp005$rho, digits = 3)
output_df2_nodupSort_nogene_cutoffp005$p_value <- formatC(output_df2_nodupSort_nogene_cutoffp005$p_value, format = "e", digits = 2)
write.table(output_df2_nodupSort_nogene_cutoffp005, file = "mecom_var_sankaran/data/hememap/output_df2_nodupSort_nogene_cutoffp005", 
            row.names = FALSE, col.names = T, sep = "\t", quote = FALSE) # large, write it by yourself
# ------------------
quantile(table(output_df2_nodupSort_nogene_cutoffp005$queryHits)) # look at how many elements link to an element
# 0%  25%  50%  75% 100%
# 1   16   28   48  251
length(unique(output_df2_nodupSort_nogene_cutoffp005$queryHits)) # 429276/450920 # # look at how many elements were covered by this analysis
# 
# processing the direct interactions
# 
output_df_p005_all <- data.frame(fread("mecom_var_sankaran/data/hememap/12102020-peakGeneCorrelation_loopFilter_p005_union.tsv", header=T, sep="\t"))
direct_candidate_pair1 <- tss250[output_df_p005_all$query_idx, c("chr", "start", "end")] # note here the tss250 is 18476
direct_candidate_pair2 <- peak_obj[output_df_p005_all$subject_idx, c("chr", "start", "end")]
# remove self interaction
non_self_idx <- direct_candidate_pair1$end == direct_candidate_pair2$end # 2,779
output_df_p005_all <- output_df_p005_all[!non_self_idx, ] # 1,216,154
direct_candidate_pair1 <- direct_candidate_pair1[!non_self_idx, ] # 1216154
direct_candidate_pair2 <- direct_candidate_pair2[!non_self_idx, ] # 1216154
# to sort
inverse_idx <- direct_candidate_pair1$start > direct_candidate_pair2$start
sum(inverse_idx)
# [1] 609894
length(inverse_idx)
# 1216154
temp <- direct_candidate_pair1[inverse_idx, ]
direct_candidate_pair1[inverse_idx, ] <- direct_candidate_pair2[inverse_idx, ]
direct_candidate_pair2[inverse_idx, ] <- temp
inverse_idx <- direct_candidate_pair1$start > direct_candidate_pair2$start
sum(inverse_idx)
# [1] 0
direct_unique_name <- paste(direct_candidate_pair1$chr, direct_candidate_pair1$start, direct_candidate_pair1$end, direct_candidate_pair2$chr, direct_candidate_pair2$start, direct_candidate_pair2$end, sep="_") # 1216154
direct_unique_name_unique <- unique(direct_unique_name) # 1,113,613
sum(direct_unique_name_unique %in% as.character(output_df2_nodupSort_gene$ATACcor_unique_name)) # 1107921 
# note that: ATAC cor mat filtered with NA and self interactions; direct interaction also contains duplications
match_idx <- match(direct_unique_name, as.character(output_df2_nodupSort_gene$ATACcor_unique_name)) # ATAC-ATAC items overlap with direct interaction # note NA contained in this idx
# any(is.na(match_idx))
# # [1] TRUE
# sum(is.na(match_idx))
# # [1] 5977  /1216154
occur_idx <- direct_unique_name %in% as.character(output_df2_nodupSort_gene$ATACcor_unique_name) # NA will be set as FALSE in this idx
direct_mergeATACinfo <- output_df2_nodupSort_gene[match_idx, ] # 1216154 # ATAC-ATAC items overlap with direct interaction 
# ----------------
direct_mergeATACinfo_rmNA <- direct_mergeATACinfo[occur_idx, ] # rm NA # 1210177 ATAC # ATAC-ATAC items overlap with direct interaction 
output_df_p005_all_rmNA <- output_df_p005_all[occur_idx, ] # rm NA # 1210177 RNA # 1216154 # overlapping direct interaction
sum(direct_mergeATACinfo_rmNA[, 1]>direct_mergeATACinfo_rmNA[, 2]) # change the weight of direct interaction to 1
direct_mergeATACinfo_rmNA_fullinfo <- data.frame(direct_mergeATACinfo_rmNA, gene=output_df_p005_all_rmNA$gene, exp_atac_rho=output_df_p005_all_rmNA$rho, exp_atac_p=output_df_p005_all_rmNA$p_value, pchic_id=output_df_p005_all_rmNA$pchic_id, weight=0) # weight set as 0
direct_mergeATACinfo_rmNA_fullinfo$rho <- round(direct_mergeATACinfo_rmNA_fullinfo$rho, digits = 3)
direct_mergeATACinfo_rmNA_fullinfo$p_value <- formatC(direct_mergeATACinfo_rmNA_fullinfo$p_value, format = "e", digits = 2)
write.table(direct_mergeATACinfo_rmNA_fullinfo, file = "mecom_var_sankaran/data/hememap/direct_mergeATACinfo_rmNA_fullinfo.tsv", 
            row.names = FALSE, col.names = T, sep = "\t", quote = FALSE) # large, write it by your self
