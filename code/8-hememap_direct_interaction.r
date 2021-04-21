#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#   HemeMap construction
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
library(annotables)
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
#   Find the target genes
#   some code inherited from 10312020-RegMapConstruction_360k.r
#####################################################################################
# -----------------------------------
# cisRE <-> gene
# -----------------------------------

# -----------------------------------
# RNA-ATAC correlation
# 
search_radius <- 500000 # search RE for each gene within sournding 500kb
# read gene filtered and promoter integrated ATAC peaks
peak_obj <- as.data.frame(data.table::fread("mecom_var_sankaran/data/hememap/10142020-atac_peak_pro.txt")) # 450920
peak_obj <- data.frame(peak_obj, paste("peak", seq(1:nrow(peak_obj)), sep="_"))
colnames(peak_obj) <- c("chr", "start", "end", "peak_type", "peak_rank") # no M
peak_sample_mat <- as.data.frame(data.table::fread("mecom_var_sankaran/data/hememap/10142020-peak_sample_mat_pro.txt", header=T))
peak_sample_mat <- subset(peak_sample_mat, select = -c(mDC, Mega)) # remove two cell types
rownames(peak_sample_mat) <- peak_obj$peak_rank

# expression data and tss
expr_mat <- data.frame(fread("mecom_var_sankaran/data/hememap/10142020-expr_mat_16bulk.txt"))[, -1] # expr profiles; 20174
tss250 <- data.frame(fread("mecom_var_sankaran/data/hememap/10142020-tss250.txt", header=F)) # tss of genes; 20174
expr_idx <- rowSums(expr_mat)!=0 # remove 1682 genes were not expressed in all the cells
tss250 <- tss250[expr_idx, ] # 18492
expr_mat <- expr_mat[expr_idx, ] # 18492
g_dup_idx <- ! duplicated(tss250[, 1:3]) # remove the dup of gene 16 genes has the same tss 
tss250 <- tss250[g_dup_idx, ] #  18476
expr_mat <- expr_mat[g_dup_idx, ] # 18476
rownames(tss250) <- paste(tss250[, 1], tss250[, 2], tss250[, 3], sep="_")
colnames(tss250) <- c("chr", "start", "end", "gene")

# normalize count tables for exp and atac
identical(colnames(expr_mat), colnames(peak_sample_mat))# make sure that the same names of exp and atac count table
expr_mat_cpm <- edgeR::cpm(expr_mat) # cpm norm for exp
peak_sample_mat_cpm <- edgeR::cpm(peak_sample_mat) # cpm norm for ATAC count

# find the gene-RE pairs using promoter and atac peaks
# convert to GRange obj
tss250_g <- makeGRangesFromDataFrame(tss250, keep.extra.columns = TRUE)
peak_obj_g <- makeGRangesFromDataFrame(peak_obj, keep.extra.columns = TRUE)

ol <- findOverlaps(tss250_g, peak_obj_g, maxgap=search_radius) # length(ol) 4114621

rho <- sapply(1:length(ol), function(x){
        rho <- cor(log2(expr_mat_cpm[queryHits(ol)[x], ]+1), log2(peak_sample_mat_cpm[subjectHits(ol)[x], ]+1))
        }  
        )
# sum(is.na((rho))) # 3371
# Define a function to get pvalue from correlation given a sample size
rhoToP <- function(rho, n = 16){
  t <- abs(rho / sqrt((1-rho)^2/(n-2)))
  2*pt(t, n-1, lower=FALSE)
}
p_value <- rhoToP(rho)
output_df <- data.frame(data.frame(ol), rho, p_value)
colnames(output_df) <- c("query_idx", "subject_idx", "rho", "p_value")
dim(output_df) # 4114621
na_idx <- !is.na(rho)
output_df <- output_df[na_idx, ] # 4111250
output_df$rho <- round(output_df$rho, digits = 3)
output_df$p_value <- formatC(output_df$p_value, format = "e", digits = 2)
output_df <- data.frame(gene=tss250$gene[output_df$query_idx], output_df)
output_df <- data.frame(output_df, peak_obj[output_df$subject_idx, ])

output_df <- data.frame(output_df, unique_id=paste(output_df$query_idx, output_df$subject_idx, sep="_")) # output_df is the data for passed bulk-ATAC_RNA links
colnames(output_df) <- c("gene", "query_idx", "subject_idx", "rho", "p_value", "chr", "start", "end", "peak_type", "peak_rank", "unique_id")

# -----------------------------------
# PHiC data 
# 
# Load protein coding annotations
pc_gene <- unique(output_df$gene) # grch38.pc <- grch38 %>%filter(biotype == "protein_coding")
# Blood PCHIC # Note that the pCHiC_Blood.tsv.gz file is too large to upload, please download it first to the destination folder
pchic_blood <- fread("zcat < mecom_var_sankaran/data/hememap/pCHiC_Blood.tsv.gz") #[1] 1037750      24
all_cells <- colnames(pchic_blood)[7:23]
# myeloid_cells <- c("Mon","Neu","MK","Ery","Mac0","Mac1","Mac2","EP")

# Only interactions with score ≥ 5 in at least 1 cell type
allcells <- TRUE
if (allcells){
  pchic_heme <- pchic_blood %>% dplyr::rename(gene="baitName") %>%
      dplyr::select(oeChr, oeStart, oeEnd, all_cells, gene)
    } else {
  pchic_heme <- pchic_blood %>% dplyr::rename(gene="baitName") %>%
    dplyr::select(oeChr, oeStart, oeEnd, all_of(myeloid_cells), gene)
}
pchic_heme <- pchic_heme[apply(pchic_heme[, 4:(ncol(pchic_heme)-1)], 1, FUN=max) > 5, ] # find score > 5 # [1] 1037750      21
colnames(pchic_heme) <- gsub("oe", "", colnames(pchic_heme)) # clean cell names
pchic_heme <- pchic_heme[pchic_heme$gene %in% pc_gene, ] # filter to retain interested gene related interactions # 688824 # sum(table(pchic_heme$gene) != 0) [1] 16332 genes
pchic_heme$Chr <- paste0("chr",pchic_heme$Chr)
pchic_heme <- data.frame(pchic_heme, pchic_id=seq(1:nrow(pchic_heme))) # pchic_heme is the data for passed pchic_RNA links
# pchic_heme.gr <- GRanges(pchic_heme)

# -----------------------------------
# integration of correlation and pchic data
# filter based on genes in the pchic data
output_df_filter <- output_df[output_df$gene %in% pchic_heme$gene, ] # 4111250 -> 3715083 # sum(table(output_df_filter$gene) != 0)  [1] 16332 genes
# loop genes to find overlap
pc_gene_filter <- unique(pchic_heme$gene)

ol_cor_list <- lapply(1:length(pc_gene_filter), function(x){
                        pchic_df <- pchic_heme[pchic_heme$gene==pc_gene_filter[x], ]
                        cor_df <- output_df_filter[output_df_filter$gene==pc_gene_filter[x], ]
                        pchic_df.gr <- GRanges(pchic_df)
                        cor_df.gr <- GRanges(cor_df)
                        ol <- findOverlaps(cor_df.gr, pchic_df.gr) # length(ol) 4114621
                        ol_cor <- data.frame(cor_df[queryHits(ol), ], pchic_id=pchic_df[subjectHits(ol), ]$pchic_id)
                        return(ol_cor)
                    })
# quantile(sapply(ol_cor_list, nrow))

ol_cor_df <- rbindlist(ol_cor_list) # 628004
sum(ol_cor_df$p_value <= 0.05)  # how many significant correlations
# [1] 130307
length(unique(ol_cor_df$gene))  # how many gene involved
# [1] 15083
length(unique(ol_cor_df$subject_idx)) # how many RE involved
# [1] 197268
length(unique(ol_cor_df$pchic_id)) # how many loop involved
# [1] 266749
table(ol_cor_df$peak_type) # loop types
#      G     RE
#  53043 574961
setwd("/broad/sankaranlab/fyu/HemeMap/MECOM-HemeMap/bulk-ATAC_RNA_pCHIC")
write.table(ol_cor_df, file = "12102020-peakGeneCorrelation_loopFilter.tsv", 
            row.names = FALSE, col.names = T, sep = "\t", quote = FALSE)
output_df_p005 <- output_df[(output_df$p_value <= 0.05), ] # note here uses output_df not output_df_filter; 4111250->702191
output_df_p005_noloop <- output_df_p005[! output_df_p005$unique_id %in% unique(ol_cor_df$unique_id), ] # 590929
output_df_p005_noloop <- data.frame(output_df_p005_noloop, pchic_id=0)
output_df_p005_all <- rbind.data.frame(ol_cor_df, output_df_p005_noloop) # 628004+590929=1218933
setwd("mecom_var_sankaran/data/hememap/")
write.table(output_df_p005_all, file = "12102020-peakGeneCorrelation_loopFilter_p005_union.tsv", 
            row.names = FALSE, col.names = T, sep = "\t", quote = FALSE) # large, write it by yourself

