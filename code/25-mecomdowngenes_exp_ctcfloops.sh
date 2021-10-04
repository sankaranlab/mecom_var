library(dplyr)
library(stringr)

CTCF_AAVS_deg <- read.table("CTCF_AAVS_GeneMat.de.txt", header=T)
sum(CTCF_AAVS_deg$RealFC>1)
# [1] 636
sum(CTCF_AAVS_deg$RealFC<1)
# [1] 557

getfpkm <- function(profiles=c("CTCF_1.genes.results", "CTCF_2.genes.results", "CTCF_3.genes.results", "AAVS1_1.genes.results", "AAVS1_2.genes.results", "AAVS1_3.genes.results"), mergedfile="../bulkrnaseq_ana/CTCF_AAVS_GeneMat.txt"){
    setwd("./bulkrnaseq")
    y <- read.table(mergedfile, header=T)
    tempmat <- data.frame(matrix(nrow=nrow(y), ncol=(length(profiles)+2)))
    rownames(tempmat) <- rownames(y)
    for (i in 1:length(profiles)){
        tempx <- read.table(profiles[i], header=T)
        tempmat[, i] <-tempx$FPKM
    }
    tempmat[, (i+1)] <- rowMeans(tempmat[, 1:3])
    tempmat[, (i+2)] <- rowMeans(tempmat[, 4:6])
    colnames(tempmat) <- c("case1", "case2", "case3", "ctrl1", "ctrl1", "ctrl1", "mean_case", "mean_ctrl")
    return(tempmat)
}
ctcf_profile <- getfpkm()
ctcfmecom_mecom_profile <- getfpkm(profiles=c("MEC_CTCF_1.genes.results", "MEC_CTCF_2.genes.results", "MEC_CTCF_3.genes.results", "MECOM_1.genes.results", "MECOM_2.genes.results", "MECOM_3.genes.results"), mergedfile="../bulkrnaseq_ana/MEC_CTCF_MECOM_GeneMat.txt")
ctcfmecom_mecom_profile_idx <- !duplicated(str_split(rownames(ctcfmecom_mecom_profile), "_", simplify = T)[, 2]) # 33514/33538 # remove duplicated genes
ctcfmecom_mecom_profile <- ctcfmecom_mecom_profile[ctcfmecom_mecom_profile_idx, ]
rownames(ctcfmecom_mecom_profile) <- str_split(rownames(ctcfmecom_mecom_profile), "_", simplify = T)[, 2]
# connect with deg
ctcf_profile_deg <- ctcf_profile[rownames(ctcf_profile) %in% rownames(CTCF_AAVS_deg), ]
ctcf_profile_deg <- ctcf_profile_deg[order(rownames(ctcf_profile_deg)), ]
CTCF_AAVS_deg <- CTCF_AAVS_deg[order(rownames(CTCF_AAVS_deg)), ]
ctcf_profile_deg <- data.frame(ctcf_profile_deg, CTCF_AAVS_deg)
rownames(ctcf_profile_deg) <- str_split(rownames(ctcf_profile_deg), "_", simplify = T)[, 2]
ctcf_profile_idx <- !duplicated(str_split(rownames(ctcf_profile), "_", simplify = T)[, 2]) # 33514/33538 # remove duplicated genes
ctcf_profile <- ctcf_profile[ctcf_profile_idx, ]
rownames(ctcf_profile) <- str_split(rownames(ctcf_profile), "_", simplify = T)[, 2]

MECOM_AAVS_deg <- read.table("MECOM_AAVS_GeneMat.de.txt", header=T)
sum(MECOM_AAVS_deg$RealFC>1)
# [1] 411
sum(MECOM_AAVS_deg$RealFC<1)
# [1] 375


RUNX1_AAVS_deg <- read.table("RUNX1_AAVS_GeneMat.de.txt", header=T)
sum(RUNX1_AAVS_deg$RealFC>1)
# [1] 436
sum(RUNX1_AAVS_deg$RealFC<1)
# [1] 361

bulkdegup_RUNX1_AAVS <- rownames(RUNX1_AAVS_deg)[RUNX1_AAVS_deg$RealFC>1] %>% 
                            str_split(., "_", simplify = T)
bulkdegup_RUNX1_AAVS <- bulkdegup_RUNX1_AAVS[, 2]
bulkdegdown_RUNX1_AAVS <- rownames(RUNX1_AAVS_deg)[RUNX1_AAVS_deg$RealFC<1] %>% 
                            str_split(., "_", simplify = T)
bulkdegdown_RUNX1_AAVS <- bulkdegdown_RUNX1_AAVS[, 2]

bulkdegup_MECOM_AAVS <- rownames(MECOM_AAVS_deg)[MECOM_AAVS_deg$RealFC>1] %>% 
                            str_split(., "_", simplify = T)
bulkdegup_MECOM_AAVS <- bulkdegup_MECOM_AAVS[, 2]
bulkdegdown_MECOM_AAVS <- rownames(MECOM_AAVS_deg)[MECOM_AAVS_deg$RealFC<1] %>% 
                            str_split(., "_", simplify = T)
bulkdegdown_MECOM_AAVS <- bulkdegdown_MECOM_AAVS[, 2]


intersect(bulkdegup_RUNX1_AAVS, bulkdegup_MECOM_AAVS) %>% length
intersect(rownames(RUNX1_AAVS_deg)[RUNX1_AAVS_deg$RealFC>1], rownames(MECOM_AAVS_deg)[MECOM_AAVS_deg$RealFC>1]) %>% length
# 38
intersect(bulkdegdown_RUNX1_AAVS, bulkdegdown_MECOM_AAVS) %>% length
intersect(rownames(RUNX1_AAVS_deg)[RUNX1_AAVS_deg$RealFC<1], rownames(MECOM_AAVS_deg)[MECOM_AAVS_deg$RealFC<1]) %>% length
# 41

intersect(bulkdegup_RUNX1_AAVS, bulkdegdown_MECOM_AAVS) %>% length

intersect(bulkdegdown_RUNX1_AAVS, bulkdegup_MECOM_AAVS) %>% length





# 10x 
degup10x <- read.table("deg_mast0035_up.txt")[, 1]
degdown10x <- read.table("deg_mast0035_down.txt")[, 1]
intersect(bulkdegup_RUNX1_AAVS, degup10x) %>% length
# [1] 11
intersect(bulkdegdown_RUNX1_AAVS, degdown10x) %>% length
# [1] 14

intersect(bulkdegup_MECOM_AAVS, degup10x) %>% length
# [1] 84
intersect(bulkdegdown_MECOM_AAVS, degdown10x) %>% length
# [1] 58

# ------------------------------------------------
# ------------------------------------------------
# CTCF editing deg analysis


# AAVS1 vs CTCF (what is overlap between MECOM down and CTCF up genes and vice versa)






# ------------------------------------------------
# ------------------------------------------------
# S5e
# 20-loop_contained_gene_analysis.r
# ----------------------------------- for ctcf_up and MECOM down
# intersect(overlap_peak$peak_gene_name, lt_down_df$Gene.Name) %>% length
setwd("/Users/fyu/Documents/GitHub/")
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(data.table)

# --------------------------------------------------------------------------------------------------------------
# look into the CTCF knockdown deg enrichment on the loop genes
setwd("/Users/fyu/Documents/GitHub/")
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
# idx2 <- rowSums(lt_hsc_exp[, 1:6]) != 0 # remove duplicates with the wrong with 0 expression all the samples
# lt_hsc_exp <- lt_hsc_exp[idx2, ] # 17206
bdg_gene_name <- rownames(lt_hsc_exp)  %>% unique 
loop_genes <- overlap_peak$peak_gene_name[overlap_peak$peak_gene_name %in% bdg_gene_name] %>% unique # 2564 unique 2200
################################### revise this part
# ctcf_deg <- read.table("mecom_var_sankaran/data/LowC_ana/ctcf_knock_deg.txt", sep="\t", header=T)
# lt_all_df <- (ctcf_deg$sample_1=="HSC_LacZ" & ctcf_deg$sample_2=="HSC_shCTCF") %>% ctcf_deg[., ] # 861
# lt_up_df <- (ctcf_deg$sample_1=="HSC_LacZ" & ctcf_deg$sample_2=="HSC_shCTCF") %>% ctcf_deg[., ] %>% subset(., log2_fold_change>0) # 296
# lt_down_df <- (ctcf_deg$sample_1=="HSC_LacZ" & ctcf_deg$sample_2=="HSC_shCTCF") %>% ctcf_deg[., ] %>% subset(., log2_fold_change<0) # 565

lt_all_df <- ctcf_profile_deg # 1193
lt_up_df <-  ctcf_profile_deg[ctcf_profile_deg$RealFC>1, ] # 636
lt_down_df <- ctcf_profile_deg[ctcf_profile_deg$RealFC<1, ] # 557


loop_genes <- overlap_peak$peak_gene_name[overlap_peak$peak_gene_name %in% bdg_gene_name] %>% unique # 2564 unique 1364
mast0035_down <- data.frame(fread("mecom_var_sankaran/data/down0035_deg_info/deg_mast0035_down.txt"))[, 1] # 321
ctcfdeg_loop_num <- intersect(mast0035_down, rownames(lt_up_df)) %>% length # 48
notctcfdeg_loop_num <- length(mast0035_down)-ctcfdeg_loop_num
ctcfdeg_notloop_num <- nrow(lt_up_df)-ctcfdeg_loop_num
notctcfdeg_notloop_num <- length(bdg_gene_name)-length(mast0035_down)-ctcfdeg_notloop_num

observed_table <- matrix(c(ctcfdeg_loop_num, notctcfdeg_loop_num, ctcfdeg_notloop_num, notctcfdeg_notloop_num),
                          nrow = 2, ncol = 2, byrow=T)
rownames(observed_table) <- c('mecom down', 'Not mecom down')
colnames(observed_table) <- c('CTCF up', 'Not CTCF up')
observed_table
#                CTCF up Not CTCF up
# mecom down          48         273
# Not mecom down     588       22376
chisq.test(observed_table)
# Pearson's Chi-squared test with Yates' continuity correction

# data:  observed_table
# X-squared = 175.04, df = 1, p-value < 2.2e-16


# ------------------------------------------------
# ------------------------------------------------
# f5j
# 202110228-mecom-Reanalyze_mast0035.r    3564-3662; 3279-3303

################ Look at the gene exp of LT-HSC and shCTCF for mecom down genes
set.seed(9527)
geneSet <- list() 
geneSet[[1]] <- overlap_peak$peak_gene_name[overlap_peak$peak_gene_name != "0"] %>% unique # 2452
geneSet[[2]] <- data.frame(fread("deg_mast0035_down.txt"))[, 1]
geneSet[[3]] <- data.frame(fread("deg_mast0035_up.txt"))[, 1]
geneSet[[4]] <- data.frame(fread("deg_mast0035_all.txt"))[, 1]
names(geneSet) <- c("loop_genes", "mast0035_down", "mast0035_up", "mast0035_all")
# lt_hsc_exp <- data.frame(fread("GSE125063_shCTCF_RNA_Seq_FPKM.tsv")) #2512
# lt_hsc_exp <- lt_hsc_exp[, 1:7] # 24331
lt_hsc_exp <- ctcf_profile # 33514
idx2 <- rowSums(lt_hsc_exp[, 1:6]) != 0 # remove duplicates with the wrong with 0 expression all the samples
lt_hsc_exp <- lt_hsc_exp[idx2, ] # 24471    
# look at mecom down gene expression
boxplot(lt_hsc_exp[rownames(lt_hsc_exp) %in% geneSet[[2]], 1:6] %>% log1p)
boxplot(lt_hsc_exp[rownames(lt_hsc_exp) %in% geneSet[[1]], 1:6] %>% log1p)
boxplot(lt_hsc_exp[rownames(lt_hsc_exp) %in% geneSet[[1]] & rownames(lt_hsc_exp) %in% geneSet[[2]], 1:6] %>% log1p)

####---------mecom down genes within ctcf loops
xx_2sample <- lt_hsc_exp[rownames(lt_hsc_exp) %in% geneSet[[1]] & rownames(lt_hsc_exp) %in% geneSet[[2]], 7:8] # 80

colnames(xx_2sample) <- c("LT_HSC_CTCFediting", "LT_HSC_Ctrl")
wilcox.test(xx_2sample[, 1], xx_2sample[, 2], paired = T)

xx_2sample_melt <- melt(xx_2sample)
# ft_sig_melt$value <- log1p(ft_sig_melt$value)
xx_2sample_melt$variable <- factor(xx_2sample_melt$variable, levels=c("LT_HSC_CTCFediting", "LT_HSC_Ctrl"))
library(ggpubr)
compare_means(value ~ variable, data = xx_2sample_melt, method="wilcox.test", paired = T)
# Visualize: Specify the comparisons you want
xx_2sample_melt$value <- log1p(xx_2sample_melt$value)
my_comparisons <- list(c("LT_HSC_CTCFediting", "LT_HSC_Ctrl"))
setwd("bulkrnaseq_ana2")
pdf("MECOM_down_genes_in_CTCF_loops.pdf")
p <- ggboxplot(xx_2sample_melt, x = "variable", y = "value",  legend = "none", color="black", palette = "jco", width=0.5, fill = "variable", title="MECOM down genes in CTCF loops")+
    # , add = "jitter", notch = TRUE)+ 
        stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label.y = c(7.2, 7), paired=T)+ # Add pairwise comparisons p-value
        stat_compare_means(label.y = 15)  +   # Add global p-value
        xlab("cell types") + ylab("log2(RPKM+1)") + ylim(1, 7.8)
print(p)
dev.off()


# ------------------------------------------------
# ------------------------------------------------
# S5f
# non-ctcf loop genes

peak_obj <- as.data.frame(data.table::fread("peak_obj_genename.tsv")) #450920
# bedpe_file <- read.table("ctcf_fp_eitherloop.bedpe")
# bedpe_file <- read.table("ctcf_fp_eitherloop.bedpe")
# bedpe_file <- read.table("ctcf_fp_eitherloop.bedpe")
bedpe_file <- read.table("no_ctcf_fp_eitherloop.bedpe") # mast0035 noCTCF loops
# bedpe_file <- bedpe_file[(bedpe_file[, 3]-bedpe_file[, 2])==25000, ]
# bedpe_file <- bedpe_file[(bedpe_file[, 5]-bedpe_file[, 2]) > 30*(25000/sqrt(2)), ] # select long-range loops to avoid peaks(anchors) too close


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


####---------mecom down genes within no-ctcf loops
xx_2sample <- lt_hsc_exp[rownames(lt_hsc_exp) %in% geneSet[[1]] & rownames(lt_hsc_exp) %in% geneSet[[2]], 7:8] # 29

# xx <- lt_hsc_exp[rownames(lt_hsc_exp) %in% geneSet[[1]] & rownames(lt_hsc_exp) %in% geneSet[[2]], 1:6] #28
# xx_2sample <- data.frame(rowMeans(xx[, 1:3]), rowMeans(xx[, 4:6]))
colnames(xx_2sample) <- c("LT_HSC_CTCFediting", "LT_HSC_Ctrl")
wilcox.test(xx_2sample[, 1], xx_2sample[, 2], paired = T)

xx_2sample_melt <- melt(xx_2sample)
# ft_sig_melt$value <- log1p(ft_sig_melt$value)
xx_2sample_melt$variable <- factor(xx_2sample_melt$variable, levels=c("LT_HSC_CTCFediting", "LT_HSC_Ctrl"))
library(ggpubr)
compare_means(value ~ variable, data = xx_2sample_melt, method="wilcox.test", paired = T)
# Visualize: Specify the comparisons you want
xx_2sample_melt$value <- log1p(xx_2sample_melt$value)
my_comparisons <- list(c("LT_HSC_CTCFediting", "LT_HSC_Ctrl"))
setwd("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/09022020-MECOM_Richard/data/20210801-cellrevision/bulkrnaseq_ana2")
pdf("MECOM_down_genes_in_noCTCF_loops.pdf")
p <- ggboxplot(xx_2sample_melt, x = "variable", y = "value",  legend = "none", color="black", palette = "jco", width=0.5, fill = "variable", title="MECOM down genes in NOT CTCF loops")+
    # , add = "jitter", notch = TRUE)+ 
        stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label.y = c(7.5, 7), paired=T)+ # Add pairwise comparisons p-value
        stat_compare_means(label.y = 15)  +   # Add global p-value
        xlab("cell types") + ylab("log2(RPKM+1)") + ylim(2, 9)
print(p)
dev.off()


################ Look at the gene exp of LT-HSC and shCTCF for mecom down genes
library(GenomicRanges)
library(data.table)
library(dplyr)
library(stringr)
set.seed(9527)
names(geneSet) <- c("loop_genes", "mast0035_down", "mast0035_up", "mast0035_all")
# lt_hsc_exp <- data.frame(fread("GSE125063_shCTCF_RNA_Seq_FPKM.tsv")) #2512
# lt_hsc_exp <- lt_hsc_exp[, 1:7] # 24331

lt_hsc_exp <- ctcfmecom_mecom_profile # 33514
idx2 <- rowSums(lt_hsc_exp[, 1:6]) != 0 # remove duplicates with the wrong with 0 expression all the samples
lt_hsc_exp <- lt_hsc_exp[idx2, ] # 24532    
# look at mecom down gene expression
boxplot(lt_hsc_exp[rownames(lt_hsc_exp) %in% geneSet[[2]], 1:6] %>% log1p)
boxplot(lt_hsc_exp[rownames(lt_hsc_exp) %in% geneSet[[1]], 1:6] %>% log1p)
boxplot(lt_hsc_exp[rownames(lt_hsc_exp) %in% geneSet[[1]] & rownames(lt_hsc_exp) %in% geneSet[[2]], 1:6] %>% log1p)

####---------mecom down genes within ctcf loops
xx_2sample <- lt_hsc_exp[rownames(lt_hsc_exp) %in% geneSet[[1]] & rownames(lt_hsc_exp) %in% geneSet[[2]], 7:8] # 80

colnames(xx_2sample) <- c("LT_HSC_MECOM+CTCFediting", "LT_HSC_MECOMediting")
wilcox.test(xx_2sample[, 1], xx_2sample[, 2], paired = T)

xx_2sample_melt <- melt(xx_2sample)
# ft_sig_melt$value <- log1p(ft_sig_melt$value)
xx_2sample_melt$variable <- factor(xx_2sample_melt$variable, levels=c("LT_HSC_MECOM+CTCFediting", "LT_HSC_MECOMediting"))
library(ggpubr)
compare_means(value ~ variable, data = xx_2sample_melt, method="wilcox.test", paired = T)
# Visualize: Specify the comparisons you want
xx_2sample_melt$value <- log1p(xx_2sample_melt$value)
my_comparisons <- list(c("LT_HSC_MECOM+CTCFediting", "LT_HSC_MECOMediting"))

pdf("MECOM_down_genes_in_CTCF_loops-mecom+ctcf_mecom.pdf")
p <- ggboxplot(xx_2sample_melt, x = "variable", y = "value",  legend = "none", color="black", palette = "jco", width=0.5, fill = "variable", title="MECOM down genes in CTCF loops")+
    # , add = "jitter", notch = TRUE)+ 
        stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label.y = c(6.2, 7), paired=T)+ # Add pairwise comparisons p-value
        stat_compare_means(label.y = 15)  +   # Add global p-value
        xlab("cell types") + ylab("log2(RPKM+1)") + ylim(1, 7.8)
print(p)
dev.off()


# ------------------------------------------------
# ------------------------------------------------
# similar to S5f
# non-ctcf loop genes

peak_obj <- as.data.frame(data.table::fread("peak_obj_genename.tsv")) #450920
bedpe_file <- read.table("no_ctcf_fp_eitherloop.bedpe") # mast0035 noCTCF loops

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
    G    RE 
  889 27760 

####---------mecom down genes within no-ctcf loops
xx_2sample <- lt_hsc_exp[rownames(lt_hsc_exp) %in% geneSet[[1]] & rownames(lt_hsc_exp) %in% geneSet[[2]], 7:8] # 29

# xx <- lt_hsc_exp[rownames(lt_hsc_exp) %in% geneSet[[1]] & rownames(lt_hsc_exp) %in% geneSet[[2]], 1:6] #28
# xx_2sample <- data.frame(rowMeans(xx[, 1:3]), rowMeans(xx[, 4:6]))
colnames(xx_2sample) <- c("LT_HSC_MECOM+CTCFediting", "LT_HSC_MECOMediting")
wilcox.test(xx_2sample[, 1], xx_2sample[, 2], paired = T)

xx_2sample_melt <- melt(xx_2sample)
# ft_sig_melt$value <- log1p(ft_sig_melt$value)
xx_2sample_melt$variable <- factor(xx_2sample_melt$variable, levels=c("LT_HSC_MECOM+CTCFediting", "LT_HSC_MECOMediting"))
library(ggpubr)
compare_means(value ~ variable, data = xx_2sample_melt, method="wilcox.test", paired = T)
# Visualize: Specify the comparisons you want
xx_2sample_melt$value <- log1p(xx_2sample_melt$value)
my_comparisons <- list(c("LT_HSC_MECOM+CTCFediting", "LT_HSC_MECOMediting"))

pdf("MECOM_down_genes_in_noCTCF_loops-mecom+ctcf_mecom.pdf")
p <- ggboxplot(xx_2sample_melt, x = "variable", y = "value",  legend = "none", color="black", palette = "jco", width=0.5, fill = "variable", title="MECOM down genes in NOT CTCF loops")+
    # , add = "jitter", notch = TRUE)+ 
        stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label.y = c(7.5, 7), paired=T)+ # Add pairwise comparisons p-value
        stat_compare_means(label.y = 15)  +   # Add global p-value
        xlab("cell types") + ylab("log2(RPKM+1)") + ylim(2, 9)
print(p)
dev.off()
