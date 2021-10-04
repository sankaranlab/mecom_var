# validation of hememap interactions with hic of hspc
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2861nnn/GSM2861708/suppl/GSM2861708_HIC1258_aligned_inter_30.hic
# conda install hicexplorer -c bioconda -c conda-forge
# 25 kb resolution
hicConvertFormat -m GSM2861708_HIC1258_aligned_inter_30.hic \
                --inputFormat hic \
                --outputFormat cool \
                -o GSM2861708_HIC1258_aligned_inter_30.hic.10k.chr20.cool \
                --resolutions 25000 \
                --chromosome 20


# ------------ in R
# BiocManager::install("HiCcompare")
library(HiCcompare)
library(data.table)
library(GenomicRanges)
# load("hic_dat.rda")
file_path <- "GSM2861708_HIC1258_aligned_inter_30.hic.10k.chr20_10000.cool"
#A list with two items. Item 1, "cis" contains the intra-chromosomal contact matrices, one per chromosome. Item 2, "trans" contains the inter-chromsomal contact matrix.
hic_dat <- cooler2bedpe(path = file_path)
save(hic_dat, file="hic_dat.rda")


d_interaction <- read.table("2102020-peakGeneCorrelation_loopFilter_p005_union.tsv", header=T, sep="\t")
peak_obj <- as.data.frame(data.table::fread("10142020-atac_peak_pro.txt")) # 450920
peak_obj <- data.frame(peak_obj, paste("peak", seq(1:nrow(peak_obj)), sep="_"))
colnames(peak_obj) <- c("chr", "start", "end", "peak_type", "peak_rank") # no M

# expression data and tss
expr_mat <- data.frame(fread("10142020-expr_mat_16bulk.txt"))[, -1] # expr profiles; 20174
tss250 <- data.frame(fread("10142020-tss250.txt", header=F)) # tss of genes; 20174
expr_idx <- rowSums(expr_mat)!=0 # remove 1682 genes were not expressed in all the cells
tss250 <- tss250[expr_idx, ] # 18492
expr_mat <- expr_mat[expr_idx, ] # 18492
g_dup_idx <- ! duplicated(tss250[, 1:3]) # remove the dup of gene 16 genes has the same tss 
tss250 <- tss250[g_dup_idx, ] #  18476
expr_mat <- expr_mat[g_dup_idx, ] # 18476
rownames(tss250) <- paste(tss250[, 1], tss250[, 2], tss250[, 3], sep="_")
colnames(tss250) <- c("chr", "start", "end", "gene")

# find the gene-RE pairs using promoter and atac peaks
# convert to GRange obj
tss250_g <- makeGRangesFromDataFrame(tss250, keep.extra.columns = TRUE)
peak_obj_g <- makeGRangesFromDataFrame(peak_obj, keep.extra.columns = TRUE)


d_interaction_pchic <- d_interaction[d_interaction$pchic_id>0, ]
d_interaction_cor <- d_interaction[d_interaction$pchic_id==0, ]
dim(d_interaction_pchic)
# [1] 628004     12
dim(d_interaction_cor)
# [1] 590929     12

##--------
##--------for pchic
peak1_pchic <- tss250[d_interaction_pchic$query_idx, ]
peak2_pchic <- peak_obj[d_interaction_pchic$subject_idx, ]
# convert to GRange obj
peak1_pchic_g <- makeGRangesFromDataFrame(peak1_pchic, keep.extra.columns = TRUE)
peak2_pchic_g <- makeGRangesFromDataFrame(peak2_pchic, keep.extra.columns = TRUE)


olap_df <- data.frame(matrix(nrow=nrow(d_interaction_pchic), ncol=22))
for (i in 1:22){
    hic_peak1 <- hic_dat$cis[[i]][, 1:3]
    hic_peak1[, 1] <- paste0("chr", hic_peak1[, 1])
    hic_peak2 <- hic_dat$cis[[i]][, 4:6]
    hic_peak2[, 1] <- paste0("chr", hic_peak2[, 1])
    colnames(hic_peak1) <- colnames(hic_peak2) <- c("chr", "start", "end")
    hic_peak1_g <- makeGRangesFromDataFrame(hic_peak1, keep.extra.columns = TRUE)
    hic_peak2_g <- makeGRangesFromDataFrame(hic_peak2, keep.extra.columns = TRUE)
    ol1 <- GenomicRanges::findOverlaps(peak1_pchic_g,
                                    hic_peak1_g,
                                    select="all")
    ol2 <- GenomicRanges::findOverlaps(peak2_pchic_g,
                                        hic_peak2_g,
                                        select="all")

    length(unique(queryHits(ol1))) # 67223
    length(unique(queryHits(ol2))) # 67245
    length(intersect(unique(queryHits(ol1)), unique(queryHits(ol2)))) # 67166

    ol1 <- as.list(ol1)
    ol2 <- as.list(ol2)
    olap <- mapply(function(o1, o2) {
        if (length(o1) == 0) out <- FALSE
        else if (length(intersect(o1,o2) > 0)) out <- TRUE
        else out <- FALSE
        return(out)
    }, ol1, ol2)
    olap_df[, i] <- olap
    print(i)    
}

# peaks are in the 25k anchors
distance_peak_pchic <- GenomicRanges::distance(peak1_pchic_g, peak2_pchic_g)
sum(distance_peak_pchic<10000)
# [1] 8455
(sum(rowSums(olap_df))+sum(distance_peak_pchic<10000))/nrow(d_interaction_pchic)
# [1] 0.5127372



##--------
##--------for cor
peak1_cor <- tss250[d_interaction_cor$query_idx, ]
peak2_cor <- peak_obj[d_interaction_cor$subject_idx, ]
# convert to GRange obj
peak1_cor_g <- makeGRangesFromDataFrame(peak1_cor, keep.extra.columns = TRUE)
peak2_cor_g <- makeGRangesFromDataFrame(peak2_cor, keep.extra.columns = TRUE)

olap_df <- data.frame(matrix(nrow=nrow(d_interaction_cor), ncol=22))
for (i in 1:22){
    hic_peak1 <- hic_dat$cis[[i]][, 1:3]
    hic_peak1[, 1] <- paste0("chr", hic_peak1[, 1])
    hic_peak2 <- hic_dat$cis[[i]][, 4:6]
    hic_peak2[, 1] <- paste0("chr", hic_peak2[, 1])
    colnames(hic_peak1) <- colnames(hic_peak2) <- c("chr", "start", "end")
    hic_peak1_g <- makeGRangesFromDataFrame(hic_peak1, keep.extra.columns = TRUE)
    hic_peak2_g <- makeGRangesFromDataFrame(hic_peak2, keep.extra.columns = TRUE)
    ol1 <- GenomicRanges::findOverlaps(peak1_cor_g,
                                    hic_peak1_g,
                                    select="all")
    ol2 <- GenomicRanges::findOverlaps(peak2_cor_g,
                                        hic_peak2_g,
                                        select="all")

    length(unique(queryHits(ol1))) # 67223
    length(unique(queryHits(ol2))) # 67245
    length(intersect(unique(queryHits(ol1)), unique(queryHits(ol2)))) # 67166

    ol1 <- as.list(ol1)
    ol2 <- as.list(ol2)
    olap <- mapply(function(o1, o2) {
        if (length(o1) == 0) out <- FALSE
        else if (length(intersect(o1,o2) > 0)) out <- TRUE
        else out <- FALSE
        return(out)
    }, ol1, ol2)
    olap_df[, i] <- olap
    print(i)    
}

# peaks are in the 25k anchors
distance_peak_cor <- GenomicRanges::distance(peak1_cor_g, peak2_cor_g)
sum(distance_peak_cor<10000)
# [1] 31608
(sum(rowSums(olap_df))+sum(distance_peak_cor<10000))/nrow(d_interaction_cor)
# [1] 0.5626801

