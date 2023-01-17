# ----------------------------------------------------
# ----------------------------------------------------
# pseudo-bulk DEG analysis
# ----------------------------------------------------
# ----------------------------------------------------

# 1. Get merged profiles of gene_by_group (2 columns) from gene_by_cell count matrix

library(Seurat)
library(dplyr)
library(ggplot2)
# library("DESeq2")

#------------------------------------------------------------------------------------
# 
# in the broad 
# 
load("aavs1_11k.Rdata")
load("mecom_16k.Rdata")

aavs1_11k[["percent.mt"]] <- PercentageFeatureSet(aavs1_11k, pattern = "^MT-")
aavs1_11k <- subset(aavs1_11k, subset = nFeature_RNA >500 & nFeature_RNA <10000 & percent.mt < 10 )
aavs1_11k <- NormalizeData(object = aavs1_11k, normalization.method = "LogNormalize", scale.factor = 1e4)

mecom_16k[["percent.mt"]] <- PercentageFeatureSet(mecom_16k, pattern = "^MT-")
mecom_16k <- subset(mecom_16k, subset = nFeature_RNA > 500 & nFeature_RNA <10000 & percent.mt < 10 )
mecom_16k <- NormalizeData(object = mecom_16k, normalization.method = "LogNormalize", scale.factor = 1e4)

mecom_count_mat <- mecom_16k@assays$RNA@counts
aavs_count_mat <- aavs1_11k@assays$RNA@counts
mecom_count_bin <- mecom_count_mat
aavs_count_bin <- aavs_count_mat
mecom_count_bin@x[mecom_count_bin@x>=1] <- 1
aavs_count_bin@x[aavs_count_bin@x>=1] <- 1
mecom_count_bin_rowsum <- rowSums(mecom_count_bin)
aavs_count_bin_rowsum <- rowSums(aavs_count_bin)


which(rownames(mecom_count_bin)=="MECOM")
# [1] 6826
mecom_count_bin_rowsum[6826]
# MECOM 
#   235 
aavs_count_bin_rowsum[6826]
# MECOM 
#   647 
235/ncol(mecom_count_bin)
# [1] 0.03964238
647/ncol(aavs_count_bin)
# [1] 0.1509214

mecom_idx <- mecom_count_bin_rowsum>=floor(ncol(mecom_count_bin)*0.039) # 9768
aavs_idx <- aavs_count_bin_rowsum>=floor(ncol(aavs_count_bin)*0.1) # 7958
mecom_count_mat <- mecom_count_mat[mecom_idx & aavs_idx, ]
aavs_count_mat <- aavs_count_mat[mecom_idx & aavs_idx, ]

mecom_pseudo <- rowSums(mecom_count_mat) # 7958
aavs_pseudo <- rowSums(aavs_count_mat)



save.image("pesudo-bulk_ana.rda")

load("pesudo-bulk_ana.rda")
mast0035down <- read.table("deg_mast0035_down.txt", header=F)[, 1]
mast0035up <- read.table("deg_mast0035_up.txt", header=F)[, 1]
# rerun 404-450
fclist <- sort((mecom_pseudo/aavs_pseudo))
log2fclist <- sort((mecom_pseudo/aavs_pseudo) %>% log2)
genename <- names(log2fclist)
deglist <- genename
length(mast0035down)
# [1] 322
length(mast0035up)
# [1] 402
sum(genename %in% mast0035down)
# [1] 313
sum(genename %in% mast0035up)
# [1] 400

deglist[genename %in% mast0035down | genename %in% mast0035up] <- ""
deg_upordown <- rep("not DEGs", length(deglist))
deg_upordown[genename %in% mast0035down] <- "MECOM down genes"
deg_upordown[genename %in% mast0035up] <- "MECOM up genes"
box_down=c('MECOM','HOPX','RBPMS','MLLT3','HEMGN','SOCS2','ALDH2','MSI2','JAML','BEX2')
box_up=c('MPO','IGFBP2','CTSG','SRGN','IGLL1','MIF','CALR','GFI1B','MYB','PRAM1')
top10gene <- c(genename[1:10], rev(genename)[1:10])
top5gene <- c(genename[1:5], rev(genename)[1:5])
xcoordinate <- seq(1, length(log2fclist))
mydata <- data.frame(xcoordinate, log2fclist, fclist, genename, deglist, deg_upordown)
mydata2 <- filter(mydata, deg_upordown!="not DEGs")
mydata3 <- filter(mydata2, genename  %in% c(box_down, box_up))
mydatatop10 <- filter(mydata2, genename  %in% top10gene)
mydatatop5 <- filter(mydata2, genename  %in% top5gene)
mycolor <- c("#810e11", "#0c0c83", "#b7b6b8") # red, blue, gray
xx <- ggplot(mydata, aes(xcoordinate, log2(fclist), color = factor(deg_upordown)))+
        geom_point() + theme_classic(base_size = 10) + scale_color_manual(breaks = c("MECOM down genes", "MECOM up genes", "not DEGs"), values=(mycolor)) +
        # scale_y_continuous(trans='log2') + scale_fill_brewer(palette="Dark2") 
        # ylim(0.5, 2)
        ylim(-2, 1) +
        ggrepel::geom_text_repel(data=mydatatop5, aes(xcoordinate, log2(fclist), label = genename), color = 'black',
                        size = 2, max.overlaps=300) + xlab("gene ranking") +ylab("log2(fold change)")
        theme(legend.position = "bottom")
xx
ggsave("pseudobulk_deganalysis3.pdf", width = 5, height = 5)

######## compare the fold change of DEGs of singlecell and pseudo-bulk data
load("pesudo-bulk_ana.rda")
mast0035down <- read.table("deg_mast0035_down.txt", header=F)[, 1]
mast0035up <- read.table("deg_mast0035_up.txt", header=F)[, 1]
# rerun 404-450
fclist <- sort((mecom_pseudo/aavs_pseudo))
log2fclist <- sort((mecom_pseudo/aavs_pseudo) %>% log2)
genename <- names(log2fclist)
deglist <- genename

mast035deg <- read.csv("mast_de_merge_norm_1mecom_2aavs_0.035.csv")
deglist[genename %in% mast0035down | genename %in% mast0035up] <- ""
deg_upordown <- rep("not DEGs", length(deglist))
deg_upordown[genename %in% mast0035down] <- "MECOM down genes"
deg_upordown[genename %in% mast0035up] <- "MECOM up genes"
box_down=c('MECOM','HOPX','RBPMS','MLLT3','HEMGN','SOCS2','ALDH2','MSI2','JAML','BEX2')
box_up=c('MPO','IGFBP2','CTSG','SRGN','IGLL1','MIF','CALR','GFI1B','MYB','PRAM1')
top10gene <- c(genename[1:10], rev(genename)[1:10])
top5gene <- c(genename[1:5], rev(genename)[1:5])
xcoordinate <- seq(1, length(log2fclist))
mydata <- data.frame(xcoordinate, log2fclist, fclist, genename, deglist, deg_upordown)
mydata2 <- filter(mydata, deg_upordown!="not DEGs")
mydata3 <- filter(mydata2, genename  %in% c(box_down, box_up))
mydatatop10 <- filter(mydata2, genename  %in% top10gene)
mydatatop5 <- filter(mydata2, genename  %in% top5gene)
colnames(mast035deg)[1] <- "genename"
combineddeg <- merge(mydata2, mast035deg, "genename")
mydatatop10 <- filter(mydata2, genename  %in% top10gene)
top5gene <- c(genename[1:5], rev(genename)[1:5])
top5gene2 <- c(mast035deg[order(mast035deg$avg_log2FC), ]$genename[1:5], rev(mast035deg[order(mast035deg$avg_log2FC), ]$genename)[1:5])
eithertop5 <- filter(combineddeg, (genename  %in% top5gene) | (genename  %in% top5gene2))

cor.test(combineddeg$avg_log2FC, combineddeg$log2fclist, method="spearman")
# 	Spearman's rank correlation rho

# data:  combineddeg$avg_log2FC and combineddeg$log2fclist
# S = 9420570, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.8498181 

# plot(combineddeg$log2fclist, combineddeg$avg_log2FC, ylim=c(-2, 1), xlim=c(-2, 1), pch=20)
# both top5gene


ggplot(combineddeg, aes(avg_log2FC, log2fclist, color=deg_upordown))+
        geom_point() + theme_classic(base_size = 10) +
        ylim(-1.8, 1) + xlim(-0.3, 0.5)

ggplot(combineddeg, aes(avg_log2FC, log2fclist, color=deg_upordown))+
        geom_point() + theme_classic(base_size = 10) + xlim(-0.5, 0.6) +  ylim(-1.7, 1) +
        ggrepel::geom_text_repel(data=eithertop5, aes(avg_log2FC, log2fclist, label = genename), color = 'black',
                        size = 2, max.overlaps=300) + xlab("log2(fold change) from single-cell data") +ylab("log2(fold change) from pseudo bulk data")

ggsave("foldchange_comparison.pdf", width = 5, height = 5)

