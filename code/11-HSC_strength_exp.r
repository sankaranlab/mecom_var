#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#   Explore the relationship between exp and hememap scores
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
library(tidyr)
library(edgeR)
library(dplyr)
library(data.table)
setwd("/Users/fyu/Documents/GitHub/")


load("mecom_var_sankaran/data/hememap/interaction_name_gene_info.rda.temp")  # too lare to upload, please contact to the aunthors
load("mecom_var_sankaran/data/hememap/01022021-heme_shortest_paths_geoActivate_info.rda.temp")  # too lare to upload, please contact to the aunthors
expr_mat <- data.frame(fread("mecom_var_sankaran/data/hememap/10142020-expr_mat_16bulk.txt"))[, -1] # expr profiles; 20174
tss250 <- data.frame(fread("mecom_var_sankaran/data/hememap/10142020-tss250.txt", header=F)) # tss of genes; 20174
expr_idx <- rowSums(expr_mat)!=0 # remove 1682 genes were not expressed in all the cells
tss250 <- tss250[expr_idx, ] # 18492
expr_mat <- expr_mat[expr_idx, ] # 18492
expr_mat_cpm <- edgeR::cpm(expr_mat)


interaction_strength_idx <- order(data.frame(interaction_geoATACmean18)$HSC, decreasing=T)
interaction_gene_name_sort <- as.character(interaction_df$gene[interaction_strength_idx])
chunk2 <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE)) 
interaction_gene_name_group <- chunk2(as.vector(t(interaction_gene_name_sort)), 200)
# interaction_gene_name_group <- sapply(interaction_gene_name_group, unique)
interaction_hsc_exp_group <- sapply(interaction_gene_name_group, function(x) expr_mat_cpm[match(x, tss250[, 4]), "HSC"] )
# sapply(interaction_hsc_exp_group, mean, na.rm=T)
sapply(rev(interaction_hsc_exp_group), mean, na.rm=T)
interaction_gene_name_group <- chunk2(as.vector(t(interaction_gene_name_sort)), 200)
# interaction_gene_name_group <- sapply(interaction_gene_name_group, unique)
interaction_hsc_exp_group <- sapply(interaction_gene_name_group, function(x) expr_mat_cpm[match(x, tss250[, 4]), "HSC"] )
sapply(interaction_hsc_exp_group, mean, na.rm=T)

# -------------- median visualization
library(RColorBrewer)
time_num <- 50
set.seed(9527)
interaction_gene_name_group <- chunk2(as.vector(t(interaction_gene_name_sort)), time_num)
# interaction_gene_name_group <- sapply(interaction_gene_name_group, unique)
interaction_hsc_exp_group <- sapply(interaction_gene_name_group, function(x) expr_mat_cpm[match(x, tss250[, 4]), "HSC"] )
# m <- sapply(interaction_hsc_exp_group, median, na.rm=T)  %>% rev() %>% barplot(., col=rev(as.character(jdb_palette("brewer_spectra", type="continuous"))[time_num:1]), ylim=c(0, 20), border=NA)
pdf("mecom_var_sankaran/data/exp_strength_cor/interaction_hsc_exp_strength_50_rev.pdf")
xx <- jdb_palette("brewer_spectra", 4, type ="discrete")
m <- sapply(interaction_hsc_exp_group, median, na.rm=T)  %>% 
        rev() %>% 
        barplot(., col=colorRampPalette(xx)(time_num)%>% rev(), ylim=c(0, 20), border=NA, ylab="Gene expression in HSCs", xlab="Increasing interaction strength", names.arg=NA)

xx <- rep(0, time_num)
for (i in 1:time_num){
        xx[i] <- expr_mat_cpm[sample(1:nrow(expr_mat_cpm), interaction_hsc_exp_group[[1]], replace=T), "HSC"] %>% median
}
points(m, xx, cex=1, lwd=2)
dev.off()

# # -------------- z-score visualization
# (sapply(interaction_hsc_exp_group, median, na.rm=T)-mean(xx))/sd(xx)
# plot((sapply(interaction_hsc_exp_group, median, na.rm=T)-mean(xx))/sd(xx))

# # -------------- boxplot visualization
# pdf("mecom_var_sankaran/data/exp_strength_cor/interaction_hsc_exp_group_box-nonrun_200_rev.pdf")
# boxplot(rev(interaction_hsc_exp_group), outline=FALSE, width=rep(0.4, 200), lwd=.5, 
#         # medcol = "indianred", medcex = 3,
#         ylab="gene expression in HSCs", col=rev(as.character(jdb_palette("brewer_spectra", type="continuous"))[200:1]))
# dev.off()

