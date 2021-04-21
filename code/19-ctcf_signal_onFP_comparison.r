# ----------------------------------------------------------------------------------------------------------------
# comparison of CTCF signal between HSC and other differentiated cells
# ----------------------------------------------------------------------------------------------------------------
setwd("/Users/fyu/Documents/GitHub/")
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(data.table)

heatmap_value <- read.table("mecom_var_sankaran/data/CTCF_enrich/deeptools_ana/heatmap_value.mat", skip=1)
plot_box_comparison <- function(bin_num, inte_fun, z_score=T, log1p_score=F){
    i=0
    CD34_signal <- heatmap_value[, ((100+i*200-bin_num)+1):(100+i*200+bin_num)]
    i=1
    Ery_signal <- heatmap_value[, ((100+i*200-bin_num)+1):(100+i*200+bin_num)]
    i=2
    T_signal <- heatmap_value[, ((100+i*200-bin_num)+1):(100+i*200+bin_num)]
    i=3
    B_signal <- heatmap_value[, ((100+i*200-bin_num)+1):(100+i*200+bin_num)]
    i=4
    Mono_signal <- heatmap_value[, ((100+i*200-bin_num)+1):(100+i*200+bin_num)]

    CD34_signal_sum <- apply(CD34_signal, 1, sum)
    Ery_signal_sum <- apply(Ery_signal, 1, sum)
    T_signal_sum <- apply(T_signal, 1, sum)
    B_signal_sum <- apply(B_signal, 1, sum)
    Mono_signal_sum <- apply(Mono_signal, 1, sum)
    if (z_score){
        CD34_signal_sum <- (CD34_signal_sum-mean(CD34_signal_sum))/sd(CD34_signal_sum)
        Ery_signal_sum <- (Ery_signal_sum-mean(Ery_signal_sum))/sd(Ery_signal_sum)
        T_signal_sum <- (T_signal_sum-mean(T_signal_sum))/sd(T_signal_sum)
        B_signal_sum <- (B_signal_sum-mean(B_signal_sum))/sd(B_signal_sum)
        Mono_signal_sum <- (Mono_signal_sum-mean(Mono_signal_sum))/sd(Mono_signal_sum)
    }
    if (log1p_score){
        CD34_signal_sum <- log1p(CD34_signal_sum)
        Ery_signal_sum <- log1p(Ery_signal_sum)
        T_signal_sum <- log1p(T_signal_sum)
        B_signal_sum <- log1p(B_signal_sum)
        Mono_signal_sum <- log1p(Mono_signal_sum)
    }

    CD34_signal_mean <- apply(CD34_signal, 1, mean)
    Ery_signal_mean <- apply(Ery_signal, 1, mean)
    T_signal_mean <- apply(T_signal, 1, mean)
    B_signal_mean <- apply(B_signal, 1, mean)
    Mono_signal_mean <- apply(Mono_signal, 1, mean)
    if (z_score){
        CD34_signal_mean <- (CD34_signal_mean-mean(CD34_signal_mean))/sd(CD34_signal_mean)
        Ery_signal_mean <- (Ery_signal_mean-mean(Ery_signal_mean))/sd(Ery_signal_mean)
        T_signal_mean <- (T_signal_mean-mean(T_signal_mean))/sd(T_signal_mean)
        B_signal_mean <- (B_signal_mean-mean(B_signal_mean))/sd(B_signal_mean)
        Mono_signal_mean <- (Mono_signal_mean-mean(Mono_signal_mean))/sd(Mono_signal_mean)
    }
    if (log1p_score){
        CD34_signal_mean <- log1p(CD34_signal_mean)
        Ery_signal_mean <- log1p(Ery_signal_mean)
        T_signal_mean <- log1p(T_signal_mean)
        B_signal_mean <- log1p(B_signal_mean)
        Mono_signal_mean <- log1p(Mono_signal_mean)
    }

    if (inte_fun=="sum"){
        ft_sig <- list(CD34_signal_sum, Ery_signal_sum, T_signal_sum, B_signal_sum, Mono_signal_sum)
    } else {
        ft_sig <- list(CD34_signal_mean, Ery_signal_mean, T_signal_mean, B_signal_mean, Mono_signal_mean)
    }
    return(ft_sig)
}

ft_sig <- plot_box_comparison(5, "sum", T)
names(ft_sig) <- c("CD34", "Ery", "T", "B", "Mono")
ft_sig_melt <- melt(ft_sig)
# ft_sig_melt$value <- log1p(ft_sig_melt$value)
ft_sig_melt$L1 <- factor(ft_sig_melt$L1, levels=c("CD34", "Ery", "T", "B", "Mono"))
library(ggpubr)
compare_means(value ~ L1, data = ft_sig_melt, method="wilcox.test", paired = T)

# Visualize: Specify the comparisons you want
my_comparisons <- list(c("CD34", "Ery"), c("CD34", "T"), c("CD34", "B"), c("CD34", "Mono"))

p <- ggboxplot(ft_sig_melt, x = "L1", y = "value", color = "L1",palette = "nejm", legend = "none")+ 
        stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label.y = c(12, 11, 10, 9), paired=T)+ # Add pairwise comparisons p-value
        stat_compare_means(label.y = 15)  +   # Add global p-value
        xlab("cell types") + ylab("z-score")
p
