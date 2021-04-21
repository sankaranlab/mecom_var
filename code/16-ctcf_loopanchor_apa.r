# ----------------------------------------------------------------------------------------------------------------
# 1. Low-C enrichment analysis 
#       (for CTCF associated loops and non-CTCF associated loops) 
# 2. CTCF footprints distribution
# ----------------------------------------------------------------------------------------------------------------

setwd("/Users/fyu/Documents/GitHub/")
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(data.table)

# ----------------------------------------------------------------------------------------------------------------
# CTCF associated loops enrichment
# ----------------------------------------------------------------------------------------------------------------
bedpe_file <- "mecom_var_sankaran/data/LowC_ana/GSM3562031_All_Chr_Merged_Loops.bedpe"
bedpe_file <- read.table(bedpe_file, header=T) # 7358
bedpe_file[, 1] <- paste0("chr", bedpe_file[, 1])
bedpe_file[, 4] <- paste0("chr", bedpe_file[, 4])
# 5000 resoluation
bedpe_file_1 <- bedpe_file[, 1:3]
bedpe_file_2 <- bedpe_file[, 4:6]
colnames(bedpe_file_1) <- colnames(bedpe_file_2) <- c("chr", "start", "end")
# all the loop data could be found in the folder

# ----------------------
##
## mast0035_down
##
ctcf_fp <- "mecom_var_sankaran/data/CTCF_enrich/ctcf_fp_re_sort_bed_info.txt"
ctcf_fp <- unique(read.table(ctcf_fp, header=F)[, 1:3]) # 4826 # note that be duplication
colnames(ctcf_fp) <- c("chr", "start", "end")

bedpe_file_1_g <- makeGRangesFromDataFrame(bedpe_file_1)
bedpe_file_2_g <- makeGRangesFromDataFrame(bedpe_file_2)
ctcf_fp_g <- makeGRangesFromDataFrame(ctcf_fp)

ol1 <- findOverlaps(bedpe_file_1_g, ctcf_fp_g, maxgap=-1L, select="first") # sum(!is.na(ol1))        188/7358 # some of them are duplicated
ol2 <- findOverlaps(bedpe_file_2_g, ctcf_fp_g, maxgap=-1L, select="first") # sum(!is.na(ol2))        184/7358

idx <- !is.na(ol1) | !is.na(ol2) # sum(idx) # 300

ctcf_fp_bedpe_file <- bedpe_file[idx, ] # 300
all(ctcf_fp_bedpe_file[, 1]==ctcf_fp_bedpe_file[, 4]) # all the interactions are intra-chromsomal interactions
# [1] TRUE
colnames(ctcf_fp_bedpe_file)[1] <- "#chr1"
write.table(ctcf_fp_bedpe_file[, 1:10], "mecom_var_sankaran/data/CTCF_enrich/ctcf_fp_eitherloop.bedpe", row.names=F, col.names=F, sep="\t", quote=F)

# # # # # # # # 
# Not in R (make sure you install the juicer and use 25kb finally )
cd /Users/fyu/Documents/GitHub/
#  three resoluation ; output all the results including a set of text files
mkdir -p mecom_var_sankaran/data/LowC_ana/mast0035_down/APA-LT-CTCFfp_either_loop
mkdir -p mecom_var_sankaran/data/LowC_ana/mast0035_down/APA-ST-CTCFfp_either_loop
## Note that the low.c data are too large to upload, please download them to the destination folder by yourself first to perform the analysis
java -Xms512m -Xmx2048m -jar /Users/fyu/Documents/software/juicer/CPU/juicer_tools.1.9.9_jcuda.0.8.jar apa -r 5000,10000,25000 -k KR -u \
    mecom_var_sankaran/data/LowC_ana/GSM4825427_LT_inter_30.hic \
    mecom_var_sankaran/data/LowC_ana/mast0035_down/ctcf_fp_eitherloop.bedpe \
    mecom_var_sankaran/data/LowC_ana/mast0035_down/APA-LT-CTCFfp_either_loop
java -Xms512m -Xmx2048m -jar /Users/fyu/Documents/software/juicer/CPU/juicer_tools.1.9.9_jcuda.0.8.jar apa -r 5000,10000,25000 -k KR -u \
    mecom_var_sankaran/data/LowC_ana/GSM4825428_ST_inter_30.hic \
    mecom_var_sankaran/data/LowC_ana/mast0035_down/ctcf_fp_eitherloop.bedpe \
    mecom_var_sankaran/data/LowC_ana/mast0035_down/APA-ST-CTCFfp_either_loop
# # # # # # # # # # 


# define a function to draw heatmap
apa_pdf <- function(mat_file, output_pdf, p2ll, color_low, color_medi, color_high, range_num=100){
    # BiocManager::install("ComplexHeatmap")
    # install.packages("circlize")
    library("stringr")
    library("circlize")
    library(ComplexHeatmap)

    mat <- read.table(mat_file, sep=",")
    mat2 <- apply(mat, 2, function(x) {str_replace_all(x, "[\\[\\]]", "")})
    storage.mode(mat2) <- "numeric"
    mat2 <- round(mat2, 3)

    rownames(mat2) <- c(paste0("-", range_num), rep("", 9), "0", rep("", 9), range_num)
    colnames(mat2) <- c(paste0("-", range_num), rep("", 9), "0", rep("", 9), range_num)
    # col_fun = colorRamp2(c(0, .2, 0.6), c("blue", "white", "red"))
    col_fun = colorRamp2(c(color_low, color_medi, color_high), c("blue", "white", "red"))
    pdf(output_pdf)
    ht <- Heatmap(mat2, 
            col=col_fun,
            column_title = paste("P2LL", p2ll), 
            row_title = "kb from the loop anchor", 
            column_title_side = "bottom", 
            row_names_side = "left", 
            cluster_rows = FALSE, 
            cluster_columns = FALSE,
            #  show_row_names = FALSE,
            show_row_names = T,
            #  show_column_names = FALSE, 
            show_column_names = T, 
            width = unit(6, "cm"), height = unit(6, "cm"), name = "Normalized Aggregate Signals")
    draw(ht)
    dev.off()
}

setwd("mecom_var_sankaran/data/LowC_ana/mast0035_down")
apa_pdf("mecom_var_sankaran/data/LowC_ana/mast0035_down/APA-LT-CTCFfp_either_loop/25000/gw/normedAPA.txt", "CTCF_loop_LT-lowC_enrichment_25k.pdf", 0.579, 0, 0.2, 0.6)
apa_pdf("mecom_var_sankaran/data/LowC_ana/mast0035_down/APA-ST-CTCFfp_either_loop/25000/gw/normedAPA.txt", "CTCF_loop_ST-lowC_enrichment_25k.pdf", 1.844, 0, 0.2, 0.6)


setwd("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/09022020-MECOM_Richard/data/02152021-Low-C_analysis/final_pdf")
apa_pdf("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/09022020-MECOM_Richard/data/02152021-Low-C_analysis/APA-lt2/25000/gw/normedAPA.txt", "CTCF_loop_LT-lowC_enrichment.pdf", 0.552)
apa_pdf("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/09022020-MECOM_Richard/data/02152021-Low-C_analysis/APA-st2/25000/gw/normedAPA.txt", "CTCF_loop_ST-lowC_enrichment.pdf", 2.754)
apa_pdf("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/09022020-MECOM_Richard/data/02152021-Low-C_analysis/APA-lt2-all/25000/gw/normedAPA.txt", "all_loop_LT-lowC_enrichment.pdf", 2.375)
apa_pdf("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/09022020-MECOM_Richard/data/02152021-Low-C_analysis/APA-st2-all/25000/gw/normedAPA.txt", "all_loop_ST-lowC_enrichment.pdf", 2.368)


# define a function to draw the boxplots
plot_box_enrich <- function(mat_file, mat_file2, dim_num=3, file_name, add_point=T){
    read_mat <- function(x){
        mat <- read.table(x, sep=",")

        mat2 <- apply(mat, 2, function(x) {str_replace_all(x, "[\\[\\]]", "")})
        storage.mode(mat2) <- "numeric"
        mat2 <- round(mat2, 3)

        ll <- c(mat2[(21-dim_num+1):21, 1:dim_num])
        print(paste("P2ll on the dim", dim_num, "is", mat2[11, 11]/mean(ll)))
        print(paste("z-score on the dim", dim_num, "is", (mat2[11, 11]-mean(ll))/sd(ll)))
        p2ll <- mat2[11, 11]/mean(ll)
        z_score <- (mat2[11, 11]-mean(ll))/sd(ll)
        ll_z <- (ll-mean(ll))/sd(ll)
        # ll <- c(c(mat2[(21-dim_num+1):21, 1:dim_num]), mat2[11, 11])
        # print(paste("z-score on the dim", dim_num, "is", (mat2[11, 11]-mean(ll))/sd(ll)))
        return(list(ll_z, p2ll, z_score))
    }
    interaction_score1 <- read_mat(mat_file)
    interaction_score2 <- read_mat(mat_file2)
    interaction_score <- c(interaction_score1[[1]], interaction_score2[[1]])
    groups <- c(rep("LT-HSC", length(interaction_score1[[1]])), rep("ST-HSC", length(interaction_score2[[1]])))
    data_df <- data.frame(interaction_score, groups)

    groups <- c("LT-HSC", "ST-HSC")
    y_point <- c(interaction_score1[[3]], interaction_score2[[3]])
    pointdata <- data.frame(groups, y_point)

    pdf(paste0(file_name, ".pdf"))
    if(add_point){
        p <- ggplot(data = data_df, aes(x = groups, y = interaction_score, fill = groups))+
            scale_fill_manual(values=c("#dee0e1", "#2596be")) +
            #  geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
            geom_boxplot(notch = F,  outlier.size = -1, color="black",lwd=1.2, show.legend = F, width=0.5) + pretty_plot() +
            geom_point(data = pointdata, 
                    mapping = aes(x = groups, y = y_point), color="darkred", show.legend = F, size = 4) 

    } else {
        p <- ggplot(data = data_df, aes(x = groups, y = interaction_score, fill = groups))+
            scale_fill_manual(values=c("#dee0e1", "#2596be")) +
            #  geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
            geom_boxplot(notch = F,  outlier.size = -1, color="black",lwd=1.2, show.legend = F, width=0.5) + pretty_plot()
    }
    print(p)
    dev.off()
}
mat_file="mecom_var_sankaran/data/LowC_ana/mast0035_down/APA-LT-CTCFfp_either_loop/25000/gw/normedAPA.txt"
mat_file2="mecom_var_sankaran/data/LowC_ana/mast0035_down/APA-ST-CTCFfp_either_loop/25000/gw/normedAPA.txt"
myfile <- "/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/09022020-MECOM_Richard/data/02262021-downMast_re-analysis/low-c_boxplot/APAinteraction_zscore_CTCFfp_either25k_mast0035_down"
plot_box_enrich(mat_file, mat_file2, dim_num=3, myfile, add_point=T)




# ----------------------------------------------------------------------------------------
# extract non-ctcf associated loop, get the loop distribution
# ----------------------------------------------------------------------------------------
bedpe_file <- "mecom_var_sankaran/data/LowC_ana/GSM3562031_All_Chr_Merged_Loops.bedpe"
bedpe_file <- read.table(bedpe_file, header=T) # 7358
bedpe_file[, 1] <- paste0("chr", bedpe_file[, 1])
bedpe_file[, 4] <- paste0("chr", bedpe_file[, 4])
# 5000 resoluation
bedpe_file_1 <- bedpe_file[, 1:3]
bedpe_file_2 <- bedpe_file[, 4:6]
colnames(bedpe_file_1) <- colnames(bedpe_file_2) <- c("chr", "start", "end")

re_data <- read.table("mecom_var_sankaran/data/down0035_deg_info/01222021-mast0035_down_interaction_re.bed")
colnames(re_data) <- c("chr", "start", "end", "type", "peak_rank", "gene_name", "peak_gene_name", "interaction_strength", "interaction_p", "interaction_type")
ctcf_fp <- "mecom_var_sankaran/data/CTCF_enrich/ctcf_fp_re_sort_bed_info.txt"
ctcf_fp <- unique(read.table(ctcf_fp, header=F)) # 2528 # note that be duplication
colnames(ctcf_fp) <- c("chr", "start", "end")
bedpe_file_1_g <- makeGRangesFromDataFrame(bedpe_file_1)
bedpe_file_2_g <- makeGRangesFromDataFrame(bedpe_file_2)
ctcf_fp_g <- makeGRangesFromDataFrame(ctcf_fp)
re_data_g <- makeGRangesFromDataFrame(re_data)

ol1 <- findOverlaps(bedpe_file_1_g, re_data_g, maxgap=-1L, select="first") # sum(!is.na(ol1))        283/7358 # some of them are duplicated
ol2 <- findOverlaps(bedpe_file_2_g, re_data_g, maxgap=-1L, select="first") # sum(!is.na(ol2))        286/7358
idx <- !is.na(ol1) | !is.na(ol2) # sum(idx) # 418
re_bedpe_file <- bedpe_file[idx, ] # 418
all(re_bedpe_file[, 1]==re_bedpe_file[, 4]) # all the interactions are intra-chromsomal interactions
# [1] TRUE

bedpe_file_1 <- re_bedpe_file[, 1:3]
bedpe_file_2 <- re_bedpe_file[, 4:6]
colnames(bedpe_file_1) <- colnames(bedpe_file_2) <- c("chr", "start", "end")
bedpe_file_1_g <- makeGRangesFromDataFrame(bedpe_file_1)
bedpe_file_2_g <- makeGRangesFromDataFrame(bedpe_file_2)
ol1 <- findOverlaps(bedpe_file_1_g, ctcf_fp_g, maxgap=-1L, select="first") # sum(!is.na(ol1))        188/7358 # some of them are duplicated
ol2 <- findOverlaps(bedpe_file_2_g, ctcf_fp_g, maxgap=-1L, select="first") # sum(!is.na(ol2))        184/7358
idx <- !is.na(ol1) | !is.na(ol2) # sum(idx) # 300 # in either anchor
idx2 <- !is.na(ol1) & !is.na(ol2) # 72 # in both anchor
ctcf_bedpe_file <- re_bedpe_file[idx, ] # 300
# which(idx2)
ol_f_1 <- findOverlaps(bedpe_file_1_g, ctcf_fp_g, maxgap=-1L) # sum(!is.na(ol1))        348/7358 # some of them are duplicated
ol_f_2 <- findOverlaps(bedpe_file_2_g, ctcf_fp_g, maxgap=-1L) # sum(!is.na(ol2))        366/7358

# 
# get the loop with no ctcf fp and perform the apa analysis
# 
no_ctcf_idx <- setdiff(1:nrow(re_bedpe_file), unique(c(queryHits(ol_f_1), queryHits(ol_f_2))))  # sum(idx) # 118
no_ctcf_fp_bedpe_file <- re_bedpe_file[no_ctcf_idx, ] # 118
colnames(no_ctcf_fp_bedpe_file)[1] <- "#chr1"
write.table(no_ctcf_fp_bedpe_file[, 1:10], "mecom_var_sankaran/data/LowC_ana/mast0035_down/no_ctcf_fp_eitherloop.bedpe", row.names=F, col.names=F, sep="\t", quote=F)

#  three resoluation ; output all the results including a set of text files
mkdir -p mecom_var_sankaran/data/LowC_ana/mast0035_down/APA-LT-noCTCFfp_either_loop
mkdir -p mecom_var_sankaran/data/LowC_ana/mast0035_down/APA-ST-noCTCFfp_either_loop
java -Xms512m -Xmx2048m -jar /Users/fyu/Documents/software/juicer/CPU/juicer_tools.1.9.9_jcuda.0.8.jar apa -r 5000,10000,25000 -k KR -u \
    mecom_var_sankaran/data/LowC_ana/GSM4825427_LT_inter_30.hic \
    mecom_var_sankaran/data/LowC_ana/mast0035_down/no_ctcf_fp_eitherloop.bedpe \
    mecom_var_sankaran/data/LowC_ana/mast0035_down/APA-LT-noCTCFfp_either_loop
java -Xms512m -Xmx2048m -jar /Users/fyu/Documents/software/juicer/CPU/juicer_tools.1.9.9_jcuda.0.8.jar apa -r 5000,10000,25000 -k KR -u \
    mecom_var_sankaran/data/LowC_ana/GSM4825428_ST_inter_30.hic \
    mecom_var_sankaran/data/LowC_ana/mast0035_down/no_ctcf_fp_eitherloop.bedpe \
    mecom_var_sankaran/data/LowC_ana/mast0035_down/APA-ST-noCTCFfp_either_loop


setwd("mecom_var_sankaran/data/LowC_ana/mast0035_down")
apa_pdf("mecom_var_sankaran/data/LowC_ana/mast0035_down/APA-LT-noCTCFfp_either_loop/25000/gw/normedAPA.txt", "noCTCF_loop_LT-lowC_enrichment_25k.pdf", 0.835, 0, 0.2, 0.6, range_num=250)
apa_pdf("mecom_var_sankaran/data/LowC_ana/mast0035_down/APA-ST-noCTCFfp_either_loop/25000/gw/normedAPA.txt", "noCTCF_loop_ST-lowC_enrichment_25k.pdf", 1.092, 0, 0.2, 0.6, range_num=250)

mat_file="mecom_var_sankaran/data/LowC_ana/mast0035_down/APA-LT-noCTCFfp_either_loop/25000/gw/normedAPA.txt"
mat_file2="mecom_var_sankaran/data/LowC_ana/mast0035_down/APA-ST-noCTCFfp_either_loop/25000/gw/normedAPA.txt"
myfile <- "/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/09022020-MECOM_Richard/data/02262021-downMast_re-analysis/low-c_boxplot/APAinteraction_zscore_noCTCFfp_either25k_mast0035_down"
plot_box_enrich(mat_file, mat_file2, dim_num=3, myfile, add_point=T)


# 
# calculate how many loops for different ctcf orientation
# 
opposite_num <- 0
same_num <- 0
for (i in which(idx2)){
    fp_1_idx <- subjectHits(ol_f_1)[queryHits(ol_f_1)==i]
    fp_2_idx <- subjectHits(ol_f_2)[queryHits(ol_f_2)==i]

    anchor1_uni <- ctcf_fp[fp_1_idx, 6] %>% as.character %>% unique
    anchor2_uni <- ctcf_fp[fp_2_idx, 6] %>% as.character %>% unique

    if((anchor1_uni  %>% length) >=2 | (anchor2_uni %>% length)>=2){
        opposite_num <- opposite_num+1
    } else if (anchor1_uni != anchor2_uni){
        opposite_num <- opposite_num+1
    } else {
        same_num <- same_num+1
    }
}
opposite_num
# [1] 60
same_num
# [1] 12


library(ggplot2)
library(dplyr)
# Create Data
data <- data.frame(
  group=factor(c("no CTCF fp", "CTCF fp in one anchor", "CTCF fp in both anchors with opposite orientation ", "CTCF fp in both anchors with same orientation"), levels=c("CTCF fp in one anchor", "CTCF fp in both anchors with opposite orientation ", "CTCF fp in both anchors with same orientation", "no CTCF fp") %>% rev),
  value=c(118, 228, 60, 12)
)

# Compute the position of labels
data <- data %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

library(scales)
blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )
ggplot(data, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none")  + scale_fill_brewer("Blues") + blank_theme +
  theme(axis.text.x=element_blank())+
  geom_text(aes(y = ypos, label = value), size=2) + ggtitle("MECOM associated loops with CTCF footprints") + guides(fill=guide_legend(title="Loop Categories")) 
# piechart_loop_CTCFfp_percentage_ori_mast0035_down.pdf
