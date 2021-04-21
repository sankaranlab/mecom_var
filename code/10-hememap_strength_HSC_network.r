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
library(igraph)
setwd("/Users/fyu/Documents/GitHub/")

###########################################################################################################
#   construction of network, 
#   find the shortestpath, 
#   calculate the hememap socre, 
#   define HSC network
###########################################################################################################
#
# merge the all the gene associated and RE associated interactions
# change the weight to 1-rho (0 for gene associated interactions)
# find all the indirect interactions
# network=direct G-RE + sig RE-RE
# search pair=all the G-RE except for direct G-RE
#
direct_mergeATACinfo_rmNA_fullinfo <- data.frame(fread("mecom_var_sankaran/data/hememap/direct_mergeATACinfo_rmNA_fullinfo.tsv")) # Note that the direct_mergeATACinfo_rmNA_fullinfo.tsv file is too large to upload, please download it first to the destination folder
output_df2_nodupSort_nogene_cutoffp005 <- data.frame(fread("mecom_var_sankaran/data/hememap/output_df2_nodupSort_nogene_cutoffp005"))
sum(output_df2_nodupSort_nogene_cutoffp005$rho < 0) #15,061,815
# [1] 0
direct_mergeATACinfo_rmNA_weight1 <- direct_mergeATACinfo_rmNA_fullinfo[, 1:9] # 1,210,177
direct_mergeATACinfo_rmNA_weight1$rho <- 0
output_df2_nodupSort_nogene_cutoffp005$rho <- 1-output_df2_nodupSort_nogene_cutoffp005$rho


# 
# create node file
# 
temp1 <- unique(data.frame(output_df2_nodupSort_nogene_cutoffp005$queryHits, output_df2_nodupSort_nogene_cutoffp005$peak1_type))
temp2 <- unique(data.frame(output_df2_nodupSort_nogene_cutoffp005$subjectHits, output_df2_nodupSort_nogene_cutoffp005$peak2_type))
temp3 <- unique(data.frame(direct_mergeATACinfo_rmNA_weight1$queryHits, direct_mergeATACinfo_rmNA_weight1$peak1_type))
temp4 <- unique(data.frame(direct_mergeATACinfo_rmNA_weight1$subjectHits, direct_mergeATACinfo_rmNA_weight1$peak2_type))
colnames(temp1) <- colnames(temp2) <- colnames(temp3) <- colnames(temp4) <- c("node_id", "type")
node_file <- unique(rbind.data.frame(temp1, temp2, temp3, temp4)) # 449869 (not cutoff filter -- 449869; almost no change)

# 
# create link file
# 
link_file <- rbind.data.frame(output_df2_nodupSort_nogene_cutoffp005[, c("queryHits", "subjectHits", "rho", "peak1_type", "peak2_type")],
                            direct_mergeATACinfo_rmNA_weight1[, c("queryHits", "subjectHits", "rho", "peak1_type", "peak2_type")]) # large, generate it by your self
colnames(link_file) <- c("from", "to", "weight", "type1", "type2") # 16,271,992
rownames(link_file) <- NULL

link_file <- unique(link_file) # remove duplicates
sum(link_file[, 1] == link_file[, 2])
# [1] 0

setwd("mecom_var_sankaran/data/hememap/")
write.table(node_file, file = "12312020-node_file.tsv", 
            row.names = FALSE, col.names = T, sep = "\t", quote = FALSE)
write.table(link_file, file = "12312020-link_file.tsv", 
            row.names = FALSE, col.names = T, sep = "\t", quote = FALSE)
# link_file <- data.frame(fread("mecom_var_sankaran/data/hememap/12312020-link_file.tsv"))


# net <- graph_from_data_frame(d=link_file, vertices=node_file, directed=F)
#
# find the search pair between G and surrounding RE
#
search_radius <- 500000 # search RE for each gene within sournding 500kb
# read gene filtered and promoter integrated ATAC peaks
peak_obj <- as.data.frame(data.table::fread("mecom_var_sankaran/data/hememap/10142020-atac_peak_pro.txt")) #450920
peak_obj <- data.frame(peak_obj, 1:nrow(peak_obj))
colnames(peak_obj) <- c("chr", "start", "end", "peak_type", "peak_rank") # no M
gene_id <- which(peak_obj$peak_type=="G")
peak_id <- 1:nrow(peak_obj)
gene_peak_obj_g <- makeGRangesFromDataFrame(peak_obj[gene_id, ], keep.extra.columns = TRUE)
peak_obj_g <- makeGRangesFromDataFrame(peak_obj, keep.extra.columns = TRUE)
ol <- findOverlaps(gene_peak_obj_g, peak_obj_g, maxgap=search_radius) # length(ol) 4394108
pair1 <- peak_obj[gene_id, ][queryHits(ol), ]
pair2 <- peak_obj[subjectHits(ol), ]
# head(subjectHits(ol))
# head(peak_obj[subjectHits(ol), "peak_rank"])

#
# find the shortest path
# 
# make sure the nodes of searching pairs in the list of network
node_exist_idx1 <- pair1$peak_rank %in% node_file$node_id # 4300384/4394108
node_exist_idx2 <- pair2$peak_rank %in% node_file$node_id # 4382911/4394108
node_exist_idx <- node_exist_idx1 & node_exist_idx2 # 4293995/4394108
pair1_input <- pair1[node_exist_idx, ]
pair2_input <- pair2[node_exist_idx, ]
# make sure the links between node and itself were removed 
self_idx <- pair1_input$peak_rank != pair2_input$peak_rank # 19266
pair1_input <- pair1_input[self_idx, ] # 4274729
pair2_input <- pair2_input[self_idx, ]
table(pair1_input$peak_type)
#       G
# 4274729
table(pair2_input$peak_type)
#       G      RE
#  306706 3968023
# you do not need to search for direct interactions, so remove such pairs in the search list
pair_input_rank1 <- paste(pair1_input$peak_rank, pair2_input$peak_rank, sep="_")
pair_input_rank2 <- paste(pair2_input$peak_rank, pair1_input$peak_rank, sep="_")
direct_rank <- paste(direct_mergeATACinfo_rmNA_fullinfo$queryHits, direct_mergeATACinfo_rmNA_fullinfo$subjectHits, sep="_") # 1210177
direct_idx <- pair_input_rank1 %in% direct_rank | pair_input_rank2 %in% direct_rank # 1169370/4274729

pair1_input_rmDirect <- pair1_input[!direct_idx, ] # 3,105,359
pair2_input_rmDirect <- pair2_input[!direct_idx, ] # 3105359
# note that the input is the node index not node itself
pair1_input_rmDirect$peak_rank_idx <- match(pair1_input_rmDirect$peak_rank, node_file$node_id)
pair2_input_rmDirect$peak_rank_idx <- match(pair2_input_rmDirect$peak_rank, node_file$node_id)
sum(is.na(pair2_input_rmDirect[, 1]))
# [1] 0

#
# for loop for all genes
# 
interested_gene_node <- unique(pair1_input$peak_rank) # 19266 G # This may cause length of some of path is 1 in the shortest paths
loop_idx <- 1:length(interested_gene_node)
cores=8
ptm <- proc.time()
heme_shortest_paths_weighted <- parallel::mclapply(loop_idx, mc.cores = cores, function(i){

      interested_node_full <- c(interested_gene_node[i], pair2_input$peak_rank[pair1_input$peak_rank == interested_gene_node[i]]) # gene and its peaks (direct peaks removed)
      interested_node <- c(interested_gene_node[i], pair2_input_rmDirect$peak_rank[pair1_input_rmDirect$peak_rank == interested_gene_node[i]]) # gene and its peaks (direct peaks removed)
      
      interested_node_full_idx <- node_file$node_id %in% interested_node_full #228
      interested_link_full_idx <- link_file$from %in% interested_node_full & link_file$to %in% interested_node_full #10834
    #   interested_node_idx <- node_file$node_id %in% interested_node #205
    #   interested_link_idx <- link_file$from %in% interested_node & link_file$to %in% interested_node #13450

      subnet <- graph_from_data_frame(d=link_file[interested_link_full_idx, ], vertices=node_file[interested_node_full_idx, ], directed=F)
      heme_shortest_paths_weighted <- get.shortest.paths(subnet, match(interested_gene_node[i], names(V(subnet))), match(pair2_input_rmDirect$peak_rank[pair1_input_rmDirect$peak_rank == interested_gene_node[i]], names(V(subnet))))
      print(i)
      return(heme_shortest_paths_weighted$vpath)
  }
  )

#
# calculate signal geometric average of path
#
peak_sample_mat <- as.data.frame(data.table::fread("/broad/sankaranlab/fyu/HemeMap/data/atacData/10142020-peak_sample_mat_pro.txt", header=T)) # dddddddddi
# peak_sample_mat <- subset(peak_sample_mat, select = -c(mDC, Mega)) # remove two cell types
# normalize count tables for atac - 18 cell types
peak_sample_mat_cpm <- edgeR::cpm(peak_sample_mat)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# did not remove 2315 self interactions
atac_signal_geo_mean_df <- rbindlist(
  lapply(
  heme_shortest_paths_weighted, function(x){
    data.frame(t(
      sapply(x, function(y){
        if(length(y)==1){
          temp <- peak_sample_mat_cpm[as.numeric(names(y)), ]
        # temp <- rep(0, ncol(peak_sample_mat_cpm))
        } else {
          temp <- apply(peak_sample_mat_cpm[as.numeric(names(y)), ], 2, gm_mean)
        }
        return(temp)
    })
    ))
  }
)
)

#
# path name (the first and last node name)
#
interaction_name <- unlist(
  parallel::mclapply(1:length(heme_shortest_paths_weighted), mc.cores = 8, function(x){
    sapply(heme_shortest_paths_weighted[[x]], function(y){
      temp <- paste(names(y)[1], names(y)[length(y)], sep="_")
      return(temp)
    })
  }
  )
) # 3105359
type_indirect <- rep("indirect", length(interaction_name))
indirect_gene <- unlist(sapply(heme_shortest_paths_weighted, function(y){
                            sapply(y, function(x) return(names(x)[1]) )  
                        }))
length(unique(indirect_gene))
# [1] 19733
#---------------- find the gene name for peak_obj ----------------#
peak_obj <- as.data.frame(data.table::fread("/broad/sankaranlab/fyu/HemeMap/data/atacData/10142020-atac_peak_pro.txt")) #450920
peak_obj <- data.frame(peak_obj, paste("peak", seq(1:nrow(peak_obj)), sep="_"))
colnames(peak_obj) <- c("chr", "start", "end", "peak_type", "peak_rank") # no M
tss250 <- data.frame(fread("/broad/sankaranlab/fyu/HemeMap/data/exprData/10142020-tss250.txt", header=F)) # tss of genes; 20174
colnames(tss250) <- c("chr", "start", "end", "gene")
tss250_2 <- tss250[!duplicated(tss250[, 1:3]),] #20130
peak_gene_name <- rep(0, length(peak_obj$peak_type))
peak_gene_name[peak_obj$peak_type=="G"] <- 1:sum(peak_obj$peak_type=="G")
peaks_g <- makeGRangesFromDataFrame(peak_obj, keep.extra.columns = TRUE)
tss_g <- makeGRangesFromDataFrame(tss250_2, keep.extra.columns = TRUE)
ol <- findOverlaps(peaks_g, tss_g, select="first") #note some genes in the peak_obj overlapping
peak_gene_name[!is.na(ol)] <- tss250_2$gene[na.omit(ol)]
peak_obj <- data.frame(peak_obj, peak_gene_name)
write.table(peak_obj, file = "/broad/sankaranlab/fyu/HemeMap/MECOM-HemeMap2/data/peak_obj_genename.tsv", 
            row.names = FALSE, col.names = T, sep = "\t", quote = FALSE)
#----------------------------------------------------------------#
indirect_gene_name <- peak_obj$peak_gene_name[as.numeric(indirect_gene)]
sum(indirect_gene_name=="0") # no gene name which were duplicate 
# [1] 652
indirect_df <- data.frame(interaction_name=interaction_name, type=type_indirect, gene=indirect_gene_name)

#
# direct interaction
#
direct_mergeATACinfo_rmNA_fullinfo <- data.frame(fread("/broad/sankaranlab/fyu/HemeMap/MECOM-HemeMap2/data/direct_mergeATACinfo_rmNA_fullinfo.tsv"))   

sum(direct_mergeATACinfo_rmNA_fullinfo$queryHits> direct_mergeATACinfo_rmNA_fullinfo$subjectHits)
# 0
interaction_name_direct <- paste(direct_mergeATACinfo_rmNA_fullinfo$queryHits, direct_mergeATACinfo_rmNA_fullinfo$subjectHits, sep="_") # 1210177
# interaction_name <- parallel::mclapply(1:nrow(direct_mergeATACinfo_rmNA_fullinfo), mc.cores = 8, function(x){
atac_signal_geo_mean_df_direct <- parallel::mclapply(1:nrow(direct_mergeATACinfo_rmNA_fullinfo), mc.cores = 8, function(x){
        temp <- apply(peak_sample_mat_cpm[c(direct_mergeATACinfo_rmNA_fullinfo$queryHits[x], direct_mergeATACinfo_rmNA_fullinfo$subjectHits[x]), ], 2, gm_mean)
        return(temp)
    })
atac_signal_geo_mean_df_direct2 <- as.data.frame(do.call(rbind, atac_signal_geo_mean_df_direct))
type_direct <- rep("direct", length(interaction_name_direct))
direct_df <- data.frame(interaction_name=interaction_name_direct, type=type_direct, gene=direct_mergeATACinfo_rmNA_fullinfo$gene)

#
# merge the direct and indirect interaction
#

interaction_df <- rbind.data.frame(indirect_df, direct_df) # 4315536
interaction_geoATACmean18 <- rbind.data.frame(atac_signal_geo_mean_df, atac_signal_geo_mean_df_direct2) # 4315536

save(interaction_df, interaction_geoATACmean18, file="mecom_var_sankaran/data/hememap/01022021-heme_shortest_paths_geoActivate_info.rda") # large, write it by yourself

# ----------------------
# test if the activation fit for chi-square distribution
ks.test(interaction_geoATACmean18$HSC, "pchisq", mean(interaction_geoATACmean18$HSC), sd(interaction_geoATACmean18$HSC))
# estimate the p value when use chi-square distribution
p_hsc_activation <-  pchisq(interaction_geoATACmean18$HSC, df=mean(interaction_geoATACmean18$HSC), lower.tail=FALSE) # 3.4 # 4,315,536
sum(p_hsc_activation<=0.05)
# [1] 372,491
sum(p_hsc_activation<=0.05)/length(p_hsc_activation)
# [1] 0.08631396
interaction_geoATACmean18$HSC[which.min(abs(p_hsc_activation-0.05))] # 0.049945 # 8.5
# HSC activation network
hsc_activation_idx <- p_hsc_activation<=0.05 # 372491  ## 
# hsc_activation_idx <- p_hsc_activation<=0.01 # 190178  ## 
interaction_df_hsc <- interaction_df[hsc_activation_idx, ] # HSC activation netwrok
length(unique(interaction_df_hsc$gene)) # how many gene involved
# [1] 12808


# ----------------------
# mecom regulated gene and cisREs
load("mecom_var_sankaran/data/signature_exp_enrichment/normalized_counts_deg.rda")
# normalized_counts_down, normalized_counts_up, normalized_counts_all, 
deg_mast0035_down <- rownames(normalized_counts_down) # 307
deg_mast0035_up <- rownames(normalized_counts_up) # 385
deg_mast0035_all <- rownames(normalized_counts_all) # 692
length(intersect(deg_mast0035_down, as.character(unique(interaction_df_hsc$gene)))) # how many genes intersect bewteen interesting and hsc network
# 228/322
length(intersect(deg_mast0035_up, as.character(unique(interaction_df_hsc$gene)))) # how many genes intersect bewteen interesting and hsc network
# 277/402
length(intersect(deg_mast0035_all, as.character(unique(interaction_df_hsc$gene)))) # how many genes intersect bewteen interesting and hsc network
# 505/724
interested_gene_list_down <- intersect(deg_mast0035_down, as.character(unique(interaction_df_hsc$gene)))
interested_gene_list_up <- intersect(deg_mast0035_up, as.character(unique(interaction_df_hsc$gene)))
interested_gene_list_all <- intersect(deg_mast0035_all, as.character(unique(interaction_df_hsc$gene)))
interested_gene_peak_idx_down <- interaction_df_hsc$gene %in% interested_gene_list_down # 7274
interested_gene_peak_idx_up <- interaction_df_hsc$gene %in% interested_gene_list_up # 8312
interested_gene_peak_idx_all <- interaction_df_hsc$gene %in% interested_gene_list_all # 15586

setwd("mecom_var_sankaran/data/down0035_deg_info")
write.table(interested_gene_list_down, "deg_mast0035_down_HemeHSC_gene.txt", row.names=F, col.names=F, sep="\t", quote=F)
write.table(interested_gene_list_up, "deg_mast0035_up_HemeHSC_gene.txt", row.names=F, col.names=F, sep="\t", quote=F)
write.table(interested_gene_list_all, "deg_mast0035_all_HemeHSC_gene.txt", row.names=F, col.names=F, sep="\t", quote=F)
# ----------------------
load("mecom_var_sankaran/data/hememap/interaction_name_gene_info.rda") # large, write it by yourself
##down
temp_interaction_hsc <- str_split(interaction_df_hsc[interested_gene_peak_idx_down, ]$interaction_name, "_", simplify = T)
temp_node_hsc <- as.numeric(unique(c(temp_interaction_hsc[, 1], temp_interaction_hsc[, 2]))) # 6932
peak_node_hsc <- peak_obj[as.numeric(temp_node_hsc), ]
table(peak_node_hsc[, 4]) # 6932 -> 6605
#    G   RE
# 1116 5816
write.table(peak_node_hsc, "down_mast0035_interactions_peak4.bed", row.names=F, col.names=F, sep="\t", quote=F)
##up
temp_interaction_hsc <- str_split(interaction_df_hsc[interested_gene_peak_idx_up, ]$interaction_name, "_", simplify = T)
temp_node_hsc <- as.numeric(unique(c(temp_interaction_hsc[, 1], temp_interaction_hsc[, 2]))) # 
peak_node_hsc <- peak_obj[as.numeric(temp_node_hsc), ]
table(peak_node_hsc[, 4]) # 7128
#    G   RE
# 1486 5642
write.table(peak_node_hsc, "up_mast0035_interactions_peak4.bed", row.names=F, col.names=F, sep="\t", quote=F)
##all
temp_interaction_hsc <- str_split(interaction_df_hsc[interested_gene_peak_idx_all, ]$interaction_name, "_", simplify = T)
temp_node_hsc <- as.numeric(unique(c(temp_interaction_hsc[, 1], temp_interaction_hsc[, 2]))) # 
peak_node_hsc <- peak_obj[as.numeric(temp_node_hsc), ]
table(peak_node_hsc[, 4]) # 13184
#     G    RE
#  2476 10708
write.table(peak_node_hsc, "all_mast0035_interactions_peak4.bed", row.names=F, col.names=F, sep="\t", quote=F)
