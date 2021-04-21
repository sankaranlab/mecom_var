require(tidyverse)
require(annotables)
require(grImport)
require(dplyr)
require(edgeR)

# get the interested gene list 
setwd("/Users/fyu/Documents/GitHub")
load("mecom_var_sankaran/data/signature_exp_enrichment/normalized_counts_deg.rda")
load("mecom_var_sankaran/data/signature_exp_enrichment/normalized_counts.rda")

# define funcitons
import_data <- function(type){
   df <- read_tsv("https://raw.githubusercontent.com/jeffverboon/blood_expressions_plots/master/data/counts.tsv")
   if (type == "cpm"){
       df[-1] <- lapply(df[-1], cpm)
   }
   if (type == "log2_cpm"){
       df[-1] <- lapply(df[-1], cpm)
       df[-1] <- lapply(df[-1], function(x) log2(x + 1))
   }
   if(type == "log2"){
    #   df[-1] <- lapply(df[-1], cpm)
         df[-1] <- lapply(df[-1], function(x) log2(x + 1))
      }
   df
}
anno_grch38 <- function(df){
   annotables::grch38 %>%
      mutate(ensg = ensgene) %>%
      dplyr::select(ensg, symbol) %>%
      right_join(., df)
}


# run this code
df <- import_data(type = "log2_cpm") # import expr count mat and convert it to log2 cpm value
df <- anno_grch38(df) # add gene symbol follows the column of ensgene in the df
df <- as.data.frame(df) 
# set.seed(9527)
# 
# down genes
# 
genes <- rownames(normalized_counts_down)
newdf <- df[df$symbol %in% genes, -(1:2) ]
background_df <- df[df$symbol %in% tss250[, 4], -(1:2) ]

obj_exp <- colMeans(newdf)
exp_exp_df <- data.frame(matrix(0, nrow=1000, ncol=length(obj_exp)))
for (i in 1:1000){
    sample_newdf <- background_df[sample(1:nrow(background_df), nrow(newdf)), ]
    exp_exp_df[i, ] <- colMeans(sample_newdf)
}
exp_exp <- colMeans(sample_newdf)
exp_sd <- apply(sample_newdf, 2, sd)

z_score <- (obj_exp-exp_exp)/exp_sd
z_score
#         B       CD4       CD8       CLP       CMP       ERY      GMPA      GMPB      GMPC      GRAN       HSC      LMPP       LSC 
# 0.8536257 0.7845451 0.7854256 0.9432778 0.9288675 0.5747249 0.9752541 1.0133017 1.0288096 0.5688033 1.128431 1.0351751 0.9282326 
#       MDC      MEGA       MEP      MONO       MPP        NK       PDC       PLT 
# 0.8645416 0.8710295 0.8353811 0.9628376 1.0700109 0.7261047 0.7867739 0.5083076 

# 
# up genes
# 
genes <- rownames(normalized_counts_up)
newdf <- df[df$symbol %in% genes, -(1:2) ]
background_df <- df[df$symbol %in% tss250[, 4], -(1:2) ]

obj_exp <- colMeans(newdf)
exp_exp_df <- data.frame(matrix(0, nrow=1000, ncol=length(obj_exp)))
for (i in 1:1000){
    sample_newdf <- background_df[sample(1:nrow(background_df), nrow(newdf)), ]
    exp_exp_df[i, ] <- colMeans(sample_newdf)
}
exp_exp <- colMeans(sample_newdf)
exp_sd <- apply(sample_newdf, 2, sd)
z_score <- (obj_exp-exp_exp)/exp_sd
z_score
#         B       CD4       CD8       CLP       CMP       ERY      GMPA      GMPB      GMPC      GRAN       HSC      LMPP       LSC 
# 0.7279153 0.7802504 0.7902271 0.6983152 1.1159091 0.7110103 1.2427163 1.2440327 1.1509442 0.7844958 0.8255038 0.9661857 0.7241834 
#       MDC      MEGA       MEP      MONO       MPP        NK       PDC       PLT 
# 0.9840710 1.1216498 1.1041458 0.8973525 0.9579781 0.7817102 0.9685910 0.4994506



