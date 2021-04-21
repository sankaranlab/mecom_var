require(tidyverse)
require(annotables)
require(grImport)
require(BuenColors)
require(cowplot)
require(edgeR)


# get the interested gene list 
setwd("/Users/fyu/Documents/GitHub")
load("mecom_var/data/signature_exp_enrichment/normalized_counts_deg.rda")
load("mecom_var/data/signature_exp_enrichment/normalized_counts.rda")

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
# --- functions to be loaded first
enrich_s <- 
names(enrich_s) <- c("B", "CD4", "CD8", "CLP", "CMP", "ERY", "GMPA", "GMPB", "GMPC", "HSC", "LMPP", "MEGA", "MEP", "MONO", "MPP", "NK", "GRAN", "MDC", "PDC", "PLT")
blood_cell_plot <- function(enrich_s, palette = "solar_rojos", output_dir = ".",
                            type = "log2_cpm",
                            pics_xml = "./data/blank_cells.eps.xml"){
   
   genes = genes
   n <- length(genes)
   df <- import_data(type = type)
   df <- anno_grch38(df)
   df <- filter_genes_long(df, genes)
   df <- generate_plot_info(df, palette, n, type = type)[[1]] %>% distinct()
   legend <- generate_plot_info(df, palette, n, type = type)[[2]]
   for(symb in unique(df$symbol)){
      curr <- df[df$symbol == symb)
      my_shape <- color_cells(myshape, curr) %>% pictureGrob(.) %>% ggdraw(.)
      p <- my_shape + annotation_custom(legend, xmin =.75, ymin = .65)
      ggsave(file.path(output_dir,
                       paste0(symb, "_bloodEnrich.pdf")),
             p, width = 8, height = 6)
   }
}

# generate color graphic information using ggplot
generate_plot_info <- function(df, palette, n_genes, type){
#    palette = jdb_palette(palette)
#  palette = jdb_palette(palette)[c(1:5,9)]
   palette = jdb_palette(palette)[c(1:7,9)]
   df <- data.frame(df)
   if(n_genes == 1){
      p1 <- ggplot(df, aes(x = cell_type, y = expression, fill = expression)) +
         geom_bar(stat = "identity") + 
         coord_flip() +
         pretty_plot() + 
         scale_fill_gradientn(colors = palette)
   }else{
      p1 <- ggplot(df, aes(x = cell_type, y = expression, fill = expression)) +
         geom_bar(stat = "identity") + 
         coord_flip() +
         facet_wrap(~factor(ensg), ncol = n_genes) +
         pretty_plot() + 
         scale_fill_gradientn(colors = palette)
   }
   p1 <- p1 + labs(fill = type)
   df <- ggplot_build(p1)[[1]] %>%
      data.frame(.) %>% 
      dplyr::select(1,3) %>%
      left_join(df, ., by = c("expression" = "y"))
   legend <- get_legend(p1)
   list(df, legend)
}



# function to change colors of cells
color_cells <- function(picture, df){
   my_shape = picture
   B <- df[df$cell_type == "B", "fill"]
   CD4 <- df[df$cell_type == "CD4", "fill"]
   CD8 <- df[df$cell_type == "CD8", "fill"]
   CLP <- df[df$cell_type == "CLP", "fill"]
   CMP <- df[df$cell_type == "CMP", "fill"]
   ERY <- df[df$cell_type == "ERY", "fill"]
   GMPA <- df[df$cell_type == "GMPA", "fill"]
   GMPB <- df[df$cell_type == "GMPB", "fill"]
   GMPC <- df[df$cell_type == "GMPC", "fill"]
   HSC <- df[df$cell_type == "HSC", "fill"]
   LMPP <- df[df$cell_type == "LMPP", "fill"]
   MEGA <- df[df$cell_type == "MEGA", "fill"]
   MEP <- df[df$cell_type == "MEP", "fill"]
   MONO <- df[df$cell_type == "MONO", "fill"]
   MPP <- df[df$cell_type == "MPP", "fill"]
   NK <- df[df$cell_type == "NK", "fill"]
   GRAN <- df[df$cell_type == "GRAN", "fill"]
   MDC <- df[df$cell_type == "MDC", "fill"]
   PDC <- df[df$cell_type == "PDC", "fill"]
   PLT <- df[df$cell_type == "PLT", "fill"]
   
   # B cell
   my_shape@paths[59]$path@rgb <- B
   my_shape@paths[181]$path@rgb <- B
   
   #CD4  53, 177
   my_shape@paths[53]$path@rgb <- CD4
   my_shape@paths[177]$path@rgb <- CD4
   
   #CD8 56, 179 
   my_shape@paths[56]$path@rgb <- CD8
   my_shape@paths[179]$path@rgb <- CD8
   
   #CLP  44, 169
   my_shape@paths[44]$path@rgb <- CLP
   my_shape@paths[169]$path@rgb <- CLP
   
   #CMP   46, 171
   my_shape@paths[46]$path@rgb <- CMP
   my_shape@paths[171]$path@rgb <- CMP
   
   #Ery   68, 70, 72, 187, 205, 207
   my_shape@paths[68]$path@rgb <- ERY
   my_shape@paths[70]$path@rgb <- ERY
   my_shape@paths[72]$path@rgb <- ERY
   my_shape@paths[187]$path@rgb <- ERY
   my_shape@paths[205]$path@rgb <- ERY
   my_shape@paths[207]$path@rgb <- ERY
   
   #GMP-A 97, 173
   my_shape@paths[97]$path@rgb <- GMPA
   my_shape@paths[173]$path@rgb <- GMPA
   
   #GMP-B 84, 197
   my_shape@paths[84]$path@rgb <- GMPB
   my_shape@paths[197]$path@rgb <- GMPB
   
   # GMP-C 48, 199
   my_shape@paths[48]$path@rgb <- GMPC
   my_shape@paths[199]$path@rgb <- GMPC
   
   # GRAN 121, 123, 125, 127, 129, 131, 133, 135, 137, 139, 141,
   #       143, 145, 147, 149, 151, 153
   my_shape@paths[121]$path@rgb <- GRAN
   my_shape@paths[123]$path@rgb <- GRAN
   my_shape@paths[125]$path@rgb <- GRAN
   my_shape@paths[127]$path@rgb <- GRAN
   my_shape@paths[129]$path@rgb <- GRAN
   my_shape@paths[131]$path@rgb <- GRAN
   my_shape@paths[133]$path@rgb <- GRAN
   my_shape@paths[135]$path@rgb <- GRAN
   my_shape@paths[137]$path@rgb <- GRAN
   my_shape@paths[139]$path@rgb <- GRAN
   my_shape@paths[141]$path@rgb <- GRAN
   my_shape@paths[145]$path@rgb <- GRAN
   my_shape@paths[147]$path@rgb <- GRAN
   my_shape@paths[149]$path@rgb <- GRAN
   my_shape@paths[143]$path@rgb <- GRAN
   my_shape@paths[151]$path@rgb <- GRAN
   my_shape@paths[153]$path@rgb <- GRAN
   
   # HSC   25, 163
   my_shape@paths[25]$path@rgb <- HSC
   my_shape@paths[163]$path@rgb <- HSC
   
   # LMPP  31, 167
   my_shape@paths[31]$path@rgb <- LMPP
   my_shape@paths[167]$path@rgb <- LMPP
   
   # mDC   115, 201
   my_shape@paths[115]$path@rgb <- MDC
   my_shape@paths[201]$path@rgb <- MDC
   
   # mega  1, 155, 157, 159, 161
   my_shape@paths[1]$path@rgb <- MEGA
   my_shape@paths[157]$path@rgb <- MEGA
   my_shape@paths[159]$path@rgb <- MEGA
   my_shape@paths[161]$path@rgb <- MEGA
   my_shape@paths[155]$path@rgb <- MEGA
   
   # MEP   50, 175
   my_shape@paths[50]$path@rgb <- MEP
   my_shape@paths[175]$path@rgb <- MEP
   
   # Mono  65, 185
   my_shape@paths[65]$path@rgb <- MONO
   my_shape@paths[185]$path@rgb <- MONO
   
   # MPP   27, 165
   my_shape@paths[27]$path@rgb <- MPP
   my_shape@paths[165]$path@rgb <- MPP
   
   # NK    62, 183
   my_shape@paths[62]$path@rgb <- NK
   my_shape@paths[183]$path@rgb <- NK
   
   # pDC   111, 203
   my_shape@paths[111]$path@rgb <- PDC
   my_shape@paths[203]$path@rgb <- PDC
   
   # Plt   190, 192, 194, 217, 219, 221, 223, 225
   my_shape@paths[190]$path@rgb <- PLT
   my_shape@paths[192]$path@rgb <- PLT
   my_shape@paths[194]$path@rgb <- PLT
   my_shape@paths[217]$path@rgb <- PLT
   my_shape@paths[219]$path@rgb <- PLT
   my_shape@paths[221]$path@rgb <- PLT
   my_shape@paths[223]$path@rgb <- PLT
   my_shape@paths[225]$path@rgb <- PLT
   return(my_shape)
}   


# -------- -------- -------- -------- -------- -------- --------
# -------- -------- -------- -------- -------- -------- --------
# run this code
df <- import_data(type = "log2_cpm") # import expr count mat and convert it to log2 cpm value
df <- anno_grch38(df) # add gene symbol follows the column of ensgene in the df
df <- as.data.frame(df) 
# set.seed(9527)
# down genes
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
#         B       CD4       CD8       CLP       CMP       ERY      GMPA      GMPB      GMPC      GRAN       HSC      LMPP       LSC 
# 0.8536257 0.7845451 0.7854256 0.9432778 0.9288675 0.5747249 0.9752541 1.0133017 1.0288096 0.5688033 1.0828431 1.0351751 0.9282326 
#       MDC      MEGA       MEP      MONO       MPP        NK       PDC       PLT 
# 0.8645416 0.8710295 0.8353811 0.9628376 1.0700109 0.7261047 0.7867739 0.5083076 

z_score[11]=1.11
generate_plot_info <- function(df, palette, n_genes, type){
#    palette = jdb_palette(palette)
    palette = jdb_palette(palette)[c(1:7,9)]

   df <- data.frame(df)
   if(n_genes == 1){
      p1 <- ggplot(df, aes(x = cell_type, y = expression, fill = expression)) +
         geom_bar(stat = "identity") + 
         coord_flip() +
         pretty_plot() + 
         scale_fill_gradientn(colors = palette)
   }else{
      p1 <- ggplot(df, aes(x = cell_type, y = expression, fill = expression)) +
         geom_bar(stat = "identity") + 
         coord_flip() +
         facet_wrap(~factor(ensg), ncol = n_genes) +
         pretty_plot() + 
         scale_fill_gradientn(colors = palette)
   }
   p1 <- p1 + labs(fill = type)
   df <- ggplot_build(p1)[[1]] %>%
      data.frame(.) %>% 
      dplyr::select(1,3) %>%
      left_join(df, ., by = c("expression" = "y"))
   legend <- get_legend(p1)
   list(df, legend)
}

p1 <- dataset20cell_expression(z_score, type="Enrichment Z-score")
setwd("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/09022020-MECOM_Richard/data/02262021-downMast_re-analysis/10selected_gene_vis")
ggsave(p1, filename = paste0("mecomDOWN_enrichment_bgnorm_update.pdf"), height = 8, width = 10)

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
#         B       CD4       CD8       CLP       CMP       ERY      GMPA      GMPB      GMPC      GRAN       HSC      LMPP       LSC 
# 0.7279153 0.7802504 0.7902271 0.6983152 1.1159091 0.7110103 1.2427163 1.2440327 1.1509442 0.7844958 0.8255038 0.9661857 0.7241834 
#       MDC      MEGA       MEP      MONO       MPP        NK       PDC       PLT 
# 0.9840710 1.1216498 1.1041458 0.8973525 0.9579781 0.7817102 0.9685910 0.4994506
generate_plot_info <- function(df, palette, n_genes, type){
#    palette = jdb_palette(palette)
    palette = jdb_palette(palette)[c(1:7,9)]

   df <- data.frame(df)
   if(n_genes == 1){
      p1 <- ggplot(df, aes(x = cell_type, y = expression, fill = expression)) +
         geom_bar(stat = "identity") + 
         coord_flip() +
         pretty_plot() + 
         scale_fill_gradientn(colors = palette)
   }else{
      p1 <- ggplot(df, aes(x = cell_type, y = expression, fill = expression)) +
         geom_bar(stat = "identity") + 
         coord_flip() +
         facet_wrap(~factor(ensg), ncol = n_genes) +
         pretty_plot() + 
         scale_fill_gradientn(colors = palette)
   }
   p1 <- p1 + labs(fill = type)
   df <- ggplot_build(p1)[[1]] %>%
      data.frame(.) %>% 
      dplyr::select(1,3) %>%
      left_join(df, ., by = c("expression" = "y"))
   legend <- get_legend(p1)
   list(df, legend)
}

p1 <- dataset20cell_expression(z_score, type="Enrichment Z-score")
setwd("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/09022020-MECOM_Richard/data/02262021-downMast_re-analysis/10selected_gene_vis")
ggsave(p1, filename = paste0("mecomUP_enrichment_bgnorm_update.pdf"), height = 8, width = 10)




