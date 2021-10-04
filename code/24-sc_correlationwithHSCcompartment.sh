######### 
# fold change correlation for HSC compartment and all cells
combineddeg <- read.csv("/Users/fyu/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vijay/project/09022020-MECOM_Richard/data/20210801-cellrevision/all_cells_vs_lthsc.csv")
# 2040   13
sum(is.na(combineddeg$avg_logFC_x) | is.na(combineddeg$avg_logFC_y)) # 1316
combineddeg <- combineddeg[!(is.na(combineddeg$avg_logFC_x) | is.na(combineddeg$avg_logFC_y)), ]
colnames(combineddeg)[2] <- "genename"
top5gene <- c(combineddeg[order(combineddeg$avg_logFC_x), 2][1:5], rev(combineddeg[order(combineddeg$avg_logFC_x), 2])[1:5])
top5gene2 <- c(combineddeg[order(combineddeg$avg_logFC_y), 2][1:5], rev(combineddeg[order(combineddeg$avg_logFC_y), 2])[1:5])
eithertop5 <- filter(combineddeg, (genename  %in% c(top5gene, "MECOM")) | (genename  %in% top5gene2))
# [1] 537  13
ggplot(combineddeg, aes(avg_logFC_x, avg_logFC_y))+
        geom_point() + theme_classic(base_size = 10) + xlim(-0.22, 0.4)  + xlab("log2(fold change) from all cells") +ylab("log2(fold change) from HSC compartment") +  #ylim(-1.7, 1) +
        ggrepel::geom_text_repel(data=eithertop5, aes(avg_logFC_x, avg_logFC_y, label = genename), color = 'black',
                        size = 2, max.overlaps=30, segment.linetype = 1)
setwd("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/09022020-MECOM_Richard/data/20210801-cellrevision")
ggsave("foldchange_comparison_allcellandHSC.pdf", width = 5, height = 5)


cor.test(combineddeg$avg_logFC_x, combineddeg$avg_logFC_y, method="spearman")
#         Spearman's rank correlation rho

# data:  combineddeg$avg_logFC_x and combineddeg$avg_logFC_y
# S = 3544866, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.8626497 
